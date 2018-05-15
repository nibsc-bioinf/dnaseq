#!/usr/bin/perl
use warnings;
use strict;
use Bio::DB::Fasta;

die ("Usage: $0 reference sample pileup\n") if ($#ARGV != 2);
my ($psReference, $psSample, $psPileup) = @ARGV[0..2];

my $lsCalls = $psPileup =~ s/pileup/calls/r;
my $loReference = Bio::DB::Fasta->new($psReference);

my $lsChromosome = "";
my $liPosition = 0;
my $lsSequence = "";
my $liSequence = -1;

print "$psPileup\t$lsCalls\n";

open I, ">$lsCalls" or die("Cannot open output file: $lsCalls\n");
if (-s $psPileup) {
  open H, "<$psPileup" or die("Cannot open input file: $psPileup\n");
  while (<H>) {
    chomp;
    my @laLine = split /\t/;
    if ($lsChromosome ne $laLine[0]) {
      if ($lsChromosome ne "" && $laLine[1] < $liSequence) {
        for (my $i = $liPosition; $i < $liSequence; $i++) {
          printLine($psSample,$lsChromosome,$i+1,substr($lsSequence,$i,1),0,"");
        }
      }
      $lsChromosome = $laLine[0];
      $liPosition = $laLine[1];
      $lsSequence = $loReference->get_Seq_by_id($lsChromosome)->seq;
      $liSequence = length($lsSequence);
      print STDERR "$lsChromosome\t$liSequence\n";
      if ($laLine[1] > 1) {
        for (my $i = 0; $i < $laLine[1] - 1; $i++) {
          printLine($psSample,$laLine[0],$i+1,substr($lsSequence,$i,1),0,"");
        }
      }
    }
    if ($liPosition < $laLine[1] - 1) {
      for (my $i = $liPosition; $i < $laLine[1] - 1; $i++) {
        printLine($psSample,$laLine[0],$i+1,substr($lsSequence,$i,1),0,"");
      }
    }
    printLine($psSample,$laLine[0],$laLine[1],$laLine[2],$laLine[3],$laLine[4],$laLine[5]);
    $lsChromosome = $laLine[0];
    $liPosition = $laLine[1];
  }
  if ($liPosition < $liSequence) {
    for (my $i = $liPosition; $i < $liSequence; $i++) {
      printLine($psSample,$lsChromosome,$i+1,substr($lsSequence,$i,1),0,"");
    }
  }
  close H;
} else {
  for $lsChromosome ($loReference->get_all_primary_ids()) {
    $lsSequence = $loReference->get_Seq_by_id($lsChromosome)->seq;
    $liSequence = length($lsSequence);
    for (my $i = 0; $i < $liSequence; $i++) {
      printLine($psSample,$lsChromosome,$i+1,substr($lsSequence,$i,1),0,"");
    }
  }
}
close I;
exit;

sub printLine {
  my ($psSample, $psChromosome, $piPosition, $pcReference, $piCoverage, $psAlignment, $psQuality) = @_;
  my $liReferenceSkip = 1 * ($psAlignment =~ s/<>//g);
  my $liPreviousDel   = 1 * ($psAlignment =~ s/\*//g);
  my $liReadStart     = 1 * ($psAlignment =~ s/\^.//g);
  my $liReadEnd       = 1 * ($psAlignment =~ s/\$//g);
  my (@laInsertion)   = ($psAlignment =~ /\+[0-9]+/g);
  my (@laDeletion)    = ($psAlignment =~ /\-[0-9]+/g);
  my $liInsertion     = $#laInsertion + 1;
  my $liDeletion      = $#laDeletion + 1;
  for (@laInsertion) {
    $_ = 1*$_;
    $psAlignment =~ s/\+$_.{$_}//;
  }
  for (@laDeletion) {
    $_ = -1*$_;
    $psAlignment =~ s/\-$_.{$_}//;
  }
  my $liForward       = () = $psAlignment =~ /[\.ACGTN]/g;
  my $liReverse       = () = $psAlignment =~ /[,acgtn]/g;
  my $liReference     = () = $psAlignment =~ /[\.,]/g;
  $psAlignment =~ s/[\.,]/$pcReference/g;
  $psAlignment = uc($psAlignment);
  my $liA             = () = $psAlignment =~ /A/g;
  my $liC             = () = $psAlignment =~ /C/g;
  my $liG             = () = $psAlignment =~ /G/g;
  my $liT             = () = $psAlignment =~ /T/g;
  my $liN             = () = $psAlignment =~ /N/g;
  if ($piCoverage != length($psAlignment) + $liPreviousDel) {
    print STDERR "ERROR: $psChromosome\t$piPosition\t$pcReference\t$piCoverage\t". length($psAlignment) . "\t$liReference\t$liA\t$liC\t$liG\t$liT\t$liForward\t$liReverse\t$liInsertion\t$liDeletion\t$liReadStart\t$liReadEnd\t$liPreviousDel\t$liReferenceSkip\n";
    $piCoverage = length($psAlignment);
  }
  $piCoverage = length($psAlignment);
  my ($lrReference, $lrA, $lrC, $lrG, $lrT, $lrN) = (0, 0, 0, 0, 0, 0);
  if ($piCoverage > 0) {
    $lrReference = sprintf("%0.3f", $liReference / $piCoverage);
    $lrA         = sprintf("%0.3f", 1.0 * $liA / $piCoverage);
    $lrC         = sprintf("%0.3f", $liC / $piCoverage);
    $lrG         = sprintf("%0.3f", $liG / $piCoverage);
    $lrT         = sprintf("%0.3f", $liT / $piCoverage);
    $lrN         = sprintf("%0.3f", $liN / $piCoverage);
  }
  print I "$psSample\t$psChromosome\t$piPosition\t$pcReference\t$piCoverage\t$liReference\t$liA\t$liC\t$liG\t$liT\t$liN\t$lrReference\t$lrA\t$lrC\t$lrG\t$lrT\t$lrN\t$liForward\t$liReverse\t$liInsertion\t$liDeletion\t$liReadStart\t$liReadEnd\t$liPreviousDel\t$liReferenceSkip\n";
}

#  # Create database from a directory of Fasta files
#  my @ids      = $db->get_all_primary_ids;
#  my $seq     = $db->get_Seq_by_id('CHROMOSOME_I');
#  my $seqstr  = $seq->seq;
#  my $subseq  = $seq->subseq(4_000_000 => 4_100_000);
#  my $trunc   = $seq->trunc(4_000_000 => 4_100_000);
#  my $length  = $seq->length;
