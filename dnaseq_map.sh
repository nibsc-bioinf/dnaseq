#!/bin/bash

date

trim="30"

dir=$1
reference=$2
sample=$3
samplename=$4
threads=$5

if [ "x$1" == "x" -o "x$2" == "x" -o "x$3" == "x" -o "x$4" == "x" -o "x$5" == "x" ]; then
  echo "Usage: $0 dir reference sample samplename threads"
  exit
fi

fasta="$dir/input/reference/$reference.fasta"
input="$dir/input/$sample"
output="$dir/processed/$samplename"
log="$output/$trim.bwa.$reference.log"
samp="$output/$trim.bwa.$reference"

mkdir -p $output
cd $output
rm -f $log

#---------------------------------------------------------------------------------------------------
echo "MDP: BWA" >> $log
bwa mem -t $threads $fasta $input.p1.fastq.gz $input.p2.fastq.gz > $samp.p.sam 2>> $log
bwa mem -t $threads $fasta $input.s1.fastq.gz > $samp.s1.sam                   2>> $log
bwa mem -t $threads $fasta $input.s2.fastq.gz > $samp.s2.sam                   2>> $log
java -jar /usr/local/share/picard/picard.jar MergeSamFiles I=$samp.p.sam I=$samp.s1.sam I=$samp.s2.sam O=$samp.sorted.bam SO=coordinate 2>> $log

#---------------------------------------------------------------------------------------------------
echo "MDP: PICARD" >> $log
java -jar /usr/local/share/picard/picard.jar AddOrReplaceReadGroups I=$samp.sorted.bam    O=$samp.readgroup.bam LB="Library" PU="PlatformUnit" SM="$samp" PL="Illumina" 2>> $log
java -jar /usr/local/share/picard/picard.jar MarkDuplicates         I=$samp.readgroup.bam O=$samp.bam           METRICS_FILE=$samp.dups ASSUME_SORTED=T                 2>> $log
java -jar /usr/local/share/picard/picard.jar BuildBamIndex          I=$samp.bam                                                                                         2>> $log
sleep 1

#---------------------------------------------------------------------------------------------------
echo "MDP: GATK" >> $log
java -jar /usr/local/share/gatk/GenomeAnalysisTK.jar -T DepthOfCoverage -I $samp.bam -R $fasta           -o $samp.gatk.depth  2>> $log
java -jar /usr/local/share/gatk/GenomeAnalysisTK.jar -T Pileup          -I $samp.bam -R $fasta           -o $samp.gatk.pileup 2>> $log
java -jar /usr/local/share/gatk/GenomeAnalysisTK.jar -T HaplotypeCaller -I $samp.bam -R $fasta -ploidy 1 -o $samp.gatk.vcf    2>> $log

#---------------------------------------------------------------------------------------------------
echo "MDP: SAMTOOLS" >> $log
samtools depth                             $samp.bam > $samp.samtools.depth           2>> $log
samtools stats -r $fasta                   $samp.bam > $samp.samtools.stats           2>> $log
samtools flagstat                          $samp.bam > $samp.samtools.flagstat        2>> $log
samtools mpileup -O -s -d 100000 -f $fasta $samp.bam > $samp.samtools.pileup          2>> $log
samtools mpileup -g    -d 100000 -f $fasta $samp.bam > $samp.samtools.bcf             2>> $log
bcftools call -c -v -Ov $samp.samtools.bcf           > $samp.samtools.vcf             2>> $log
/home/mark/cluster/software/pileup2calls.pl $fasta $samp $samp.samtools.pileup        2>> $log
#  /home/mark/cluster/software/calls2stats.r $directory $reference $samp $trim   2>> $log

rm -f $samp.p.sam $samp.s1.sam $samp.s2.sam $samp.mapped.bam $samp.sorted.bam $samp.readgroup.bam 2>> $log

date
