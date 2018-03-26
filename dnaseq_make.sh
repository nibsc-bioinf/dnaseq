#!/bin/bash

if [ -z "$2" ]; then
  echo "Usage: $0 <type> <projectid>"
  exit
fi

type=$1
project=$2
dir="/sequencing/projects/$project/"
executable="/usr/local/mark/"
#executable="$dir/scripts/"
threads=8

cd $dir
echo `pwd`
rm -f dnaseq.$type.*

#  Write command to output files
#    Assume dir, project, type already set.
function writeCommand {
  command=$1
  date=$2
  r=$3
  t1=$4
  t2=$5
  case "$command" in
    index)
      hold="0"
      name="dnaseq.index.$project.$r"
      shell="/bin/bash"
      commandLine="$executable/dnaseq_index.sh $dir $r"
      log="$r"
      ;;
    map)
      hold="dnaseq.index.$project.$r"
      name="dnaseq.map.$project"
      shell="/bin/bash"
      commandLine="$executable/dnaseq_map.sh $dir $r $t1 ${date}_$t2 $threads"
      log="$date.$r.$t2"
      ;;
    report)
      hold="dnaseq.map.$project"
      name="dnaseq.report.$project"
      shell="/usr/bin/Rscript"
      commandLine="$executable/dnaseq_report.r $dir $project $date"
      log="$date"
      ;;
    *)
      echo "Bad: $command"
      exit 1
  esac
  err="$dir/log/qsub.err.$command.$log"
  log="$dir/log/qsub.log.$command.$log"
  case "$type" in
    qsub)
      echo qsub -N "$name" -hold_jid "$hold" -pe parallel "$threads" -cwd -S "$shell" -o "$log" -e "$err" "$commandLine" >> dnaseq.$type.$project.sh
      ;;
    bash)
      echo "$commandLine > $log 2> $err" >> dnaseq.$type.$project.$command.sh
      ;;
    *)
      echo "Bad: $type"
      exit 1
  esac
}

#  Indices
for r in `cut -f1 input/references | grep -v "^#"`
do
  writeCommand "index" "X" $r # X = no date
done

#  Alignments (mapping) and reports
for date in `ls input/samples.* | cut -b15-`
do
  for r in `cut -f1 input/references | grep -v "^#"`
  do
    for s in `cut -f2 input/samples.$date | grep -v "^#" | sed 's/,/\n/g'`
    do
      t1=`echo $s  | sed 's/.p1.fastq.gz//'`
      t2=`echo $t1 | cut -f4 -d'/'`
      writeCommand "map" $date $r $t1 $t2
    done
  done
  writeCommand "report" $date
done

chmod +x dnaseq.$type.*
