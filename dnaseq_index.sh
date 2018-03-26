date

cd $1
mkdir -p log

reference=$2
echo $reference
bowtie2-build input/reference/$reference.fasta input/reference/$reference.fasta                         &>> log/index.$reference
bwa index input/reference/$reference.fasta                                                              &>> log/index.$reference
samtools faidx input/reference/$reference.fasta                                                         &>> log/index.$reference
java -jar /usr/local/share/picard/picard.jar CreateSequenceDictionary R=input/reference/$reference.fasta O=input/reference/$reference.dict    &>> log/index.$reference

date
