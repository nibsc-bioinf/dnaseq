#This is a script to look at a given raw directory with some fastq.gz files in it which all start with a given project number
#It then generates a samples.date file with lines of the form:
#181010_200-R0104-day5-A_S45_L001        181010/fastq/30/200-R0104-day5-A_S45_L001
#The rawdir has filenames of the form:
#213-R0100-day13-A_S20_L001_R1_001.fastq.gz

import sys
import os

rawdir = sys.argv[1] #/sequencing/miseq/output/181012_M01745_0220_000000000-BJGCW/Data/Intensities/BaseCalls
seqdate = sys.argv[2] #181012
projectnumber = sys.argv[3] #213

forwardfiles = []
for filename in os.listdir(rawdir):
    if "_R1_001.fastq.gz" in filename and filename[:(len(projectnumber)+1)] == (projectnumber+"-"):
        forwardfiles.append(filename)
for filename in forwardfiles:
    towrite = seqdate+"_"+(filename.replace("_R1_001.fastq.gz",""))+"\t"+seqdate+"/fastq/30/"+(filename.replace("_R1_001.fastq.gz",""))
    print(towrite)
