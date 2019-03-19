#Modular bioinformatics pipeline
#NIBSC - NGS/bioinf - Thomas Bleazard
#Given location of raw fastq sequencing data, will write a bash script to perform processing from trimming and alignment to analysis

#Usage:
#mobipi.py -rawdir /sequencing/miseq/output/180928_M01745_0217_000000000-C4JYW/Data/Intensities/BaseCalls -working /sequencing/projects/207 -project 207

import os
import sys

#Here are the input command line argument options
def parseoptions(commargs):
    print("Parsing command-line arguments into script settings")
    options = {}
    #options["parallel"] = "bash" #bash default or qsub for Sun Grid Engine submission commands written
    options["rawdir"] = None #/sequencing/miseq/output/181005_M01745_0218_000000000-BJ78Y/Data/Intensities/BaseCalls or /sequencing/nextseq/processed/180904/bcl2fastq/182
    options["working"] = None #/sequencing/projects/207 and note that the bash script created by mobipi will be written here
    options["project"] = None #207 which should match the first characters of the fastq files for this project
    options["reference"] = None #J01917-1.fasta and it is a requirement that $working/reference/J01917-1.fasta be present
    options["sequencer"] = "miseq" #or nextseq
    options["miseq-adapter"] = "/usr/local/bin/adapters/NexteraPE-PE.fa" #change this if doing scriptseq
    options["nextseq-adapter"] = "/usr/local/bin/adapters/TruSeq-adapters-recommended.fa" #change this if using other nextseq adapters or to add poly(A) and poly(G) if doing transcriptome assembly
    options["genegtf"] = None #/sequencing/projects/202/reference/hsa.ensembl.92.chr.gtf
    #options["filetable"] = None #this allows for an input raw file list, with associated different references for each
    for i in range(len(commargs)):
        if commargs[i][0] == "-" and commargs[i][1:] in options:
            options[commargs[i][1:]] = commargs[i+1]
    return options

#Requirements for the input arguments
def passinginput(options):
    print("Checking arguments acceptable")
    if options["working"] == None or not os.path.isdir(options["working"]):
        print("Bad argument working: please create a working directory")
        return False
    #more tests to be added
    return True

#Create directory structure
def setupdirectory(working, shellscript):
    print("Creating subdirectories")
    shellscript.write("mkdir -p "+working+"/log\n") #logs of stderr and stdout
    shellscript.write("mkdir -p "+working+"/trimmed\n") #trimmed fastq data now goes here
    shellscript.write("mkdir -p "+working+"/aligned\n") #will put subdirectories for each sample in here with alignment and other processing files
    shellscript.write("mkdir -p "+working+"/analysis\n") #results and figures go here
    shellscript.write("mkdir -p "+working+"/reference\n") #This should already be here, with a reference fasta inside
    shellscript.write("mkdir -p "+working+"/custom\n") #Custom analysis and code goes here
    shellscript.write("mkdir -p "+working+"/qc\n") #Fastqc and multiqc output go here
    shellscript.write("mkdir -p "+working+"/denovo\n") #de novo assemblies go here
    shellscript.write("mkdir -p "+working+"/tmp\n")
    shellscript.write("mkdir -p "+working+"/deseq\n") #in case we do differential expression analysis

#Collect filenames for the raw fastq.gz data, taking only forward files for miseq, and only forward lane one files for nextseq
def getfilelist(rawdir, projectnumber):
    print("Fetching list of forward lane one filenames from the raw directory")
    forlones = []
    for filename in os.listdir(rawdir):
        if "_L001_R1_001.fastq.gz" in filename and filename[:len(projectnumber)] == projectnumber:
            forlones.append(filename)
    print("Found "+str(len(forlones))+" samples")
    return forlones

#Run fastqc
def runfastqc(working, shellscript, rawdir, forwardfiles, sequencer):
    print("Running fastqc")
    if sequencer == "miseq":
        for filename in forwardfiles:
            shellscript.write("fastqc -o "+working+"/qc "+rawdir+"/"+filename+"\n")
            shellscript.write("fastqc -o "+working+"/qc "+rawdir+"/"+(filename.replace("_L001_R1_001.fastq.gz","_L001_R2_001.fastq.gz"))+"\n")
    if sequencer == "nextseq":
        for filename in forwardfiles:
            for laneid in ["_L001_","_L002_","_L003_","_L004_"]:
                shellscript.write("fastqc -o "+working+"/qc "+rawdir+"/"+(filename.replace("_L001_",laneid))+"\n")
                shellscript.write("fastqc -o "+working+"/qc "+rawdir+"/"+((filename.replace("_L001_",laneid)).replace("R1_001.fastq.gz","R2_001.fastq.gz"))+"\n")

#Directly make trimmer auxiliary script and write it out now ready for use
def maketrimmerscript(working):
    fileout = open(working+"/runcutadapt.sh", "w")
    fileout.write("#!/bin/bash\n")
    fileout.write("cd "+working+"\n")
    fileout.write("adapter=$1\n")
    fileout.write("forwardfile=$2\n")
    fileout.write("reversefile=$3\n")
    fileout.write("outfor=$4\n")
    fileout.write("outrev=$5\n")
    fileout.write("cutadapt -a file:$adapter -A file:$adapter -g file:$adapter -G file:$adapter -o $outfor -p $outrev $forwardfile $reversefile -q 30,30 --minimum-length 50 --times 40 -e 0.1 --max-n 0\n")
    fileout.close()

#Run cutadapt then merge if nextseq
def trimreads(working, shellscript, rawdir, forwardfiles, sequencer, miseqadapter, nextseqadapter):
    print("Trimming reads using auxiliary script "+working+"/runcutadapt.sh")
    shellscript.write("chmod u+x "+working+"/runcutadapt.sh\n")
    if sequencer == "miseq":
        adapter = miseqadapter
        print("Using adapter "+adapter+" to trim. Please change script if this is incorrect.")
        for forwardfile in forwardfiles:
            shellscript.write(working+"/runcutadapt.sh "+adapter+" "+rawdir+"/"+forwardfile+" "+rawdir+"/"+(forwardfile.replace("_L001_R1_001.fastq.gz","_L001_R2_001.fastq.gz"))+" ")
            shellscript.write(working+"/trimmed/"+forwardfile+" "+working+"/trimmed/"+(forwardfile.replace("_L001_R1_001.fastq.gz","_L001_R2_001.fastq.gz"))+"\n")
    if sequencer == "nextseq":
        adapter = nextseqadapter
        print("Using adapter "+adapter+" to trim. Please change script if this is incorrect.")
        for forwardfile in forwardfiles:
            for laneid in ["_L001_","_L002_","_L003_","_L004_"]:
                forraw = forwardfile.replace("_L001_",laneid)
                revraw = forwardfile.replace("_L001_",laneid).replace("_R1_001.fastq.gz","_R2_001.fastq.gz")
                shellscript.write(working+"/runcutadapt.sh "+adapter+" "+rawdir+"/"+forraw+" "+rawdir+"/"+revraw+" ")
                shellscript.write(working+"/trimmed/"+forraw+" "+working+"/trimmed/"+revraw+" > "+working+"/log/"+forraw+".cutadapt.log\n")
            shellscript.write("cat")
            for laneid in ["_L001_","_L002_","_L003_","_L004_"]:
                forraw = forwardfile.replace("_L001_",laneid)
                shellscript.write(" "+working+"/trimmed/"+forraw)
            shellscript.write(" > "+working+"/trimmed/merging.fastq.gz\n")
            for laneid in ["_L001_","_L002_","_L003_","_L004_"]:
                forraw = forwardfile.replace("_L001_",laneid)
                shellscript.write("rm "+working+"/trimmed/"+forraw+" \n")
            shellscript.write("mv "+working+"/trimmed/merging.fastq.gz "+working+"/trimmed/"+forwardfile+"\n")
            shellscript.write("cat")
            for laneid in ["_L001_","_L002_","_L003_","_L004_"]:
                revraw = forwardfile.replace("_L001_",laneid).replace("_R1_001.fastq.gz","_R2_001.fastq.gz")
                shellscript.write(" "+working+"/trimmed/"+revraw)
            shellscript.write(" > "+working+"/trimmed/merging.fastq.gz\n")
            for laneid in ["_L001_","_L002_","_L003_","_L004_"]:
                revraw = forwardfile.replace("_L001_",laneid).replace("_R1_001.fastq.gz","_R2_001.fastq.gz")
                shellscript.write("rm "+working+"/trimmed/"+revraw+" \n")
            shellscript.write("mv "+working+"/trimmed/merging.fastq.gz "+working+"/trimmed/"+(forwardfile.replace("_R1_001.fastq.gz","_R2_001.fastq.gz"))+"\n")



#Directly make auxiliary de novo command script ready for use now
def makedenovoscript(working):
    fileout = open(working+"/rundenovo.sh", "w")
    fileout.write("#!/bin/bash\n")
    fileout.write("date\n")
    fileout.write("velvet=$1\n")
    fileout.write("rootdir=$2\n")
    fileout.write("trimp1=$3\n")
    fileout.write("trimp2=$4\n")
    fileout.write("threads=2\n")
    fileout.write("readlength=255\n")
    fileout.write("rm -r $velvet\n")
    fileout.write("mkdir -p $rootdir\n")
    fileout.write("cd $rootdir\n")
    velveter = 'velvetCommand="VelvetOptimiser.pl -d $velvet --s 51 --e $readlength -m "-0.001" --threads $threads -f '
    velveter = velveter+"'-fastq.gz -shortPaired $trimp1 $trimp2'"
    velveter = velveter+'"\n'
    fileout.write(velveter)
    fileout.write("eval $velvetCommand\n")
    fileout.close()
     
#Perform de novo assembly with Velvet Optimiser using rundenovo.sh auxiliary script
def dodenovo(working, shellscript, forwardfiles):
    print("Performing de novo assembly with auxiliary script "+working+"/rundenovo.sh")
    shellscript.write("chmod u+x "+working+"/rundenovo.sh\n")
    for forwardfile in forwardfiles:
        velvet = working+"/denovo/"+(forwardfile.replace("_L001_R1_001.fastq.gz",""))
        rootdir = working+"/denovo"
        trimp1 = working+"/trimmed/"+forwardfile
        trimp2 = working+"/trimmed/"+(forwardfile.replace("_L001_R1_001.fastq.gz","_L001_R2_001.fastq.gz"))
        shellscript.write(working+"/rundenovo.sh "+velvet+" "+rootdir+" "+trimp1+" "+trimp2+"\n")

#Directly make auxiliary reference indexer command script ready for use now
def makeindexerscript(working):
    fileout = open(working+"/runindexer.sh", "w")
    fileout.write("#!/bin/bash\n")
    fileout.write("cd "+working+"/reference\n")
    fileout.write("reference=$1\n")
    fileout.write("bowtie2-build $reference $reference\n")
    #adding for large index for tophat2 use
    fileout.write("refbase=`echo $reference | sed 's/.fasta//g' | sed 's/.fa//g'`\n")
    fileout.write("bowtie2-build --large-index $reference $refbase\n")
    fileout.write("bwa index $reference\n")
    fileout.write("samtools faidx $reference\n")
    fileout.write("java -jar /usr/local/share/picard/picard.jar CreateSequenceDictionary R=$reference O=$reference.dict\n")
    fileout.close()

#Index the provided reference fasta file
def doreferenceindex(working, shellscript, reference):
    print("Indexing reference provided at: "+reference+" using auxiliary script "+working+"/runindexer.sh")
    shellscript.write("chmod u+x "+working+"/runindexer.sh\n")
    shellscript.write(working+"/runindexer.sh "+reference+"\n")

#Directly make auxiliary script for running bwa for a single sample
#Will put the aligned files all in the aligned folder together
#Also build bam index and run flagstat
def makealignmentscript(working):
    fileout = open(working+"/runalignment.sh", "w")
    fileout.write("#!/bin/bash\n")
    fileout.write("cd "+working+"/aligned\n")
    fileout.write("reffasta=$1\n") #full path to the reference
    fileout.write("fortrimmed=$2\n") #full path to the forward trimmed fastq
    fileout.write("revtrimmed=$3\n")
    fileout.write("shortname=$4\n") #the forward filename with no path, and with _L001_R1_001.fastq.gz removed from the name
    fileout.write("bwa mem -t 4 $reffasta $fortrimmed $revtrimmed > "+working+"/aligned/$shortname.unsorted.sam\n")
    fileout.write("java -jar /usr/local/share/picard/picard.jar AddOrReplaceReadGroups I="+working+"/aligned/$shortname.unsorted.sam O="+working+"/aligned/$shortname.readgrouped.bam")
    fileout.write(' LB="Library" PU="PlatformUnit" SM="$shortname" PL="Illumina" TMP_DIR='+working+'/tmp SO=coordinate\n')
    fileout.write("java -jar /usr/local/share/picard/picard.jar MarkDuplicates I="+working+"/aligned/$shortname.readgrouped.bam O="+working+"/aligned/$shortname.readgrouped.deduped.bam")
    fileout.write(" METRICS_FILE="+working+"/aligned/$shortname.dups TMP_DIR="+working+"/tmp ASSUME_SORTED=T\n")
    fileout.write("java -jar /usr/local/share/picard/picard.jar BuildBamIndex TMP_DIR="+working+"/tmp I="+working+"/aligned/$shortname.readgrouped.deduped.bam\n")
    fileout.close()

#Now run through each sample with that alignment script
def doalignment(working, shellscript, reference, forwardfiles):
    print("Writing commands to use runalignment.sh to align the trimmed fastq files")
    shellscript.write("chmod u+x "+working+"/runalignment.sh\n")
    reffasta = working+"/reference/"+reference
    for forwardfile in forwardfiles:
        fortrimmed = working+"/trimmed/"+forwardfile
        revtrimmed = working+"/trimmed/"+(forwardfile.replace("_L001_R1_001.fastq.gz","_L001_R2_001.fastq.gz"))
        shortname = forwardfile.replace("_L001_R1_001.fastq.gz","")
        shellscript.write(working+"/runalignment.sh "+reffasta+" "+fortrimmed+" "+revtrimmed+" "+shortname+"\n")

#Make auxiliary script to run lofreq for a single sample
def makelofreqscript(working):
    fileout = open(working+"/runlofreq.sh", "w")
    fileout.write("#!/bin/bash\n")
    fileout.write("cd "+working+"/analysis\n")
    fileout.write("reffasta=$1\n") #full path to the reference
    fileout.write("shortname=$2\n") #no path, and _L001_R1_001.fastq.gz removed from the filename
    fileout.write("lofreq indelqual --dindel -f $reffasta -o "+working+"/aligned/$shortname.readgrouped.deduped.indelqual.bam "+working+"/aligned/$shortname.readgrouped.deduped.bam\n")
    fileout.write("samtools index "+working+"/aligned/$shortname.readgrouped.deduped.indelqual.bam\n")
    fileout.write("lofreq call -f $reffasta -o "+working+"/analysis/lofreq.$shortname.vcf --call-indels "+working+"/aligned/$shortname.readgrouped.deduped.indelqual.bam\n")

#Now use this script to run lofreq for each aligned sample
def dolofreq(working, shellscript, reference, forwardfiles):
    print("Writing commands to use runlofreq.sh to call variants")
    shellscript.write("chmod u+x "+working+"/runlofreq.sh\n")
    reffasta = working+"/reference/"+reference
    for forwardfile in forwardfiles:
        shortname = forwardfile.replace("_L001_R1_001.fastq.gz","")
        shellscript.write(working+"/runlofreq.sh "+reffasta+" "+shortname+"\n")

#Write out commands to run tophat2 - this is hardcoded for project 202 for now
def preparetophat(working, shellscript, forwardfiles):
    print("Writing commands to run tophat2")
    threads = "4"
    innerdistance = "50"
    gtfindex = genegtf+".indexed"
    referencebase = ((working+"/reference/"+reference).replace(".fasta","")).replace(".fa","")
    shellscript.write("tophat2 -G "+genegtf+" --transcriptome-index "+gtfindex+" "+referencebase+"\n")
    for forwardfile in forwardfiles:
        shortname = forwardfile.replace("_L001_R1_001.fastq.gz","")
        fortrimmed = working+"/trimmed/"+forwardfile
        revtrimmed = working+"/trimmed/"+(forwardfile.replace("_L001_R1_001.fastq.gz","_L001_R2_001.fastq.gz"))
        shellscript.write("mkdir -p "+working+"/aligned/"+shortname+"\n")
        shellscript.write("tophat2 -p "+threads+" -o "+working+"/aligned/"+shortname+" --mate-inner-dist "+innerdistance+" --coverage-search -G "+genegtf+" --transcriptome-index "+gtfindex)
        shellscript.write(" "+referencebase+" "+fortrimmed+" "+revtrimmed+" > "+working+"/log/"+shortname+".tophat2.log\n")
        shellscript.write("java -jar /usr/local/share/picard/picard.jar AddOrReplaceReadGroups I="+working+"/aligned/"+shortname+"/accepted_hits.bam O="+working+"/aligned/"+shortname+"/aligned.readgroup.bam")
        shellscript.write(' LB="Library" PU="PlatformUnit" SM="'+shortname+'" PL="Illumina"\n')
        shellscript.write("java -jar /usr/local/share/picard/picard.jar BuildBamIndex I="+working+"/aligned/"+shortname+"/aligned.readgroup.bam\n")

#Write commands to run htseq-count - need source activate py2.7 from the tbleazar conda and also hardcode the gtf for project 202 for now
def runhtseqcount(working, shellscript, forwardfiles, genegtf):
    print("Writing commands to run htseq-count")
    shellscript.write("source activate py2.7\n")
    for forwardfile in forwardfiles:
        shortname = forwardfile.replace("_L001_R1_001.fastq.gz","")
        shellscript.write("htseq-count -f bam -r pos -a 5 -s reverse -t exon -i gene_name -m union --nonunique=all "+working+"/aligned/"+shortname+"/aligned.readgroup.bam "+genegtf)
        shellscript.write(" > "+working+"/deseq/"+shortname+".htseq.a5.txt\n")

#Main function
def main():
    print("")
    print("Modular bioinformatics pipeline")
    print("")
    #Parse options
    options = parseoptions(sys.argv)
    for option in sorted(options.keys()):
        print(option+"\t"+str(options[option]))
    print("")
    #Check for input argument errors
    if not passinginput(options):
        print("Usage:")
        print("mobipi.py -rawdir /sequencing/miseq/output/180928_M01745_0217_000000000-C4JYW/Data/Intensities/BaseCalls -working /sequencing/projects/207 -project 207")
        sys.exit()
    #Now beginning writing shell script in working directory
    if os.path.isfile(options["working"]+"/processing."+options["project"]+".sh"):
        print(options["working"]+"/processing."+options["project"]+".sh already exists. Exiting.")
        sys.exit()
    print("Writing commands to run pipeline to file:")
    print(options["working"]+"/processing."+options["project"]+".sh")
    shellscript = open(options["working"]+"/processing."+options["project"]+".sh", "w")
    #Create directory structure
    setupdirectory(options["working"], shellscript)
    #Get list of forward files
    forwardfiles = getfilelist(options["rawdir"], options["project"])
    #Run fastqc to check quality (will run multiqc on this output manually if needed)
    runfastqc(options["working"], shellscript, options["rawdir"], forwardfiles, options["sequencer"])
    #Make auxiliary trimmer script
    maketrimmerscript(options["working"])
    #Trim with cutadapt (and then merge if nextseq data)
    trimreads(options["working"], shellscript, options["rawdir"], forwardfiles, options["sequencer"], options["miseq-adapter"], options["nextseq-adapter"])
    #Make auxiliary de novo script
    makedenovoscript(options["working"])
    #Now perform de novo assembly
    dodenovo(options["working"], shellscript, forwardfiles)
    #Make auxiliary script to index a reference fasta
    makeindexerscript(options["working"])
    #Now index the reference provided   
    if options["reference"] == None or not os.path.isfile(options["working"]+"/reference/"+options["reference"]):
        print("Could not find reference file: "+options["working"]+"/reference/"+options["reference"])
        print("Terminating script generation here")
        sys.exit()
    doreferenceindex(options["working"], shellscript, options["reference"])
    #Make auxiliary script to perform alignment and tidy processing
    makealignmentscript(options["working"])
    #Now use this script with the forward filenames
    doalignment(options["working"], shellscript, options["reference"], forwardfiles)
    #Make auxiliary script to run lofreq
    makelofreqscript(options["working"])
    #Now use this script to run lofreq for each aligned sample
    dolofreq(options["working"], shellscript, options["reference"], forwardfiles)
    #If nextseq then run tophat2
    preparetophat(options["working"], shellscript, forwardfiles, options["genegtf"])
    #If nextseq also now run htseq-count
    runhtseqcount(options["working"], shellscript, forwardfiles, options["genegtf"])

    shellscript.close() 

#Execute main
if __name__ == "__main__":
    main()





