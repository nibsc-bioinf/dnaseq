#!/usr/bin/Rscript
options(stringsAsFactors=FALSE,width=200)

library("data.table")
library("bit64")
library("xlsx")
library("xtable")
library("stringr")
library("seqinr")

argv <- commandArgs(T)
dir <- argv[1] #eg /home/AD/tbleazar/165
project <- argv[2] #eg 165
date <- argv[3] #eg 180430
#will look for an annotationfile of form paste0(dir,"/input/reference/",reference,".gff3")

setwd(paste0(dir,"/analysis"))

cnames <- c("sample","reference","chr","pos","ref","cov","R","A","C","G","T","N","pR","pA","pC","pG","pT","pN",
            "forward","reverse","insertions","deletions","readstart","readend","prevdel","refskip")

#
#  REFERENCES
#
references <- fread("../input/references",header=F,sep="\t")$V1
#references <- as.character(references[1,2,with=F])
#references <- unlist(strsplit(references,","))
print(references)
baseCalls <- data.frame()
for (reference in references) {
  fasta <- read.fasta(sprintf("../input/reference/%s.fasta",reference))
  for (n in names(fasta)) {
    f <- toupper(unlist(fasta[n]))
    nfasta <- length(f)
    baseCalls <- rbind(baseCalls,cbind(rep(reference,nfasta),rep(n,nfasta),1:nfasta,f))
  }
}
baseCalls <- data.table(baseCalls)
baseCalls <- cbind(rep("sample",nrow(baseCalls)),baseCalls,matrix(0,nrow=nrow(baseCalls),ncol=21))
colnames(baseCalls) <- cnames
baseCalls$pos <- as.numeric(baseCalls$pos)
baseCalls$plotPos <- 1:nrow(baseCalls)   #### HACK
baseCalls$chrEnd <- baseCalls$pos == 1
baseCalls$chrEnd[1] <- F
#print(baseCalls)

#
#  HERE CREATING GENOME ANNOTATIONS IF THE APPROPRIATE ANNOTATION FILES EXIST
#
annotatetable = FALSE

#initialise an empty data frame with our appropriate types
allredgreen = data.frame(REFERENCE=character(), CHROMOSOME=character(), POSITION=integer(), REPEATREGION=logical(), GENEREGION=logical(), GENENAMES=character())

for (reference in references) {
  print("Starting looking for annotation for reference:")
  annotationfile = paste0(dir,"/input/reference/",reference,".gff3")
  print(annotationfile)
  if (!file.exists(annotationfile)) {
    print("Failed to find annotation, will not use for the report")
    next
  } else {
    print("Found file, creating redgreen reference table")
    annotatetable = TRUE
    foundchromosome = FALSE
  }
  gfflines = readLines(con=annotationfile)
  for (gffline in gfflines) {
    if ((substr(gffline,1,1) == "#") || (nchar(gffline) < 2)) { #ignore comment or empty lines
      next
    }
    collect = unlist(strsplit(gffline, split="\t"))
    if ((collect[3] == "region") && (!foundchromosome)) {
      print("Creating redgreen dataframe for this reference, of size:")
      genomesize = strtoi(collect[5])
      print(genomesize)
      redgreen = data.frame(REFERENCE=rep(reference, genomesize), CHROMOSOME=rep(collect[1], genomesize), POSITION=1:genomesize, REPEATREGION=rep(FALSE, genomesize), GENEREGION=rep(FALSE, genomesize), GENENAMES=rep("", genomesize))
      print(summary(redgreen))
      print(head(redgreen))
      foundchromosome = TRUE
      print("Now going through the gff3 annotation lines slowly")
    }
    if (grepl("repeat", collect[3], fixed=TRUE)) {
      #print("Found a repeat")
      startpoint = strtoi(collect[4])
      stoppoint = strtoi(collect[5])
      if (stoppoint <= genomesize) {
        redgreen$REPEATREGION[startpoint:stoppoint] = TRUE
      } else {
        print("Error, gff3 file records a repeat that goes outside the genome size")
      }
    }
    if (grepl("gene", collect[3], fixed=TRUE)) {
      #print("Found a gene")
      startpoint = strtoi(collect[4])
      stoppoint = strtoi(collect[5])
      if (stoppoint <= genomesize) {
        redgreen$GENEREGION[startpoint:stoppoint] = TRUE
        genename = (unlist(strsplit((unlist(strsplit(collect[9], split="Name=")))[2], split=";")))[1]
        for (i in startpoint:stoppoint) {
          if (nchar(redgreen[i,'GENENAMES'])<1) {
            redgreen[i,'GENENAMES'] = genename
          } else {
            redgreen[i,'GENENAMES'] = paste(redgreen[i,'GENENAMES'], genename, sep=",")
          }
        }
      } else {
        print("Error, gff3 file records a gene that goes outside the genome size")
      }
    }
  }
  allredgreen = rbind(allredgreen, redgreen)
}
redgreen = allredgreen
#redgreen is of the form
#      REFERENCE  CHROMOSOME POSITION REPEATREGION GENEREGION GENENAMES
#1 HSV1-strain17 NC_001798.2        1         TRUE       TRUE       LAT
#2 HSV1-strain17 NC_001798.2        2         TRUE       TRUE       LAT
#3 HSV1-strain17 NC_001798.2        3         TRUE       TRUE       LAT
 
#
#  Now merge in the redgreen annotations to the baseCalls, if we have found the gff3 file
#
if (annotatetable) {
  print("Now merging the redgreen annotations with the baseCalls")
  baseCalls = merge(x=baseCalls, y=redgreen, by.x=c("reference","chr","pos"), by.y=c("REFERENCE","CHROMOSOME","POSITION"), all.x=TRUE)
  print(head(baseCalls))
}



#
#  SAMPLES
#
sampleFiles <- list.files("../input/","samples.*")
samples <- data.frame()
for (s in sampleFiles) {
  samples <- rbind(samples,read.delim(paste("../input/",s,sep=""),header=F))
}
colnames(samples) <- c("name","file")
samples <- data.table(samples)
samples <- samples[!grepl("#",samples$name)]
samples$sample <- gsub("/fastq/30/","_",samples$file)
print(samples)


#
#  CALLS
#
makeCalls <- function() {
  cat("makeCalls\n")
  allCalls <- c()
  for (sample in unique(samples$name)) {
    for (reference in references) {
      filename <- sprintf("../processed/%s/30.bwa.%s.samtools.calls",sample,reference)
      if (file.exists(filename) & file.size(filename) != 0) {
        calls <- fread(filename,header=F)
        calls[,1] <- sample
        calls[,reference:=reference]
        setcolorder(calls,c(1,ncol(calls),2:(ncol(calls)-1)))
        colnames(calls) <- cnames
      }
      allCalls <- rbind(allCalls,calls)
    }
  }
  colnames(allCalls) <- cnames
  cols <- 6:26
  allCalls[,(cols):=lapply(.SD,as.numeric),.SDcols=cols]
  #tom addition: also insert the columns from baseCalls table for repeat and gene region info if we have set logical annotatetable
  if (annotatetable) {
    print("Merging the allCalls table with additional annotations")
    allCalls = merge(x=allCalls, y=baseCalls[,c("chr","pos","REPEATREGION","GENEREGION","GENENAMES")], by.x=c("chr","pos"), by.y=c("chr","pos"), all.x=TRUE)
  }
  print("Now writing to tsv files")
  print(head(allCalls))
  write.table(allCalls,file=sprintf("%s.calls",project),row.names=F,quote=F,sep="\t")
  bA <- (allCalls$A   >= 100 & allCalls$pA >= 0.01)
  bC <- (allCalls$C   >= 100 & allCalls$pC >= 0.01)
  bG <- (allCalls$G   >= 100 & allCalls$pG >= 0.01)
  bT <- (allCalls$T   >= 100 & allCalls$pT >= 0.01)
  bP <- (allCalls$cov >= 100 & allCalls$pR <= 0.01)
  write.table(file=sprintf("%s.mutations.bi.tsv",project)    ,allCalls[(bA + bC + bG + bT) == 2,],row.names=F,col.names=T,quote=F,sep="\t")
  write.table(file=sprintf("%s.mutations.tri.tsv",project)   ,allCalls[(bA + bC + bG + bT) == 3,],row.names=F,col.names=T,quote=F,sep="\t")
  write.table(file=sprintf("%s.mutations.quad.tsv",project)  ,allCalls[(bA + bC + bG + bT) == 4,],row.names=F,col.names=T,quote=F,sep="\t")
  write.table(file=sprintf("%s.mutations.nonref.tsv",project),allCalls[bP,],row.names=F,col.names=T,quote=F,sep="\t")
  write.table(file=sprintf("%s.insertions.tsv",project),allCalls[allCalls$insertions / allCalls$cov >= 0.01 & allCalls$cov >= 100,],row.names=F,col.names=T,quote=F,sep="\t")
  write.table(file=sprintf("%s.deletions.tsv",project) ,allCalls[allCalls$deletions  / allCalls$cov >= 0.01 & allCalls$cov >= 100,],row.names=F,col.names=T,quote=F,sep="\t")
  write.table(file=sprintf("%s.prevdel.tsv",project)   ,allCalls[allCalls$prevdel    / allCalls$cov >= 0.01 & allCalls$cov >= 100,],row.names=F,col.names=T,quote=F,sep="\t")
}

#
#  ENHANCED SAMPLES
#
makeSamples <- function() {
  cat("makeSamples\n")
  stats <- c()
  for (sample in unique(samples$name)) {
    sSmall <- substring(sample,8)
    stats <- rbind(stats,c(sample,
               "raw",unlist(read.delim(sprintf("../input/%s/fastq/raw/%s_R1_001.stats",date,sSmall),header=F)[1,2:4]),
               "raw",unlist(read.delim(sprintf("../input/%s/fastq/raw/%s_R2_001.stats",date,sSmall),header=F)[1,2:4]),
               "30",unlist(read.delim(sprintf("../input/%s/fastq/30/%s.p1.stats",date,sSmall),header=F)[1,2:4]),
               "30",unlist(read.delim(sprintf("../input/%s/fastq/30/%s.p2.stats",date,sSmall),header=F)[1,2:4]),
               "30",unlist(read.delim(sprintf("../input/%s/fastq/30/%s.s1.stats",date,sSmall),header=F)[1,2:4]),
               "30",unlist(read.delim(sprintf("../input/%s/fastq/30/%s.s2.stats",date,sSmall),header=F)[1,2:4])))
  }
  colnames(stats) <- c("sample","r1type","r1reads","r1bases","r1mean","r2type","r2reads","r2bases","r2mean",
                                "p1type","p1reads","p1bases","p1mean","p2type","p2reads","p2bases","p2mean",
                                "s1type","s1reads","s1bases","s1mean","s2type","s2reads","s2bases","s2mean")
  stats <- data.table(stats)
  cols <- c(3:5,7:9,11:13,15:17,19:22,23:25)
  stats[,(cols):=lapply(.SD,as.numeric),.SDcols=cols]

  flagstats <- c()
  for (sample in unique(samples$name)) {
    for (reference in references) {
      fs <- read.delim(sprintf("../processed/%s/30.bwa.%s.samtools.flagstat",sample,reference),header=F,sep=" ")
      fs[,1] <- as.numeric(fs[,1])
      flagstats <- rbind(flagstats,c(sample,reference,fs[1,1],fs[5,1],fs[5,1]/fs[1,1],fs[4,1],fs[3,1]))
    }
  }
  colnames(flagstats) <- c("sample","reference","reads","mapped","percent","dups","supp")
  flagstats <- data.table(flagstats)
  cols <- c(3:7)
  flagstats[,(cols):=lapply(.SD,as.numeric),.SDcols=cols]
	
  mean   <- calls[,mean(cov),by=.(sample,reference)]
  m.bi   <- m.bi[,.N,by=.(sample,reference)]
  m.nr   <- m.nr[,.N,by=.(sample,reference)]
  colnames(mean)[3] <- "mean"
  colnames(m.bi)[3] <- "m.bi"
  colnames(m.nr)[3] <- "m.nr"

  samples <- merge(samples,stats,by=c("sample"),all=F)
  samples <- merge(samples,flagstats,by=c("sample"),all=T)

  samples[,readPercent:=50*(reads-supp)/r1reads]
  samples[,readsUseful:=mapped-dups-supp]

  samples <- merge(samples,mean,by=c("sample","reference"),all=T)
  if (nrow(m.nr) > 0) {
    samples <- merge(samples,m.nr,by=c("sample","reference"),all=T)
  } else {
    samples[,m.nr:=0]
  }
  samples$m.nr[is.na(samples$m.nr)] <- 0

  if (nrow(m.bi) > 0) {
    samples <- merge(samples,m.bi,by=c("sample","reference"),all=T)
  } else {
    samples[,m.bi:=0]
  }
  samples$m.bi[is.na(samples$m.bi)] <- 0

  write.table(samples,file=sprintf("%s.samples.tsv",project),row.names=F,col.names=T,quote=F,sep="\t")
  write.xlsx(samples,file=sprintf("%s.samples.xlsx",project),sheetName="Mapping Statistics")
  samples
}

#
#  RESULTS
#
makeResults <- function() {
  cat("makeResults\n")
  results <- c()
  for (r in references) {
    for (s in unique(samples$name)) {
      coverage <- calls[reference == r & sample == s]
      nr <- m.nr[sample==s & reference == r,pos,]
      bi <- m.bi[sample==s & reference == r,pos,]
      lc <- coverage$cov < 100
      results <- rbind(results,c(s,r,nrow(coverage),sum(lc),length(nr),length(bi)))
    }
  }
  results <- data.table(results)
  colnames(results) <- c("sample","reference","length","low","nonref","snps")
  results[,length:=as.numeric(length)]
  results[,low:=as.numeric(low)]
  results[,nonref:=as.numeric(nonref)]
  results[,snps:=as.numeric(snps)]
  write.table(results,file=sprintf("%s.results.tsv",project),row.names=F,col.names=T,quote=F,sep="\t")
  results
}

#
#  FIGURES
#
makeFigures <- function() {
  cat("makeFigures\n")
  smallBaseCalls <- baseCalls[,c("reference","chr","pos","plotPos","chrEnd")]
  for (r in references) {
    for (s in unique(samples$name)) {
      #cat("On the first reference/sample combo\n")
      cat(r,s,"\n")
      coverage <- merge(smallBaseCalls,calls[reference==r & sample==s],by=c("reference","chr","pos"))
      print("Just merged the coverage table with the calls:")
      print(head(coverage))
      if (nrow(m.nr) > 0) {
      nr       <- merge(smallBaseCalls, m.nr[reference==r & sample==s],by=c("reference","chr","pos"))
      } else {
      nr <- data.table(reference=character(0), chr=character(0), pos=integer(0)) #clause added by tom to handle empty m.nr
      }      
      bi       <- merge(smallBaseCalls, m.bi[reference==r & sample==s],by=c("reference","chr","pos"))
      print("Made bi merge")
      print(head(bi))
      col <- rep("black",nrow(baseCalls))
      lc <- coverage$plotPos[coverage$cov < 100]
      col[lc] <- "grey"
      if (annotatetable) {
        print("Using repeat region annotations in plot")
        repeatregions = coverage$plotPos[coverage$REPEATREGION]
        col[repeatregions] = "purple" #a long vector of colors for all the points along the coverage plot, with those in repeat regions set to purple
      }
      png(sprintf("figures/coverage_%s_%s.png",r,s),width=1400,height=700)
      par(cex=2)
      plot(coverage$plotPos,log10(coverage$cov),col=col,pch=".",xlim=c(0,maxpos),ylim=c(0,log10(maxcov)),xlab=sprintf("%s (%dbp)",r,nrow(coverage)),ylab="log10 Coverage")
      if (nrow(nr) > 0) { abline(v=nr$plotPos,col="red") }
      if (nrow(bi) > 0) { abline(v=bi$plotPos,col="green") }
      abline(h=2,col="grey")
      abline(v=smallBaseCalls$plotPos[smallBaseCalls$chrEnd],col="grey")
      if (annotatetable) {
        legend("bottomright",legend=c(sprintf("Non-reference : %d",nrow(nr)),sprintf("SNVs : %d",nrow(bi)),sprintf("Repeat region : %d",length(repeatregions)),sprintf("Low coverage : %d",length(lc))),col=c("red","green","purple","grey"),pch=19)
      } else {
        legend("bottomright",legend=c(sprintf("Non-reference : %d",nrow(nr)),sprintf("SNVs : %d",nrow(bi)),sprintf("Low coverage : %d",length(lc))),col=c("red","green","grey"),pch=19)
      }
      dev.off()
    }
  }
}

#
#  MUTATIONS
#
makeMutations <- function () {
  cat("makeMutations\n")
  mutations <- rbind(m.bi,m.tr,m.qu,m.nr)
#  mutations <- merge(results,mutations,by=c("sample","reference"),all.x=F,all.y=F,sort=F)
  mutations[,pA:=100*pA]
  mutations[,pC:=100*pC]
  mutations[,pG:=100*pG]
  mutations[,pT:=100*pT]
  mutations <- mutations[order(pos)]
  write.table(mutations,file=sprintf("%s.mutations.tsv",project),row.names=F,col.names=T,quote=F,sep="\t")
  mutations
}

#
#  REPORT
#
makeReport <- function() {
  cat("makeReport\n")
  sink(sprintf("%s.tex",project))
  cat("\\documentclass[a4paper]{article}\n")
  cat("\\usepackage{graphicx}\n")
  cat("\\usepackage{longtable}\n")
  cat("\\usepackage[margin=1.5cm]{geometry}\n")
  cat("\\begin{document}\n")
  for (r in references) {
    for (s in unique(samples$name)) {
        if (annotatetable) {
          m <- mutations[reference==r & sample==s,c("chr","pos","ref","cov","pA","pC","pG","pT","REPEATREGION","GENEREGION","GENENAMES")]
          m$genes = m$GENENAMES
          m$repeats = rep("",nrow(m))
          for (geneind in 1:nrow(m)) {
            if (m$REPEATREGION[geneind]) {
              m$repeats[geneind] = "repeat region"
            }
          }
          m = m[,c("chr","pos","ref","cov","pA","pC","pG","pT","repeats","genes")]
        } else {
          m <- mutations[reference==r & sample==s,c("chr","pos","ref","cov","pA","pC","pG","pT")]
        }
      #sorting the mutations list by chromosome name
      m <- m[order(m$chr),]
      cat("\\begin{center}\\Large \\verb|",r,"| : \\verb|",s,"|\\end{center}\n",sep="")
      cat("\\includegraphics[width=\\textwidth]{figures/coverage_",r,"_",s,"}\n",sep="")
      if (nrow(m) > 0) {
        if (annotatetable) {
          print(xtable(m[1:min(35,nrow(m))],align="|lrrrr|rrrr|cc|"),tabular.environment="longtable",floating=F)
        } else {
          print(xtable(m[1:min(35,nrow(m))],align="|lrrrr|rrrr|"),tabular.environment="longtable",floating=F)
        }
        if (nrow(m) > 35) { cat("Only top 35 rows shown \\\\ \n") }
      } else {
        cat("No SNPs \\\\ \n")
      }
      cat("\\newpage\n")
    }
  }
  cat("\\end{document}\n")
  sink()
  system(sprintf("pdflatex %s.tex",project))
}

#  ------------------------------------------------------------------
#
#  MAIN
#
if (! file.exists(sprintf("%s.calls",project))) {
  makeCalls()
}
calls <- fread(sprintf("%s.calls",project))
m.bi <- fread(sprintf("%s.mutations.bi.tsv",project))
m.tr <- fread(sprintf("%s.mutations.tri.tsv",project))
m.qu <- fread(sprintf("%s.mutations.quad.tsv",project))
m.nr <- fread(sprintf("%s.mutations.nonref.tsv",project))

print(head(m.bi))

if (! file.exists(sprintf("%s.samples.tsv",project))) {
  makeSamples()
}
samples <- fread(sprintf("%s.samples.tsv",project))

#if (! file.exists(sprintf("%s.results.tsv",project))) {
  dir.create("figures",showWarnings=F,recursive=T)
  maxcov <- max(calls$cov)
  maxpos <- nrow(baseCalls)
  makeResults()
  makeFigures()
#}
results <- fread(sprintf("%s.results.tsv",project))

if (! file.exists(sprintf("%s.mutations.tsv",project))) {
  makeMutations()
}
mutations <- fread(sprintf("%s.mutations.tsv",project))

makeReport()
