setwd("")
load("ROHall.RData")
load("inputForXROH.RData")
options   <- commandArgs(trailingOnly = TRUE)
chrom     <- as.numeric(options[1])
cat(paste0("# Chromosome ",chrom,"\n"))
ROHhere   <- ROHall[which(ROHall[,"CHR"]==chrom),]
nROH      <- nrow(ROHhere)

aa    <- as.matrix(listAlleleAge[[chrom]][,-1])
aa    <- aa[order(aa[,"BP"]),]

tmrca <- as.matrix(listTMRCA[[chrom]][,-1])
tmrca <- tmrca[order(tmrca[,"BP"]),]

ibp   <- listBP[[chrom]]
ival  <- listValues[[chrom]]

annotROH <- function(iROH){
  chrom <- ROHhere[iROH,"CHR"]
  pos1  <- ROHhere[iROH,"POS1"]
  pos2  <- ROHhere[iROH,"POS2"]
  len   <- (pos2-pos1)/1e6
  
  xroh         <- rep(NA,ncol(ival)+7)
  names(xroh)  <- c("IID","CHR","POS1","POS2","NSNP",cnvalues,"Allele_Age","TMRCA")
  xroh["IID"]  <- ROHhere[iROH,"IID"]
  xroh["CHR"]  <- ROHhere[iROH,"CHR"]
  xroh["POS1"] <- ROHhere[iROH,"POS1"]
  xroh["POS2"] <- ROHhere[iROH,"POS2"]
  xroh["NSNP"] <- sum(ibp>=pos1 & ibp<=pos2)
  
  ## Allele Age
  who <- which(aa[,"BP"]>=pos1 & aa[,"BP"]<=pos2)
  if(length(who)>0){
    xroh["Allele_Age"] <- sum(aa[who,"AlleleAge"])
  }else{
    xroh["Allele_Age"] <- 0
  }
  
  ## TMRCA
  who   <- which(tmrca[,"BP"]>=pos1 & tmrca[,"BP"]<=pos2)
  if(length(who)>0){
    xroh["TMRCA"] <- sum(tmrca[who,"TMRCA"])
  }else{
    xroh["TMRCA"] <- 0
  }
  
  ## Other annotations
  for(annot in cnvalues){
    who <- which(ibp>=pos1 & ibp<=pos2)
    if(length(who)>0){
      xroh[annot] <- sum(ival[who,annot])
    }else{
      xroh[annot] <- 0
    }
  }
  return(xroh)
}

library(parallel)
nROH <- nrow(ROHhere)
cat(paste0("# Analysing ",nROH," ROHs...\n"))
XROH <- do.call("rbind",mclapply(1:nROH,annotROH,mc.cores=16))
save(list="XROH",file=paste0("XROH/chrom",chrom,".RData"))

