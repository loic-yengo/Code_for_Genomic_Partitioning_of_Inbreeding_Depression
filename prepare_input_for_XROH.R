setwd("")
cc  <- c("numeric","numeric","character",rep("numeric",191-3))
ccb <- c("numeric","character","numeric","numeric","character","character")
i1  <- "annot-files"
i2  <- "ukb-geno"

QC_annot <- function(chrom){
  cat(paste0("\tChromosome ",chrom,"...\n"))
  d   <- read.table(paste0(i1,"/baselineLF2.2.UKB.",chrom,".annot.gz"),header=TRUE,stringsAsFactors=FALSE,colClasses=cc)
  cnd <- colnames(d)
  out <- c(grep("500",cnd),grep("MAFbin",cnd),grep("Predicted_Allele",cnd),grep("AFR",cnd),grep("ASMC",cnd),grep("LindbladToh",cnd))
  d   <- d[,-out]
  cnd <- colnames(d)
  X   <- as.matrix(d[,-(1:4)])
  an  <- colnames(X)
  ## filter annotations
  ss  <- gsub("_lowfreq","",an)
  ss  <- gsub("_common","",ss)
  
  dhs           <- ss[grep("DHS_",ss)]
  #conserved     <- ss[grep("Conserved_",ss)]
  H3K27ac       <- ss[grep("H3K27ac",ss)]
  H3K27ac       <- H3K27ac[-grep("BLUEPRINT_",H3K27ac)]
  H3K4me1       <- ss[grep("H3K4me1",ss)]
  H3K4me1       <- H3K4me1[-grep("BLUEPRINT_",H3K4me1)]
  H3K4me3       <- ss[grep("H3K4me3",ss)]
  H3K9ac        <- ss[grep("H3K9ac",ss)]
  Enhancer      <- ss[ss%in%c("Enhancer_Andersson","Enhancer_Hoffman")]
  
  ss[which(ss%in%dhs)]       <- "DHS"
  #ss[which(ss%in%conserved)] <- "Conserved"
  ss[which(ss%in%H3K27ac)]   <- "H3K27ac"
  ss[which(ss%in%H3K4me1)]   <- "H3K4me1"
  ss[which(ss%in%H3K4me3)]   <- "H3K4me3"
  ss[which(ss%in%H3K9ac)]    <- "H3K9ac"
  ss[which(ss%in%Enhancer)]  <- "Enhancer"
  
  nval         <- apply(X,2,function(u) length(unique(u)))
  binary_annot <- unique(ss[which(nval==2)])
  contin_annot <- unique(ss[which(nval>2)])
  
  si <- sort(unique(ss))
  na <- length(si)
  x  <- matrix(0,nrow=nrow(d),ncol=na)
  for(i in 1:na){
    k <- which(ss==si[i])
    if(length(k)>1){
      x[,i] <- rowSums(X[,k])
    }else{
      x[,i] <- X[,k]
    }
    nval <- length(unique(x[,i]))
    if(si[i]%in%binary_annot){
      x[,i] = as.numeric(x[,i]>0)
    }
  }
  colnames(x)     <- si
  cnx             <- colnames(x)
  nval_x          <- apply(x,2,function(u) length(unique(u)))
  binary_annot_x  <- cnx[which(nval_x==2)]
  contin_annot_x  <- cnx[which(nval_x>2)]
  whoCONT         <- which(nval_x>2) 
  colnames(x)[whoCONT] <- paste0("CONT_",cnx[whoCONT])
  chrbp <- as.matrix(d[,1:2])
  return(cbind(chrbp,x))
}

dtannot   <- do.call("rbind",mclapply(1:22,QC_annot,mc.cores=10))
chrbp     <- dtannot[,1:2]
values    <- as.matrix(dtannot[,-(1:2)])
cnvalues  <- colnames(values)
kqt_annot <- grep("CONT_",colnames(values))
for(k in kqt_annot){
  values[,k] <- values[,k] / max(values[,k])
}
values[,"CONT_Nucleotide_Diversity_10kb"] <- 1-values[,"CONT_Nucleotide_Diversity_10kb"]
propAnnot <- colMeans(values)


### Deal with ROH now
load("ROHall.RData")

## Age of the allele
## ASMC - TMRCA: Time to Most Recent Common Ancestor (Haplotype/Allele Age) - r~0.29
alleleAge     <- read.table("download/continuous/alleleage.bed",stringsAsFactors = F,colClasses=c("character",rep("numeric",3)))
alleleAge[,1] <- as.numeric(gsub("chr","",alleleAge[,1],fixed=TRUE))
alleleAge     <- rbind(as.matrix(alleleAge[,c(1,2,4)]),as.matrix(alleleAge[,c(1,3,4)])); colnames(alleleAge) <- c("CHR","BP","AlleleAge")
alleleAge[,3] <- 1-alleleAge[,3] / max(alleleAge[,3]) # changed on Nov 4th

asmc          <- read.table("download/continuous/ASMC.bed",stringsAsFactors = F,colClasses=rep("numeric",5))[,-5]
asmc          <- rbind(as.matrix(asmc[,c(1,2,4)]),as.matrix(asmc[,c(1,3,4)])); colnames(asmc) <- c("CHR","BP","TMRCA")
asmc[,3]      <- 1-asmc[,3] / max(asmc[,3]) # changed on Nov 4th

propAnnot     <- c(propAnnot,Allele_Age=mean(alleleAge[,3]),TMRCA=mean(asmc[,3]))
save(list="propAnnot",file="propAnnotForROH.RData")


listBP     <- lapply(1:22, function(chrom){
  chrbp[which(chrbp[,"CHR"]==chrom),"BP"]
})
listValues <- lapply(1:22, function(chrom){
  values[which(chrbp[,"CHR"]==chrom),]
})
listAlleleAge <- lapply(1:22, function(chrom){
  alleleAge[which(alleleAge[,"CHR"]==chrom),]
})
listTMRCA <- lapply(1:22, function(chrom){
  asmc[which(asmc[,"CHR"]==chrom),]
})

## Needed as input
## listBP, listValues, listAlleleAge, listTMRCA
## cnvalues

save(list=c("listBP","listValues","listAlleleAge","listTMRCA","propAnnot","cnvalues"),
     file="inputForXROH.RData")

