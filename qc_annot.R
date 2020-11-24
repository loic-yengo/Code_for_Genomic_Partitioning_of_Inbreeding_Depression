setwd("annots/")
cc  <- c("numeric","numeric","character",rep("numeric",191-3))
ccb <- c("numeric","character","numeric","numeric","character","character")
i1  <- "annot-files"
i2  <- "ukb-geno"

## Age of the allele
alleleAge     <- read.table("download/continuous/alleleage.bed",stringsAsFactors = F,colClasses=c("character",rep("numeric",3)))
alleleAge[,1] <- as.numeric(gsub("chr","",alleleAge[,1],fixed=TRUE))

## ASMC ??? - TMRCA: Time to Most Recent Common Ancestor (Haplotype/Allele Age) - r~0.29
asmc <- read.table("download/continuous/ASMC.bed",stringsAsFactors = F,colClasses=rep("numeric",5))

load("ukb-geno/rsq/info.RData")

QC_annot <- function(chrom){
  cat(paste0("\tChromosome ",chrom,"...\n"))
  b    <- read.table(paste0(i2,"/chrom",chrom,".bim"),stringsAsFactors=FALSE,colClasses=ccb)
  d    <- read.table(paste0(i1,"/baselineLF2.2.UKB.",chrom,".annot.gz"),header=TRUE,stringsAsFactors=FALSE,colClasses=cc)
  pos  <- intersect(b[,4],d[,"BP"])
  pos  <- intersect(pos,info[which(info[,"CHR"]==chrom & info[,"INFO"]>=0.3),"BP"])
  mlt  <- unique(c(b[which(duplicated(b[,4])),4],d[which(duplicated(d[,"BP"])),"BP"]))
  outp <- which(pos%in%mlt)
  if(length(outp)>0){
    pos <- pos[-outp]
  }
  ## Freq
  freq <- read.table(paste0("freq/chrom",chrom,".afreq"),header=TRUE,stringsAsFactors=FALSE,
                     comment.char = "!",colClasses=c("numeric",rep("character",3),"numeric","numeric"))[,-c(1,2)]
  b    <- b[which(b[,4]%in%pos),]
  freq <- freq[which(b[,4]%in%pos),]
  d    <- d[which(d[,"BP"]%in%pos),]

  maf  <- ifelse(freq[,"ALT_FREQS"]>0.5,1-freq[,"ALT_FREQS"],freq[,"ALT_FREQS"])
  grp  <- cut(maf,breaks = c(0,0.05,0.1,0.2,0.3,0.4,0.5),include.lowest = TRUE)
  het  <- 2 * maf * (1-maf)
  xgrp <- model.matrix(~grp-1); colnames(xgrp) <- c("MAF_1_to_5","MAF_5_to_10","MAF_10_to_20","MAF_20_to_30","MAF_30_to_40","MAF_40_to_50")
  
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
  H3K27ac       <- ss[grep("H3K27ac",ss)]
  H3K27ac       <- H3K27ac[-grep("BLUEPRINT_",H3K27ac)]
  H3K4me1       <- ss[grep("H3K4me1",ss)]
  H3K4me1       <- H3K4me1[-grep("BLUEPRINT_",H3K4me1)]
  H3K4me3       <- ss[grep("H3K4me3",ss)]
  H3K9ac        <- ss[grep("H3K9ac",ss)]
  Enhancer      <- ss[ss%in%c("Enhancer_Andersson","Enhancer_Hoffman")]
  
  ss[which(ss%in%dhs)]       <- "DHS"
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
  colnames(x) <- si
  
  ## Add a few more annotations
  ## allele age, ASMC, LD score
  aa <- as.matrix(alleleAge[which(alleleAge[,1]==chrom),-1])
  aa <- aa[which(aa[,1]%in%d[,2] | aa[,2]%in%d[,2]),]
  aa <- rbind(aa[,c(1,3)],aa[,c(2,3)]); colnames(aa) <- c("BP","AlleleAge")
  aa <- aa[order(aa[,"BP"]),]
  aa <- merge(d[,c("BP","SNP")],aa,by="BP",all.x=TRUE)
  aa <- aa[-which(duplicated(aa[,"BP"])),]
  
  whoNA <- which(is.na(aa[,"AlleleAge"]))
  if(length(whoNA)>0){
    aa[whoNA,"AlleleAge"] <- 0
  }
  
  check <- mean(aa[,"BP"]==d[,"BP"])
  if(check!=1){
    cat(paste0("\tWell there is an issue with chromsome ",chrom," and allele age!\n"))
  }
  
  ## Time to most recent common ancestor - background selection
  tmrca <- as.matrix(asmc[which(asmc[,1]==chrom),-c(1,5)])
  tmrca <- tmrca[which(tmrca[,1]%in%d[,2] | tmrca[,2]%in%d[,2]),]
  tmrca <- rbind(tmrca[,c(1,3)],tmrca[,c(2,3)]); colnames(tmrca) <- c("BP","TMRCA")
  tmrca <- tmrca[order(tmrca[,"BP"]),]
  tmrca <- merge(d[,c("BP","SNP")],tmrca,by="BP",all.x=TRUE)
  tmrca <- tmrca[-which(duplicated(tmrca[,"BP"])),]
  whoNA <- which(is.na(tmrca[,"TMRCA"]))
  if(length(whoNA)>0){
    tmrca[whoNA,"TMRCA"] <- 0
  }
  
  check <- mean(tmrca[,"BP"]==d[,"BP"])
  if(check!=1){
    cat(paste0("\tWell there is an issue with chromsome ",chrom," and TMRCA!\n"))
  }
  
  ## LD scores
  ldsc <- read.table(paste0(i1,"/weights.UKB.",chrom,".l2.ldscore.gz"),
                     header=TRUE,stringsAsFactors=FALSE,colClasses=c("numeric","character","numeric","numeric"))[,-c(1,2)]
  ldsc  <- ldsc[which(ldsc[,"BP"]%in%d[,"BP"]),]
  ldsc  <- merge(d[,c("BP","SNP")],ldsc,by="BP",all.x=TRUE)
  whoNA <- which(is.na(ldsc[,"L2"]))
  if(length(whoNA)>0){
    ldsc[whoNA,"L2"] <- 0
  }
  check <- mean(ldsc[,"BP"]==d[,"BP"])
  if(check!=1){
    cat(paste0("\tWell there is an issue with chromsome ",chrom," and LDSC!\n"))
  }
  
  ## Combine all
  x               <- cbind(x,Allele_Age=aa[,"AlleleAge"],TMRCA=tmrca[,"TMRCA"],LDSC=ldsc[,"L2"],xgrp,HET=het)
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

## LD stratification
grp       <- cut(values[,"CONT_LDSC"],breaks = quantile(values[,"CONT_LDSC"]),include.lowest = TRUE)
xld       <- model.matrix(~grp-1); colnames(xld) <- paste0("LD_",1:4)
values    <- cbind(values,xld)
snpid     <- paste0(chrbp[,1],":",chrbp[,2])

## GWAS -- latest addition
traits <- read.table("gwas/Input.txt",stringsAsFactors = FALSE)[,1]
ccgwas <- c("numeric","character","numeric","character","character",rep("numeric",5))
getCHISQ <- function(trait){
  gwas   <- read.table(paste0("gwas/results/",trait,".fastGWA"),h=T,stringsAsFactors = F,colClasses = ccgwas)
  insnp  <- paste0(gwas[,"CHR"],":",gwas[,"POS"])
  whogw  <- which(insnp%in%snpid)
  CHISQ  <- (gwas[whogw,"BETA"]/gwas[whogw,"SE"])^2
  CHISQ
}
tt <- system.time( chisqGWAS <- do.call("cbind",mclapply(traits,getCHISQ,mc.cores=10)) )
colnames(chisqGWAS) <- paste0("CONT_",traits)
values    <- cbind(values,chisqGWAS)

cnvalues  <- colnames(values)
kqt_annot <- grep("CONT_",colnames(values))

for(k in kqt_annot){
  values[,k] <- values[,k] / max(values[,k])
}

## Change a few annotations
## TMRCA, Allele_Age and Nucleotide_Diversity_10kb
values[,"CONT_Nucleotide_Diversity_10kb"] <- 1-values[,"CONT_Nucleotide_Diversity_10kb"]
values[,"CONT_Allele_Age"]                <- 1-values[,"CONT_Allele_Age"]
values[,"CONT_TMRCA"]                     <- 1-values[,"CONT_TMRCA"]



write(colnames(values),"baselineLF2.2.UKB.annot-names.txt")
dt <- cbind(BP=dtannot[,2],values)
A  <- mclapply(1:22,function(chrom){
  print(chrom)
  write.table(dt[which(chrbp[,1]==chrom),],paste0("annot-files/baselineLF2.2.UKB.",chrom,".qced.annot"),quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
},mc.cores=10)

propAnnot <- colMeans(values)
save(list="propAnnot",file="propAnnot.RData")
Rannot <- cor(values)
save(list="Rannot",file="correlation-between-annotations.RData")

