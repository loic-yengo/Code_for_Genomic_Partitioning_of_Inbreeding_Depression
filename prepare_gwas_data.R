setwd("")

## Bim file
ccb <- c("numeric","character","numeric","numeric","character","character")
bim <- do.call("rbind",mclapply(1:22,function(chrom){
  read.table(paste0("ukb-geno/chrom",chrom,".bim"),stringsAsFactors=FALSE,colClasses=ccb)[,c(1,4)]
},mc.cores=10))
colnames(bim) <- c("CHR","BP")
snpid <- paste0(bim[,"CHR"],":",bim[,"BP"])

## Get all the Chisq
ccgwas <- c("numeric","character","numeric","character","character",rep("numeric",5))
getCHISQ <- function(trait){
  gwas   <- read.table(paste0("gwas/results/",trait,".fastGWA"),h=T,stringsAsFactors = F,colClasses = ccgwas)
  CHISQ  <- (gwas[,"BETA"]/gwas[,"SE"])^2
  CHISQ  <- (CHISQ-1)/gwas[,"N"]
  return(CHISQ)
}
traits <- read.table("gwas/Input.txt",stringsAsFactors = FALSE)[,1]
CHISQ  <- do.call("cbind",mclapply(traits,getCHISQ,mc.cores=10))
colnames(CHISQ) <- paste0("CHISQ_",traits)

## Get all Z_F
ccadddom <- c(rep("numeric",2),rep("character",5),rep("numeric",5),"character")
getadgwas <- function(trait){
  print(trait)
  Zf_t <- do.call("c",mclapply(1:22,function(chrom){
    gwas   <- read.table(paste0("gwas/add-dom/chrom",chrom,".",trait,".glm.linear"),h=T,stringsAsFactors = F,colClasses = ccadddom,comment.char = "!")
    gwas   <- gwas[which(gwas[,"TEST"]=="DOMDEV"),]
    Zf     <- -gwas[,"T_OR_F_STAT"] / sqrt(gwas[,"OBS_CT"])
    return(Zf)
  },mc.cores=10))
  return(Zf_t)
}
Zfm <- do.call("cbind",lapply(traits,getadgwas))
colnames(Zfm) <- paste0("Zf_",traits)

GWAS       <- cbind.data.frame(SNPID=snpid,CHISQ,Zfm)
save(list='GWAS',file="GWAS.RData")



