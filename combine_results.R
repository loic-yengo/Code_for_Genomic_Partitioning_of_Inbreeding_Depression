#load("datAll.RData")
#an <- readLines("baselineLF2.2.UKB.annot-names.txt")
na <- length(an)
traits <- colnames(dat)[-(1:2)]
N <- 348501
g <- matrix(0,nrow=N,ncol=2)
f <- n <- matrix(0,nrow=N,ncol=na)
for(chrom in 1:22){
  print(chrom)
  g <- g + as.matrix(read.table(paste0("ibc_combined/g.chrom",chrom),stringsAsFactors=F,colClasses="numeric"))
  f <- f + as.matrix(read.table(paste0("ibc_combined/m.chrom",chrom,".f"),stringsAsFactors=F,colClasses="numeric"))
  n <- n + as.matrix(read.table(paste0("ibc_combined/m.chrom",chrom,".n"),stringsAsFactors=F,colClasses="numeric"))
}
fx <- f/n
fg <- g[,1]/g[,2]
summary(fg)
iid <- read.table("ukb-geno/chrom22.fam")[,2]

save(list=c("f","n","g","fg","iid","N","traits","dat","an"),file=paste0("combined_common_",na,"annots.RData"))

