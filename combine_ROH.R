library(parallel)
setwd("")
XROH <- do.call("rbind",mclapply(1:22,function(k){
  load(paste0("XROH/chrom",k,".RData"))
  return(XROH)
},mc.cores=10))

load("inputForXROH.RData")
load("propAnnotForROH.RData")
thingsToAggregate <- c("NSNP",names(propAnnot))
ROHwgt  <- aggregate(XROH[,thingsToAggregate],by=list(XROH[,"IID"]),FUN=sum)

### scaling
scaling <- c(
  NSNP=sum(sapply(listBP,length)),
  colSums(do.call("rbind",lapply(listValues,colSums))),
  Allele_Age=sum(sapply(listAlleleAge,function(u) sum(u[,"AlleleAge"]))),
  TMRCA=sum(sapply(listTMRCA,function(u) sum(u[,"TMRCA"])))
)
for(u in names(scaling)){
  ROHwgt[,u] <- ROHwgt[,u] / scaling[u]
}
ROHunw <- ROHwgt
for(u in names(propAnnot)){
  ROHunw[,u] <- ROHunw[,u] - ROHunw[,"NSNP"]
  ROHwgt[,u] <- propAnnot[u] * ROHunw[,u] / (1-propAnnot[u])
}
colnames(ROHwgt)[1:2] <- colnames(ROHunw)[1:2] <- c("IID","FNROH")
save(list=c("propAnnot","ROHwgt","ROHunw"),file="ROHwgt.RData")

