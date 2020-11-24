N <- 348501
a <- seq(0,348500,by=20500)
r <- NULL
for(i in 1:17){
  r <- rbind(r,c(a[i],a[i+1]-1))
}
r[17,2] <- N-1
input   <- cbind(rep(1:22,each=17),do.call("rbind",lapply(1:22,function(k) r)))
write.table(input,"input_parameters_sorted.txt",quote=F,row.names=F,col.names=F,sep="\t")

set.seed(1234)
input   <- input[sample(1:nrow(input)),]
write.table(input,"input_parameters_for_fannot.txt",quote=F,row.names=F,col.names=F,sep="\t")

