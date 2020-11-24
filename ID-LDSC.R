# library(parallel)
# library(MASS)
# setwd("")
# load("propAnnot.RData")
# load("GWAS.RData")
# load("ukb-geno/rsq/info.RData")
# irsq   <- which(info[,"INFO"]>=0.3)
# rsq3   <- paste0(info[irsq,"CHR"],":",info[irsq,"BP"])
# rsq3   <- rsq3[which(!duplicated(rsq3))]
# 
# traits  <- read.table("gwas/Input.txt",stringsAsFactors = FALSE)[,1]
# nt      <- length(traits)
# snpid   <- rsq3
# who     <- which(GWAS[,"SNPID"]%in%snpid)
# GWAS    <- na.omit(GWAS[who,c("SNPID",paste0("Zf_",traits))])
# GWAS    <- GWAS[-which(duplicated(GWAS[,"SNPID"])),]
# 
# ## Get LD scores
# ccldsc <- c("character",rep("numeric",7))
# ldsc <- do.call("rbind",mclapply(1:22,function(chrom){
#   read.table(paste0("ldscores/chrom",chrom,".score.ld"),colClasses = ccldsc,header=TRUE)
# },mc.cores=10))
# 
# ldsc$SNPID <- paste0(ldsc[,"chr"],":",ldsc[,"bp"])
# dat        <- merge(GWAS,ldsc,by="SNPID")

load("data_for_ID-LDSC.RData")
dat$L2     <- dat[,"ldscore"] / nrow(dat)

grp <- NULL
for(k in 1:22){
  tmp <- ldsc[which(dat[,"chr"]==k),"bp"] 
  grp <- c(grp,paste0(k,":",cut(tmp,breaks = seq(min(tmp),max(tmp),by=1e7),include.lowest=TRUE)))
}
grp  <- as.numeric(as.factor(grp))
B    <- max(grp)
o    <- lapply(1:B,function(k) which(grp==k))



## Last QC
x  <- dat[,"L2"]
w  <- 1/dat[,"ldscore"]
X  <- cbind(1,x)

Y      <- as.matrix(dat[,paste0("Zf_",traits)]); colnames(Y) <- traits
X      <- cbind(I=1,L2=x)
w      <- w / mean(w)
V_X    <- X
for(j in 1:ncol(X)){
  V_X[,j] <- X[,j]*w
}
V_Y    <- Y
for(j in 1:ncol(Y)){
  V_Y[,j] <- Y[,j]*w
}

XTV_X    <- crossprod(X,V_X)
XTV_Y    <- crossprod(X,V_Y)
OLS      <- solve(XTV_X,XTV_Y) ## Singular!!!!!

JK <- lapply(1:B,function(k){
  who    <- o[[k]]
  xTV_x  <- crossprod(X[who,],V_X[who,])
  xTV_y  <- crossprod(X[who,],V_Y[who,])
  ols    <- solve(XTV_X-xTV_x,XTV_Y-xTV_y)
  return(ols)
})

## Main effect
Bg    <- OLS["L2",]
Bgs   <- do.call("rbind",lapply(1:B,function(k) JK[[k]]["L2",]))
Bg_se <- sapply(1:nt,function(j) sqrt( sum((1-1/B)*(Bgs[,j]-Bg[j])^2)  ))
Bg_pv <- pchisq(q=(Bg/Bg_se)^2,df=1,lower.tail = FALSE)

Bg_r  <- cbind.data.frame(Trait=traits,B=Bg,SE=Bg_se,Pval=Bg_pv)



## Constrained analysis
X  <- cbind(x)
Y      <- as.matrix(dat[,paste0("Zf_",traits)]); colnames(Y) <- traits
X      <- cbind(L2=x)
w      <- w / mean(w)
V_X    <- X
for(j in 1:ncol(X)){
  V_X[,j] <- X[,j]*w
}
V_Y    <- Y
for(j in 1:ncol(Y)){
  V_Y[,j] <- Y[,j]*w
}

XTV_X    <- crossprod(X,V_X)
XTV_Y    <- crossprod(X,V_Y)
OLS      <- solve(XTV_X,XTV_Y) ## Singular!!!!!

JK <- lapply(1:B,function(k){
  who    <- o[[k]]
  xTV_x  <- crossprod(X[who,],V_X[who,])
  xTV_y  <- crossprod(X[who,],V_Y[who,])
  ols    <- solve(XTV_X-xTV_x,XTV_Y-xTV_y)
  return(ols)
})

## Main effect
Bg    <- OLS["L2",]
Bgs   <- do.call("rbind",lapply(1:B,function(k) JK[[k]]["L2",]))
Bg_se <- sapply(1:nt,function(j) sqrt( sum((1-1/B)*(Bgs[,j]-Bg[j])^2)  ))
Bg_pv <- pchisq(q=(Bg/Bg_se)^2,df=1,lower.tail = FALSE)

Bg_r  <- cbind.data.frame(Trait=traits,B=Bg,SE=Bg_se,Pval=Bg_pv)


