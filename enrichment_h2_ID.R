library(parallel)
library(MASS)
load("propAnnot.RData")
load("GWAS.RData")
annots <- do.call("rbind",mclapply(1:22,function(chrom){
  print(chrom)
  tmp <- read.table(paste0("annot-files/baselineLF2.2.UKB.",chrom,".qced.annot"),sep="\t",colClasses = "numeric")
  colnames(tmp) <- c("BP",names(propAnnot))
  tmp <- cbind.data.frame(CHR=chrom,tmp)
  return(tmp)
},mc.cores=10))

traits  <- read.table("gwas/Input.txt",stringsAsFactors = FALSE)[,1]
nt      <- length(traits)
an      <- colnames(annots)
out     <- c(1,2,grep("MAF",an),grep("HET",an),grep("LDSC",an),grep("LD_",an),
             sapply(traits,function(trait) grep(trait,an)))
indexes <- which(!an%in%an[out])
cn      <- an[indexes]

snpid <- paste0(annots[,"CHR"],":",annots[,"BP"])
who   <- which(GWAS[,"SNPID"]%in%snpid)
table(GWAS[who,"SNPID"]==snpid)
GWAS <- GWAS[who,]

## Get LD scores
ccldsc <- c("character",rep("numeric",7))
ldsc <- do.call("rbind",mclapply(1:22,function(chrom){
  read.table(paste0("ldscores/chrom",chrom,".score.ld"),colClasses = ccldsc,header=TRUE)
},mc.cores=10))
snpid_ldsc <- paste0(ldsc[,"chr"],":",ldsc[,"bp"])
who_ldsc   <- which(snpid_ldsc%in%GWAS[,"SNPID"])
table(GWAS[,"SNPID"]==snpid_ldsc[who_ldsc])
ldsc <- ldsc[who_ldsc,]

## Final cleaning?
naDom  <- which( rowSums(is.na(GWAS[,paste0("Zf_",traits)])) > 0 )
annots <- annots[-naDom,]
GWAS   <- GWAS[-naDom,]
ldsc   <- ldsc[-naDom,]
M      <- nrow(ldsc)

mldsc_g  <- mean(ldsc[,"ldscore"])
mchisq_g <- colMeans(GWAS[,paste0("CHISQ_",traits)]); names(mchisq_g) <- traits
m2pqd_g  <- colMeans(GWAS[,paste0("Zf_",traits)]); names(m2pqd_g) <- traits
mldsc_a  <- as.numeric(crossprod(as.matrix(annots[,cn]),ldsc[,"ldscore"])) / M; names(mldsc_a) <- cn

mchisq_t <- crossprod(as.matrix(annots[,cn]),as.matrix(GWAS[,paste0("CHISQ_",traits)])) / M
colnames(mchisq_t) <- traits
rownames(mchisq_t) <- cn

m2pqd_t <- crossprod(as.matrix(annots[,cn]),as.matrix(GWAS[,paste0("Zf_",traits)])) / M
colnames(m2pqd_t) <- traits
rownames(m2pqd_t) <- cn

enr_h2 <- do.call("cbind",lapply(traits,function(trait) (mchisq_t[,trait]/mldsc_a)/(mchisq_g[trait]/mldsc_g)))
colnames(enr_h2) <- traits
rownames(enr_h2) <- cn

enr_ID <- do.call("cbind",lapply(traits,function(trait) (m2pqd_t[,trait]/mldsc_a)/(m2pqd_g[trait]/mldsc_g)))
colnames(enr_ID) <- traits
rownames(enr_ID) <- cn

ID <- colSums(GWAS[,paste0("Zf_",traits)],na.rm=TRUE) / mldsc_g
save(list=c("ID","enr_h2","enr_ID"),file="enrichment_h2_ID.RData")

#### ---
## Should we get s.e.?
## Block Jackknife
grp <- NULL
for(k in 1:22){
  tmp <- ldsc[which(ldsc[,"chr"]==k),"bp"] 
  grp <- c(grp,paste0(k,":",cut(tmp,breaks = seq(min(tmp),max(tmp),by=1e7),include.lowest=TRUE)))
}
grp  <- as.numeric(as.factor(grp))
B    <- max(grp)
tab  <- as.numeric(table(grp))
mem  <- model.matrix(~factor(grp)-1)
Moth <- M-tab

indexGrp <- lapply(1:B,function(k) which(grp==k))
set_of_mldsc_g  <- (sum(ldsc[,"ldscore"])-as.numeric(crossprod(mem,ldsc[,"ldscore"])))/Moth

## Set of mchisq
set_of_mchisq_g <- crossprod(mem,as.matrix(GWAS[,paste0("CHISQ_",traits)]))
colnames(set_of_mchisq_g) <- traits
for(u in traits){
  set_of_mchisq_g[,u] <- (sum(GWAS[,paste0("CHISQ_",u)])-set_of_mchisq_g[,u])/Moth
}

## Set of m2pqd
set_of_m2pqd_g  <- crossprod(mem,as.matrix(GWAS[,paste0("Zf_",traits)])); colnames(set_of_m2pqd_g) <- traits
for(u in traits){
  set_of_m2pqd_g[,u] <- (sum(GWAS[,paste0("Zf_",u)])-set_of_m2pqd_g[,u])/Moth
}

## Set of mldsc
A <- as.matrix(annots[,cn])
for(u in cn){
  A[,u] <- A[,u] * ldsc[,"ldscore"]
}
set_of_mldsc_a  <-  crossprod(mem,A)
colnames(set_of_mldsc_a) <- cn
for(u in cn){
  set_of_mldsc_a[,u] <- (sum(A[,u])-set_of_mldsc_a[,u]) / Moth
}

## Set of mchisq_t
C <- crossprod(as.matrix(annots[,cn]),as.matrix(GWAS[,paste0("CHISQ_",traits)]))
set_of_mchisq_t <- lapply(1:B,function(k){
  c_k <- crossprod(as.matrix(annots[indexGrp[[k]],cn]),
                   as.matrix(GWAS[indexGrp[[k]],paste0("CHISQ_",traits)]))
  A <- (C-c_k)/Moth[k]
  colnames(A) <- traits
  rownames(A) <- cn
  return(A)
})

C <- crossprod(as.matrix(annots[,cn]),as.matrix(GWAS[,paste0("Zf_",traits)]))
set_of_m2pqd_t <- lapply(1:B,function(k){
  c_k <- crossprod(as.matrix(annots[indexGrp[[k]],cn]),
                   as.matrix(GWAS[indexGrp[[k]],paste0("Zf_",traits)]))
  A <- (C-c_k)/Moth[k]
  colnames(A) <- traits
  rownames(A) <- cn
  return(A)
})

enr_h2_jk <- function(k){
  Ain <- do.call("cbind",lapply(traits,function(trait){
    (set_of_mchisq_t[[k]][,trait]/set_of_mldsc_a[k,])/(set_of_mchisq_g[k,trait]/set_of_mldsc_g[k])
  }))
  colnames(Ain) <- traits
  rownames(Ain) <- cn
  Ain
}
enr_h2_se <- sqrt( (1-1/B)*Reduce("+",lapply(1:B,function(k) (enr_h2_jk(k)-enr_h2)^2)) )
enr_h2_pv <- pchisq(q=((enr_h2-1)/enr_h2_se)^2,df=1,lower.tail = FALSE)

enr_ID_jk <- function(k){
  Ain <- do.call("cbind",lapply(traits,function(trait){
    (set_of_m2pqd_t[[k]][,trait]/set_of_mldsc_a[k,])/(set_of_m2pqd_g[k,trait]/set_of_mldsc_g[k])
  }))
  colnames(Ain) <- traits
  rownames(Ain) <- cn
  Ain
}
enr_ID_se <- sqrt( (1-1/B)*Reduce("+",lapply(1:B,function(k) (enr_ID_jk(k)-enr_ID)^2)) )
enr_ID_pv <- pchisq(q=((enr_ID-1)/enr_ID_se)^2,df=1,lower.tail = FALSE)


## SE for the average h2
avg_enr_h2 <- rowMeans(enr_h2)
enr_avg_h2_jk <- function(k){
  Ain <- do.call("cbind",lapply(traits,function(trait){
    (set_of_mchisq_t[[k]][,trait]/set_of_mldsc_a[k,])/(set_of_mchisq_g[k,trait]/set_of_mldsc_g[k])
  }))
  colnames(Ain) <- traits
  rownames(Ain) <- cn
  rowMeans(Ain)
}
enr_avg_h2_se <- sqrt( (1-1/B)*Reduce("+",lapply(1:B,function(k) (enr_avg_h2_jk(k)-avg_enr_h2)^2)) )
enr_avg_h2_pv <- pchisq(q=((avg_enr_h2-1)/enr_avg_h2_se)^2,df=1,lower.tail = FALSE)
avg_h2        <- cbind(Enr=avg_enr_h2,SE=enr_avg_h2_se,Pval=enr_avg_h2_pv); rownames(avg_h2) <- names(avg_enr_h2)


## SE for the average ID
avg_enr_ID <- rowMeans(enr_ID)
enr_avg_ID_jk <- function(k){
  Ain <- do.call("cbind",lapply(traits,function(trait){
    (set_of_m2pqd_t[[k]][,trait]/set_of_mldsc_a[k,])/(set_of_m2pqd_g[k,trait]/set_of_mldsc_g[k])
  }))
  colnames(Ain) <- traits
  rownames(Ain) <- cn
  rowMeans(Ain)
}
enr_avg_ID_se <- sqrt( (1-1/B)*Reduce("+",lapply(1:B,function(k) (enr_avg_ID_jk(k)-avg_enr_ID)^2)) )
enr_avg_ID_pv <- pchisq(q=((avg_enr_ID-1)/enr_avg_ID_se)^2,df=1,lower.tail = FALSE)
avg_ID        <- cbind(Enr=avg_enr_ID,SE=enr_avg_ID_se,Pval=enr_avg_ID_pv); rownames(avg_ID) <- names(avg_enr_ID)

## Correlation for each trait
cor_enr_h2_ID <- diag(cor(enr_h2,enr_ID))
cor_enr_h2_ID_jk <- function(k){
  ## h2
  Ain <- do.call("cbind",lapply(traits,function(trait){
    (set_of_mchisq_t[[k]][,trait]/set_of_mldsc_a[k,])/(set_of_mchisq_g[k,trait]/set_of_mldsc_g[k])
  }))
  colnames(Ain) <- traits
  rownames(Ain) <- cn
  X <- Ain
  
  ## ID
  Ain <- do.call("cbind",lapply(traits,function(trait){
    (set_of_m2pqd_t[[k]][,trait]/set_of_mldsc_a[k,])/(set_of_m2pqd_g[k,trait]/set_of_mldsc_g[k])
  }))
  colnames(Ain) <- traits
  rownames(Ain) <- cn
  Y <- Ain
  return(diag(cor(X,Y)))
}
cor_enr_h2_ID_se <- sqrt( (1-1/B)*Reduce("+",lapply(1:B,function(k) (cor_enr_h2_ID_jk(k)-cor_enr_h2_ID)^2)) )
cor_enr_h2_ID_pv <- pchisq(q=(cor_enr_h2_ID/cor_enr_h2_ID_se)^2,df=1,lower.tail = FALSE)
cor_enr_h2_ID    <- cbind(Enr=cor_enr_h2_ID,SE=cor_enr_h2_ID_se,Pval=cor_enr_h2_ID_pv)

## Correlation of average
cor_avg_enr_h2_ID <- cor(avg_enr_ID,avg_enr_h2)
cor_avg_enr_h2_ID_jk <- function(k){
  ## h2
  Ain <- do.call("cbind",lapply(traits,function(trait){
    (set_of_mchisq_t[[k]][,trait]/set_of_mldsc_a[k,])/(set_of_mchisq_g[k,trait]/set_of_mldsc_g[k])
  }))
  colnames(Ain) <- traits
  rownames(Ain) <- cn
  X <- rowMeans(Ain)
  
  ## ID
  Ain <- do.call("cbind",lapply(traits,function(trait){
    (set_of_m2pqd_t[[k]][,trait]/set_of_mldsc_a[k,])/(set_of_m2pqd_g[k,trait]/set_of_mldsc_g[k])
  }))
  colnames(Ain) <- traits
  rownames(Ain) <- cn
  Y <- rowMeans(Ain)
  return(cor(X,Y))
}
cor_avg_enr_h2_ID_se <- sqrt( (1-1/B)*sum(sapply(1:B,function(k) (cor_avg_enr_h2_ID_jk(k)-cor_avg_enr_h2_ID)^2)) )
cor_avg_enr_h2_ID_pv <- pchisq(q=(cor_avg_enr_h2_ID/cor_avg_enr_h2_ID_se)^2,df=1,lower.tail = FALSE)
cor_avg_enr_h2_ID    <- cbind(Enr=cor_avg_enr_h2_ID,SE=cor_avg_enr_h2_ID_se,Pval=cor_avg_enr_h2_ID_pv)

cor_enr_h2_ID        <- rbind.data.frame(AVERAGE=cor_avg_enr_h2_ID,cor_enr_h2_ID)

save(list=c("ID",
            "enr_h2","enr_h2_se","enr_h2_pv",
            "enr_ID","enr_ID_se","enr_ID_pv",
            "avg_ID","avg_h2","cor_enr_h2_ID"
            ),file="enrichment_h2_ID.RData")


