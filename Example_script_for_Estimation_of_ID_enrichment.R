## This scripts provide an example to estimate the enrichment of inbreeding depression (ID)
## using a SNP-based measure F. This example is based on the inbreeding measure FUNI.
## For simplicity, we assume no genetic correlation between traits.

library(MASS)
## Expected enrichemnt is M/m = 10.
## Simulations parameters 
N  <- 10000 # sample size
M  <-  5000 # Total number of SNPs
m  <-   500 # Number of SNPs inside the effect annotation
p  <- rep(0.2,M) # allele frequency
h2 <-   0.5 # Heritability 
nt <-    10 # Numner of traits
nt_sq <- nt^2
B  <- runif(nt,min=-10,max=-3) # Inbreeding depression parameter for each trait
R  <- outer(1:nt,1:nt,function(i,j) 0.5**abs(i-j)) ## Phenotypic correlation betwen traits

## Data simulation
h  <- 2*p*(1-p)
s  <- sqrt(h)
X  <- do.call("cbind",lapply(1:M,function(j) rbinom(N,2,p[j])))
Fx3<- do.call("cbind",lapply(1:M,function(j) (X[,j]*X[,j]-(1+2*p[j])*X[,j]+2*p[j]^2)/h[j] ))
Fx2<- do.call("cbind",lapply(1:M,function(j) 1-X[,j]*(2-X[,j])/h[j]))
fh3<- rowMeans(Fx3)
fh2<- c(Fx2[,1:m]%*%h[1:m])/sum(h[1:m]) ## Inbreeding coefficients at causal SNps (in the annotation)
Z  <- do.call("cbind",lapply(1:m,function(j) (X[,j]-2*p[j])/s[j]))
Y  <- mvrnorm(n=N,mu=rep(0,nt),Sigma = R)
for(k in 1:nt){
  b     <- rnorm(m,mean=0,sd=sqrt(h2/m))
  g     <- c(Z%*%b)
  Y[,k] <- B[k] * fh2 + g + Y[,k]*sqrt(1-B[k]*B[k]/m-h2)
}

## Analysis
Fg  <- rowMeans(Fx3)       ## genome inbreeding - FUNI
Fk  <- rowMeans(Fx3[,1:m]) ## Annotation specific inbreeding
Pik <- m / M               ## proportion of SNPs in annotations
Dk  <- Pik * (Fk - Fg) / (1-Pik)

traits      <- paste0("Trait",1:nt)
colnames(Y) <- traits
dataFinal   <- cbind.data.frame(Y,Fg=Fg,Dk=Dk)

analyseSingleTrait <- function(trait){
  frm    <- as.formula(paste0(trait,"~Fg + Dk"))
  mod    <- summary(lm(frm,data=dataFinal))$coefficients[-1,]
  Bt     <- mod[,1]; SE <- mod[,2]; Pval <- mod[,4]
  Ek     <- 1 + Bt[2]/Bt[1]
  if(Ek>1){
    Direction="Enrichment"
    Fold = Ek
  }else{
    Direction="Depletion"
    Fold = 1-(Ek-1) * propAnnot[an[k]] / (1-propAnnot[an[k]])
  }
  
  Result <- cbind.data.frame(ID_genome=Bt[1],SE_ID_genome=SE[1],P_ID_genome=Pval[1],
                             Tau_k=Bt[2],SE_Tau_k=SE[2],P_Tau_k=Pval[2],
                             Delta_k=Fold,Trait=trait,Direction=Direction)
  return(Result)
}

## Trait-specific enrichment
Rpm  <- do.call("rbind",lapply(traits,analyseSingleTrait))

## Average across traits
Ry  <- cor(Y) # sample correlation matrix between traits
avg <- function(Rp){
  mBg <- mean(Rp[,"ID_genome"]); mBk <- mean(Rp[,"Tau_k"])
  vBk <- sum(tcrossprod(Rp[,"SE_Tau_k"]) * Ry) / (nt_sq); zBk <- mBk/sqrt(vBk)
  pBk <- pchisq(q=mBk*mBk/vBk,df=1,lower.tail=FALSE); mEk <- 1 + mBk/mBg; sEk = mEk/abs(zBk)
  if(mEk > 1){
    Direction="Enrichment"
    Fold = mEk
  }else{
    Direction="Depletion"
    Fold = 1-(mEk-1) * propAnnot[an] / (1-propAnnot[an])
  }
  cbind.data.frame(Delta_bar_k=Fold,mean_Tau_k=mBk,se_mean_Tau_k=sqrt(vBk),Delta_bar_k=sEk,P_Delta_bar_k=pBk,mean_ID_genome=mBg,Direction=Direction)
}
avgRpm <- avg(Rpm)
print(avgRpm)






