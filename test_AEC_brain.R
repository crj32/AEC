load("C:/Users/chris/Dropbox/Spectral_code_and_figures/19032019/spectralpaperfigures/suppdata/Supplementary_file1_V2.RData")

#load("C:/Users/chris/Dropbox/Spectral_manuscript/bioinformatics/Supplementary_file2.RData")

library(survival)
source('C:/Users/chris/Desktop/M3C2V4.R')
library(doSNOW)
library(ggplot2)
library(corpcor)

## single omic
mydata <- breast[[1]]
clinicali <- clinical_bladder

## synthetic
x <- clusterlab::clusterlab(centers=3,r=1.25)
mydata <- x$synthetic_data

## 
set.seed(123)
mydata <- read.csv('unifiedScaledFiltered.tsv',sep='\t')
mydata <- Spectrum::brain[[2]]

#
r <- M3C(mydata,method=2,iters=20,cores=4)

assignments <- r$realdataresults[[6]]$assignments

## test clinical
clinicali$Death <- as.numeric(as.character(clinicali$Death))
coxFit <- coxph(Surv(time = Time, event = Death) ~ as.factor(assignments), data = clinicali, ties = "exact")
coxresults <- summary(coxFit)
print(coxresults$logtest[3])

### code for tuning lambda

##
likelihood <- function(m,m2){
  m <- m[upper.tri(m)]
  m2 <- m2[upper.tri(m2)]
  L <- sum(log(m2*m+(1-m)*(1-m2))) # m = ground truth
  #L <- sum((m-m2)^2)
  return(L)
}

getl <- function(r,krange){
  for (k in krange){
    ## get ground truth
    clust <- cluster::pam(t(r$realdataresults[[k]]$ordered_data),k=k)
    clus <- clust$clustering
    rowNum <- nrow(t(mydata))
    S <- matrix(0, rowNum, rowNum)
    for (j in 1:k) {
      X <- rep(0, rowNum)
      X[which(clus == j)] <- 1
      S <- S + X %*% t(X)
    }
    ## get perturbed probs
    rm <- r$realdataresults[[k]]$consensus_matrix
    rm <- as.matrix(rm)
    ##
    message(paste('log-likelihood',likelihood(S,rm)))
  }
}

## diff matrix
for (lambda in seq(0.05,0.1,by=0.01)){
  message(paste('tuning lambda:',lambda))
  ent <- log(r$scores$ENTROPY)+lambda*r$scores$K
  print(ent)
  min <- which.min(ent)+1
  getl(r,min)
}

for (k in seq(2,10)){
  ## get ground truth
  clust <- cluster::pam(t(r$realdataresults[[k]]$ordered_data),k=k)
  clus <- clust$clustering
  rowNum <- nrow(t(r$realdataresults[[k]]$ordered_data))
  S <- matrix(0, rowNum, rowNum)
  for (j in 1:k) {
    X <- rep(0, rowNum)
    X[which(clus == j)] <- 1
    S <- S + X %*% t(X)
  }
  ## get perturbed probs
  rm <- r$realdataresults[[k]]$consensus_matrix
  rm <- as.matrix(rm)
  ##
  print(likelihood(S,rm))
}
