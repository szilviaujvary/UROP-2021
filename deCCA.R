deCCA <- function(x, y, zeroToOneRatioX, zeroToOneRatioY, it=100, reps=1){
  K=ncol(x)
  L=ncol(y)
  
  control <- function(ids){
    idsX = ids[1:K]
    idsY = ids[(K+1):(L+K)]
    #print(ids)
    
    if(length(which(idsX>=0.9))>10){
      inds <- which(order(idsX, decreasing=TRUE)>=0.9)
      idsX <- rep(0, K)
      idsX[inds]=1
    }
    
    if(length(which(idsY>=0.9))>10){
      inds <- which(order(idsY, decreasing=TRUE)>=0.9)
      idsY <- rep(0, L)
      idsY[inds]=1
    }
    ids <- c(idsX, idsY)
    return(ids)
    
    #cat("genes", which(idsX>=0.9), "\n")
    #cat("phens", which(idsY>=0.9), "\n")
  }
  fitness3 <- function(ids){
    
    idsX = ids[1:K]
    idsY = ids[(K+1):(L+K)]
    #print(ids)
    
    #if(length(which(idsX>=0.9))>10){
    #  inds <- sample(which(idsX>=0.9), 10)
    #  idsX <- rep(0, 10)
    #  idsX[inds] = 1
      
    #}
    
    #if(length(which(idsY>=0.9))>10){
    #  inds <- sample(which(idsY>=0.9), 10)
    #  idsY <- rep(0, 10)
    #  idsY[inds] = 1
    #}
    
    #cat("genes", which(idsX>=0.9), "\n")
    #cat("phens", which(idsY>=0.9), "\n")
    if(all(idsX<=0.9) | all(idsY<=0.9)==TRUE){
      return(30)}
    else{
      return(gene.assocGenSel(as.matrix(x[,(which(idsX>=0.9))]),as.matrix(y[,(which(idsY>=0.9))]), verbose==T))
    }  
  }
  popX = sample(c(rep(0, zeroToOneRatioX), 1), K, rep = TRUE)
  popY = sample(c(rep(0, zeroToOneRatioY),1), L, rep = TRUE)
  par1=c(popX, popY)
  #cat("\n start:", which(par1==1), "\n")
  

  library(foreach)
  library(Rmpi)
  library(doMPI)
  
  cl <- startMPIcluster(count=8)
  registerDoMPI(cl)

output <-  foreach(j = 1:reps)%dopar%{
    source("C:/Users/szilv/OneDrive/Dokumentumok/Imperial/UROP/Week 6/JDDEoptim.R")
    source("C:/Users/szilv/OneDrive/Dokumentumok/Imperial/UROP/Week 1-2/UROP week 2 gaCCA/DIGLOM_gene-based.cca3.2.R")
    result = JDDEoptim(rep(0, (L+K)), rep(1, (L+K)), fitness3, maxiter=it, add_to_init_pop =NULL)
    
    cols <- which(result$par>=0.9)
    genes <- cols[which(cols <= K)]
    phens_ind <- which(cols > K)
    phens = cols[phens_ind]-K
    out = list(
      p_value = result$value,
      genes = genes,
      phens = phens
    )
    #output[[j]]=out
    out
  }
  
  closeCluster(cl)
  return(output)
}  



