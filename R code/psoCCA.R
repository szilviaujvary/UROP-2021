psoCCA <- function(x, y, zeroToOneRatioX, zeroToOneRatioY, it=100, reps=1){
  K=ncol(x)
  L=ncol(y)
  fitness3 = function(ids){
  
    idsX = ids[1:K]
    idsY = ids[(K+1):(L+K)]
  #print(ids)
  #print(which(ids>=0.9))
  
    if(all(idsX<=0.9) | all(idsY<=0.9)==TRUE){
      return(30)}
    else{
    
      p <- gene.assocGenSel(as.matrix(x[,(which(idsX>=0.9))]),as.matrix(y[,(which(idsY>=0.9))]), verbose==T)
      return(p)
    }  
  }
  #result <- genalg3(sizeX=K, sizeY=L, evalFunc=fitness3)
  popX = sample(c(rep(0, zeroToOneRatioX), 1), K, rep = TRUE)
  popY = sample(c(rep(0, zeroToOneRatioY),1), L, rep = TRUE)
  par1=c(popX, popY)
  #print(which(par1==1))

  library(foreach)
  library(Rmpi)
  library(doMPI)
  
  cl <- startMPIcluster(count=8)
  registerDoMPI(cl)
  

 output <- foreach(j = 1:reps)%dopar%{
    source("C:/Users/szilv/OneDrive/Dokumentumok/Imperial/UROP/Week 6/PSO/R/psoptim.R")
    source("C:/Users/szilv/OneDrive/Dokumentumok/Imperial/UROP/Week 1-2/UROP week 2 gaCCA/DIGLOM_gene-based.cca3.2.R")
    
    result <- psoptim(par=par1, fn=fitness3, lower=rep(0, (K+L)), upper=rep(1, (K+L)), it=it)
    
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
  #if(out$p_value==0){
    
  #  x=as.matrix(x[, genes])
  #  colnames=genes
    
  #  y=as.matrix(x[, phens])
  #  colnames=phens
    
  #  closeCluster(cl)         
  #  return(psoCCA(x, y, zeroToOneRatioX=length(genes), zeroToOneRatioY=length(phens), it=100, reps=reps))
#  }
#  else{
    closeCluster(cl)
    return(output)
    #}
}  
  
#set.seed(1234567)





