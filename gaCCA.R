### genetic algorithm selection of best CCA genes/phen subsets

library(genalg)

#source("C:/Users/Asus/OneDrive/Dokumentumok/Imperial/UROP/Week 2/UROP week 2 gaCCA/genealgTwoCromosom.R")
#source("C:/Users/Asus/OneDrive/Dokumentumok/Imperial/UROP/Week 2/UROP week 2 gaCCA/DIGLOM_gene-based.cca3.2.R")

gaCCA <- function(X, Y, zeroToOneRatioX, zeroToOneRatioY, reps=2, it=10, verbose=TRUE)
  {
  l.genes = ncol(X)
  phentype.l = ncol(Y)
  
  x = scale(X)
  Y = scale(Y)


reps = reps # number of rules to be extracted

# define fitness function 
# calculate the gene based association value for gene /phenotype selection
  fitness3 = function(ids){
    dimX = l.genes
    
    dimY = phentype.l

    idsX = ids[1:dimX]
    idsY = ids[(dimX+1):(dimX+dimY)]
      
    if(all(idsX==0) | all(idsY==0)==TRUE){
      return(30)}
    else{
      return(gene.assocGenSel(as.matrix(X[,(which(idsX==1))]),as.matrix(Y[,(which(idsY==1))]), verbose==T))
    }  
 }

 it = it  # total number of iterations of the GA

 genes <- colnames(X)
 phenotypes <- colnames(Y)
 
 #monitor function to control the population homogeneity using hamming distance
  monitor <- function(obj) {
    minEval = min(obj$evaluations);
    pop = obj$population
    if(verbose==TRUE){
      cat(sum(e1071::hamming.distance(pop)))
      cat("\n")
      cat(minEval)
      cat("\n")
      bestpop = pop[which.min(obj$evaluations),]
      cat(genes[as.logical(bestpop[1:l.genes])])
      cat("\n")
      cat(phenotypes[as.logical(bestpop[(l.genes+1):(l.genes+phentype.l)])])
      cat("\n")
    }

  }


  
# population size for the gene/phenotype selection approach
populationSize = 1000
mutationChanceX = 1/l.genes  #number of total genes dependent mutation chance for gene selection

mutationChanceY= 1/phentype.l #number of total phenotype dependent mutation chance for phen selection

library(foreach)
library(Rmpi)
library(doMPI)

cl <- startMPIcluster(count=8)
registerDoMPI(cl)

output <- foreach(j =1:reps)%dopar%{
  source("C:/Users/szilv/OneDrive/Dokumentumok/Imperial/UROP/Week 1-2/UROP week 2 gaCCA/genealgTwoCromosom.R")
  source("C:/Users/szilv/OneDrive/Dokumentumok/Imperial/UROP/Week 1-2/UROP week 2 gaCCA/DIGLOM_gene-based.cca3.2.R")
  
  	rbga.results = genalg3(sizeX=l.genes, sizeY=phentype.l,  zeroToOneRatioX=zeroToOneRatioX,zeroToOneRatioY=zeroToOneRatioY,iters=it, elitism= 100, evalFunc=fitness3,popSize =populationSize , verbose=F,monitorFunc=NULL) 

    
    pvalue <- rbga.results$best[it]
 
    bestInd = rbga.results$population[which.min(rbga.results$evaluations),]
    
    geneInd = bestInd[1:l.genes]
    #pheneInd = bestInd[(l.genes+1):(l.genes+phentype.l)]
    #list1[["genesInd"]] = geneInd
    genes = genes[as.logical(geneInd)]

    phenInd = bestInd[(l.genes+1):(l.genes+phentype.l)]
    #list1[["phenInd"]] = phenInd
    phens = phenotypes[as.logical(phenInd)]

    
    
    out = list(
      p_value = pvalue,
      genes = which(geneInd==1),
      phens = which(phenInd==1)
    )
    #output[[j]]=out
    out
}

closeCluster(cl)
return(output)

}

#function to check p-value for a particular selection manually

gaCCA_manual <- function(X, Y, ids)
{
  l.genes = ncol(X)
  phentype.l = ncol(Y)
  
  x = scale(X)
  Y = scale(Y)
  
  # define fitness function 
  # calculate the gene based association value for gene /phenotype selection
  fitness3 = function(ids){
    dimX = l.genes
    
    dimY = phentype.l

    idsX = ids[1:dimX]
    idsY = ids[(dimX+1):(dimX+dimY)]
    
    if(all(idsX==0) | all(idsY==0)==TRUE){
      return(30)}
    else{
      return(gene.assocGenSel(as.matrix(X[,(which(idsX==1))]),as.matrix(Y[,(which(idsY==1))]), verbose==T))
    }  
  }
  
  return(fitness3(ids))
}




