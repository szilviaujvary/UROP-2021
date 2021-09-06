library("doMC")
library(genalg)

source("C:/Users/Asus/OneDrive/Dokumentumok/Imperial/UROP/Week 2/UROP week 2 gaCCA/genealgTwoCromosom.R")
source("C:/Users/Asus/OneDrive/Dokumentumok/Imperial/UROP/Week 2/UROP week 2 gaCCA/DIGLOM_gene-based.cca3.2.R")

genalg_calibrated = function (sizeX = 10,sizeY = 10, suggestions = NULL, popSize = 200, iters = 100, 
                    mutationChanceX = NA, elitism = NA, zeroToOneRatioX = 10,
                    mutationChanceY = NA,  zeroToOneRatioY = 10, monitorFunc = NULL, checkerFunc = NULL,
                    evalFunc = NULL, showSettings = FALSE, verbose = FALSE, effness=50) 
{
  
  
  if (is.null(evalFunc)) {
    stop("A evaluation function must be provided. See the evalFunc parameter.")
  }
  varsX = sizeX
  varsY = sizeY
  if (is.na(mutationChanceX)) {
    mutationChanceX = 1/(varsX + 1)
  }
  if (is.na(mutationChanceY)) {
    mutationChanceY = 1/(varsY + 1)
  }
  
  if (is.na(elitism)) {
    elitism = floor(popSize/5)
  }
  
  if (verbose) 
    cat("Testing the sanity of parameters...\n")
  if (popSize < 5) {
    stop("The population size must be at least 5.")
  }
  if (iters < 1) {
    stop("The number of iterations must be at least 1.")
  }
  if (!(elitism < popSize )) {
    stop("The population size must be greater than the elitism.")
  }
  if (showSettings) {
    if (verbose) 
      cat("The start conditions:\n")
    result = list(sizeX = sizeX, sizeY = sizeY, suggestions = suggestions, 
                  popSize = popSize, iters = iters, elitism = elitism, 
                  mutationChanceX = mutationChanceX, mutationChanceY=mutationChanceY)
    class(result) = "rbga"
    cat(summary(result))
  }
  else {
    if (verbose) 
      cat("Not showing GA settings...\n")
  }
  if (varsX > 0 & varsY>0) {
    if (!is.null(suggestions)) {
      if (verbose) 
        cat("Adding suggestions to first population...\n")
      population = matrix(nrow = popSize, ncol = varsX+varsY)
      suggestionCount = dim(suggestions)[1]
      for (i in 1:suggestionCount) {
        population[i, ] = suggestions[i, ]
      }
      if (verbose){ 
        cat("Filling others with random values in the given domains...\n")}
      for (child in (suggestionCount + 1):popSize) {
        
        popX = sample(c(rep(0, zeroToOneRatioX), 
                        1), varsX, rep = TRUE)
        popY = sample(c(rep(0, zeroToOneRatioY), 
                        1), varsY, rep = TRUE)
        population[child, ] = c(popX,popY)
        while (sum(population[child, 1:varsX]) == 0 | sum(population[child, (varsX+1):(varsX+varsY)]) == 0 ) {
          popX = sample(c(rep(0, zeroToOneRatioX), 
                          1), varsX, rep = TRUE)
          popY = sample(c(rep(0, zeroToOneRatioY), 
                          1), varsY, rep = TRUE)
          population[child, ] = c(popX,popY)
          
        }
      }
    }
    else {
      if (verbose) 
        cat("Starting with random values in the given domains...\n")
      population = matrix(nrow = popSize, ncol = varsX+varsY)
      for (child in 1:popSize) {
        popX = sample(c(rep(0, zeroToOneRatioX), 
                        1), varsX, rep = TRUE)
        popY = sample(c(rep(0, zeroToOneRatioY), 
                        1), varsY, rep = TRUE)
        population[child, ] = c(popX,popY)
        while (sum(population[child, 1:varsX]) == 0 | sum(population[child, (varsX+1):(varsX+varsY)]) == 0 ) {
          popX = sample(c(rep(0, zeroToOneRatioX), 
                          1), varsX, rep = TRUE)
          popY = sample(c(rep(0, zeroToOneRatioY), 
                          1), varsY, rep = TRUE)
          population[child, ] = c(popX,popY)
          
        }
      }
    }
    bestEvals = rep(NA, iters)
    meanEvals = rep(NA, iters)
    evalVals = rep(NA, popSize)
    it = 1
    percent_overlap=0
    genes_wrong = 50000000000000
    phens_wrong = 50000000000000
    while (percent_overlap < effness[1] | genes_wrong > effness[2] | phens_wrong > effness[3]) {
      #cat("effness not high enough...", "\n")
      if (verbose) 
        cat(paste("Starting iteration", it, "\n"))
      if (verbose) 
        cat("Calculating evaluation values... ")
      
      tmpVal = foreach(ii = 1:popSize,.combine=c)%do%{
        if (is.na(evalVals[ii])) {
          evalFunc(population[ii,]);
        }
        
      }
      idtmp = which(is.na(evalVals))
      evalVals[idtmp]=unlist(tmpVal)
      
      
      bestEvals[it] = min(evalVals)
      
      meanEvals[it] = mean(evalVals)
      if (verbose) 
        cat(" done.\n")
      if (!is.null(monitorFunc)) {
        if (verbose) 
          cat("Sending current state to rgba.monitor()...\n")
        result = list(type = "binary chromosome", sizeX = sizeX, sizeY,sizeY, 
                      popSize = popSize, iter = iter, iters = it, 
                      population = population, elitism = elitism,
                      mutationChanceX = mutationChanceX,mutationChanceY = mutationChanceY, evaluations = evalVals, 
                      best = bestEvals, mean = meanEvals)
        class(result) = "rbga"
        monitorFunc(result)
        check =checkerFunc(result)
        percent_overlap = check[1]
        genes_wrong = check[2]
        phens_wrong = check[3]
      
      }

        
      
      if (percent_overlap < effness[1] | genes_wrong > effness[2] | phens_wrong > effness[3]) {
        if (verbose) 
          cat("Creating next generation...\n")
        it= it+1
        cat("\n")
        cat("Starting it number", it)
        cat("\n")
        newPopulation = matrix(nrow = popSize, ncol = varsX+varsY)
        newEvalVals = rep(NA, popSize)
        if (verbose) 
          cat("  sorting results...\n")
        sortedEvaluations = sort(evalVals, index = TRUE)
        sortedPopulation = matrix(population[sortedEvaluations$ix, 
        ], ncol = varsX+varsY)
        if (elitism > 0) {
          if (verbose) 
            cat("  applying elitism...\n")
          newPopulation[1:elitism, ] = sortedPopulation[1:elitism, 
          ]
          newEvalVals[1:elitism] = sortedEvaluations$x[1:elitism]
        }
        if (varsX > 1) {  #varsX>1 and varsY >1??
          if (verbose) 
            cat("  applying crossover...\n")
          for (child in (elitism + 1):popSize) {
            parentProb = dnorm(1:popSize, mean = 0, sd = (popSize/3))
            parentIDs = sample(1:popSize, 2, prob = parentProb)
            parents = sortedPopulation[parentIDs, ]
            crossOverPointX = sample(0:varsX, 1)
            crossOverPointY = sample(0:varsY, 1)
            if (crossOverPointX == 0) {
              childX = parents[2, 1:varsX]
              #newEvalVals[child] = sortedEvaluations$x[parentIDs[2]]
            }
            else if (crossOverPointX == varsX) {
              childX = parents[1, 1:varsX]              
              #              newEvalVals[child] = sortedEvaluations$x[parentIDs[1]]
            }
            else {
              childX = c(parents[1,1:varsX ][1:crossOverPointX], 
                         parents[2,1:varsX ][(crossOverPointX + 1):varsX])
              
            }
            
            if (crossOverPointY == 0) {
              childY = parents[2, (varsX+1):(varsX+varsY)]
              #newEvalVals[child] = sortedEvaluations$x[parentIDs[2]]
            }
            else if (crossOverPointY == varsY) {
              childY = parents[1, (varsX+1):(varsX+varsY)]
              #newEvalVals[child] = sortedEvaluations$x[parentIDs[1]]
            }
            else {
              childY = c(parents[1,(varsX+1):(varsX+varsY) ][1:crossOverPointY], 
                         parents[2,(varsX+1):(varsX+varsY) ][(crossOverPointY + 1):varsY])
              
            }
            
            while(sum(childX)==0){
              childX = sample(c(rep(0, zeroToOneRatioX), 
                                1), varsX, rep = TRUE)
            }
            while(sum(childY)==0){
              childY = sample(c(rep(0, zeroToOneRatioY), 
                                1), varsY, rep = TRUE)
            }
            
            newPopulation[child,] = c(childX,childY)
            
            
          } #close for
        }  #close if vars >1
        else {
          if (verbose) 
            cat("  cannot crossover (#vars=1), using new randoms...\n")
          newPopulation[(elitism + 1):popSize, ] = sortedPopulation[sample(1:popSize, 
                                                                           popSize - elitism), ]
        }
        population = newPopulation
        evalVals = newEvalVals
        
        if (mutationChanceX > 0) {
          if (verbose) 
            cat("  applying mutations in X... ")
          mutationCount = 0
          for (object in (elitism + 1):popSize) {
            for (var in 1:varsX) {
              if (runif(1) < mutationChanceX) {
                population[object, var] = sample(c(rep(0, 
                                                       zeroToOneRatioX), 1), 1)
                mutationCount = mutationCount + 1
              }
            }
          }
          if (verbose) 
            cat(paste(mutationCount, "mutations applied in X\n"))
        }
        
        if (mutationChanceY > 0) {
          if (verbose) 
            cat("  applying mutations in Y... ")
          mutationCount = 0
          for (object in (elitism + 1):popSize) {
            for (var in (varsX+1):(varsX+varsY)) {
              if (runif(1) < mutationChanceY) {
                population[object, var] = sample(c(rep(0, 
                                                       zeroToOneRatioY), 1), 1)
                mutationCount = mutationCount + 1
              }
            }
          }
          if (verbose) 
            cat(paste(mutationCount, "mutations applied in Y\n"))
        }
        
        
        
      }
    }
  }
  result = list(type = "binary chromosome", sizeX = sizeX, sizeY,sizeY,popSize = popSize, 
                iters = it, suggestions = suggestions, population = population, 
                elitism = elitism, mutationChanceX = mutationChanceX, mutationChanceY =mutationChanceY,
                evaluations = evalVals,    best = bestEvals, mean = meanEvals, iters = it)
  class(result) = "rbga"
  return(result)
}







########################################################################################################

univariate_sign_check <- function(x, y, inds_x, inds_y){
  require(broom)
  list_x=c()
  list_y=c()
  
  for (i in inds_y){
    #list2=list()
    for (j in inds_x){
      
      if(summary(lm(y[, i] ~ x[, j]))$coefficients[8] <= 0.025){
        list_x <- c(list_x, j)
        list_y <- c(list_y, i)
      }
    }
  }
  list1 = c(unique(list_x), 0, unique(list_y))
  
  return(list1)
}

gaCCA_calibrated <- function(N, l.genes, phentype.l, zeroToOneRatioX, zeroToOneRatioY, x_prop, y_prop, effness=60, verbose=TRUE, Nutri=FALSE){ 
  x <- matrix( rnorm(N*l.genes,mean=0,sd=1), N, l.genes) 
  rownames(x) <- seq(1:N)
  colnames(x) <- seq(1:l.genes)
  
  inds_x <- sample(1:l.genes,l.genes*x_prop,replace=FALSE)
  beta <- rep(0,l.genes)
  beta[inds_x]=1+ rnorm(length(inds_x),0,2)
  
  
  y <- matrix( rnorm(N*phentype.l,mean=0,sd=1), N, phentype.l) 
  rownames(y) <- seq(1:N)
  colnames(y) <- seq(1:phentype.l)
  
  inds_y <- sample(1:phentype.l,phentype.l*y_prop,replace=FALSE)
  
  
  for(i in inds_y){
    k <- sample(1:l.genes, 1)
    y_1 <- x%*%beta +rnorm(N,0, 1)+ 2*x[, k]
    inds_x = c(inds_x, k)
    y[, i]= y_1
  }
  
  

  X = scale(x)
  Y = scale(y)
  sign_list <- univariate_sign_check(x, y, inds_x=inds_x, inds_y = inds_y)
  
  inds_x_new <- unique(sign_list[1: ((which(sign_list==0))-1)])
  inds_y_new <- unique(sign_list[((which(sign_list==0))+1):(length(sign_list))])
  
  cat("ground truth gene associations:", inds_x_new)
  cat("\n")
  cat("ground truth phen associations:", inds_y_new)
  cat("\n")
  ## generate filtered phenotype matrix
  
  
  # prepare parallel job
  registerDoMC(4)
  
  options(cores=4)
  getDoParWorkers()
  
  
  
  #reps = reps # number of rules to be extracted
  
  # define fitness function 
  # calculate the gene based association value for gene /phenotype selection
  fitness3 = function(ids){

    dimX = l.genes
    
    dimY = phentype.l
    
    
    #verbose = F
    
    
    idsX = ids[1:dimX]
    idsY = ids[(dimX+1):(dimX+dimY)]
    
    if(all(idsX==0) | all(idsY==0)==TRUE){
      return(30)}
    else{
      #selGenes = X[idsX]
      #snpsIdx = which(snps_locusID[,2]%in%selGenes)    
      
      #ix = as.logical(ids)
      #print("p value:")
      return(gene.assocGenSel(as.matrix(X[,(which(idsX==1))]),as.matrix(Y[,(which(idsY==1))]), verbose==T))
    }  
  }   

  
   # total number of iterations of the GA

  
  genes <- colnames(X)
  phenotypes <- colnames(Y)
  
  #monitor function to control the population homogeneity using hamming distance
  monitor <- function(obj) {
    minEval = min(obj$evaluations);
    #plot(obj, type="hist");
    pop = obj$population
    if(verbose==TRUE){
      cat("\n")
      cat("Hamming distance:", sum(e1071::hamming.distance(pop)))
      cat("\n")
      cat("Min. p-value:", minEval)
      cat("\n")
      bestpop = pop[which.min(obj$evaluations),]
      cat("Assoc gene numbers:", genes[as.logical(bestpop[1:l.genes])])
      cat("\n")
      cat("Assoc phenotype numbers:",phenotypes[as.logical(bestpop[(l.genes+1):(l.genes+phentype.l)])])
      cat("\n")
    }
    
    
  }
  
  list2 = list()
  
  # population size for the gene/phenotype selection approach
  populationSize = 1000
  mutationChanceX = 1/l.genes  #number of total genes dependent mutation chance for gene selection
  
  mutationChanceY= 1/phentype.l #number of total phenotype dependent mutation chance for phen selection
  
 
  checkerFunc <- function(obj){
    bestInd = obj$population[which.min(obj$evaluations),]
    
    geneInd = bestInd[1:l.genes]
    pheneInd = bestInd[(l.genes+1):(l.genes+phentype.l)]
    
    
    phenInd = bestInd[(l.genes+1):(l.genes+phentype.l)]
    
    
    genes_overlap = intersect(which(geneInd=='1'), inds_x_new)
    percent_true_pos_g = length(genes_overlap)/length(inds_x_new)*100
    
    false_pos_g = setdiff(which(geneInd=='1'), inds_x_new)
    percent_false_pos_g = length(false_pos_g)/length(inds_x_new)*100
    false_neg_g = setdiff(inds_x_new, which(geneInd=='1'))
    percent_false_neg_g = length(false_neg_g)/length(inds_x_new)*100
    
    phens_overlap = intersect(which(phenInd=='1'), inds_y_new)
    percent_true_pos_p = length(phens_overlap)/length(inds_y_new)*100
    
    false_pos_p = setdiff(which(phenInd=='1'), inds_y_new)
    percent_false_pos_p = length(false_pos_p)/length(inds_y_new)*100
    false_neg_p = setdiff(inds_y_new, which(phenInd=='1'))
    percent_false_neg_p = length(false_neg_p)/length(inds_y_new)*100
    
    effness=rep(0, 3)
    effness[1]= min(percent_true_pos_g, percent_true_pos_p)
    effness[2]= percent_false_pos_g+percent_false_neg_g
    effness[3]= percent_false_pos_p+percent_false_neg_p
    
    cat("Percent overlap genes:", percent_true_pos_g)
    cat("\n")
    cat("Percent overlap phens:", percent_true_pos_p)
    cat("\n")
    return(effness)
  }
  


  rbga.results = genalg_calibrated(sizeX=l.genes, sizeY=phentype.l,  zeroToOneRatioX=zeroToOneRatioX,zeroToOneRatioY=zeroToOneRatioY, elitism= 100, evalFunc=fitness3,popSize =populationSize , verbose=F,monitorFunc=monitor, checkerFunc = checkerFunc, effness=effness) 
    
  return(rbga.results$iters)
  
}






gaCCA_calibrated(500, 500, 500, 250, 250, 0.02, 0.025, effness=c(70, 40, 40), verbose=TRUE, Nutri=FALSE)





