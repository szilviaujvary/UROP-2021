JDDEoptim <-function (lower, upper, fn, constr = NULL, meq = 0, eps = 1e-05, 
          NP =1000
            #5*floor(10+2*sqrt(length(d)))
            , Fl = 0.1, Fu = 1, tau_F = 0.1, tau_CR = 0.1, 
          tau_pF = 0.1, jitter_factor = 0.001, tol = 0, maxiter = it, fnscale = 1e-15, compare_to = c("median", "max"), add_to_init_pop = NULL, 
          trace = FALSE, triter = 1, details = FALSE, ...) 
{
  handle.bounds <- function(x, u) {
    bad <- x > upper
    x[bad] <- 0.5 * (upper[bad] + u[bad])
    bad <- x < lower
    x[bad] <- 0.5 * (lower[bad] + u[bad])
    x
  }
  performReproduction <- function() {
    ignore <- runif(d) > CRtrial
    if (all(ignore)) 
      ignore[sample(d, 1)] <- FALSE
    trial <- if (runif(1) <= pFtrial) 
      X.base + Ftrial * (X.r1 - X.r2)
    else X.base + 0.5 * (Ftrial + 1) * (X.r1 + X.r2 - 2 * 
                                          X.base)
    trial[ignore] <- X.i[ignore]
    trial
  }
  which.best <- if (!is.null(constr)) 
    function(x) {
      ind <- TAVpop <= mu
      if (all(ind)) 
        which.min(x)
      else if (any(ind)) 
        which(ind)[which.min(x[ind])]
      else which.min(TAVpop)
    }
  else which.min
  compare_to <- match.arg(compare_to)
  d <- length(lower)
  if (length(upper) != d) 
    stop("'lower' must have same length as 'upper'")
  stopifnot(is.numeric(lower), is.numeric(upper), is.finite(lower), 
            is.finite(upper), lower <= upper, length(fnscale) == 
              1, is.finite(fnscale), fnscale > 0, is.function(fn))
  if (!is.null(constr)) {
    stopifnot(is.function(constr))
    stopifnot(length(meq) == 1, meq == as.integer(meq), 
              meq >= 0, is.numeric(eps), is.finite(eps), eps > 
                0)
    if (length(eps) == 1) 
      eps <- rep.int(eps, meq)
    else if (length(eps) != meq) 
      stop("eps must be either of length meq, or length 1")
  }
  stopifnot(length(NP) == 1, NP == as.integer(NP), length(Fl) == 
              1, is.numeric(Fl), length(Fu) == 1, is.numeric(Fu), 
            Fl <= Fu)
  stopifnot(length(tau_F) == 1, is.numeric(tau_F), 0 <= tau_F, 
            tau_F <= 1, length(tau_CR) == 1, is.numeric(tau_CR), 
            0 <= tau_CR, tau_CR <= 1, length(tau_pF) == 1, is.numeric(tau_pF), 
            0 <= tau_pF, tau_pF <= 1)
  if (!is.null(jitter_factor)) 
    stopifnot(length(jitter_factor) == 1, is.numeric(jitter_factor))
  stopifnot(length(tol) == 1, is.numeric(tol), length(maxiter) == 
              1, maxiter == as.integer(maxiter), length(triter) == 
              1, triter == as.integer(triter))
  if (!is.null(add_to_init_pop)) 
    stopifnot(NROW(add_to_init_pop) == d, is.numeric(add_to_init_pop), 
              add_to_init_pop >= lower, add_to_init_pop <= upper)
  
  child <- if (is.null(constr)) {
    expression({
      ftrial <- fn1(trial)
      if (ftrial <= fpop[i]) {
        pop[, i] <- trial
        fpop[i] <- ftrial
        F[, i] <- Ftrial
        CR[i] <- CRtrial
        pF[i] <- pFtrial
      }
    })
  }
  else if (meq > 0) {
    expression({
      htrial <- constr1(trial)
      TAVtrial <- sum(pmax(htrial, 0))
      if (TAVtrial > mu) {
        if (TAVtrial <= TAVpop[i]) {
          pop[, i] <- trial
          hpop[, i] <- htrial
          F[, i] <- Ftrial
          CR[i] <- CRtrial
          pF[i] <- pFtrial
          TAVpop[i] <- TAVtrial
        }
      } else if (TAVpop[i] > mu) {
        pop[, i] <- trial
        fpop[i] <- fn1(trial)
        hpop[, i] <- htrial
        F[, i] <- Ftrial
        CR[i] <- CRtrial
        pF[i] <- pFtrial
        TAVpop[i] <- TAVtrial
      } else {
        ftrial <- fn1(trial)
        if (ftrial <= fpop[i]) {
          pop[, i] <- trial
          fpop[i] <- ftrial
          hpop[, i] <- htrial
          F[, i] <- Ftrial
          CR[i] <- CRtrial
          pF[i] <- pFtrial
          TAVpop[i] <- TAVtrial
          FF <- sum(TAVpop <= mu)/NP
          mu <- mu * (1 - FF/NP)
        }
      }
    })
  }
  else {
    expression({
      htrial <- constr1(trial)
      TAVtrial <- sum(pmax(htrial, 0))
      if (TAVtrial > mu) {
        if (TAVtrial <= TAVpop[i]) {
          pop[, i] <- trial
          hpop[, i] <- htrial
          F[, i] <- Ftrial
          CR[i] <- CRtrial
          pF[i] <- pFtrial
          TAVpop[i] <- TAVtrial
        }
      } else if (TAVpop[i] > mu) {
        pop[, i] <- trial
        fpop[i] <- fn1(trial)
        hpop[, i] <- htrial
        F[, i] <- Ftrial
        CR[i] <- CRtrial
        pF[i] <- pFtrial
        TAVpop[i] <- TAVtrial
        FF <- sum(TAVpop <= mu)/NP
        mu <- mu * (1 - FF/NP)
      } else {
        ftrial <- fn1(trial)
        if (ftrial <= fpop[i]) {
          pop[, i] <- trial
          fpop[i] <- ftrial
          hpop[, i] <- htrial
          F[, i] <- Ftrial
          CR[i] <- CRtrial
          pF[i] <- pFtrial
          TAVpop[i] <- TAVtrial
          FF <- sum(TAVpop <= mu)/NP
          mu <- mu * (1 - FF/NP)
        }
      }
    })
  }
  fn1 <- function(par) fn(par, ...)
  if (!is.null(constr)) 
    constr1 <- if (meq > 0) {
      eqI <- 1:meq
      function(par) {
        h <- constr(par, ...)
        h[eqI] <- abs(h[eqI]) - eps
        h
      }
    }
  else function(par) constr(par, ...)
  use.jitter <- !is.null(jitter_factor)
  conv <- expression((do.call(compare_to, list(fpop)) - fpop[x.best.ind])/fnscale)
  pop <- matrix(runif(NP * d, lower, upper), nrow = d)
  if (!is.null(add_to_init_pop)) {
    pop <- unname(cbind(pop, add_to_init_pop))
    NP <- ncol(pop)
  }
  stopifnot(NP >= 4)
  F <- if (use.jitter) 
    (1 + jitter_factor * runif(d, -0.5, 0.5)) %o% runif(NP, 
                                                        Fl, Fu)
  else matrix(runif(NP, Fl, Fu), nrow = 1)
  CR <- runif(NP)
  pF <- runif(NP)
 
  fpop <- apply(pop, 2, fn1)
  #print(which(fpop>=0.9))
  if (!is.null(constr)) {
    hpop <- apply(pop, 2, constr1)
    if (any(is.na(hpop))) 
      stop("value of meq is invalid")
    if (is.vector(hpop)) 
      dim(hpop) <- c(1, length(hpop))
    TAVpop <- apply(hpop, 2, function(x) sum(pmax(x, 0)))
    mu <- median(TAVpop)
  }
  popIndex <- 1:NP
  x.best.ind <- which.best(fpop)
  converge <- eval(conv)
  rule <- if (!is.null(constr)) 
    expression(converge >= tol || any(hpop[, x.best.ind] > 0))
  else {expression(converge >= tol)}

  convergence <- 0
  iteration <- 0
  while (eval(rule)) {
    if (iteration >= maxiter) {
      warning("maximum number of iterations reached without convergence")
      convergence <- 1
      break
    }
    iteration <- iteration + 1
    for (i in popIndex) {
      i <- ((iteration + i)%%NP) + 1
      Ftrial <- if (runif(1) <= tau_F) {
        if (use.jitter) 
          runif(1, Fl, Fu) * (1 + jitter_factor * runif(d, 
                                                        -0.5, 0.5))
        else runif(1, Fl, Fu)
      }
      else F[, i]
      CRtrial <- if (runif(1) <= tau_CR) 
        runif(1)
      else CR[i]
      pFtrial <- if (runif(1) <= tau_pF) 
        runif(1)
      else pF[i]
      X.i <- pop[, i]
      r <- sample(popIndex[-i], 3)
      X.base <- pop[, r[1L]]
      X.r1 <- pop[, r[2L]]
      X.r2 <- pop[, r[3L]]
      trial <- handle.bounds(performReproduction(), X.base)
      eval(child)
      x.best.ind <- which.best(fpop)
    }
    converge <- eval(conv)
    if (trace && (iteration%%triter == 0)) 
      cat(iteration, ":", "<", converge, ">", "(", fpop[x.best.ind], 
          ")", pop[, x.best.ind], if (!is.null(constr)) 
            paste("{", which(hpop[, x.best.ind] > 0), 
                  "}"), fill = TRUE)
  }
  res <- list(
    par = pop[, x.best.ind], 
    value = fpop[x.best.ind], 
              iter = iteration
    #, convergence = convergence
    )
  if (details) {
    res$poppar <- pop
    res$popcost <- fpop
  }
  res
}
