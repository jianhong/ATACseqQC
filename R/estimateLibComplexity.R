#' Library complexity estimation
#' @description Estimating the library complexity.
#' @param histFile A two-column matrix of integers. The 1st column is the frequency 
#' j = 1,2,3,.... The 2nd column is the number of genomic regions with the same 
#' fequency (j) of duplication. This file should be sorted by the first column 
#' in ascending order. For example, one row of a histogram file:
#' 10    20
#' means there are 10 genomic regions, each of which is covered by 20 identical 
#' fragments at a given sequencing depth of a sequencing library.
#' @param times An positive integer representing the minimum required number of 
#' successful estimation. Default is 100.
#' @param interpolate.sample.sizes A numeric vector with values between (0, 1].
#' @param extrapolate.sample.sizes A numeric vector with values greater than 1.
#' @export
#' @importFrom polynom polynomial polylist
#' @author Haibo Liu, Feng Yan, Chao Deng
#' @return invisible estimates, a data frame of 3 columns: relative sequence depth, 
#' number of distinct fragments, number of putative sequenced reads.
#' @seealso \link{readsDupFreq}
#' @examples
#' FisherButterfly <- readRDS(system.file('extdata', 'FisherButterfly.rds',
#'                            package='ATACseqQC'))
#' estimateLibComplexity(histFile=FisherButterfly, times=100)

estimateLibComplexity <- function(histFile, times=100, 
                                  interpolate.sample.sizes=seq(0.1, 1, by=0.1),
                                  extrapolate.sample.sizes=seq(5, 20, by=5)){
    total <- histFile[,1] %*% histFile[,2]  ## added
    suppressWarnings({result = ds.rSAC.bootstrap(histFile, r=1, times=times)})
    sequences <- c(interpolate.sample.sizes, extrapolate.sample.sizes)
    estimates <- data.frame(relative.size=sequences, values=rep(NA, length(sequences)))
    for ( i in seq_along(sequences))
    {
        estimates$values[i] <- result$f(sequences[i])
    }
    suppressWarnings(estimates$reads <- estimates$relative.size * total) ## added
    plot(x=estimates$reads/10^6, y=estimates$values/10^6,  
         type="o", xlab =expression(Putative~sequenced~fragments~x~10^6), 
         ylab=expression(Distinct~fragments~x~10^6),
         main="Estimation of ATAC-seq\nlibrary complexity")
    return(invisible(estimates))  
}

## copied from preseqR::rSAC.R, since preseqR is archived.
# Copyright (C) 2016 University of Southern California and
#          Chao Deng and Andrew D. Smith and Timothy Daley
#
# Authors: Chao Deng
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>

## continued fraction approximant to a power series based on
## QD algorithm
## coefs, coefficients of the power series;
## mt, the number of terms in the power series used for constructing
## the continued fraction approximation
## ref pp. 131, 147 and 148 in the book Pad\'{e} Approximants 2ed
ps2cfa <- function(coefs, mt) {
  ## use nonzero terms which are required by QD algorithm
  index <- which(coefs == 0)
  if (length(index) == 0) {
    mt <- min(mt, length(coefs))
  } else {
    mt <- min(mt, index[1] - 1)
  }
  if (mt == 1) {
    return(coefs[1])
  }
  ## QD algorithm
  qd.table <- matrix(data=0, nrow=mt, ncol=mt)
  ## initialize the table
  ## the first column is 0
  qd.table[1:(mt-1), 2] <- coefs[2:mt] / coefs[1:(mt-1)]
  if (mt == 2) {
    return(c(coefs[1], -qd.table[1, 2]))
  }
  ## two types of columns e or q
  for (i in 3:mt) {
    n <- mt + 1 - i
    if (i %% 2 == 1) {
      ## number of entries in the column
      qd.table[1:n, i] <- qd.table[2:(n+1), i-1] - qd.table[1:n, i-1] +
        qd.table[2:(n+1), i-2]
      if (!is.finite(qd.table[1, i]) || qd.table[1, i] == 0)
        return(c(coefs[1], -qd.table[1, 2:(i-1)]))
    } else {
      qd.table[1:n, i] <- qd.table[2:(n+1), i-1] / qd.table[1:n, i-1] *
        qd.table[2:(n+1), i-2]
      if (!is.finite(qd.table[1, i]) || qd.table[1, i] == 0)
        return(c(coefs[1], -qd.table[1, 2:(i-1)]))
    }
  }
  return(c(coefs[1], -qd.table[1, 2:mt]))
}

## convert truncated continued fraction to a series of rational functions
## numerators are stored in set A and denumerators are stored in set B
## equation (2.14a), (2.14b), (2.15) in the book Pad\'{e} Approximants 2ed
cfa2rf <- function(CF) {
  ## A, B are sets of polynomials based on recursive formula
  A <- list()
  B <- list()
  if (length(CF) < 2) {
    return(polynomial(CF))
  }
  A[[1]] <- polynomial(CF[[1]])
  A[[2]] <- polynomial(CF[[1]])
  B[[1]] <- polynomial(1)
  B[[2]] <- polynomial(c(1, CF[[2]]))
  if (length(CF) == 2) {
    return(list(A=A, B=B))
  }
  for (i in 3:length(CF)) {
    A[[i]] <- A[[i-1]] + polynomial(c(0, CF[[i]])) * A[[i-2]]
    B[[i]] <- B[[i-1]] + polynomial(c(0, CF[[i]])) * B[[i-2]]
  }
  return(list(A=A, B=B))
}

## Pad\'{e} approximant by picking out the numerator and the denominator
## input: two sets of polynomials for numerators and denominators
##        the degree m
## output: Pad\'{e} approximant
rf2rfa <- function(RF, m) {
  return(polylist(RF$A[[m]], RF$B[[m]]))
}

## simplify the rational function, eliminate defects and partial-fraction
## decompoistion
rfa.simplify <- function(rfa) {
  ## solving roots
  numer.roots <- solve(rfa[[1]])
  denom.roots <- solve(rfa[[2]])
  
  ## finite
  if (any(!is.finite(c(numer.roots, denom.roots))))
    return(NULL)
  
  ## identify defects
  ## the root and the pole is a defect if the difference is less than
  ## the predefined precision, which is defined by the variable PRECISION
  PRECISION <- 1e-3
  tmp.roots <- c()
  for (i in 1:length(denom.roots)) {
    if (length(numer.roots) > 0) {
      d <- Mod(denom.roots[i] - numer.roots)
      ind <- which.min(d)
      if (d[ind] < PRECISION) {
        numer.roots <- numer.roots[-ind]
        tmp.roots <- c(tmp.roots, denom.roots[i])
      }
    }
  }
  
  ## eliminate defects
  denom.roots <- denom.roots[!denom.roots %in% tmp.roots]
  ## convert roots from t - 1 to t
  poles <- denom.roots + 1
  
  ## treat both numerator and denuminator in the rational function as
  ## monic polynomials
  ## the difference from the original rational function is up to a factor
  if (length(numer.roots) == 0) {
    poly.numer <- as.function(polynomial(1))
  } else {
    ## construct polynomials using all the roots
    p <- 1
    for (x in numer.roots) {
      p <- c(0, p) - c(x * p, 0)
    }
    ## in theory coefficients p of the polynomial should be real numbers
    ## Re(p) == p
    poly.numer <- as.function(polynomial(Re(p)))
  }
  l <- length(denom.roots)
  
  ## coefficients in the partial fraction
  coefs <- sapply(1:l, function(x) {
    poly.numer(denom.roots[x]) / prod(denom.roots[x] - denom.roots[-x])})
  ## calculate the factor
  C <- coef(rfa[[1]])[length(coef(rfa[[1]]))] /
    coef(rfa[[2]])[length(coef(rfa[[2]]))]
  coefs <- coefs * C
  return(list(coefs=coefs, poles=poles))
}

### check the input histogram in an appropriate format
checking.hist <- function(n)
{
  if (ncol(n)!=2 || is.numeric(n[,1])==FALSE || is.numeric(n[,2])==FALSE) {
    stop("Input must be a two-column matrix")
  }
  ## the first column is the frequency i
  ## the second column is the number of species represented i times
  ## in the sample
  freq <- n[, 1]
  num <- n[, 2]
  
  ## check whether frequencies are at least one and the histogram is sorted
  ## based on frequencies
  for (i in 1:length(freq))
    if (freq[i] <= 0 || freq[i] != floor(freq[i])) {
      stop("The first column must be positive integers!")
    } else if (num[i] < 0) {
      stop("The second column must be non negative")
    }
  else {
    if (i > 1 && freq[i - 1] >= freq[i])
      stop("The first column is not sorted in the ascending order")
  }
  
  return(n)
}

## coefficients for the power series of E(S_1(t)) / t
## return the first mt terms
discoveryrate.ps <- function(n, mt)
{
  ## transform a histogram into a vector of frequencies
  hist.count <- vector(length=max(n[, 1]), mode="numeric")
  hist.count[n[, 1]] <- n[, 2]
  
  PS.coeffs <- sum(hist.count)
  if (mt == 1) {
    return(PS.coeffs)
  }
  
  change.sign <- 0
  for (i in 1:(min(mt-1, length(hist.count)))) {
    PS.coeffs <- c(
      PS.coeffs,
      (-1)^change.sign * hist.count[i] - PS.coeffs[length(PS.coeffs)])
    change.sign <- change.sign + 1
  }
  
  ## truncate at coefficients where it is zero
  zero.index <- which(PS.coeffs == 0)
  if (length(zero.index) > 0) {
    PS.coeffs[1:(min(zero.index) - 1)]
  } else {
    PS.coeffs
  }
}


ds.rSAC <- function(n, r=1, mt=20)
{
  checking.hist(n)
  
  n[, 2] <- as.numeric(n[, 2])
  
  ## coefficients of average discovery rate for the first mt terms
  PS.coeffs <- discoveryrate.ps(n, mt=mt)
  
  if (is.null(PS.coeffs)) {
    write("the size of the initial experiment is insufficient", stderr())
    return(NULL)
  }
  
  ## use nonzero coefficients
  mt <- min(mt, length(PS.coeffs))
  PS.coeffs <- PS.coeffs[ 1:mt ]
  
  ## check whether sample size is sufficient
  if (mt < 2)
  {
    m <- paste("max count before zero is less than min required count (2)",
               " sample not sufficiently deep or duplicates removed", sep = ',')
    write(m, stderr())
    return(NULL)
  }
  
  ## construct the continued fraction approximation to the power seies
  cf <- ps2cfa(coefs=PS.coeffs, mt=mt)
  rf <- cfa2rf(CF=cf)
  ## the length of cf could be less than mt
  ## even if ps do not have zero terms, coefficients of cf may have
  mt <- length(cf)
  mt <- mt - (mt %% 2)
  valid.estimator <- FALSE
  m <- mt
  while (valid.estimator == FALSE && m >= 2) {
    
    ## rational function approximants [m / 2 - 1,  m / 2]
    rfa <- rf2rfa(RF=rf, m=m)
    rfa <- rfa.simplify(rfa)
    if (is.null(rfa)) {
      m <- m - 2
      next
    }
    ## check stability
    if (any(Re(rfa$poles) >= 0)) {
      m <- m - 2
      next
    }
    
    coefs <- rfa$coefs
    poles <- rfa$poles
    ## check whether the estimator is non-decreased
    ## NOTE: it only checks for t >= 1 !!!
    deriv.f <- function(t) {
      Re(sapply(t, function(x) {-(coefs*poles) %*% ( 1 / ((x-poles)^2))}))}
    if (any( deriv.f(seq(1, 100, by=0.05)) < 0 )) {
      m <- m - 2
      next
    } else {
      f.rSAC <- function(t) {
        sapply(t, function(x) {
          Re(coefs %*% (x / (x - poles))^r)})}
      valid.estimator <- TRUE
    }
  }
  
  if (valid.estimator == TRUE) {
    return(f.rSAC)
  } else {
    ## the case S1 = S2 where the numbe of species represented exactly once
    ## is 0
    return(function(t) {sapply(t, function(x) return(sum(n[, 2])))})
  }
}

## the bootstrap version of ds.rSAC
## with confidence interval
#' @importFrom stats rmultinom qnorm coef 
ds.rSAC.bootstrap <- function(n, r=1, mt=20, times=30, conf=0.95)
{
  n[, 2] <- as.numeric(n[, 2])
  ## individuals in the sample
  N <- n[, 1] %*% n[, 2]
  
  ## returned function
  f.rSACs <- vector(length=times, mode="list")
  
  f.bootstrap <- function(n, r, mt) {
    n.bootstrap <- matrix(c(n[, 1], rmultinom(1, sum(n[, 2]), n[, 2])), ncol=2)
    N.bootstrap <- n.bootstrap[, 1] %*% n.bootstrap[, 2]
    N <- n[, 1] %*% n[, 2]
    t.scale <- N / N.bootstrap
    f <- ds.rSAC(n.bootstrap, r=r, mt=mt)
    return(function(t) {f(t * as.vector(t.scale))})
  }
  
  while (times > 0) {
    f.rSACs[[times]] <- f.bootstrap(n=n, r=r, mt=mt)
    ## prevent later binding!!!
    f.rSACs[[times]](1)
    times <- times - 1
  }
  
  estimator <- function(t) {
    result <- sapply(f.rSACs, function(f) f(t))
    if (length(t) == 1) {
      return(median(result))
    } else {
      return(apply(result, FUN=median, MARGIN=1))
    }
  }
  
  variance <- function(t) {
    result <- sapply(f.rSACs, function(f) f(t))
    if (length(t) == 1) {
      return(var(result))
    } else {
      return(apply(result, FUN=var, MARGIN=1))
    }
  }
  
  se <- function(x) sqrt(variance(x))
  
  ## prevent later binding!!!
  estimator(1); estimator(1:2)
  variance(1); variance(1:2)
  ## confidence interval using lognormal
  q <- (1 + conf) / 2
  lb <- function(t) {
    C <- exp(qnorm(q) * sqrt(log( 1 + variance(t) / (estimator(t)^2) )))
    return(estimator(t) / C)
  }
  ub <- function(t) {
    C <- exp(qnorm(q) * sqrt(log( 1 + variance(t) / (estimator(t)^2) )))
    return(estimator(t) * C)
  }
  return(list(f=estimator, se=se, lb=lb, ub=ub))
}

