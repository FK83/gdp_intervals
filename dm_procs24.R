# Implementation of the EWC estimator proposed by Lazarus et al (JBES, 2018, see e.g. Table 1, line 3)
# Loosely based on Matlab code by Lazarus et al, available at https://ebenlazarus.github.io/#research

# Helper functions for ewc variance estimator 
cosine_weight <- function(T, q){
  # T is length of time series, q is nr of cosine terms
  fvec <- pi*(1:q)
  tvec <- seq(from = .5, by = 1, length.out = T)/T
  psi <- sqrt(2/T)*cos(outer(tvec, fvec))
  # Returns T*q matrix of weights
  psi
}
cutoff <- function(T){
  # See Equation (4) and Footnote (3) in Lazarus et al (2018)
  floor(.4*(T^(2/3)))
}

# input: fit, linear model object (output of lm function)
# output: list with three objects 
# (covariance matrix of parameters, cutoff value for lags, vector of t statistics for each parameter in fit)
ewc <- function(fit, nu = NULL){
  # Get regressor matrix and dimensions
  x <- model.matrix(fit)
  T <- nrow(x)
  k <- ncol(x)
  # Get residuals, construct x*u
  u <- residuals(fit)
  z <- x*u
  # Truncation parameter (called either nu or B in paper?)
  if (is.null(nu)){
    nu <- cutoff(T)
  }
  # Omega matrix (Equation 10 in paper)
  w <- cosine_weight(T, nu)
  Omega <- (t(z) %*% w %*% t(w) %*% z)/nu
  # Get inverse of Q matrix (see e.g. Hansen, Ch. 14.33)
  Q_inv <- summary(fit)$cov.unscaled*T
  # Multiply by Q^{-1} (both sides)
  vcv <- Q_inv %*% Omega %*% Q_inv
  # t ratios
  t <- sqrt(T)*coefficients(fit)/sqrt(diag(vcv))
  # p-values (two-sided)
  p <- 2*pt(-abs(t), df = nu)
  list(Omega = Omega, 
       vcv = vcv, nu = nu, 
       t = t, 
       p = p)
}

# Benchmark implementation of DM for Online Supplement
# (called `CM13' there)
# Inputs: d, series of loss differences
# h, forecast horizon
# fsa, whether or not to compute finite sample adjustment of variance
# rect.kern, whether or not to use rectangular kernel
dmtest <- function(d,h, fsa=TRUE, rect.kern = TRUE){
  if (all(d == 0)){ list(stat=NA,p1=NA,p2=NA) } else {
  
    T <- length(d)
    
    if (rect.kern == TRUE){
      v <- mean( (d-mean(d))^2 )
      ind <- 0
      if (h > 1){
        for (hh in 1:(h-1)){
          v <- v + 2*sum( (d-mean(d))[(hh+1):T] * (d-mean(d))[1:(T-hh)] )/T
        }
      }
      
      if (v <= 0){
        rect.kern <- FALSE # Switch type in case type 1 produced negative variance
        warning("Rectangular kernel gave negative variance - other kernel used")
      } 
    }
    
    if (rect.kern == FALSE){
      aux <- data.frame(d=d)
      aux <- lm(d~1,data=aux)
      v <- unname(NeweyWest(aux)*T)
    }
    t <- mean(d)/sqrt(v/T)
    if (fsa == TRUE) t <- t*sqrt( (T+1-2*h+h*(h-1)/T)/T )
    list(stat=t,p1=2*pnorm(-abs(t)),p2 = pnorm(t))  
  }
}