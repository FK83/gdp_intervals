# Define helper function for printing
print_helper <- function(x, nd = 2){
  format(round(x, nd), nsmall = nd)
}


# postprocessing via Gaussian heteroskedastic regression
gaussreg_pp <- function(y_training, h_training, tau_gauss = c(.1, .9),
                        h_test = h_training, use_crps = TRUE,
                        theta0 = NULL, zero_mean = FALSE){
  if (is.null(theta0)){
    theta0 <- c(rep(0, 2-zero_mean), 40, 0)
  }
  if (use_crps){
    scoring_rule <- crps_norm
  } else {
    scoring_rule <- logs_norm
  }
  opt <- optim(theta0, ff, y = y_training, x = h_training,
               scoring_rule = scoring_rule, zero_mean = zero_mean)
  theta <- opt$par
  # parameter location depends on whether mean is zero
  if (zero_mean){
    inds <- 1:3
    m <- 0
  } else {
    inds <- 2:4
    m <- theta[1]
  }
  s_pred <- exp(theta[inds[1]])*plogis(h_test,
                                       location = theta[inds[2]],
                                       scale = exp(theta[inds[3]]))
  pred <- matrix(NA, length(h_test), length(tau_gauss))
  for (jj in 1:length(tau_gauss)){
    pred[,jj] <- m + qnorm(tau_gauss[jj])*s_pred
  }
  return(list(opt = opt, pred = pred))
}

# objective function
ff <- function(theta, y, x, scoring_rule, zero_mean){
  if (zero_mean){
    inds <- 1:3
    a <- 0
  } else {
    inds <- 2:4
    a <- theta[1]
  }
  b <- exp(theta[inds[1]])
  c <- theta[inds[2]]
  d <- exp(theta[inds[3]])
  s <- b*plogis(x, location = c, scale = d)
  mean(scoring_rule(y, mean = rep(a, length(y)), sd = s))
}

# Postprocessing via isotonic distributional regression
# (assumption about first-order dominance of absolute forecast errors)
isodistrreg_pp <- function(y_training, h_training,
                           tau_idr = c(.1, .9), h_test = h_training,
                           simple = FALSE){
  if (length(tau_idr) != 2){
    stop("Implemented for two quantile levels only")
  }
  if (simple){
    # Assume that distribution of e is symmetric about zero
    train_df <- data.frame(h = h_training)
    test_df <- data.frame(h = h_test)
    fit <- idr(y = abs(y_training), X = train_df)
    predictions <- predict(fit, data = test_df)
    level_l <- 1 - 2*tau_idr[1]
    level_u <- 2*(tau_idr[2] - .5)
    pred <- matrix(NA, length(h_test), 2)
    pred[,1] <- -qpred(predictions, quantiles = level_l)
    pred[,2] <- qpred(predictions, quantiles = level_u)
    list(fit = fit, pred = pred)
  } else {
    # Separate analysis for negative and positive errors
    sel_neg <- y_training <= 0
    sel_pos <- y_training > 0
    p_neg <- mean(sel_neg)
    p_pos <- mean(sel_pos)
    fit_neg <- idr(y = abs(y_training[sel_neg]),
                   X = data.frame(h = h_training[sel_neg]))
    fit_pos <- idr(y = abs(y_training[sel_pos]),
                   X = data.frame(h = h_training[sel_pos]))
    predictions_neg <- predict(fit_neg,
                               data = data.frame(h = h_test))
    predictions_pos <- predict(fit_pos,
                               data = data.frame(h = h_test))
    level_l <- 1 - tau_idr[1]/p_neg
    level_u <- (tau_idr[2] - p_neg)/p_pos
    pred <- matrix(NA, length(h_test), 2)
    pred[,1] <- -qpred(predictions_neg, quantiles = level_l)
    pred[,2] <- qpred(predictions_pos, quantiles = level_u)
    list(fit_neg = fit_neg, fit_pos = fit_pos, pred = pred)
  }
}

# Function to get RMSE values from structural parameters (see function ar1_model)
get_rmse <- function(hh, P_t_t, F_power_list, F_prime_power_list, 
                     Q, gamma0){
  if (hh == 0){
    tmp1 <- P_t_t
  } else {
    tmp1 <- F_power_list[[hh]] %*% P_t_t %*% F_prime_power_list[[hh]]  
  }
  tmp2 <- Q
  if (hh > 1){
    for (jj in 1:(hh-1)){
      tmp2 <- tmp2 + F_power_list[[jj]] %*% Q %*% 
        F_prime_power_list[[jj]]
    }    
  }
  tmp3 <- tmp1 + tmp2
  sqrt(t(gamma0) %*% tmp3 %*% gamma0)
}

# function to simulate from ar(1) model (outcomes + forecasts)
# forecasts are based on imperfect observations (measurement error)
ar1_model <- function(freq = "w", type_agg = "annual_average", 
                      phi = .3, s2e = .09, s2w = .003, 
                      n_years = 30, n_mc = 1e3){
  # collect inputs (returned for reference below)
  inputs <- as.list(environment())
  # nr of time units per year
  n_units <- ifelse(freq == "q", 4, 
                    ifelse(freq == "m", 12, 52))
  # determine dimension of state vector
  n <- ifelse(type_agg == "annual_average", 2*n_units, n_units)
  # construct vector of aggregation weights
  if (type_agg == "annual_average"){
    gamma0 <- matrix(c(1:n_units, (n_units-1):0)/n_units, ncol = 1)
  } else {
    gamma0 <- matrix(1, n_units, 1)
  }
  # nr of time periods
  n_periods <- n_years*n_units
  
  # Draw realization data
  y <- matrix(0, n_periods, n_mc)
  for (rr in 1:n_mc){
    y[, rr] <- arima.sim(list(ar = phi), n = n_periods,
                         sd = sqrt(s2e))
  }
  y_tilde <- y + matrix(rnorm(n_periods*n_mc, sd = sqrt(s2w)), 
                        n_periods, n_mc)
  
  # Construct Kalman filter matrices
  F <- rbind(c(phi, rep(0, n-1)), cbind(diag(n-1), 0))
  F_prime <- t(F)
  Q <- matrix(0, n, n)
  Q[1,1] <- s2e
  #H <- matrix(c(1, rep(0, n-1)), n, 1) // not actually needed below (H = 1 is hardcoded)
  R <- s2w
  
  # Make list of matrices of matrix powers (needed below)
  # see https://stackoverflow.com/questions/14284676/initialize-a-list-of-matrices-in-r
  F_power_list <- F_prime_power_list <- 
    rep(list(matrix(0, n, n)), n)
  for (hh in 1:(2*n_units)){
    F_power_list[[hh]] <- bvarsv:::matmult(F, hh)
    F_prime_power_list[[hh]] <- bvarsv:::matmult(F_prime, hh)  
  }
  
  # Initialize matrices for MC results
  rlz_all <- fcst_all <- matrix(NA, n_periods*3, n_mc + 2)
  s_all <- matrix(NA, n_periods*3, 3)
  
  # initial mean and VCV of state forecast
  eta_tp1_t <- matrix(0, n, n_mc) # separate for each MC iteration
  P_tp1_t <- s2e/(1-phi^2) * (phi^as.matrix(dist(1:n))) 
  
  # Auxiliary matrices (to be modified during loop)
  est <- matrix(NA, n_periods, n_mc)
  
  # iterations over time
  ct_mat_rows <- 0
  for (tt in 1:n_periods){
    # print progress
    if ((tt %% 100) == 0){
      print(paste("Now running iteration", tt))
    }
    # loop over MC iterations
    eta_t_t <- matrix(NA, n, n_mc)
    for (rr in 1:n_mc){
      # Compute state nowcast
      eta_t_t[, rr] <- eta_tp1_t[, rr] + 
        (P_tp1_t[,1,drop = FALSE]/(P_tp1_t[1,1] + R))*
        (y_tilde[tt, rr]-eta_tp1_t[1, rr])
      # Save nowcast of most recent obs for checking
      est[tt,rr] <- eta_t_t[1,rr]  
      # Make forecast for next period's state
      eta_tp1_t[, rr] <- F %*% eta_t_t[,rr]
    }
    # Compute updates of Kalman filter matrices
    # (these are the same for each MC iteration)
    P_t_t <- P_tp1_t - (P_tp1_t[,1,drop = FALSE] %*% 
                          P_tp1_t[1,,drop = FALSE])/(P_tp1_t[1,1] + R)
    P_tp1_t <- F %*% P_t_t %*% F_prime + Q
    # Make forecasts
    # nr of years elapsed
    years_elapsed <- floor(tt/n_units)
    # nr of time units elapsed
    units_elapsed <- tt %% n_units
    # determine horizons
    hors <- c(0, n_units, 2*n_units) - units_elapsed
    hors <- hors[hors >= 0]
    # determine target years
    target_years <- (years_elapsed*n_units + units_elapsed + hors)/n_units
    if (any(target_years != as.integer(target_years))){
      stop("Smth wrong with forecast date calculations")
    }
    # vector to store forecast std. dev. (same for all MC iterations)
    s <- rep(NA, length(hors))
    # matrix to store forecasts, outcomes and errors
    fcst <- rlz <- matrix(NA, length(hors), n_mc)
    for (hh in hors){
      hh_ind <- which(hors == hh)
      s[hh_ind] <- get_rmse(hh, P_t_t, F_power_list,
                            F_prime_power_list, Q, gamma0)
      # check whether outcome data is available
      ind1 <- tt+hh-n+1
      ind2 <- tt+hh
      available <- (ind1 >= 1) & (ind2 <= n_periods)
      for (rr in 1:n_mc){
        # Make forecast
        if (hh == 0){
          fcst_tmp <- eta_t_t[,rr]  
        } else {
          fcst_tmp <- F_power_list[[hh]] %*% eta_t_t[,rr]
        }
        fcst[hh_ind,rr] <- sum(gamma0 * fcst_tmp)
        # Make outcome
        if (available){
          rlz[hh_ind,rr] <- sum(y[(tt+hh-n+1):(tt+hh),rr]*rev(gamma0))
        } 
      }
    }
    # Append results to matrix
    rlz_all[(ct_mat_rows+1):(ct_mat_rows+nrow(rlz)), ] <- 
      cbind(target_years, hors, rlz)
    fcst_all[(ct_mat_rows+1):(ct_mat_rows+nrow(fcst)), ] <- 
      cbind(target_years, hors, fcst)
    s_all[(ct_mat_rows+1):(ct_mat_rows+length(s)),] <- 
      cbind(target_years, hors, s)
    ct_mat_rows <- ct_mat_rows + nrow(rlz)
  }
  
  # Remove superfluous rows
  rlz_all <- rlz_all[1:ct_mat_rows,]
  fcst_all <- fcst_all[1:ct_mat_rows,]
  s_all <- s_all[1:ct_mat_rows,]
  
  # most recent (~steady state) version of P_tp1_t
  P1 <- P_tp1_t
  
  # most recent version of P_t_t
  P0 <- P1 - (P1[,1,drop = FALSE]) %*% P1[1,,drop = FALSE] / (P1[1,1] + R)
  
  # Compute rmse for each horizon 
  rmse <- rep(NA, 2*n_units)
  for (hh in 1:(2*n_units)){
    rmse[hh] <- get_rmse(hh, P0, F_power_list, F_prime_power_list, 
                         Q, gamma0)
  }
  
  # Get forecast errors
  e_all <- cbind(rlz_all[,1:2], rlz_all[,-(1:2)]-fcst_all[,-(1:2)])
  
  # Return list of outputs
  list(inputs = inputs, y = y, y_tilde = y_tilde, 
       rlz = rlz_all, fcst = fcst_all, e = e_all, s = s_all, 
       rmse = rmse)
  
}

#' Quantile and interval scores
#' @param y vector of observations
#' @param x vector of quantile predictions
#' @param x_lower,x_upper vector of quantile predictions (lower and upper endpoints of prediction intervals)
#' @param alpha quantile level of interest 
#' @param target_coverage target (i.e., nominal) coverage level of prediction interval
#' @return A vector of score values. Smaller values indicate better forecasts. Note that 
#' the interval score refers to the central prediction interval at level \code{target_coverage}.
#' @references 
#' Quantile score:
#' Koenker, R. and G. Bassett (1978): `Regression quantiles', Econometrica 46, 33-50. \doi{https://doi.org/10.2307/1913643}
#' Interval score:
#' Gneiting, T. and A.E. Raftery (2007):
#' `Strictly proper scoring rules, prediction and estimation',
#' Journal of the American Statistical Association 102, 359-378. \doi{10.1198/016214506000001437}
#' @examples
#' # Interval score is proportional to sum of two quantile scores
#' 10*(qs(y = 1, x = qnorm(.1), alpha = .1) + qs(y = 1, x = qnorm(.9), alpha = .9))
#' ints(y = 1, x_lower = qnorm(.1), x_upper = qnorm(.9), target_coverage = .8)
#' @name scores_quantiles
#' @export
qs <- function(y, x, alpha){
  if (!is.numeric(alpha) || length(alpha) != 1 ){
    stop("alpha must be numeric of length 1")
  } else if (alpha <= 0 || alpha >= 1) {
    stop("alpha must satisfy 0 < alpha < 1")
  }
  # score calculation
  ((y < x)-alpha)*(x-y)
}

#' interval score
#' @export
#' @rdname scores_quantiles
ints <- function(y, x_lower, x_upper, target_coverage){
  if (!is.numeric(target_coverage) || length(target_coverage) != 1 ){
    stop("target_coverage must be numeric of length 1")
  } else if (target_coverage <= 0 || target_coverage >= 1) {
    stop("target_coverage must satisfy 0 < target_coverage < 1")
  }
  if (any(x_lower > x_upper)){
    stop("'x_lower' contains values greater than corresponding values in 'x_upper'.")  
  }
  # get interval score as sum of two quantile scores
  aux1 <- qs(y, x_lower, .5*(1-target_coverage)) 
  aux2 <- qs(y, x_upper, .5*(1+target_coverage)) 
  (2/(1-target_coverage))*(aux1 + aux2)
}

# SPF upper bounds for real gdp (from P41 of SPF documentation pdf file)
upper_bounds_gdp <- function(origin_date){
  if (origin_date >= 1981.5 & origin_date <= 1991.75){
    ub <- c(-2, c(0, 2, 4, 6) - .05, Inf)
    n_bins <- 6
  } else if (origin_date >= 1992 & origin_date <= 2009){
    ub <- c(-2, seq(from = -1, to = 6, by = 1) - .05, Inf)
    n_bins <- 10
  } else if (origin_date >= 2009.25 & origin_date <= 2020){
    ub <- c(-3, seq(from = -2, to = 6, by = 1) - .05, Inf)
    n_bins <- 11
  } else if (origin_date > 2020){
    ub <- c(-12, c(-6, -3, 0, 1.5, 2.5, 4, 7, 10, 16) - .05, Inf)
    n_bins <- 11
  }
  list(ub = ub, n_bins = n_bins)
}

# SPF upper bounds for inflation
upper_bounds_inflation <- function(origin_date){
  if (origin_date >= 1981.5 & origin_date <= 1985){
    ub <- c(4, c(6, 8, 10, 12) - .05, Inf)
    n_bins <- 6
  } else if (origin_date >= 1985.25 & origin_date <= 1991.75){
    ub <- c(2, c(4, 6, 8, 10) - .05, Inf)
    n_bins <- 6
  } else if (origin_date >= 1992 & origin_date <= 2013.75){
    ub <- c(0, (1:8) - .05, Inf)
    n_bins <- 10
  } else if (origin_date >= 2014){
    ub <- c(0, seq(from = .5, to = 4, by = .5) - .05, Inf)
    n_bins <- 10
  }
  list(ub = ub, n_bins = n_bins)
}

# SPF upper bounds for either variable
upper_bounds <- function(origin_date, var){
  if (var == "gdp"){
    out <- upper_bounds_gdp(origin_date)
  } else if (var == "inf"){
    out <- upper_bounds_inflation(origin_date)
  }
  if (length(out$ub) != out$n_bins) stop("Wrong info on nr of bins")
  out
}

# helper function to rescale probabilities
rescale_p <- function(p){
  if (sum(p) != 1){
    sum_old <- sum(p)
    p <- p * (1/sum(p))
    print(paste("Sum (old):", sum_old))
    print(paste("Sum (new):", sum(p)))
  }
  p
}
