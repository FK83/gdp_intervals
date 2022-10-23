# predict quantiles based on a fitted crch model
pred_helper_crch <- function(fit, new_x, prob){
  pred_location <- predict(fit, 
                           newdata = new_x, 
                           type = "location")
  pred_scale <- predict(fit, 
                        newdata = new_x, 
                        type = "scale")
  if (any(pred_scale < 0)){
    stop("Negative predicted std. dev.")
  }
  pred_location_mat <- matrix(pred_location, ncol = length(prob), 
                              nrow = length(pred_scale),
                              byrow = FALSE)
  pred_scale_mat <- matrix(pred_scale, ncol = length(prob), 
                           nrow = length(pred_scale),
                           byrow = FALSE)
  out <- pred_location_mat +
    (pred_scale_mat * matrix(qnorm(prob), nrow = nrow(pred_scale_mat),
                                  ncol = length(prob), byrow = TRUE))
  out
}

# quantile evaluation (sum of quantile scores at various levels) 
quantile_eval <- function(fc, y, tau){
  fc <- as.matrix(fc)
  out <- rep(0, nrow(fc)) # return sum of quantile scores (across all levels) for each case
  for (jj in 1:ncol(fc)){
    out <- out + apl_score(x = fc[, jj], y = y, alpha = tau[jj])
  }
  out
}

# find best threshold (theta) for quantile prediction of forecast errors
# can be applied to quantile regression ("rq") or crch ("crch") model
find_best_theta <- function(df_train, quantile_levels = c(.1, .9), 
                            theta_grid = seq(from = 5, to = 20, by = .5), 
                            model = "rq"){
  loss <- rep(NA, length(theta_grid))
  for (tt in theta_grid){
    df_train <- df_train %>% mutate(h_kink = if_else(h >= tt, tt, h))
    if (model == "rq"){
      fit <- rq(e~h_kink, data = df_train, tau = quantile_levels, 
                method = "fn")
      pred <- predict(fit)
    } else if (model == "crch"){
      fit <- crch(e~1|h_kink, data = df_train, 
                  dist = "gaussian", link = "identity", type = "crps")
      pred <- pred_helper_crch(fit, new_x = df_train, prob = quantile_levels)    
    }
    # compute quantile loss for present choice of threshold (tt)
    loss[which(theta_grid == tt)] <- 
      quantile_eval(fc = pred, y = df_train$e, tau = quantile_levels) %>% mean
  }
  list(loss = loss, best_theta = theta_grid[which.min(loss)])
}

# fit quantile prediction model for forecast errors, using given fixed threshold theta
# can be applied to quantile regression ("rq") or crch ("crch") model
fit_fixed_theta <- function(df_train, df_test, theta = 12, quantile_levels = c(.1, .9), 
              		     model = "rq"){
  # define training sample data frame
  df_train <- df_train %>% mutate(h_kink = if_else(h >= theta, theta, h))
  # define test sample data frame
  df_test <- df_test %>% mutate(h_kink = if_else(h >= theta, theta, h))
  if (model == "rq"){
    # fit model using training sample data
    fit <- rq(e~h_kink, data = df_train, tau = quantile_levels, method = "fn")
    # predictions for test sample
    pred <- predict(fit, newdata = df_test)
  } else if (model == "crch"){
    fit <- crch(e~1|h_kink, data = df_train, 
                dist = "gaussian", link = "identity", type = "crps")
    pred <- pred_helper_crch(fit, new_x = df_test, prob = quantile_levels)    
  }
  # compute test sample loss
  loss <- quantile_eval(fc = pred, y = df_test$e, tau = quantile_levels)
  list(pred = pred, loss = loss, theta = theta)
}

make_load_mat <- function(h_max){
  # make forecast covariance matrix of fixed-event forecast (= sums over log growth rates)
  load_mat <- matrix(0, h_max, h_max)
  aux2 <- matrix(1, 12, 12)
  aux2[upper.tri(aux2)] <- 0
  load_mat[1:12, 1:12] <- aux2
  for (jj in 13:(h_max)){
    load_mat[jj,(jj-12+1):(jj)] <- rep(1, 12)
  }    
  load_mat
}

make_fe_vcv <- function(h_max, a, s2, load_mat = NULL){
  # make covariance matrix of forecasts at different horizons
  aux <- as.matrix(dist(1:h_max))
  aux <- a^aux
  aux[upper.tri(aux)] <- 0
  fc_v <- s2 * aux %*% t(aux)
  
  if (is.null(load_mat)){
    load_mat <- make_load_mat(h_max)
  }

  fc_v_fe <- load_mat %*% fc_v %*% t(load_mat)
  fc_v_fe
}

# forecast standard deviation for fixed-event forecast (assuming AR1 model)
makesd_alphasigma <- function(h, a, s2, load_mat = NULL){
  # longest horizon considered
  h_max <- max(c(ceiling(max(h)), 13))
  
  # make fixed-event forecast error covariance matrix
  fc_v_fe <- make_fe_vcv(h_max, a, s2, load_mat = load_mat)
  
  # forecast standard deviations
  fc_s_fe <- data.frame(h = 0:h_max, s = c(0, sqrt(diag(fc_v_fe)))) 
  # add linear interpolation for .5 horizons
  tmp <- approx(x = fc_s_fe$h, y = fc_s_fe$s, xout = .5:(h_max - 0.5)) 
  fc_s_fe <- rbind(fc_s_fe, data.frame(h = tmp$x, s = tmp$y)) %>% arrange(h)
  
  # return standard deviations at the horizons in h
  fc_s_fe$s[match(h, fc_s_fe$h)]
}

# objective function of alphasigma model
fn_alphasigma <- function(theta, dat, load_mat = NULL){
  a <- theta[1]
  s2 <- theta[2]^2
  if (abs(a) > 1){
    out <- 1e5
  } else {
    fc_s_fe_out <- makesd_alphasigma(dat$h, a, s2, load_mat = load_mat)
    
    # compute crps
    out <- crps_norm(dat$e, mean = 0, sd = fc_s_fe_out) %>% mean
  }
  return(out)
}

# fit alphasigma model
fit_alphasigma <- function(df_train, df_test, quantile_levels = c(.1, .9), 
                           params = NULL, starting_vals = NULL){
  # if params are not provided: fit model based on training data
  if (is.null(params)){
    if (is.null(starting_vals)){
      starting_vals <- c(.4, .4)
    }
    params <- optim(starting_vals, fn_alphasigma, dat = df_train, 
                    load_mat = make_load_mat(ceiling(max(df_train$h))))$par  
  } 
  # compute predictions
  pred_s <- makesd_alphasigma(df_test$h, a = params[1], s2 = params[2]^2)
  pred <- matrix(NA, nrow(df_test), length(quantile_levels))
  for (jj in 1:length(quantile_levels)){
    pred[,jj] <- qnorm(quantile_levels[jj], mean = 0, sd = pred_s)
  }
  # compute test sample loss
  loss <- quantile_eval(fc = pred, y = df_test$e, tau = quantile_levels)
  list(pred = pred, loss = loss, theta = params)
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

upper_bounds <- function(origin_date, var){
  if (var == "gdp"){
    out <- upper_bounds_gdp(origin_date)
  } else if (var == "inf"){
    out <- upper_bounds_inflation(origin_date)
  }
  if (length(out$ub) != out$n_bins) stop("Wrong info on nr of bins")
  out
}

rescale_p <- function(p){
  if (sum(p) != 1){
    sum_old <- sum(p)
    p <- p * (1/sum(p))
    print(paste("Sum (old):", sum_old))
    print(paste("Sum (new):", sum(p)))
  }
  p
}

# Compute clustered standard errors
# fit is a fitted model (via lm)
# cluster_var is a vector containing the categorical variable used for clustering
# function follows example code in Hansen ("Econometrics", 2022, Section 4.21)
compute_cluster_VCV <- function(fit, cluster_var){
  # Get data and dimensions
  X <- model.matrix(fit)
  k <- ncol(X)
  n <- nrow(X)
  e <- residuals(fit)
  y <- predict(fit) + e
  XXi <- solve(t(X) %*% X)
  # X times e (n*k matrix)
  Xe <- X * (e %*% matrix(1, 1, k))
  # Sum across clusters
  Xe_sum <- rowsum(Xe, cluster_var)
  # Number of clusters
  G <- nrow(Xe_sum)
  # double-check
  if (G != length(unique(cluster_var))) stop("cluster_var seems wrong")
  omega <- t(Xe_sum) %*% Xe_sum
  # Scaling factor
  scl <- (G/(G-1))*((n-1)/(n-k))
  V_clustered <- scl* XXi %*% omega %*% XXi
  se_clustered <- sqrt(diag(V_clustered))
  # Return list of coefficients, standard errors and VCV matrix
  list(beta = coefficients(fit), se = se_clustered, V = V_clustered)
}

# Diebold-Mariano test using clustered standard errors
dm_cluster <- function(loss1, loss2, cluster_var, type = "base"){
  d <- data.frame(ld = loss1 - loss2, cluster_var = cluster_var)
  fit <- lm(ld~1, data = d)
  if (type == "base"){
    aux <- compute_cluster_VCV(fit, cluster_var)
    t <- unname(aux$beta/aux$se)
    # based on sandwich package (robustness check, should yield same result)
  } else if (type == "sandwich"){
    aux <- vcovCL(fit, cluster = ~cluster_var)
    t <- (coefficients(fit)/sqrt(aux)) %>% unname
  }
  list(t = t, p = 2*pnorm(-abs(t))) 
}
