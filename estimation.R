#-------------------------------------------------------------------------------
# This file contains functions for the estimation of conditional distributions
# under the increasing concave and increasing convex stochastic order. It 
# requires the additional file "iso_icv_icx.cpp" that is available from the
# same source as this file.
#
# Last modified: May 2022

#-------------------------------------------------------------------------------
# Required packages and files
require(isodistrreg)
require(Iso)
require(Rcpp)
Rcpp::sourceCpp("iso_icv_icx.cpp")

#-------------------------------------------------------------------------------
#' Estimate distribution functions under increasing convex/concave order
#' (R version)
#'
#' @param y response variable (numeric vector)
#' @param X covariate (data.frame with one column, same length as \code{y})
#' @param concave estimate under increasing concave (\code{TRUE}) or convex
#'     (\code{FALSE}) order.
#' @param thresholds thresholds at which the conditional distribution functions
#'     are estimated. Equals \code{sort(unique(y))} if omitted.
#' 
#' @details 
#' The vector \code{thresholds}, if specified, should be sorted increasingly. No
#' checks are performed to verify this!
#' 
#' @note 
#' This R implementation is purely for illustrative purposes. It requires
#' \code{\link[Iso]{pava}} for monotone regression. In applications, use the
#' more efficient C++ function \code{iso_icx_icv}
#' 
#' @return 
#' Object of class \code{"idrfit"}, to be used with the functions of the 
#' \code{isodistrreg} package for prediction and forecast evaluation.
iso_icv_icx_r <- function(y, X, concave = TRUE, thresholds = NULL) {
  ## Check input
  if (!is.data.frame(X) || ncol(X) != 1)
    stop("X should be a data.frame with one column")
  if (!is.vector(y, "numeric") || length(y) != nrow(X))
    stop("'y' should be a numeric vector with the same length as X")
  if (anyNA(y) || anyNA(X))
    stop("'X' and 'y' should not contain NA")
  if (!(isTRUE(concave) || isFALSE(concave)))
    stop("'concave' should be TRUE or FALSE")
  if (!is.null(thresholds)) {
    if (!is.vector(thresholds, "numeric"))
      stop("'thresholds' should be a numeric vector")
    y_range <- range(y)
    if (thresholds[1] > y_range[1] ||
        thresholds[length(thresholds)] < y_range[2])
      stop("endpoints of 'thresholds' should span range of 'y'")
  } else {
    thresholds <- sort(unique(y))
  }
  
  ## Aggregate data
  x <- X[, 1, drop = TRUE]
  data <- data.frame(y = y, ind = seq_along(x))
  data <- aggregate(data, by = list(x = x), FUN = identity, simplify = FALSE)
  data$w <- lengths(data$y)
  w <- data$w
  y <- data$y
  n <- nrow(data)
  m <- length(thresholds)
  cdf <- matrix(nrow = n, ncol = m)
  
  ## Estimate CDFs
  if (concave) {
    for (j in 1:m) {
      # Estimate \tilde{M}_{x_i}(y_j)
      thr <- thresholds[j]
      vec <- sapply(y, function(z) mean(pmax(0, thr - z)))
      cdf[, j] <- Iso::pava(y = vec, w = w, decreasing = TRUE)
    }

    dy <- diff(thresholds)
    for (i in 1:n) {
      # Compute \tilde{F}_{x_i}(y_j)
      cdf[i, ] <- c(Iso::pava(y = diff(cdf[i, ]) / dy, w = dy), 1)
    }
  } else {
    for (j in 1:m) {
      # Weighted monotone regression to obtain \hat{F}_{x_i}(y_j)
      thr <- thresholds[j]
      vec <- sapply(y, function(z) mean(pmax(0, z - thr)))
      cdf[, j] <- Iso::pava(y = vec, w = w)
    }

    dy <- diff(thresholds)
    for (i in 1:n) {
      cdf[i, ] <- 1 + c(Iso::pava(y = diff(cdf[i, ]) / dy, w = dy), 0)
    }
  }
  
  ## Export in 'idrfit' class
  name <- colnames(X)
  X <- data.frame(x = data$x)
  colnames(X) <- name
  out <- list(
    thresholds = thresholds,
    cdf = cdf,
    orders = c("comp", 1),
    groups = setNames(1, "x"),
    constraints = NULL,
    indices = data$ind,
    stoch = ifelse(concave, "icv", "icx"),
    X = X,
    diagnostic = list(precision = 0, convergence = 1)
  )
  structure(out, class = "idrfit")
}

#-------------------------------------------------------------------------------
#' Estimate distribution functions under increasing convex/concave order
#' (more efficient C++ version)
#'
#' @param y response variable (numeric vector)
#' @param X covariate (data.frame with one column, same length as \code{y})
#' @param concave estimate under increasing concave (\code{TRUE}) or convex
#'     \code{FALSE} order.
#' @param thresholds thresholds at which the conditional distribution functions
#'     are estimated. Equals \code{sort(unique(y))} if omitted.
#' 
#' @return 
#' Object of class \code{"idrfit"}, to be used with the functions of the 
#' \code{isodistrreg} package for prediction and forecast evaluation.
iso_icv_icx <- function(y, X, concave = TRUE, thresholds = NULL) {
  ## Check input
  if (!is.data.frame(X) || ncol(X) != 1)
    stop("X should be a data.frame with one column")
  if (!is.vector(y, "numeric") || length(y) != nrow(X))
    stop("'y' should be a numeric vector with the same length as X")
  if (anyNA(y) || anyNA(X))
    stop("'X' and 'y' should not contain NA")
  if (!(isTRUE(concave) || isFALSE(concave)))
    stop("'concave' should be TRUE or FALSE")
  if (!is.null(thresholds)) {
    if (!is.vector(thresholds, "numeric"))
      stop("'thresholds' should be a numeric vector")
    y_range <- range(y)
    if (thresholds[1] > y_range[1] ||
        thresholds[length(thresholds)] < y_range[2])
      stop("endpoints of 'thresholds' should be equal to range of 'y'")
  } else {
    thresholds <- sort(unique(y))
  }
  
  ## Prepare data
  x <- X[, 1, drop = TRUE]
  data <- data.frame(y = y, ind = seq_along(x))
  data <- aggregate(data, by = list(x = x), FUN = identity, simplify = FALSE)
  w <- lengths(data$y)
  Y <- unlist(data$y)
  n <- nrow(data)
  pos_Y <- rep(seq_len(n) - 1L, times = w)
  
  ## Estimate CDFs
  isofun <- if (concave) iso_icv else iso_icx
  cdf <- isofun(Y = Y, pos_Y = pos_Y, w = w, thr = thresholds)$CDF
  cdf_range <- range(cdf)
  if (cdf_range[1] < -1e-5 | cdf_range[2] > 1 + 1e-5)
    stop("CDF values outside [0,1] produced!")
  
  ## Export in 'idrfit' class
  name <- colnames(X)
  X <- data[, "x", drop = FALSE]
  colnames(X) <- name
  out <- list(
    thresholds = thresholds,
    cdf = cdf,
    orders = c("comp", 1),
    groups = setNames(1, name),
    constraints = NULL,
    indices = data$ind,
    stoch = ifelse(concave, "icv", "icx"),
    X = X,
    diagnostic = list(precision = 0, convergence = 1)
  )
  structure(out, class = "idrfit")
}

#-------------------------------------------------------------------------------
#' Estimate distribution functions under increasing convex/concave order with
#' subsample aggregating (subagging)
#'
#' @param y response variable (numeric vector)
#' @param X covariate (data.frame with one column, same length as \code{y})
#' @param concave estimate under increasing concave (\code{TRUE}) or convex
#'     \code{FALSE} order.
#' @param thresholds thresholds at which the conditional distribution functions
#'     are estimated. Equals \code{sort(unique(y))} if omitted.
#' @param p fraction of the data in each subsample (number between 0 and 1)
#' @param b number of subagging samples (integer).
#' @param x_new new covariate values at which the conditional distributions
#'     are to be estimated.
#' @param digits precision of the estimated conditional distribution functions.
#'     The conditional distribution functions are rounded to \code{digits}
#'     digits after zero. Making this number larger reduces the required
#'     storage.
#'     
#' @return 
#' Object of class \code{"idr"}, to be used with the functions of the 
#' \code{isodistrreg} package for prediction and forecast evaluation.
iso_icv_icx_bag <- function(
    y,
    X,
    concave = TRUE,
    thresholds = NULL,
    p,
    b,
    x_new,
    digits = 3
  ) {
  if (is.null(thresholds)) thresholds <- sort(unique(y))
  
  m <- length(thresholds)
  n <- length(x_new)
  N <- length(x)
  Np <- ceiling(N * p)
  
  cdf <- matrix(nrow = n, ncol = m, 0)
  for (i in seq_len(b)) {
    sel <- sample(N, Np, replace = FALSE)
    fit <- iso_icv_icx(
      y = y[sel],
      X = data.frame(x = x[sel]),
      concave = concave,
      thresholds = thresholds
    )
    cdf <- cdf + cdf(predict(fit, data.frame(x = x_new)), thresholds)
  }
  cdf <- asplit(pmin(pmax(cdf / b, 0), 1), 1)
  cdf <- lapply(
    cdf,
    function(prb) {
      prb <- round(prb, digits)
      ss <- c(prb[1] > 0, diff(prb) > 0)
      data.frame(points = thresholds[ss], cdf = sort(prb[ss]))
    }
  )
  structure(cdf, class = "idr")
}