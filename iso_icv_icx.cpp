/* 
 * This file contains functions for the estimation of conditional distributions
 * under the increasing concave and increasing convex stochastic order. It 
 * requires the additional file "iso_icv_icx.R" that is available from the
 * same source as this file.
 * 
 * Last modified: May 2021
 * 
 */

#include <Rcpp.h>
using namespace Rcpp;

//[[Rcpp::export]]
List iso_icv(
    NumericVector Y,
    NumericVector pos_Y,
    NumericVector w,
    NumericVector thr) {
  /* Estimation of conditional distribution F_x(y) under increasing concave
   * order constraints.
   * 
   * Parameters:
   * - Y: Y-values
   * - w: Number of observations for each of the distinct covariate values x
   * - pos_Y: Position of X[i] in the vector of sorted, distinct covariate 
   *     values x
   * - thr: Thresholds for y at which F_x(y) is to be estimated. It is
   *     assumed that thr[0] = min(Y), thr[thr.length() - 1] = max(Y)
  */
  
  // Preparation
  int k = thr.length();
  int d = w.length();
  int n = Y.length();
  
  NumericMatrix M(d, k); // container for \tilde{M}
  NumericMatrix CDF(d, k); // container for \hat{F}
  
  int n_slopes = k - 1;
  NumericVector dt(n_slopes); // differences between subsequent thresholds
  for (int j = 0; j < n_slopes; j++) dt[j] = thr[j + 1] - thr[j];
  
  double current_thr = thr[n_slopes]; // current threshold t_j for estimation of \tilde{M}_{x_i}(t_j)
  NumericVector trunc_Y(d, 0.0); // corresponds to h_i(t_j), i = 1, \dots, d.

  // Containers for PAV-Algorithm
  int N = k;
  if (d > k) N = d;
  IntegerVector PP(N + 1); // upper boundary of PAVA index partitions
  NumericVector WW(N + 1); // aggregated weights in PAVA
  NumericVector HH(N + 1); // averages in PAVA
  PP[0] = -1;
  WW[0] = 0;
  HH[0] = R_PosInf; // to avoid checking that the running index dd is not zero
  int dd = 1;
  int lwr = 0; // indices for coyping values to M and CDF
  int upr = 0;
  double val = 0;
  
  
  // Estimate \tilde{M}_{x_i}(t_j), i = 1, \dots, d, j = 1, \dots, k
  for (int i = 0; i < d; i++) {
    // Known: \tilde{M}_{x_i}(t_1) = 0 if t_1 = \min(Y_1, \dots, Y_n), so start 
    // at j = 1. \hat{F}_{x_i}(t_k) = 1 if t_k = \max(Y_1, \dots, Y_n), so set
    // CDF(_, n_slopes) to 1.
    M(i, 0) = 0.0;
    CDF(i, n_slopes) = 1.0;
  }
  
  for (int j = n_slopes; j > 0; j--) {
    // compute trunc_Y for current threshold
    for (int s = 0; s < n; s++) {
      if (Y[s] < current_thr) {
        trunc_Y[pos_Y[s]] += (current_thr - Y[s]);
      }
    }
    
    // PAV for given threshold
    PP[1] = 0;
    WW[1] = w[0];
    HH[1] = trunc_Y[0] / w[0];

    for (int i = 1; i < d; i++) {
      dd = dd + 1;
      PP[dd] = i;
      WW[dd] = w[i];
      HH[dd] = trunc_Y[i] / w[i];

      while (HH[dd] >= HH[dd - 1]) {
        dd = dd - 1;
        HH[dd] = HH[dd] * WW[dd] + HH[dd + 1] * WW[dd + 1];
        WW[dd] = WW[dd] + WW[dd + 1];
        HH[dd] = HH[dd] / WW[dd];
        PP[dd] = i;
      }
    }

    // Copy results to matrix MM, set trunc_Y back to 0
    lwr = 0;
    for (int r = 1; r <= dd; r++) {
      upr = PP[r] + 1;
      val = HH[r];
      for (int q = lwr; q < upr; q++) {
        M(q, j) = val;
        trunc_Y[q] = 0.0;
      }
      lwr = upr;
    }

    // Next threshold, reset running index dd to 1
    current_thr = thr[j - 1];
    dd = 1;
  }
  
  // Estimation of CDFs
  HH[0] = R_NegInf; // to avoid checking that dd is positive

  // Compute \tilde{F}_{x_i}(t_j)
  for (int i = 0; i < d; i++) {
    // PAVA
    PP[1] = 0;
    WW[1] = dt[0];
    HH[1] = (M(i, 1) - M(i, 0)) / dt[0];

    for (int j = 1; j < n_slopes; j++) {
      dd = dd + 1;
      PP[dd] = j;
      WW[dd] = dt[j];
      HH[dd] = (M(i, j + 1) - M(i, j)) / dt[j];

      while (HH[dd] <= HH[dd - 1]) {
        dd = dd - 1;
        HH[dd] = HH[dd] * WW[dd] + HH[dd + 1] * WW[dd + 1];
        WW[dd] = WW[dd] + WW[dd + 1];
        HH[dd] = HH[dd] / WW[dd];
        PP[dd] = j;
      }
    }

    // Copy results to matrix CDF
    lwr = 0;
    for (int r = 1; r <= dd; r++) {
      upr = PP[r] + 1;
      val = HH[r];
      for (int q = lwr; q < upr; q++) {
        CDF(i, q) = val;
      }
      lwr = upr;
    }
    
    // Reset running index dd to 1
    dd = 1;
  }
  
  return List::create(_["CDF"] = CDF, _["M"] = M);
}

//[[Rcpp::export]]
List iso_icx(
    NumericVector Y,
    NumericVector pos_Y,
    NumericVector w,
    NumericVector thr) {
  int k = thr.length();
  int d = w.length();
  int n = Y.length();
  
  NumericMatrix M(d, k);
  NumericMatrix CDF(d, k);
  
  double current_thr = thr[0];
  NumericVector trunc_Y(d, 0.0);
  
  int n_slopes = k - 1;
  NumericVector dt(n_slopes);
  for (int j = 0; j < n_slopes; j++) dt[j] = thr[j + 1] - thr[j];
  
  int N = k;
  if (d > k) N = d;
  IntegerVector PP(N + 1);
  NumericVector WW(N + 1);
  NumericVector HH(N + 1);
  PP[0] = -1;
  WW[0] = 0;
  HH[0] = R_NegInf;
  int dd = 1;
  int lwr = 0;
  int upr = 0;
  double val = 0;
  
    for (int i = 0; i < d; i++) {
    M(i, n_slopes) = 0.0;
    CDF(i, n_slopes) = 1.0;
  }
  
  thr.push_back(0.0);
  
  for (int j = 0; j < k; j++) {
    for (int s = 0; s < n; s++) {
      if (Y[s] > current_thr) {
        trunc_Y[pos_Y[s]] += (Y[s] - current_thr);
      }
    }
    
    PP[1] = 0;
    WW[1] = w[0];
    HH[1] = trunc_Y[0] / w[0];
    
    for (int i = 1; i < d; i++) {
      dd = dd + 1;
      PP[dd] = i;
      WW[dd] = w[i];
      HH[dd] = trunc_Y[i] / w[i];
      
      while (HH[dd] <= HH[dd - 1]) {
        dd = dd - 1;
        HH[dd] = HH[dd] * WW[dd] + HH[dd + 1] * WW[dd + 1];
        WW[dd] = WW[dd] + WW[dd + 1];
        HH[dd] = HH[dd] / WW[dd];
        PP[dd] = i;
      }
    }
    
    lwr = 0;
    for (int r = 1; r <= dd; r++) {
      upr = PP[r] + 1;
      val = HH[r];
      for (int q = lwr; q < upr; q++) {
        M(q, j) = val;
        trunc_Y[q] = 0.0;
      }
      lwr = upr;
    }
    
    current_thr = thr[j + 1];
    dd = 1;
  }
  
  for (int i = 0; i < d; i++) {
    PP[1] = 0;
    WW[1] = dt[0];
    HH[1] = (M(i, 1) - M(i, 0)) / dt[0];
    
    for (int j = 1; j < n_slopes; j++) {
      dd = dd + 1;
      PP[dd] = j;
      WW[dd] = dt[j];
      HH[dd] = (M(i, j + 1) - M(i, j)) / dt[j];
      
      while (HH[dd] <= HH[dd - 1]) {
        dd = dd - 1;
        HH[dd] = HH[dd] * WW[dd] + HH[dd + 1] * WW[dd + 1];
        WW[dd] = WW[dd] + WW[dd + 1];
        HH[dd] = HH[dd] / WW[dd];
        PP[dd] = j;
      }
    }
    
    lwr = 0;
    for (int r = 1; r <= dd; r++) {
      upr = PP[r] + 1;
      val = HH[r] + 1.0;
      for (int q = lwr; q < upr; q++) {
        CDF(i, q) = val;
      }
      lwr = upr;
    }
    
    dd = 1;
  }
  
  return List::create(_["CDF"] = CDF, _["M"] = M);
}