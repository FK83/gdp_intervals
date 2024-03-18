This repository contains replication files for the working paper "Prediction intervals for economic fixed-event forecasts" (Krüger and Plett, 2024):

- "Section3_plot_rmsfe.R" plots the root mean squared forecast error (RMSFE) as a function of the forecast horizon

- "Section6_simulation.R" runs the simulation study of the paper, and saves the results in the subfolder "fe_sim/". "Section6_collect_results.R" collects the results in "fe_sim/" and makes the summary table in the paper.

- "Section7_main_results.R" runs the paper's main empirical analysis. The evaluation results are saved as "gdp_us_eval.csv", "inf_us_eval.csv" and "gdp_de_eval.csv". "Section7_hist_comparison.R" runs comparisons to the histogram-based prediction intervals (for US data). "Supplement_zero_mean.R" runs a robustness check that sets the mean of the Gaussian model to zero. "Supplement_tests.R" runs additional analyses reported in the Online Supplement for the paper.

- "gdp_procs23.R" and "dm_procs24.R" collect several R functions that are used in the project. 

- "estimation.R" and "iso_icv_icx.cpp" contain code written by Alexander Henzi, for implementing the method by Henzi (Journal of Business and Economic Statistics 41, 2023). Posted with kind permission. 

- Subfolder "data/": "gdp_us.csv" and  "inf_us.csv" contain point forecasts and realizations of US GDP and inflation. "histograms_gdp.csv" and "individual_histograms_gdp.csv" contain average and individual-level prediction intervals (based on survey 'histograms') for US GDP; analogous files (with "gdp" replaced by "inf") are available for inflation. The US data have been obtained from the Federal Reserve Bank of Philadelphia (Survey of Professional Forecasters and Real-Time Data Research Center). Accessed on 2023-09-21 (individual-level forecast histograms) and 2023-09-06 (other data). "gdp_de.csv" contains forecast and realization data for German GDP kindly provided by the Halle Institute for Economic Research: IWH Forecasting Dashboard, https://www.iwh-halle.de/ForDas. Accessed on 2023-08-14. "rmse_gdp_de.csv" contains empirical root mean squared errors computed from these data (and reported in Figure 1 of our paper). *Disclaimer note*: Of course, the original data providers cannot take responsibility for the accuracy of the data posted here. The original sources should be consulted for the official and most recent version of the data.  
