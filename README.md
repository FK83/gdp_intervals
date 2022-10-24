This repository contains replication files for the working paper "Prediction intervals for economic fixed-event forecasts" (Krüger and Plett, 2022):

- "quantile_procs.R" collects several R functions that are used in the project.

- "Section6_run_simulation.R" runs the simulation study of the paper, and saves the results in the subfolder "fe_sim/". "Section6_collect_results.R" collects the results in "fe_sim/" and makes the summary plot in the paper.

- "gdp_forecasts_DE.csv" contains point forecasts and realizations of German GDP. "gdp_forecasts_US.csv" and  "inf_forecasts_US.csv" contain point forecasts and realizations of US GDP and inflation. "forecast_histograms_gdp.csv" and "forecast_histograms_inf.csv" contain predictions intervals (based on survey 'histograms') for the US. *Disclaimer note*: These forecast and realization data files have been constructed from public sources as described in the paper. Of course, the original data providers cannot take responsibility for the accuracy of the data posted here. The original sources should be consulted for the official and most recent version of the data. 

- "Section7_stats.R" prints summary information on the data sets. "Section7_main_analysis.R" runs the main empirical analysis. The evaluation results are saved as "evaluation_DE.csv", "evaluation_US_gdp.csv" and "evaluation_US_inf.csv". "Section7_hist_comparison.R" runs comparisons to the histogram-based prediction intervals (for US data).
