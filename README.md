# Weibull Survival Analysis

This repository contains a course project completed for the subject **Statistical Analysis of Biomedical Data**.

## Project Description

A hospital conducted a study to evaluate the effectiveness of a surgical procedure for a certain type of cancer. The study lasted 10 years. During the first 4 years, new patients were enrolled. Some of them underwent surgery, while others did not. Patients' survival times were observed until the end of the study. The collected data are in file *Dane 8.1*.  

Using the parametric Weibull model, plot the estimated survival curves for both the treatment group (patients who underwent surgery) and the control group (patients who did not) on the same graph. Estimate the median survival time in both groups along with confidence intervals. Test the hypothesis of equality of medians.

## Methods Used

- **Parametric modeling** using the **Weibull distribution**
- **Maximum Likelihood Estimation (MLE)** to estimate shape and scale parameters
- **Delta method** and **bootstrap** to compute confidence intervals for survival medians
- **Likelihood ratio test** and **bootstrap test** to test the equality of medians
- **Visualization** of survival functions with `ggplot2`

The entire analysis was conducted in R.

## Full Analysis

The full implementation with code, results, and plots is available in the Jupyter Notebook:  
📓 [Notebook](project_en.ipynb)


## Authors

This project was prepared by:

- *Krzysztof Madej*
- *Mateusz Mglej*  
