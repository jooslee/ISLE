# ISLE
Code to identify clinically relevant synthetic lethality, ISLE

This repository provides the code to Identify clinically relevant Synthetic LEthality (ISLE) assciated with our manuscript J.S. Lee et al. Harnessing synthetic lethality to predict response to cancer treatment. Nature Communications 9:2546 (2018). The code was tested using R version 3.3.1 (2016-06-21), on a x86 64-pc-linux-gnu (64-bit) platform, using libraries data.table v1.9.6, Rcpp v0.12.7, RcppArmadillo v0.5.500.2.0, survival v2.39-5, ROCR v1.0-7, and caTools v1.17.1. This repository contains the code to identify the main SL networks of our manuscript - genome-wide SL network (the main ISLE pipeline, isle.r) and drug-cSL network for drug response analysis (isle.drug.r).

1. Obtain the code
```
get clone https://github/jooslee/ISLE.git
```
2. Download the data and put them into proper data folders
```
cd ISLE/data
wget ftp://ftp.umiacs.umd.edu/pub/jooslee/prob.TCGA.RData
```
3. Move to the root folder and run R
```
cd ISLE
R
```
4. Install required R libraries
```
> install.packages("data.table")
> install.packages("Rcpp")
> install.packages("RcppArmadillo")
> install.packages("survival")
> install.packages("ROCR")
> install.packages("caTools")
```
5. Launch the code of interest
```
> setwd("ISLE")
> source("./R/isle.r")
> source("./R/isle.drug.r")
