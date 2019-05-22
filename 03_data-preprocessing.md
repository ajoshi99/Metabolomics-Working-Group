---
title: "3. Data pre-processing"
author: "Anu Joshi"
date: "5/20/2019"
output:
  html_document:
    keep_md: yes
    theme: united
---




```r
library(dplyr)                  # for data manipulation
library(magrittr)               # for pipe operator
library(tibble)                 # for row and data manipulation
```
<br>

The data loaded here has been filtered for metabolites with more than 30% missing values. It has samples in rows and metabolites in columns. 

```r
df = read.csv("data/2_metabolites-clean.csv") %>%
    column_to_rownames(var = "Assay")
```
<br> 

**1. Data imputation**     

Missing values are imputed with the minimum value of the metabolite. 

```r
df_impute <- df

# replace all missing values with column minimum (metabolite minimum value)
for(i in 1:ncol(df_impute)){
  df_impute[is.na(df_impute[,i]), i] <- min(df_impute[,i], na.rm = TRUE)
}

# check for missing values to verify that all NAs have been replaced
sum(is.na(df_impute))
```
<br>


**2. Scale data**     

Rescale the data, then setting the median = 1          

Theoretically, the sample mean is generally a more consistent estimator than the sample median, but in skewed distributions and low sample sizes the efficiency of the mean can be impaired. Due to the propensity for extreme outliers in metabolomic data, which could adversely affect the sample mean, the median is used instead. This normalization procedure will be referred to as “MED”, henceforth. [Wuff, J.E. and Mitchell, M.W., 2018](10.4236/abb.2018.98022)   

```r
df_impute_scale <- df_impute

# For each metabolite, divide its value by median across the experimental samples (MED)
for(i in 1:ncol(df_impute_scale)){
  df_impute_scale[, i] <- df_impute_scale[, i]/median(df_impute_scale[,i])
}
```
<br>

**3. Volume Normalization**       

We have one sample (RW10351 M) that had their sample extracted at a volume below the rest. This can affect the LC-MS analysis.         
All normalization methods assume that the signal intensities scale linearly with the metabolite concentration. We are extending this idea and assuming that the metabolite signal intensities scale linearly with the sample volume. So, in order to have the same volume

```r
df_impute_scale["RW10351 M", ] = (df_impute_scale["RW10351 M", ]/60)*85
```
<br>


```r
df_final = df_impute_scale %>% rownames_to_column(var = "Assay")
write.csv(df_final, "data/3_metabolite-VolNormImp.csv")
```
