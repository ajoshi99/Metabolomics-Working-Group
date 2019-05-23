---
title: "3. Data pre-processing"
author: "Anu Joshi"
date: "5/20/2019"
output:
  html_document:
    keep_md: yes
    theme: united
    highlight: tango
---




```r
library(dplyr)                  # for data manipulation
library(magrittr)               # for pipe operator
library(tibble)                 # for row and data manipulation
```
<br>

This dataset has samples in rows and metabolites in columns. <br>   
_Note_: If you would like to remove data that has >30% missing values, then use the filtered dataset. 

```r
df = read.csv("data/2_metabolites.csv") %>%
    column_to_rownames(var = "Assay")
```
<br> 

**1. Data imputation**    
Missing values in the dataset can be due to two major reasons - either the metabolites are absent or they are below the detection limit. The missing data is Missing Not At Random (MNAR) so we should impute the missing values instead of completely removing it from the analysis. <br>

```r
sum(is.na(df))
```
The dataset has 120,292 missing datapoints and it is therefore not a statistically viable option to remove the missing data without introducing censoring bias. 
<br> 


Using the minimum value of the metabolite divided by the square root of 2 for data imputation. <br>

```r
# Step 1: create a dummy copy of the dataset
df_impute <- df

# Step 2: impute missing values with minimum/sqrt(2) 
for(i in 1:ncol(df_impute)){
  df_impute[is.na(df_impute[,i]), i] <- min(df_impute[, i]/sqrt(2), na.rm = TRUE)
}

# Step 3: check for missing values
# verify that all NAs have been replaced
sum(is.na(df_impute))
```
<br>
Other methods to consider for data imputation:      
* Use the minimum value/minimum metabolite abundance value                  
* kNN (k-Nearest Neighbors)         
* No-Skip kNN algorithm to estimate missing metabolite abundances ([link](https://link.springer.com/article/10.1007/s11306-018-1451-8))   
<br>


**2. Scale data**      
Theoretically, the sample mean is generally a more consistent estimator than the sample median, but in skewed distributions and low sample sizes the efficiency of the mean can be impaired. Due to the propensity for extreme outliers in metabolomic data, which could adversely affect the sample mean, the median is used instead. This normalization procedure will be referred to as “MED”, henceforth. [Wuff, JE and Mitchell, MW, 2018](10.4236/abb.2018.98022)   <br>

```r
# step 1: create a dummy copy of the dataset
df_impute_scale <- df_impute

# For each metabolite, divide its value by median across the experimental samples (MED)
# This will set the median value as 1
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
