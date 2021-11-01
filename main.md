---
title: "PRISM data"
date: "2019-May-20"
output:
  html_document:
    keep_md: yes
    theme: united
    highlight: tango
---


```r
knitr::opts_chunk$set(echo = TRUE, 
                      message = FALSE, 
                      warning = FALSE, 
                      comment = "", 
                      tidy = TRUE)

# load libraries
library(here)                   # for universal pathways
library(janitor)                # for clean names
library(readxl)                 # read Excel files
library(tidyverse)              # data manipulation
library(magrittr)               # data manipulation
```

Data pre-processing pipeline for metabolite data from the **Programming of Inter-generational Stress Mechanisms (PRISM)** study

### **Data Extraction**
















<br>


**1. Extract metabolite identifiers** 

This metabolite dataset has metadata in the first 13 columns. Extract that and save it separately.


```r
metabolite_IDs <- metabolite[14:1321, 1:13] %>%
    set_colnames(metabolite[13, 1:13]) %>%
    janitor::clean_names()

write_csv(metabolite_IDs, here("data", "1_metabolite-details.csv"))
```
<br>


**2. Extract metabolite concentrations**     

In genetic dataset, the general convention is to have the participants IDs in rows and the metabolites in columns when the number of variables are greater than the number of participants (n > p).  
<br> 

Using `Assay ID` as participant identifier (in rows) and `COMP ID` as metabolite identifiers (in columns). Of all the metabolite identifiers, `COMP ID` has the most complete data so we are going to use it to extract our metabolite information. It is an arbitrary value assigned for identifying metabolites. 


```r
# transpose the information because this is going to be used as column names
compID <- metabolite[c(14:1321), 5] %>%
    t() %>%
    paste0("C", .)


df <- metabolite[c(3, 14:1321), 14:513] %>%
    t() %>%
    as_tibble() %>%
    set_colnames(compID)
colnames(df)[1] = "Assay"
```




We can match `Assay ID` with the `studyID` from the covariate data to merge and create a complete participant history file.  
<br>

### Data Cleaning















<br>

The extracted data set has 500 observations of 1314 variables



1.  **Remove faulty sample points**

    -   Removed 30 individuals who had their blood drawn incorrectly or have withdrawn/been removed from the study

    -   470 observations of 1314 variables


```r
# step 1: load file of individuals that have been marked for second review
remove = read_csv(here("data", "raw", "outliers.csv"))

# step 2: select individuals that have been removed from the main data
remove0 = subset(remove, remove$Exclude == 1) %>%
    select(Assay)

# step 3: remove the erroneous Assay ID
df_clean = df[!df$Assay %in% remove0$Assay, ]

# step 4: remove unnecessary files
rm(remove, remove0)
```

2.  **Add household ID and individual ID information for easy identification**

    -   Removed 11 individuals with duplicate ID

    -   459 observations of 1312 variables


```r
# step 1: extract necessary variables from the ID file
id <- read_csv(here("data", "raw", "PRISM linking IDs_clean.csv")) %>%
    select("Assay", "type", "HID", "IID_2", "Cohort")

# step 2: change the nomenclature to reflect the file to match
id <- id %>%
    mutate(Assay = if_else(id$type == "original", paste0(id$Assay, " M"), paste0(id$Assay,
        " A2")))

# step 3: merge the ID file with the full data file, keep all the metabolite
# data 470 observations of 1313 variables
df_clean_id = merge(id, df_clean, by = "Assay", all.y = TRUE)

# step 4: remove duplicate ID 459 observations of 1312 variables Removed 11
# duplicated IDs
df_clean_id0 = subset(df_clean_id, df_clean_id$type == "original") %>%
    select(-type)

# step 4: remove unnecessary files
rm(id, df_clean_id)
```

3.  **Create a new `studyID` variable**

    \
    There are a few households (HID) that have more than one child in the study, so HID cannot be used as an unique identifier. Concatenate HID-IID to create a new variable: `studyID_2`

    -   459 observations of 1313 variables


```r
# create a duplicate file
clean <- df_clean_id0

# add a new variable, studyID_2
clean0 = clean %>%
    mutate(studyID_2 = paste0(clean$HID, "-", clean$IID_2)) %>%
    select(Assay, studyID_2, everything())
```

4.  **Add serum information**

    -   459 observations of 1314 variables


```r
# step 1: load serum data
serum = read_csv(here("data", "raw", "serum_wk.csv"))

# step 2: create a new unique variable, studyID_2 = HID-IID
serum$iid[which(serum$hhid > "7000")] = serum$iid[which(serum$hhid > "7000")] * 10
serum0 = serum %>%
    mutate(studyID_2 = paste0(serum$hhid, "-", serum$iid), .) %>%
    select(-hhid, -iid)


# step 3: add serum information to the main file
df0 <- merge(serum0, clean0, by = "studyID_2", all.y = TRUE) %>%
    select(Assay, studyID_2, serum_wk, everything())
```



### Pre-processing














<br>


Dataset has samples (459 observations) in rows and covariates (1308 metabolites + 6 variables) in columns. 

_Note_: If you would like to remove data that has >30% missing values, use the filtered dataset. 


```r
# extract all covariates to merge later 459 observations of 6 variables
df_var <- df0[, c(1:6)]


# remove covariates for the preprocessing steps use Assay as sample identifier
# in rownames 459 observation of 1308 metabolites
df <- df0[, -c(2:6)] %>%
    column_to_rownames(var = "Assay")
```
<br> 


**1. Data imputation**    

Missing values in the dataset can be due to two major reasons - either the metabolites are absent or they are below the limit of detection (LOD). The missing data is Missing Not At Random (MNAR) so we should impute the missing values instead of completely removing it from the analysis. <br>

The dataset has `r `sum(is.na(df))` missing datapoints and it is therefore not a statistically viable option to remove the missing data without introducing censoring bias. 
<br> 


Using the minimum value of the metabolite divided by the square root of 2 for data imputation. <br>


```r
# Step 1: create a dummy copy of the dataset
df_impute <- df

# Step 2: impute missing values with minimum/sqrt(2)
for (i in 1:ncol(df_impute)) {
    df_impute[is.na(df_impute[, i]), i] <- min(df_impute[, i]/sqrt(2), na.rm = TRUE)
}

# Step 3: check for missing values verify that all NAs have been replaced
sum(is.na(df_impute))

# Step 4: Save file
df_impute0 = df_impute %>%
    rownames_to_column(var = "Assay") %>%
    merge(df_var, ., by = "Assay")
write_csv(df_impute0, here("data", "3_metabolite-Imp.csv"))
```

<br>
Other methods to consider for data imputation:      
* Use the minimum value/minimum metabolite abundance value                  
* kNN (k-Nearest Neighbors)         
* [No-Skip kNN algorithm](https://link.springer.com/article/10.1007/s11306-018-1451-8) to estimate missing metabolite abundances  
<br>


**2. Volume Normalization**       

We have one sample (RW10351 M) that had their sample extracted at a $60 \mu l$ while the rest of the samples were measured at $85 \mu l$. This can affect the overall abundance value produced by the LC-MS analysis. All normalization methods assume that the signal intensities scale linearly with the metabolite concentration. We are extending this idea and assuming that the metabolite signal intensities scale linearly with the sample volume. 


```r
# Step 1: create a dummy copy of the dataset
df_impute_scale <- df_impute

# Step 2: Volume Normalizatin
df_impute_scale["RW10351 M", ] <- (df_impute_scale["RW10351 M", ]/60) * 85

# Step 3: Save file
df_impute_scale0 = df_impute_scale %>%
    rownames_to_column(var = "Assay") %>%
    merge(df_var, ., by = "Assay")
```



