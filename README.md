# cRegulomedb   
## Overview   
The scripts in this repository are used to build and update the database file used with the R pacakge cRegulome. 
cRegulome is R language interface to access pre-caculated expression correlation of regulome (transcription factors & microRNA)
and genes in cancers.  

## Data sources  
The transcription factors correlation data are obtained from [Cistrome Cancer](http://cistrome.org/CistromeCancer/).
This data is based on integrative analysis of public tarnscription factor ChIP-seq data and The Cancer Genome Atlas
[(TCGA)](https://tcga-data.nci.nih.gov) data. The microRNA correlation data are obtained from [miRCancerdb](https://mahshaaban.shinyapps.io/miRCancerdb/) 
database. [miRCancerdb](https://mahshaaban.shinyapps.io/miRCancerdb/) is based on a similar analysis using [(TCGA)](https://tcga-data.nci.nih.gov)
data and [TargetScan](http://www.targetscan.org) annotations.  

## Build/update the database file  
This repository contains two R scripts; `build_functions.R` and `build_script.R` and a simple `makefile`. By running `make` in
the a UNIX shell, the `build_script.R` loads the required libraries and sources the `build_functions.R`. To do that locally:  

```
git clone https://github.com/MahShaaban/cRegulomedb
cd cRegulomedb
make
```

The make command proceeds to compress the database file and upload it to (). To only buil/update the database for local use,
the `build_db` section of the `makefile` suffices.  

```
make build_db
```

## More    
+ For details the database building steps, [here]().  
+ For more about the cRegulome R package, [here]().  
+ For Citations, [here]().  


