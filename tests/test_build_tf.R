library(readr)
library(dplyr)
library(purrr)
library(sqlome)
library(tidyr)

test_that("Downloading transcription factors-gene correlations", {
    
    ## get a list of transcription factors
    if (!file.exists('~/cRegulomedb/tmp/tf_list.txt')) {
        tf <- read_lines('https://www.dropbox.com/s/dprzegcahgx6pjv/tf_list.txt?raw=1')
    } else {
        tf <- read_lines('~/cRegulomedb/tmp/tf_list.txt')
    }
    tf <- tf[1:3]
    
    ## read correlations from cistrome cancer
    ## round and transform the data
    df <- map(tf, function(x) {
        try({
            fl <- paste('~/cRegulomedb/tmp/tf/', x, '.all.cor.csv', sep = '')
            if (!file.exists(fl)) {
                url <- tf_url(x, all = TRUE)
            } else {
                url <- fl
            }
            
            tf_cor <- read_csv(url, n_max = 10) # for testing
            names(tf_cor)[1] <- 'feature'
            
            tf_cor <- mutate_if(tf_cor, is.numeric, function(x) round(x, 2) * 100)
            
            cat("Correlation data for", x, "downloaded.\n")
            
            return(tf_cor)
        }, {
            cat("Correlation data for", x, "failed.\n")
            return(NA)
        })
    })
    
    expect_true(is.list(df))
    expect_equal(dim(df[[1]]), c(10, 30))
    
    ## bind rows of indifidual file to one data.frame
    names(df) <- tf
    df <- df[!is.na(df)]
    df <- bind_rows(df, .id = 'tf')
    
    expect_identical(class(df), c("tbl_df", "tbl", "data.frame"))
    expect_equal(unique(df$tf), tf)
})


test_that("Downloading transcription factors-gene targets", {
    
    ## get a list of transcription factors
    if (!file.exists('~/cRegulomedb/tmp/tf_list.txt')) {
        tf <- read_lines('https://www.dropbox.com/s/dprzegcahgx6pjv/tf_list.txt?raw=1')
    } else {
        tf <- read_lines('~/cRegulomedb/tmp/tf_list.txt')
    }
    tf <- tf[1:3]
    
    ## read targets data from cistrome cancer
    ## remove zero correlations and drop cor
    df <- map(tf, function(x) {
        try({
            fl <- paste('~/cRegulomedb/tmp/targets/', 'AFF4', '.cor.csv', sep = '')
            if (!file.exists(fl)) {
                url <- tf_url(x, all = FALSE)
            } else {
                url <- fl
            }
            
            tf_targets <- read_csv(url, n_max = 10) # for testing
            names(tf_targets)[1] <- 'feature'
            
            tf_targets <- gather(tf_targets, study, cor, -feature) %>%
                filter(cor > 0) %>%
                select(-cor)
            
            cat("Gene targets for" , x, "downloaded.\n")
            
            return(tf_targets)
        }, {
            cat("Gene targets for", x, "failed.\n")
            return(NA)
        })
    })
    
    expect_true(is.list(df))
    expect_equal(dim(df[[1]]), c(198, 2))
    
    ## merge targets in one data.frame
    names(df) <- tf
    df <- df[!is.na(df)]
    df <- bind_rows(df, .id = 'tf')
    
    expect_identical(class(df), c("tbl_df", "tbl", "data.frame"))
    expect_equal(names(df), c('tf', 'feature',  'study'))
    expect_equal(unique(df$tf), tf)
})