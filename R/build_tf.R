tictoc::tic()

library(sqlome)

library(readr)
library(purrr)
library(tidyr)
library(dplyr)
library(RSQLite)

## get a list of transcription factors
## get a list of transcription factors
if (!file.exists('./tmp/tf_list.txt')) {
    tf <- read_lines('https://www.dropbox.com/s/dprzegcahgx6pjv/tf_list.txt?raw=1')
} else {
    tf <- read_lines('./tmp/tf_list.txt')
}

## read correlations from cistrome cancer
## round and transform the data
cat("1: Downloading transcription factors-gene correlations.\n")

df <- map(tf, function(x) {
    try({
        fl <- paste('./tmp/tf/', x, '.all.cor.csv', sep = '')
        if (!file.exists(fl)) {
            url <- tf_url(x, all = TRUE)
        } else {
            url <- fl
        }
        tf_cor <- read_csv(url) # for testing
        names(tf_cor)[1] <- 'feature'

        tf_cor <- mutate_if(tf_cor, is.numeric, function(x) round(x, 2) * 100)

        cat("Correlation data for", x, "downloaded.\n")

        return(tf_cor)
    }, {
        cat("Correlation data for", x, "failed.\n")
        return(NA)
    })
})

## bind rows of indifidual file to one data.frame
cat("2: Merging correlation data in one table.\n")

names(df) <- tf
df <- df[!is.na(df)]
df <- bind_rows(df, .id = 'tf')

## write table to db
cat("3: Writing correlations in a table.\n")

db <- dbConnect(SQLite(), 'cRegulome.db')
dbWriteTable(db,
             name = 'cor_tf',
             value =  df,
             row.names = FALSE,
             overwrite = TRUE)
dbSendQuery(db,
            statement = 'create index idx3 on cor_tf (tf);')

# targets
## read targets data from cistrome cancer
## remove zero correlations and drop cor
cat("4: Downloading transcription factors-gene targets.\n")

df <- map(tf, function(x) {
    try({
        fl <- paste('./tmp/targets/', x, '.cor.csv', sep = '')
        if (!file.exists(fl)) {
            url <- tf_url(x, all = FALSE)
        } else {
            url <- fl
        }
        
        tf_targets <- read_csv(url) # for testing
        names(tf_targets)[1] <- 'feature'

        tf_targets <- gather(tf_targets, study, cor, -feature) %>%
            filter(abs(cor) > 0) %>%
            select(-cor)

        cat("Gene targets for" , x, "downloaded.\n")

        return(tf_targets)
    }, {
        cat("Gene targets for", x, "failed.\n")
        return(NA)
    })
})

## merge targets in one data.frame
cat("5: Merging targets data in one table.\n")

names(df) <- tf
df <- df[!is.na(df)]
df <- bind_rows(df, .id = 'tf')

## write targets table
cat("6: Writing targets in a table.\n")

dbWriteTable(db,
             name = 'targets_tf',
             value =  df,
             row.names = FALSE,
             overwrite = TRUE)
dbSendQuery(db,
            statement = 'create index idx4 on targets_tf (tf);')

dbDisconnect(db)

tictoc::toc()
