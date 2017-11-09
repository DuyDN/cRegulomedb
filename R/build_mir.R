tictoc::tic()

library(sqlome)

library(RSQLite)
library(dplyr)
library(tidyr)
library(stringr)
library(data.table)
library(purrr)
library(readr)

library(RTCGA)
library(RTCGA.rnaseq)
library(RTCGA.miRNASeq)

# get the the TCGA cohorts with all three assays
cohorts <- as.character(sqlome_info()$Cohort)

# open a connection to a db
db <- dbConnect(SQLite(), 'cRegulome.db')

# microRNA gene correlaitons
## make cor_mir table
### calculate correlations
cat('1: Calculating microRNA-gene correlations.\n')

df <- map(cohorts[-20], function(x) {  # for testing
    # make names of RTCGA data.frames
    mi <- paste(x, 'miRNASeq', sep = '.')
    m <- paste(x, 'rnaseq', sep = '.')

    # read data.frames
    mi <- get(mi)
    m <- get(m)

    # tidy data.frames
    mi <- mirna_tidy(mi)  # for testing
    m <- mrna_tidy(m)  # for testing

    # calcualte correlations in a tidy data.table
    corr <- cor_make(mi, m, x, tidy = TRUE)
    corr <- as.data.table(corr)

    # print progress
    cat(paste('microRNA-gene correlation for', x, 'is done.\n'))

    # return tidy data.table
    return(corr)
})

### use reduce to merge data.tables
cat('2: Merging microRNA-gene correlations.\n')
df <- Reduce(function(x, y) merge(x, y, all=TRUE), df)

### write cor_mir table to connection db
cat('3: Writing microRNA-gene correlations.\n')
dbWriteTable(db,
             name = 'cor_mir',
             df,
             row.names = FALSE,
             overwrite = TRUE)

### making index on cor_mir
dbSendQuery(db,
            statement = 'create index idx1 on cor_mir (mirna_base);',
            overwrite = TRUE)

## make targets for/microRNA-gene mapping
cat('4: Extracting microRNA-gene targets.\n')

df <- get_targets(unique(df$mirna_base), 'gene')


# write targets table
cat('5: Writing microRNA gene and protein targets.\n')

dbWriteTable(db,
             name = 'targets_mir',
             df,
             overwrite = TRUE)

### making index on targets_mir
dbSendQuery(db,
            statement = 'create index idx2 on targets_mir (mirna_base);',
            overwrite = TRUE)

# disconnect from the db file
dbDisconnect(db)
tictoc::toc()
