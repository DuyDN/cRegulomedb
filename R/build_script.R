################################################################################################################
# File Name: Build_script.R
# Description: This script calls custom functions to extract, tidy and process TCGA and Cistrome Cancer data
# Created by: Mahmoud Ahmed <mahmoud.s.fahmy@students.kasralainy.edu.eg>
# Usage: Used with build_functions.R to build the cRegulome.db
################################################################################################################

# load/require libraries
# load/require RTCGA and data pacakges
require(RTCGA)
require(RTCGA.miRNASeq)
require(RTCGA.rnaseq)
require(RTCGA.RPPA)

# load/require annotation packages
require(targetscan.Hs.eg.db)
require(org.Hs.eg.db)

# load/require packages for data wrangling
require(readr)
require(dplyr)
require(purrr)
require(reshape2)
require(stringr)
require(data.table)

# load/require packages for creating the db
require(RSQLite)

# sourece the custom functions
source('R/build_functions.R')

# It's sellf-explainatory if you ask me. But just in case I need to use it a week latter :)
# First, a call to mirna_info to get only TCGA studies with all three assays performed.
# Second, a call to cor_write. This is where the magic happen, but it takes a while so hold on.
#         Meanwhile you can read what the function does in build_functions.R. cor_write cheats
#         and write the microRNA expression profiles to files on the way.
# Third, read the cor_write output and write it to a db, a sqlite db.
# Finlly, get the targets data and write them to the db, the same sqlite db.

# build mir tables
## get cohorts
cohorts <- mirna_info()$Cohort[-20] %>%
  as.character

## make correlations
dir.create('tmp')
dir.create('tmp/mir/')
map(cohorts, cor_write, write_profiles = FALSE, cor_rppa = FALSE)

## tidy correlations
cor_mir <- cor_tidy('rnaseq')
cor_mir <- cor_mir %>%
  mutate_at(vars(3:ncol(cor_mir)), 
            function(x) x * 100)
mirna <- unique(cor_mir$mirna_base)

## tidy and write tables to the db
db <- dbConnect(SQLite(), 'cRegulome.db')
dbWriteTable(db, 
             name = 'cor_mir',
             value =  cor_mir, row.names = FALSE)
dbSendQuery(db,
            statement = 'create index idx1 on cor_mir (mirna_base);')
rm(cor_mir)

## targets
targets_mir <- get_targets('gene') %>%
  filter(mirna_base %in% mirna)

db <- dbConnect(SQLite(), 'cRegulome.db')
dbWriteTable(db,
             name = 'targets_mir',
             value = targets_mir,
             row.names = FALSE)
dbSendQuery(db, 
            statement = 'create index idx2 on targets_mir (mirna_base);')
dbDisconnect(db)
rm(targets_mir)

# build tf tables
## get correlations
tf <- read_lines('https://www.dropbox.com/s/dprzegcahgx6pjv/tf_list.txt?raw=1')

dir.create('tmp/tf')
map(tf, tf_get, dir = 'tmp/tf/')

## tidy correlations
cor_tf <- tf_read() %>%
  bind_rows(.id = 'tf')
cor_tf <- cor_tf %>%
  mutate_at(vars(3:ncol(cor_tf)), 
            function(x) round(x, 2) * 100)

names(cor_tf)[2] <- 'feature'

## write table to db
db <- dbConnect(SQLite(), 'cRegulome.db')
dbWriteTable(db, 
             name = 'cor_tf',
             value =  cor_tf, row.names = FALSE)
dbSendQuery(db,
            statement = 'create index idx3 on cor_tf (tf);')
dbDisconnect(db)
rm(cor_tf)

## targets
## download target files
dir.create('tmp/targets')
map(tf, tf_get, dir = 'tmp/targets/', all = FALSE)

## read and tidy targets
targets_tf <- tf_read('tmp/targets/')
targets_tf <- lapply(targets_tf, '[', 1)
targets_tf <- melt(targets_tf) %>%
  setNames(c('featrue', 'tf'))

## write targets_tf to db
db <- dbConnect(SQLite(), 'cRegulome.db')
dbWriteTable(db, 
             name = 'targets_tf',
             value =  targets_tf, row.names = FALSE)
dbSendQuery(db,
            statement = 'create index idx4 on targets_tf (tf);')
dbDisconnect(db)
rm(targets_tf)
