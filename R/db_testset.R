require(RSQLite)
require(dplyr)


db <- dbConnect(SQLite(), './cRegulome.db')

mirna <-  c('hsa-let-7b', 'hsa-mir-134')
cor_mir <- db %>%
    tbl('cor_mir') %>%
    dplyr::select(mirna_base, feature, ACC, BLCA) %>%
    filter(mirna_base %in% mirna) %>%
    collect

targets_mir <- db %>%
  tbl('targets_mir') %>%
  filter(mirna_base %in% mirna) %>%
  collect

tf_id <- c('AFF4', 'AR')

cor_tf <- db %>%
    tbl('cor_tf') %>%
    dplyr::select(tf, feature, ACC, BLCA) %>%
    filter(tf %in% tf_id) %>%
    collect()

targets_tf <- db %>%
  tbl('targets_tf') %>%
  filter(tf %in% tf_id) %>%
  collect()

dbDisconnect(db)

db <- dbConnect(SQLite(), 'test.db')
dbWriteTable(db, 'cor_mir', cor_mir)
dbWriteTable(db, 'targets_mir', targets_mir)
dbWriteTable(db, 'cor_tf', cor_tf)
dbWriteTable(db, 'targets_tf', targets_tf)
dbDisconnect(db)
