require(RSQLite)
require(dplyr)

db <- dbConnect(SQLite(), '../cRegulome.db')

mirna <-  c("hsa-mir-206", "hsa-mir-613", "hsa-mir-10a", "hsa-let-7b", "hsa-let-7c", "hsa-let-7d")
cor_mir <- db %>%
  tbl('cor_mir') %>%
  filter(mirna_base %in% mirna) %>%
  collect
targets_mir <- db %>%
  tbl('targets_mir') %>%
  filter(mirna_base %in% mirna) %>%
  collect

tf_id <- c("AFF4", "AR", "ARID3A", "ARNT", "ARNTL", "ASCL1")
  
cor_tf <- db %>%
  tbl('cor_tf') %>%
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