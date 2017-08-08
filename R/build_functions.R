############################################################################################################
# File Name: Build_functions.R
# Description: These are custom functions used get and tidy data to build an image of miRCancer and Cistrome Cancer databases.
# Created by: Mahmoud Ahmed <mahmoud.s.fahmy@students.kasralainy.edu.eg>
# Usage: Used with build_script.R to build the cRegulome.db
############################################################################################################

#' Extract TCGA info 
#'
#' Extract TCGA info for mRNASeq, miRSeq and RPPA. The function calls the RTCGA::infoTCGA and
#' removes entries that doesn't have all three assays performed; mRNASeq, miRSeq and RPPA.
#'
#' @return A data.frame with Cohort, BCR and a binary column for each of the assays.
#'
#' @examples
#' \dontrun{
#' mirna_info()
#' }
#' @export
mirna_info <- function() {
  infoTCGA() %>%
    dplyr::select(Cohort, BCR, mRNASeq, miRSeq, RPPA) %>%
    dplyr::filter(mRNASeq != 0, miRSeq != 0, RPPA != 0)
}

#' Tidy RTCGA.rnaseq data.frame
#' 
#' This is a custom function to tidy up the RTCGA.rnaseq data.frames. The aim is to strip out unnessecary data
#' and get to a typical bioconductor expression matrix. To achieve that, the function drops the first column,
#' extract the patient bcr and gene symbols and add them to the transposed matrix as column and row names.
#' 
#' @param mrna A RTCGA.rnaseq data.frame
#'
#' @return A tidy expression matrix with genes in the rows and patient samples in the columns.
#' 
#' @examples 
#' \dontrun{
#' rtcga.df <- ACC.rnaseq
#' ACC.mat <- mrna_tidy(rtcga.df)
#' }
#' @export
mrna_tidy <- function(mrna) {
  # drop the first column of the data frame and transpose
  mrna_mat <- t(mrna[, -1])
  
  # extract bcr
  mrna_bcr <- mrna$bcr_patient_barcode
  mrna_bcr <- stringr::str_split(mrna_bcr, 
                                 pattern = '\\.',
                                 simplify = TRUE)[ , 1] %>% 
    unique
  
  # extract gene ids
  mrna_gene <- colnames(mrna)[-1]
  mrna_gene <- stringr::str_split(mrna_gene,
                                  pattern = '\\|',
                                  simplify = TRUE)[, 1]
  
  # confirm dimensions match the bcr and gene ids
  if(dim(mrna_mat)[2] != length(mrna_bcr))
    stop("BCR are not unique.")
  
  if(dim(mrna_mat)[1] != length(mrna_gene))
    stop("Gene IDs are not unique.")
  
  # make columns names with bcr
  colnames(mrna_mat) <- mrna_bcr
  
  # make rownames with gene ids
  ind <- mrna_gene != '?' & !duplicated(mrna_gene)
  mrna_mat <- mrna_mat[ind, ]
  rownames(mrna_mat) <- mrna_gene[ind]
  
  # return tidy matrix
  return(mrna_mat)
}

#' Tidy RTCGA.miRNASeq data.frame
#' 
#' This is a custom function to tidy up the RTCGA.miRNASeq data.frames. The aim is to strip out unnessecary data
#' and get to a typical bioconductor expression matrix. To achieve that, the function drops unneeded columns and rows,
#' extract the patient bcr and microRNA miRBase IDs and add them to the transposed matrix as column and row names.
#' 
#' @param mirna A RTCGA.miRNASeq data.frame
#'
#' @return A tidy expression matrix with microRNAs in the rows and patient samples in the columns.
#' 
#' @examples 
#' \dontrun{
#' rtcga.df <- ACC.miRNASeq
#' ACC.mat <- mirna_tidy(rtcga.df)
#' }
#' @export
mirna_tidy <- function(mirna) {
  # drop unneeded rows from data frame and transpose
  mirna_mat <- mirna %>%
    dplyr::filter(miRNA_ID == 'read_count') %>%
    dplyr::select(-1, -2) %>%
    t
  
  # convert to numeric
  mirna_mat <- apply(mirna_mat, 2, as.numeric)

  # extract bcr
  mirna_bcr <- rownames(mirna)
  mirna_bcr <- stringr::str_split(mirna_bcr, 
                                  pattern = '\\.', 
                                  simplify = TRUE)[, 1] %>%
    unique
  
  # extract miRBase ids
  mirna_mibase <- colnames(mirna)[-1:-2]
  
  # confirm dimensions match the bcr and miRBase ids
  if(dim(mirna_mat)[2] != length(mirna_bcr))
    stop("BCR are not unique.")
  
  if(dim(mirna_mat)[1] != length(mirna_mibase))
    stop("miRBase IDs are not unique.")
  
  # make col and row names
  colnames(mirna_mat) <- mirna_bcr
  rownames(mirna_mat) <- mirna_mibase
  
  # return tidy matrix
  return(mirna_mat)
}

#' Tidy RTCGA.RPPA data.frame
#'
#' This is a custom function to tidy up the RTCGA.RPPA data.frames. The aim is to strip out unnessecary data
#' and get to a typical bioconductor expression matrix. To achieve that, the function drops the first column,
#' extract the patient bcr and antibody ID and add them to the transposed matrix as column and row names.
#' 
#' @param ptn A RTCGA.RPPA data.frame
#'
#' @return A tidy expression matrix with proteins in the rows and patient samples in the columns.
#' 
#' @examples 
#' \dontrun{
#' rtcga.df <- ACC.RPPA
#' ACC.mat <- mirna_tidy(rtcga.df)
#' }
#' @export
ptn_tidy <- function(ptn) {
  # drop unneeded first column of data fram
  ptn_mat <- t(ptn[, -1])
  
  # extract bcr
  ptn_bcr <- ptn$bcr_patient_barcode
  ptn_bcr <- stringr::str_split(ptn_bcr,
                                pattern = '\\.',
                                simplify = TRUE)[ , 1] %>%
    unique
  
  # confirm dimensions match the bcr and antibody ids
  if(dim(ptn_mat)[2] != length(ptn_bcr))
    stop("BCR are not unique.")

  # make colnames
  colnames(ptn_mat) <- ptn_bcr
  
  # return tidy matrix
  return(ptn_mat)
}

#' Make mircoRNA feature correlations
#'
#' Calculate the pearson's correlation between microRNA and feature (gene/protein) matrices.
#' The function calls the R base cor function. More importantly the function prepare the input
#' matrices before the call. First, if matches the col.names/bcr. Second, remove duplicated entries.
#' Finally, it removes the microRNAs and features with expression values less than one. In addation,
#' The function returns a tidy data.frame rather than the raw matrix output of cor.
#' 
#' @param mi A tidy matrix of microRNA expression data.
#' @param m A tidy matrix of gene or protein expression data.
#' @param cohort A charachter string of the TCGA ID.
#'
#' @return A tidy data.frame of four columns; miRBase ID, feature id, correlation and TCGA study name.
#' \dontrun{
#' # calculate correlation for genes
#' ACC.mi <- mirna_tidy(ACC.miRNASeq)
#' ACC.m <- mrna_tidy(ACC.rnaseq)
#' corr <- cor_make(ACC.mi, ACC.m)
#' 
#' # calculate correlation for proteins
#' ACC.mi <- mirna_tidy(ACC.miRNASeq)
#' ACC.m <- mrna_tidy(ACC.rppa)
#' corr <- cor_make(ACC.mi, ACC.m)
#' }
#' @export
cor_make <- function(mi, m, cohort) {
  # get patient bcr for microRNA and remove duplicates
  bcr_mi <- stringr::str_split(colnames(mi),
                               pattern = '\\-',
                               simplify = TRUE)[, 3]
  mi <- mi[, !duplicated(bcr_mi)]
  colnames(mi) <- unique(bcr_mi)
  mi <- as.matrix(mi)
  
  # get patient bcr for feature matrix and remove duplicates
  bcr_m <- stringr::str_split(colnames(m), 
                              pattern = '\\-',
                              simplify = TRUE)[, 3]
  m <- m[, !duplicated(bcr_m)]
  colnames(m) <- unique(bcr_m)
  
  # get intersection of bcr and subset matrices
  ind <- intersect(colnames(m), colnames(mi))
  mi <- mi[rowSums(mi) > 0, ind]
  m <- m[rowSums(m) > 0, ind]
  
  # calculate correlation
  c <- cor(t(m), t(mi))
  
  # tidy correlation
  nms <- c('feature', 'mirna_base', paste(cohort))
  d <- reshape2::melt(c) %>%
    dplyr::mutate(value = round(value, 2)) %>%
    dplyr::filter(abs(value) > .1) %>%
    setNames(nms)
  
  # return tidy correlation data.frame
  return(d)
}

#' Write correlation files
#' 
#' This is a wrapper fuction to call the previous functions. It starts by makeing names of RTCGA data.frames
#' for three assays. Then, calls get, *_tidy on each of the corresponding data.frames. Finally, it writes the
#' expression profile for microRNAs and calls cor_make and writes its ouput to a file.
#' 
#' @param cohort A TCGA study identifier.
#' @param write_profiles A logical to decide whether or not the function exports the microRNA expression profiles.
#' @param cor_rppa A logical to decide whether or not the function calculates the microRNA-proetin correlations.
#'
#' @return Write correlation and profile data to corresponding files. It writes three kind of tables in a 
#'    temporary directory for each study; *_rnaseq.csv, *_rppa.csv and *_profile.csv
#' 
#' @examples 
#' \dontrun{
#' write_cor('ACC')
#' }
#' @export
cor_write <- function(cohort, write_profiles = TRUE, cor_rppa = TRUE) {
  l <- list()
  
  # make TCGA data.frame names
  mi <- paste(cohort, 'miRNASeq', sep = '.')
  m.rnaseq <- paste(cohort, 'rnaseq', sep = '.')
  m.rppa <- paste(cohort, 'RPPA', sep = '.')
  
  # get TCGA data in a list
  l$mi <- mirna_tidy(get(mi))
  l$m.rnaseq <- mrna_tidy(get(m.rnaseq))
  l$m.rppa <- ptn_tidy(get(m.rppa))
  
  # write mirna_profiles
  if(write_profiles == TRUE) {
    melt(l$mi) %>%
      dplyr::select(-2) %>%
      setNames(c('mirna_base', 'count')) %>%
      write.table(file = paste('tmp/mir/', cohort, '_profile.csv', sep = ''),
                  sep = ',',
                  quote = FALSE,
                  row.names = FALSE)
  }
  
  # calculte mrna correlation and write to file
  cor_make(l$mi, l$m.rnaseq, cohort) %>%
    write.table(file = paste('tmp/mir/', cohort, '_rnaseq.csv', sep = ''),
                sep = ',',
                quote = FALSE,
                row.names = FALSE)
  
  # calculate ptn correlation and write to file
  if(cor_rppa == TRUE) {
    cor_make(l$mi, l$m.rppa, cohort) %>%
      write.table(file = paste('tmp/mir/', cohort, '_rppa.csv', sep = ''),
                  sep = ',',
                  quote = FALSE,
                  row.names = FALSE)
  }
}

#' Tidy microRNA-feature correlation tables
#' 
#' This function reads the individual correlation csv files and join them to a single data.frame.
#' 
#' @param dat_type A character string. One of rnaseq or rppa.
#'
#' @return A tidy data.frame of a column for feature ID, microRNA ID and a single column fore each TCGA study
#' @examples 
#' \dontrun{
#' cor_tidy('rnaseq')
#' }
#' @export
cor_tidy <- function(dat_type) {
  # get individual correlation tables
  fls <- list.files('tmp/mir/', 
                    pattern = dat_type,
                    full.names = TRUE)
  
  # read first table and subset it from file names
  frst <- read.csv(fls[1],
                   stringsAsFactors = FALSE) %>%
    data.table::data.table()
  
  rst <- fls[-1]
  
  # loop over the rest of the file names to read and full join them to the first table
  for(i in seq_along(rst)) {
    # read the next correlation file
    nxt <- read.csv(rst[i],
                    stringsAsFactors = FALSE) %>%
      data.table::data.table()
    
    # exclude damaged files; with no entries
    if(nrow(nxt) > 1) {
      # join the next on the first
      frst <- dplyr::full_join(frst, nxt)
    }
    
    # remove the next for effeciecy before exit the iteration
    rm(nxt)
  }
  
  # return the first table with the rest full joined
  return(frst)
}

#' Get microRNA targets
#' 
#' This function uses the targetscan database annotations to get the microRNA gene and protien tragets. For microRNA-gene targets
#' it is a straight forward subsetting the database and converting to ta tidy data.frame with the appropriate IDs. For microRNA-protein
#' targets, it first obtain the gene targets as explained and then change transform the IDs to the corresponding protein/antibody IDs
#' provided by the manufacturers of the used assays.
#'
#' @param feature_type A character string to choose the desired type of annotations. Takes one of 'gene' or 'protein'.
#'
#' @return A tidy data.frame of microRNA in one column 'mirna_base' and targets column 'feature'. When used with the 'protein' another
#'   column is added to distinguish the feature_types genes or proteins.
#'
#' @examples 
#' \dontrun{
#' targets <- list()
#' targets$gene <- get_targets('gene')
#' targets$ptn <- get_targets('protein')
#' targets <- bind_rows(targets, .id = 'feature_type')
#' }
#' 
#' @export
get_targets <- function(feature_type) {
  gene_targets <- AnnotationDbi::toTable(targetscan.Hs.eg.db::targetscan.Hs.egMIRBASE2FAMILY) %>%
    dplyr::left_join(AnnotationDbi::toTable(targetscan.Hs.eg.db::targetscan.Hs.egTARGETS)) %>%
    dplyr::select(-name) %>%
    setNames(c('mirna_base', 'feature')) %>%
    dplyr::mutate(mirna_base = tolower(mirna_base),
                  feature = AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
                                                  feature, 
                                                  'SYMBOL',
                                                  'ENTREZID')$SYMBOL)
  switch(feature_type,
         'gene' = {
           gene_targets
         },
         'protein' = {
           readr::read_csv('https://www.dropbox.com/s/pvcf4y8x4hbkm89/ab_info.csv?raw=1') %>%
             dplyr::left_join(gene_targets,
                              by = c('gene_id' = 'feature')) %>%
             dplyr::select(mirna_base, ab_id) %>%
             setNames(c('mirna_base', 'feature'))
         })
  
}

#' Get transcription factor correlation data
#' 
#' This function downlaods the tarnscription factors correlation data from Cistrome Cancer. Briefly, the function constructs a url
#' for the file using a predefined handle and a suffix. Then, constructs a file name and calls the base function downlaod.file
#'
#' @param tf A character string of the official symbol of the transcription factor
#' @param dir A character string for the temporary directory.
#' @param all A logicall of whether to get the correlation files of all genes (default) or only targets.
#'
#' @return NULL. A csv file is downlaoded to dir
#' \dontrun{
#' # obtain a single file from Cistrome Cancer
#' tf_get('AFF4', dir = 'tmp/tf/')
#' 
#' # obtain multiple files
#' tf <- c('AFF4', 'AR', 'ARID3A')
#' map(tf, tf_get, dir = 'tmp/tf')
#' }
#' @export
tf_get <- function(tf, dir, all = TRUE) {
  # construct a suffix depending of type of file
  # all refers to files containing correlations of all genes not only targets
  if(all == TRUE) {
    suf <- '.all.cor.csv'
  } else {
    suf <- '.cor.csv'
  }
  
  # construct url by adding cistrome handle and file suffix
  tf_url <- paste('http://cistrome.org/CistromeCancer/CancerTarget/examples/', tf, suf, sep = '')
  
  # construct a file name 
  tf_fname <- paste(dir, tf, suf, sep = '')
  
  # downlaod file
  download.file(tf_url, destfile = tf_fname)
}

#' Read transcription factor correlation data
#'
#' This is a custom function to read the local csv files obtained from Cistrome Cancer into a list of tidy data.frames.
#' First, the function constructs a vector of available files in the given directory and reads each of them to a data.frame.
#' Then, it drops the damaged files and names the items of the list with the corresponding transcription factor name. 
#'
#' @param dir A character string of the directory containing the Cistrome Cancer files
#'
#' @return A named list of tidy data.frames. Each item of the list correspond to a transcription factor. The data.frame consist of 
#'   a column containing the names of the features and a single column for each of the TCGA studies containing the correspoinding
#'   correlation value.
#' @examples 
#' \dontrun{
#' # obtain multiple files from Cistrome Cancer
#' tf <- c('AFF4', 'AR', 'ARID3A')
#' map(tf, tf_get, dir = 'tmp/tf')
#' 
#' # read files 
#' cor_tf <- tf_read('tmp/tf/')
#' }
#' @export
tf_read <- function(dir = 'tmp/tf/') {
  # get file names and reading intro a list
  fls <- list.files(dir,
                    pattern = '*.csv',
                    full.names = TRUE)
  tf <- lapply(fls,
               readr::read_csv)
  
  # subset damaged files
  ind <- unlist(lapply(tf, nrow)) > 1
  tf <- tf[ind]
  
  # get names of okay transcription factors' files
  nms <- stringr::str_split(list.files(dir,
                                       pattern = '*.csv'),
                            pattern = '/|\\.',
                            simplify = TRUE)[, 1]
  nms <- nms[ind]
  names(tf) <- nms
  
  return(tf)
}
