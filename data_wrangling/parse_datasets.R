## libraries
require(GEOquery)
require(SummarizedExperiment)

## datasets we want parse
datasets_codes <- c('GSE100550')

## dont judge this mega hardcoded stuff, they should have given a table with
## the labels in the paper

## parse GSE100550 -------------------------------------------------------------

d_file <- list.files('outdated/rawdata', pattern = datasets_codes[1], full.names = TRUE)

gds <- getGEO(filename = d_file)

expression_matrix <- gds@assayData$exprs
sample_data <- gds@phenoData@data
gene_data <- gds@featureData@data

nrow(gene_data) == nrow(expression_matrix)
nrow(sample_data) == ncol(expression_matrix)

## remove the extra sample data that we do not need
colnames(sample_data)
sample_data <- sample_data[,c('title', 'geo_accession', 'source_name_ch1', 
                              'characteristics_ch1.3','characteristics_ch1.4')]

## add the CMS classification as described in the paper
sample_data$CMS <- NA

## based on figure S4
cms_class <- c('CMS1', 'CMS1', 'CMS1', 'CMS1', 'CMS2', 'CMS2', 'CMS2', 'CMS2', 
               'CMS2', 'CMS3', 'CMS3', 'CMS3', 'CMS4', 'CMS4', 'CMS4', 'CMS4',
               'CMS4', 'CMS4')
sample_data$CMS[1:18] <- cms_class

## based on figure 4
## organoids
cms_class <- c('CMS2', 'CMS2', 'CMS1', 'CMS2', 'CMS4', 'CMS2', 'CMS2', 'CMS2', 
               'CMS2', 'CMS1', 'CMS2', 'CMS1', 'CMS2', 'CMS2', 'CMS2', 'CMS4',
               'CMS2', 'CMS2')
sample_data$CMS[19:36] <- cms_class

## spheroids
cms_class <- c('CMS2', 'CMS2', 'CMS2', 'CMS1', 'CMS2', 'CMS4', 'CMS2',
               'CMS4', 'CMS4', 'CMS2', 'CMS2', 'CMS1', 'CMS3', 'CMS2', 'CMS4')
sample_data$CMS[86:100] <- cms_class

## based on figure 3
cms_class <- c('CMS4','CMS2', 'CMS2', 'CMS1', 'CMS2', 'CMS4', 'CMS4', 'CMS1',
               'CMS2','CMS4', 'CMS4', 'CMS2', 'CMS4', 'CMS1', 'CMS4')
sample_data$CMS[71:85] <- cms_class

## based on figure 4
cms_class <- c('CMS2', 'CMS2', 'CMS2', 'CMS4', 'CMS2','CMS4')
sample_data$CMS[101:106] <- cms_class

## based on the names
cms_class <- c('CMS1', 'CMS1', 'CMS1', # co123
               'CMS2', #co128
               'CMS4', #co144
               'CMS4', 'CMS4', 'CMS4', #co149
               'CMS2', 'CMS2', 'CMS2', #co150
               'CMS4', 'CMS4', 'CMS4', #co152
               'CMS1', 'CMS1', 'CMS1', #co153
               'CMS4', 'CMS4', 'CMS4', #co164
               'CMS4', 'CMS4', 'CMS4', #lm27
               'CMS2', 'CMS2', 'CMS2', #lm28
               'CMS2', #lm29
               'CMS1', #lm37
               'CMS2', 'CMS2', #lm38
               'CMS4', 'CMS4', #lm39
               'CMS4', 'CMS4') #lm40
sample_data$CMS[37:70] <- cms_class

sample_data$dataset <- 'GSE100550'

## remove the extra gene info that we do not need
colnames(gene_data)
gene_data <- gene_data[,c('ID', 'GB_ACC', 'Gene Title', 
                          'Gene Symbol', 'ENTREZ_GENE_ID')]


ds_100550 <- SummarizedExperiment(assays = list(rawdata = expression_matrix),
                                  rowData = gene_data, 
                                  colData = sample_data)

save(ds_100550, file = 'rdata/ds_100550.RData')

## parse GSE35566 --------------------------------------------------------------

d_file <- list.files('rawdata', pattern = datasets_codes[2], full.names = TRUE)

gds <- getGEO(filename = d_file)

expression_matrix <- gds@assayData$exprs
sample_data <- gds@phenoData@data
gene_data <- gds@featureData@data

nrow(gene_data) == nrow(expression_matrix)
nrow(sample_data) == ncol(expression_matrix)

## remove the extra sample data that we do not need
colnames(sample_data)
sample_data <- sample_data[,c('title', 'geo_accession', 'source_name_ch1', 
                              'replicate:ch1','mss/msi status:ch1')]

## add the CMS classification as described in the paper
sample_data$CMS <- NA

## based on figure S2
cms_class <- c('CMS1', 'CMS1', #rko
               'CMS1', #sw48
               'CMS?', #ht29
               'CMS2', #t84
               'CMS?', #caco2
               'CMS1', #hct116
               'CMS?', 'CMS?', #skco1
               'CMS?', 'CMS?', #lim2405
               'CMS?', 'CMS?', #hcc2998
               'CMS?', 'CMS?', #ls174t 
               'CMS?', #caco2
               'CMS?', #ht29
               'CMS1', #hct116
               'CMS1') #sw48
sample_data$CMS <- cms_class

sample_data$dataset <- 'GSE35566'

## remove the extra gene info that we do not need
colnames(gene_data)
gene_data <- gene_data[,c('ID', 'GB_ACC', 'Gene Title', 
                          'Gene Symbol', 'ENTREZ_GENE_ID')]


ds_35566 <- SummarizedExperiment(assays = list(rawdata = expression_matrix),
                                 rowData = gene_data, 
                                 colData = sample_data)

save(ds_35566, file = 'rdata/ds_35566.RData')


## now the best part, mix them all ---------------------------------------------

colData(ds_100550)
colData(ds_35566)

## column names that we want
# - cell line
# - CMS
# - geo
# - dataset
# - needs batch
# - batch

sample_data <- as.data.frame(colData(ds_100550))
colnames(sample_data)[1] <- 'cell_line'
sample_data <- sample_data[,c(1, 6, 2, 7)]
sample_data$needs_batch <- FALSE
sample_data$batch <- NA
sample_data$needs_batch[colData(ds_100550)$characteristics_ch1.3 != ''] <- TRUE
sample_data$batch <- sapply(colData(ds_100550)[,4], function(x) unlist(strsplit(as.character(x), ': '))[2])
sample_data$batch[86:100] <- as.numeric(sapply(colData(ds_100550)[,5], function(x) unlist(strsplit(as.character(x), ': '))[2])[86:100]) + 4

colData(ds_100550) <- DataFrame(sample_data)


sample_data <- as.data.frame(colData(ds_35566))
colnames(sample_data)[3] <- 'cell_line'
sample_data <- sample_data[,c(3, 6, 2, 7)]
sample_data$needs_batch <- FALSE
sample_data$batch <- NA

colData(ds_35566) <- DataFrame(sample_data)

dim(ds_100550)
dim(ds_35566)

save(ds_100550, file = 'rdata/ds_100550.RData')
save(ds_35566, file = 'rdata/ds_35566.RData')
