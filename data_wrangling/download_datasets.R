## libraries
require(GEOquery)

## datasets we want to download
datasets_codes <- c('GSE100550')

data_list <- list()
for (dataset in datasets_codes) {
  data_list[[dataset]] <- getGEO(dataset, destdir = 'outdated/rawdata')
}




