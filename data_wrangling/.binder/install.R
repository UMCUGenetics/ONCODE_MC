install.packages("e1071")
install.packages("rlang")
install.packages("ggplot2")
install.packages("gplots")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("preprocessCore")
BiocManager::install("sva")
BiocManager::install("SummarizedExperiment")
BiocManager::install("pcaMethods")
BiocManager::install("affy")
BiocManager::install("GEOquery")
