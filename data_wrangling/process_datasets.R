## libraries
require(SummarizedExperiment)
require(affy)
require(sva)
require(preprocessCore)
require(pcaMethods)
require(ggplot2)

# basic ideas
# - have two datasets from only one, the main one to play with, and then
# the second one that needs to be later integrated
# - train a model with the first one, and see if the second dataset agrees with
# it.

# fake stuff to add
# - fake missing data
# - fake bad samples
# - normalization
# - scaling
# - batch correction

###

load('rdata/ds_100550.RData') 

colData(ds_100550)$batch <- 1
colData(ds_100550)$batch[19:36] <- 2
colData(ds_100550)$batch[37:70] <- 3
colData(ds_100550)$batch[71:77] <- 4
colData(ds_100550)$batch[78:85] <- 5
colData(ds_100550)$batch[86:100] <- 6
colData(ds_100550)$batch[101] <- 4
colData(ds_100550)$batch[102:106] <- 5


require(ggplot2)

boxplot(assays(ds_100550)[[1]])

pca_res <- pca(t(assays(ds_100550)[[1]]), nPcs = 5)

plotDf <- as.data.frame(pca_res@scores)
plotDf$CMS <- colData(ds_100550)$CMS
plotDf$batch <- colData(ds_100550)$batch

ggplot(data = plotDf) + 
  geom_point(aes(x = PC1, y = PC2, color = CMS))

ggplot(data = plotDf) + 
  geom_point(aes(x = PC1, y = PC2, color = as.factor(batch)))


pca_res <- pca(t(assays(ds_100550)[[1]][,colData(ds_100550)$batch == 6]), nPcs = 5)

plotDf <- as.data.frame(pca_res@scores)
plotDf$CMS <- colData(ds_100550)$CMS[colData(ds_100550)$batch == 6]
plotDf$batch <- colData(ds_100550)$batch[colData(ds_100550)$batch == 6]
plotDf$batch <- colData(ds_100550)[,5][colData(ds_100550)$batch == 6]



pca_res <- pca(t(assays(ds_100550)[[1]][,colData(ds_100550)$batch == 1]), nPcs = 5)

plotDf <- as.data.frame(pca_res@scores)
plotDf$CMS <- colData(ds_100550)$CMS[colData(ds_100550)$batch == 1]

pca_res <- pca(t(assays(ds_100550)[[1]][,colData(ds_100550)$batch == 3]), nPcs = 5)

plotDf <- as.data.frame(pca_res@scores)
plotDf$CMS <- colData(ds_100550)$CMS[colData(ds_100550)$batch == 3]
plotDf$batch <- colData(ds_100550)[,4][colData(ds_100550)$batch == 3]

## use only the cell line data and create three batches

reduced_dataset <- ds_100550[,1:18]

expr_mat <- assays(reduced_dataset)[[1]]

pca_res <- pca(t(expr_mat), nPcs = 5)


expr_mat1 <- expr_mat
expr_mat2 <- expr_mat
expr_mat3 <- expr_mat

imp_genes <- which(pca_res@loadings[,1] > 0.005)

for (i in 1:ncol(expr_mat)) {
  set.seed(i)
  noise_vec <- abs(rnorm(nrow(expr_mat), mean = 4, sd = 0)[imp_genes])
  expr_mat2[imp_genes, i] <- expr_mat2[imp_genes, i] +  noise_vec
}

imp_genes <- which(pca_res@loadings[,2] > 0.005)
for (i in 1:ncol(expr_mat)) {
  set.seed(i)
  noise_vec <- abs(rnorm(nrow(expr_mat), mean = 4, sd = 0)[imp_genes])
  expr_mat3[imp_genes, i] <- expr_mat3[imp_genes, i] +  noise_vec
}

expr_mat2 <- expr_mat2 + 0.5
expr_mat3 <- expr_mat3 + 0.2


final_mat <- cbind(expr_mat1, expr_mat2, expr_mat3)


pca_res <- pca(t(final_mat), nPcs = 5)

plotDf <- as.data.frame(pca_res@scores)
plotDf$cms <- rep(reduced_dataset$CMS, times = 3)
plotDf$batch <- rep(1:3, each = 18)

ggplot(data = plotDf) + 
  geom_point(aes(x = PC1, y = PC2, color = cms))

ggplot(data = plotDf) + 
  geom_point(aes(x = PC1, y = PC2, color = as.factor(batch)), alpha = .4)

boxplot(final_mat)

fake_coldata <- rbind(colData(ds_100550)[1:18,],
                      colData(ds_100550)[1:18,],
                      colData(ds_100550)[1:18,])
fake_coldata$batch[19:36] <- 2
fake_coldata$batch[37:54] <- 3
fake_coldata <- fake_coldata[,c(1, 6, 7)]

fake_dataset <- SummarizedExperiment(assays = list(expression = final_mat),
                                     rowData = rowData(ds_100550),
                                     colData = fake_coldata)

## remove probes without a gene
fake_dataset <- fake_dataset[which(!is.na(rowData(fake_dataset)[,'Gene Symbol'])),]

expr_data <- as.data.frame(assays(fake_dataset)[['expression']])
expr_data$gene <- rowData(fake_dataset)[,'Gene Symbol']

## summarize the data per gene

expr_data2 <- aggregate(expr_data[, 1:54], list(expr_data$gene), mean)

row_data <- as.data.frame(rowData(fake_dataset))
row_data <- row_data[-which(duplicated(row_data[,'Gene.Symbol'])),]

row_data <- row_data[order(row_data$Gene.Symbol),]
all(row_data$gene == expr_data2$Group.1) 

expr_data2$Group.1 <- NULL

se <- SummarizedExperiment(assays = list(rma = expr_data2), 
                           rowData = row_data,
                           colData = fake_coldata)

## add fake noise
se1 <- se

mat <- assays(se1)[['rma']]
for (i in 1:ncol(mat)) {
  mat[,i] <- mat[,i] + rnorm(1, sd = 0.2)
}
boxplot(mat)

assays(se1)[['rma_fake']] <- as.matrix(mat)


# add fake bad samples

se_fake_samples <- se1[,1]
se_fake_samples$title <- 'K56'
assays(se_fake_samples)[['rma_fake']] <- assays(se_fake_samples)[['rma_fake']] - 4
assays(se_fake_samples)[['rma_fake']][assays(se_fake_samples)[['rma_fake']] < min(assays(se1)[[1]])] <- NA
se1 <- cbind(se1[,c(1:8)], se_fake_samples, se1[,c(9:54)])

se_fake_samples <- se1[,30]
se_fake_samples$title <- 'K7'
assays(se_fake_samples)[['rma_fake']] <- assays(se_fake_samples)[['rma_fake']] - 4
assays(se_fake_samples)[['rma_fake']][assays(se_fake_samples)[['rma_fake']] < min(assays(se1)[[1]])] <- NA
se1 <- cbind(se1[,c(1:25)], se_fake_samples, se1[,c(26:55)])

boxplot(assays(se1)[[2]])

# add missing data

## missing completely at random

set.seed(1)
nas <- sample(seq_len(length(assays(se1)[[2]])), size = 500)

assays(se1)[[2]][nas] <- NA

apply(assays(se1)[[2]], 2, function(x) sum(is.na(x)))

table(apply(assays(se1)[[2]], 1, function(x) sum(is.na(x))))

## missing because of detection limit

low_expressed_genes <- order(apply(assays(se1)[[2]], 1, mean, na.rm = T))[1:50]

for (i in seq_along(low_expressed_genes)) {
  
  set.seed(2 + i)
  nas <- sample(seq_len(ncol(se1)), size = 20)
  assays(se1)[[2]][low_expressed_genes[i], nas] <- NA
  
}

## missing because of really missing (0 expression)

set.seed(15)
nas <- sample(seq_len(nrow(se1)), size = 1)

assays(se1)[[2]][nas, which(colData(se1)$CMS == 'CMS2')] <- NA

## prepare data as fake

fake_se1_coldata <- colData(se1)

fake_se1_coldata$SampleId <- seq_along(fake_se1_coldata$title)
fake_se1_coldata$sample_name <- paste0('sample_', as.character(as.numeric(as.factor(fake_se1_coldata$title))))
fake_se1_coldata$batch <- c(rep(1, 19), rep(2, 19), rep(3, 18))

fake_se1_coldata <- fake_se1_coldata[,c('SampleId', 'sample_name', 'batch', 'CMS')]
rownames(fake_se1_coldata) <- seq_along(fake_se1_coldata$SampleId)

write.table(as.data.frame(fake_se1_coldata), file = 'rawdata/samples_data.txt', sep = '\t')

expr_data <- assays(se1)[[2]]
colnames(expr_data) <- rownames(fake_se1_coldata)

write.table(as.data.frame(expr_data), file = 'rawdata/expr_data.txt', sep = '\t')

gene_data <- as.data.frame(rowData(se1))[,c(1:4)]

write.table(as.data.frame(gene_data), file = 'rawdata/gene_data.txt', sep = '\t')


#### private dataset



# load the datasets

load('rdata/data.cellline.media.RData') 

se <- SummarizedExperiment(assays = list(rma = data$data$rma), 
                           rowData = data$row.annot,
                           colData = data$col.annot)

## remove probes without a gene
se <- se[which(!is.na(rowData(se)$gene)),]

expr_data <- as.data.frame(assays(se)[['rma']])
expr_data$gene <- rowData(se)$gene

## summarize the data per gene

expr_data2 <- aggregate(expr_data[, 1:52], list(expr_data$gene), mean)

row_data <- as.data.frame(rowData(se))
row_data <- row_data[-which(duplicated(row_data$gene)),]

row_data <- row_data[order(row_data$gene),]
all(row_data$gene == expr_data2$Group.1) 

expr_data2$Group.1 <- NULL

se <- SummarizedExperiment(assays = list(rma = expr_data2), 
                           rowData = row_data,
                           colData = data$col.annot)

table(colData(se)$CMS)

# second dataset
se2 <- se[, c(2, 8, 35, 13, 38, 42, 18, 27, 46, 52)]
se1 <- se[, -c(2, 8, 35, 13, 38, 42, 18, 27, 46, 52)]


# add fake batch effect on the first dataset

assays(se1)[['rma_fake']] <- assays(se1)[['rma']]

assays(se1)[['rma_fake']][,se1$Replicate == 2] <- assays(se1)[['rma_fake']][,se1$Replicate == 2] + 1


# add sample noise to normalize for

mat <- assays(se1)[['rma_fake']]
for (i in 1:ncol(mat)) {
  mat[,i] <- mat[,i] + rnorm(1, sd = 0.2)
}
boxplot(mat)

assays(se1)[['rma_fake']] <- as.matrix(mat)

# add fake bad samples

se_fake_samples <- se1[,1]
se_fake_samples$CellLine <- 'K56'
assays(se_fake_samples)[['rma_fake']] <- assays(se_fake_samples)[['rma_fake']] - 4
assays(se_fake_samples)[['rma_fake']][assays(se_fake_samples)[['rma_fake']] < min(assays(se1)[[1]])] <- NA
se1 <- cbind(se1[,c(1:8)], se_fake_samples, se1[,c(9:42)])

se_fake_samples <- se1[,30]
se_fake_samples$CellLine <- 'K7'
assays(se_fake_samples)[['rma_fake']] <- assays(se_fake_samples)[['rma_fake']] - 4
assays(se_fake_samples)[['rma_fake']][assays(se_fake_samples)[['rma_fake']] < min(assays(se1)[[1]])] <- NA
se1 <- cbind(se1[,c(1:25)], se_fake_samples, se1[,c(26:43)])

boxplot(assays(se1)[[2]])

# add missing data

## missing completely at random

set.seed(1)
nas <- sample(seq_len(length(assays(se1)[[2]])), size = 500)

assays(se1)[[2]][nas] <- NA

apply(assays(se1)[[2]], 2, function(x) sum(is.na(x)))

table(apply(assays(se1)[[2]], 1, function(x) sum(is.na(x))))

## missing because of detection limit

low_expressed_genes <- order(apply(assays(se1)[[2]], 1, mean, na.rm = T))[1:50]

for (i in seq_along(low_expressed_genes)) {
  
  set.seed(2 + i)
  nas <- sample(seq_len(ncol(se1)), size = 20)
  assays(se1)[[2]][low_expressed_genes[i], nas] <- NA
  
}

## missing because of really missing (0 expression)

set.seed(15)
nas <- sample(seq_len(nrow(se1)), size = 1)

assays(se1)[[2]][nas, which(colData(se1)$CMS == 'CMS2')] <- NA



## prepare the data as fake data

fake_se1_coldata <- colData(se1)

fake_se1_coldata$SampleId <- seq_along(fake_se1_coldata$SampleId)
fake_se1_coldata$sample_name <- paste0('sample_', as.character(as.numeric(as.factor(fake_se1_coldata$CellLine))))
fake_se1_coldata$batch <- fake_se1_coldata$Replicate

fake_se1_coldata <- fake_se1_coldata[,c('SampleId', 'sample_name', 'batch', 'CMS')]
rownames(fake_se1_coldata) <- seq_along(fake_se1_coldata$SampleId)

write.table(as.data.frame(fake_se1_coldata), file = 'rawdata/samples_data.txt', sep = '\t')

expr_data <- assays(se1)[[2]]
colnames(expr_data) <- rownames(fake_se1_coldata)

write.table(as.data.frame(expr_data), file = 'rawdata/expr_data.txt', sep = '\t')

gene_data <- as.data.frame(rowData(se1))[,c(1:4)]

write.table(as.data.frame(gene_data), file = 'rawdata/gene_data.txt', sep = '\t')



### how to solve the notebook

bad_samples <- which(apply(expr_data, 2, function(x) sum(is.na(x))) > 10000)

expr_data <- expr_data[,-bad_samples]
fake_se1_coldata <- fake_se1_coldata[-bad_samples,]

expr_data<- preprocessCore::normalize.quantiles(as.matrix(expr_data))

pca_res <- pca(t(expr_data), nPcs = 5)


plotDf <- as.data.frame(pca_res@scores)
plotDf <- cbind(plotDf, fake_se1_coldata)
plotDf <- as.data.frame(plotDf)

ggplot() + 
  geom_point(data = plotDf, aes(x = PC1, y = PC4, color = CMS))


p_vals_list <- list()

for (i in seq_len(nrow(expr_data))) {
  print(i)
  tempDf <- data.frame(expr = expr_data[i,],
                       cms = fake_se1_coldata$CMS)
  
  aov_res <- aov(expr~cms, data = tempDf)
  p_vals_list[[i]] <- unlist(summary(aov_res)[[1]][5])[1]
  
}

hist(p.adjust(unlist(p_vals_list)))



cand <- which(p.adjust(unlist(p_vals_list)) < 0.05)

reduced_mat <- expr_data[cand, ]
reduced_mat <- reduced_mat[which(apply(reduced_mat, 1, function(x) sum(is.na(x))) == 0),]

reduced_mat <- as.data.frame(t(reduced_mat))

svm_model <- svm(x = reduced_mat, y = factor(fake_se1_coldata$CMS), kernel = 'linear')

pred <- predict(svm_model, reduced_mat)
table(pred, fake_se1_coldata$CMS[as.numeric(names(pred))])


