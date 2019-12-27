x <- read.delim(file = 'counts.out', header = FALSE, comment.char = '#')
samples <- c('B14.5', 'B15.5', 'B17.5', 'B20', 'B34',
             'C14.5', 'C15.5', 'C17.5', 'C20', 'C34')
colnames(x) <- c('gene', samples)
corr_data <- cor(x[2:11], method = 'spearman')
heatmap(x = corr_data, symm = TRUE,
        distfun = function(x) {as.dist(1-x)},
        revC = TRUE, main = 'Correlation')

# age <- as.numeric(substring(colnames(x)[-1], 2))
# btissue <- rgb(0, 0, age[1:5], max = 34)
# ctissue <- rgb(age[6:10], 0, 0, max = 34)
# col <- c(btissue, ctissue)
col <- c(rep('blue', 5), rep('red', 5))
pca <- prcomp(t(x[2:11]))
dev <- round(pca$sdev / sum(pca$sdev) * 100, 1)
plot(-pca$x[, 1], -pca$x[, 2], pch = 19, cex = 2,
     col = col, main = 'PCA (normalized)',
     xlab = paste0('PC1 (', dev[1], '%)'), ylab = paste0('PC2 (', dev[2], '%)'))
text(-pca$x[, 1], -pca$x[, 2], labels = row.names(pca$x), pos = 1)

# normalization
keep <- apply(x[2:11], 1, function(row) sum(row >= 1) >= 1)
x <- x[keep, ]
x[2:11] <- sweep(x[2:11], 2, 1e-6 * colSums(x[2:11]), '/')

# retake
library(Rsubread)
files <- list.files(path = 'ME/HSE/NGS/1-3/rebams/', 
                    pattern = '*.bam$',
                    full.names = TRUE)
x <- featureCounts(
  files = files,
  annot.ext = 'ME/HSE/NGS/1-3/merged.gtf',
  isGTFAnnotationFile = TRUE,
  GTF.featureType = 'exon',
  GTF.attrType = 'gene_id',
  nthreads = 8
)
x <- x$counts
colnames(x) <- c('B14.5', 'B15.5', 'B17.5', 'B20', 'B34',
                 'C14.5', 'C15.5', 'C17.5', 'C20', 'C34')
keep <- apply(x, 1, function(row) sum(row >= 1) >= 1)
x <- x[keep, ]
x <- sweep(x, 2, 1e-6 * colSums(x), '/')

corr_data <- cor(x = x, method = 'pearson')
heatmap(x = corr_data, symm = TRUE,
        distfun = function(x) {as.dist(1-x)},
        revC = TRUE, main = 'Correlation')

col <- c(rep('blue', 5), rep('red', 5))
pca <- prcomp(t(x))
dev <- round(pca$sdev / sum(pca$sdev) * 100, 1)
plot(-pca$x[, 1], -pca$x[, 2], pch = 19, cex = 2,
     col = col, main = 'PCA (normalized)',
     xlab = paste0('PC1 (', dev[1], '%)'), ylab = paste0('PC2 (', dev[2], '%)'))
text(-pca$x[, 1], -pca$x[, 2], labels = row.names(pca$x), pos = 1)
