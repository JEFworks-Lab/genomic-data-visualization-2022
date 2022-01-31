## see for ourselves what the real data looks like
data <- read.csv('~/Desktop/class/genomic-data-visualization/data/Visium_Cortex_varnorm.csv.gz')
dim(data)
data[1:5,1:5]

pos <- data[, c('x', 'y')]
rownames(pos) <- data[,1]
head(pos)
plot(pos, pch=".")

gexp <- data[, 4:ncol(data)]
rownames(gexp) <- data[,1]
gexp[1:5,1:5]

# how do I quantiy the number of total genes per spot?
?rowSums
numgenes <- rowSums(gexp)
numgenes
barplot(numgenes)
barplot(sort(numgenes))
hist(numgenes)
# adding pseudocounts when we do logs
hist(log10(numgenes+1), breaks=20)

# let's CPM normalize
normgexp <- gexp/rowSums(gexp)*1e6
normgexp[1:5,1:5]
## sanity check
rowSums(normgexp)

# look at a gene
testgene <- normgexp[, 'Gad1']
hist(testgene, breaks=20)
hist(log10(testgene+1), breaks=20)

# plot 2 genes
library(ggplot2)
c('Gad1', 'Slc32a1') %in% colnames(data)
df <- data.frame(log10(normgexp[, c('Gad1', 'Slc32a1')]+1))
head(df)
ggplot(data = df,
       mapping = aes(x = Gad1, y = Slc32a1)) +
  geom_point()

# If I wanted to make this plot for every pair of genes
# 1600 choose 2 -> 1,279,200
# Use dimensionality reduction -> PCA, tSNE


