data <- read.csv('~/Desktop/class/genomic-data-visualization/data/MERFISH_Slice2Replicate2_halfcortex.csv.gz')
dim(data)

pos <- data[, c('x', 'y')]
rownames(pos) <- data[,1]

gexp <- data[, 4:ncol(data)]
rownames(gexp) <- data[,1]

## subsaample cells
vi <- sample(rownames(pos), 10000)
head(vi)
length(vi)
pos <- pos[vi,]
gexp <- gexp[vi, ]
plot(pos, pch=".")

## CPM normalize
numgenes <- rowSums(gexp)
normgexp <- gexp/numgenes*1e6

## PCA (SVD)
mat <- log10(normgexp + 1)
pcs <- prcomp(mat)
names(pcs)
dim(pcs$x)
plot(pcs$sdev[1:30], type="l")

library(ggplot2)
df <- data.frame(pc1 = pcs$x[,1], 
                 pc2 = pcs$x[,2], 
                 col = pcs$x[,3])
head(df)
ggplot(data = df,
       mapping = aes(x = pc1, y = pc2)) +
  geom_point(mapping = aes(col=col))

plot(pcs$x[,1], pcs$x[,2])

library(Rtsne)
## install.packages('Rtsne')
?Rtsne
emb <- Rtsne(pcs$x[,1:10], dims=2, perplexity=30)$Y
dim(emb)
rownames(emb) <- rownames(mat)
plot(emb, pch=".")

df <- data.frame(x = emb[,1], 
                 y = emb[,2], 
                 col = mat[,'Gad1'])
head(df)
library(scattermore) ## faster scatterplot
ggplot(data = df,
       mapping = aes(x = x, y = y)) +
  geom_scattermore(mapping = aes(col=col), pointsize=4)

## tinker around with is
## perplexity parameter
## number of PCs
## why even bother PCs? tSNE on the gene expression

## save for next time
save.image()
