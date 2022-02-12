load('mystery_data.RData')
## what is this?

## first look at the data
class(gexp)
class(pos)
head(pos)
head(gexp)

## Should we normalize?
## because the data is counts -> perhaps we should normalize
## if every cell has the same total number of genes, 
## would normalization be needed?

## recall CPM normalization was
numgenes <- rowSums(gexp)
hist(log10(numgenes+1), breaks=20)
normgexp <- gexp/numgenes*1e6
head(rowSums(normgexp))
mat <- log10(normgexp+1)

## look at how many cells are expressing each gene
gexp[1:5,1:5]
hist(log10(colSums(gexp)+1), breaks=20)

## plot position
plot(pos, pch=".")
'Gad1' %in% colnames(gexp)
'Aqp4' %in% colnames(gexp)
## since ggplot not working, try base plotting
MERINGUE::plotEmbedding(pos, col=mat[, 'Aqp4'])

## Real data analysis, you'd know more about 
## what platform was used to generate this data
## experimental collaborators will tell you what it is

## PCA
pcs <- prcomp(mat)
names(pcs)
dim(pcs$x)
par(mfrow=c(1,2), mar=rep(2,4))
plot(pcs$sdev, type="l")
plot(pcs$sdev[1:30], type="l")
## scree plot looks not quite log decrease, more linear

## plot first two PCs
par(mfrow=c(1,1))
plot(pcs$x[,1:2], pch=".")
plot(pcs$x[,2:3], pch=".")
plot(pcs$x[,3:4], pch=".")

## look at PCS on the spatial coordinates
## maybe PCs capture some type of spatial pattern
MERINGUE::plotEmbedding(pos, col=pcs$x[,1])
MERINGUE::plotEmbedding(pos, col=pcs$x[,2])
MERINGUE::plotEmbedding(pos, col=pcs$x[,3])

## kmeans clustering
set.seed(0)
com <- kmeans(pcs$x[, 1:30], centers=100)
cluster <- as.factor(com$cluster)
par(mfrow=c(1,1))
MERINGUE::plotEmbedding(pos, groups=cluster)

## marker genes in these clusters?
## highly expressed genes
hist(log10(colSums(gexp)+1), breaks=20)
which(log10(colSums(gexp)+1) > 5)

hist(mat[, 'Ttr'])
col = mat[, 'Ttr'] > 3
head(col)
MERINGUE::plotEmbedding(pcs$x[,1:2], groups= col)
MERINGUE::plotEmbedding(pcs$x[,1:2], col=mat[, 'Ttr'])

## guess liver based on markers and lobule structures?
## marker genes Ubb and Ttr
## so we should be able to find a 
## cluster of cells that are hepatocytes

## fast differential expression for now
dg <- MUDAN::getDifferentialGenes(t(mat), cluster)
## z = 1.96 -> p.value 0.05
ubbdg <- sapply(dg, function(x) x['Ttr',])
t(ubbdg)

## testing 2000 genes, in 100 clusters
## by chance see something significant
table(unlist(ubbdg['Z',]) > 1.96)

## tSNE
library(Rtsne)
set.seed(0) ## reproducibility
emb <- Rtsne(pcs$x[,1:30], dims=2, perplexity=30)$Y
dim(emb)
rownames(emb) <- rownames(mat)

plot(emb, pch=".")
MERINGUE::plotEmbedding(emb, col=scale(mat[, 'Ttr'])[,1])
MERINGUE::plotEmbedding(emb, col=scale(mat[, 'Aqp4'])[,1])
MERINGUE::plotEmbedding(emb, groups=cluster)

com2 <- kmeans(emb, centers=8)
cluster2 <- as.factor(com2$cluster)
head(cluster2)
MERINGUE::plotEmbedding(emb, groups=cluster2)

com3 <- kmeans(mat, centers=8)
cluster3 <- as.factor(com3$cluster)
head(cluster3)
MERINGUE::plotEmbedding(emb, groups=cluster3, mark.clusters = TRUE)
dg <- MUDAN::getDifferentialGenes(t(mat), cluster3)
head(dg[[8]][dg[[8]][, 'Z'] < -1.96,])

foo = cluster3==8
names(foo) <- rownames(mat)
MERINGUE::plotEmbedding(pos, groups=foo)


library(Rtsne)
set.seed(10) ## reproducibility
emb2 <- Rtsne(pcs$x[,1:10], dims=2, perplexity=5)$Y
dim(emb2)
rownames(emb2) <- rownames(mat)

MERINGUE::plotEmbedding(emb2, groups=cluster3)

gin <- rownames(dg[[8]][dg[[8]][, 'Z'] < -1.96,])
length(gin)
g = gin[1]

MERINGUE::plotEmbedding(emb2, col=scale(mat[,g])[,1])
MERINGUE::plotEmbedding(pos, col=scale(mat[,g])[,1])
MERINGUE::plotEmbedding(pos, groups=cluster3)

MERINGUE::plotEmbedding(emb2, groups=cluster3, show.legend=TRUE)

## final guesses
## brain
## liver
## artistic portrait
## noisy dataset
## real answer: noise 



                        