data <- read.csv('~/Desktop/class/genomic-data-visualization/data/Visium_Cortex_varnorm.csv.gz')
dim(data)

pos <- data[, c('x', 'y')]
rownames(pos) <- data[,1]

gexp <- data[, 4:ncol(data)]
rownames(gexp) <- data[,1]

## CPM normalize
numgenes <- rowSums(gexp)
normgexp <- gexp/numgenes*1e6
dim(normgexp)

## PCA 
mat <- log10(normgexp + 1)
pcs <- prcomp(mat)
names(pcs)
dim(pcs$x)
save(pcs, file='Visium_pcs.RData')
## load from saved
#load('Visium_pcs.RData')

plot(pcs$sdev[1:100], type="l")
## how many pcs should we use?
## why even use pcs? why not just go from gene expression to tsne?
library(Rtsne)
## what happens if we change perplexity?
emb <- Rtsne(pcs$x[, 1:8], dims=2, perplexity = 30)$Y
rownames(emb) <- rownames(mat)
head(emb)
save(emb, file="Visium_emb.RData")
plot(emb, pch=".")

# Cx3cr1 is not detected in this Visium dataset
# If we're interested in microglia, there are other markers
# Fth1 (ferritin)
df <- data.frame(x = emb[,1],
                 y = emb[,2],
                 col = mat[, 'Fth1'])
library(ggplot2)
library(scattermore)
p <- ggplot(data = df,
            mapping = aes(x = x, y = y)) +
  geom_scattermore(mapping = aes(col = col), pointsize=3)
p

## themes
p0 <- p + theme_classic() + 
  scale_color_gradientn("col", 
                        colours = colorRampPalette(
                          c("darkblue", "white", "red"))(50)
                        )


#######
?kmeans
set.seed(0) ## recommend you all set seeds
#com <- kmeans(emb, centers=4)
#com <- kmeans(pcs$x[,1:8], centers=10)
#com <- kmeans(mat, centers=10)
com <- kmeans(mat[, 'Fth1'], centers=2)
names(com)
head(com$cluster)

df <- data.frame(x = emb[,1],
                 y = emb[,2],
                 col = as.factor(com$cluster))
library(ggplot2)
library(scattermore)
p <- ggplot(data = df,
            mapping = aes(x = x, y = y)) +
  geom_scattermore(mapping = aes(col = col), pointsize=3)
p1 <- p + theme_classic()
p1

df <- data.frame(x = pos[,1],
                 y = pos[,2],
                 col = as.factor(com$cluster))
p <- ggplot(data = df,
            mapping = aes(x = x, y = y)) +
  geom_scattermore(mapping = aes(col = col), pointsize=3)
p2 <- p + theme_classic()

library(gridExtra) ## for arranging plots side by side
grid.arrange(p1, p2, ncol=2)
grid.arrange(p0, p1, ncol=2)

# what happens if I rerun tSNE, will these clusters stay together?
## test: rerun tSNE with more PCs
set.seed(0)
emb2 <- Rtsne(pcs$x[, 1:20], dims=2, perplexity = 30)$Y
rownames(emb2) <- rownames(mat)
head(emb2)
df <- data.frame(x = emb2[,1],
                 y = emb2[,2],
                 col = as.factor(com$cluster))
p <- ggplot(data = df,
            mapping = aes(x = x, y = y)) +
  geom_scattermore(mapping = aes(col = col), pointsize=3)
p2 <- p + theme_classic()
grid.arrange(p1, p2, ncol=2)

## test rerun tsne on random gene expression?
emb3 <- Rtsne(mat[,1:30], dims=2, perplexity = 30)$Y
rownames(emb3) <- rownames(mat)
head(emb3)
df <- data.frame(x = emb3[,1],
                 y = emb3[,2],
                 col = as.factor(com$cluster))
p <- ggplot(data = df,
            mapping = aes(x = x, y = y)) +
  geom_scattermore(mapping = aes(col = col), pointsize=3)
p + theme_classic()
