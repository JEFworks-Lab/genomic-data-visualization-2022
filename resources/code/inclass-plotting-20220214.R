data <- read.csv('~/Desktop/class/genomic-data-visualization/data/MERFISH_Slice2Replicate2_halfcortex.csv.gz')
dim(data)

pos <- data[, c('x', 'y')]
rownames(pos) <- data[,1]

gexp <- data[, 4:ncol(data)]
rownames(gexp) <- data[,1]

## subsample cells
vi <- sample(rownames(pos), 2000)
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

set.seed(0) ## recommend you all set seeds
com <- kmeans(pcs$x[,1:30], centers=20)
names(com)
head(com$cluster)
cluster <- as.factor(com$cluster)
names(cluster) <- rownames(pos)

library(Rtsne)
set.seed(0) ## reproducibility
emb <- Rtsne(pcs$x[,1:30], dims=2, perplexity=30)$Y
dim(emb)
rownames(emb) <- rownames(mat)
MERINGUE::plotEmbedding(emb, groups=cluster)
MERINGUE::plotEmbedding(pos, groups=cluster, cex=2)

set.seed(1) ## reproducibility
emb2 <- Rtsne(mat, dims=2, perplexity=10)$Y
dim(emb2)
rownames(emb2) <- rownames(mat)

## animation
library(gganimate)

a <- data.frame(emb, cluster, 'emb1')
b <- data.frame(emb2, cluster, 'emb2')
colnames(a) <- colnames(b) <- c('x', 'y', 'cluster', 'type')
df.trans <- rbind(a, b)

p <- ggplot(df.trans, aes(x = x, y = y, col=cluster)) +
  geom_point(show.legend = FALSE) +
  theme_void() 

anim <- p +
  transition_states(type,
                    transition_length = 5,
                    state_length = 1) +
  labs(title = '{closest_state}') +
  theme(plot.title = element_text(size = 28)) +
  enter_fade()

anim

