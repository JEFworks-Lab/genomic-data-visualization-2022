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

library(Rtsne)
set.seed(0) ## reproducibility
emb <- Rtsne(pcs$x[,1:10], dims=2, perplexity=30)$Y
dim(emb)
rownames(emb) <- rownames(mat)
plot(emb, pch=".")

set.seed(0) ## recommend you all set seeds
#com <- kmeans(emb, centers=8)
#com <- kmeans(pcs$x[,1:10], centers=8)
com <- kmeans(mat, centers=8)
names(com)
head(com$cluster)

df <- data.frame(x = emb[,1],
                 y = emb[,2],
                 col = as.factor(com$cluster))
p <- ggplot(data = df,
            mapping = aes(x = x, y = y)) +
  geom_scattermore(mapping = aes(col = col), 
                   pointsize=3)
p1 <- p + theme_classic()
p1

df <- data.frame(x = pos[,1],
                 y = pos[,2],
                 col = as.factor(com$cluster))
p <- ggplot(data = df,
            mapping = aes(x = x, y = y)) +
  geom_scattermore(mapping = aes(col = col), 
                   pointsize=3)
p2 <- p + theme_classic()
p2

library(gridExtra) 
## for arranging plots side by side
grid.arrange(p1, p2, ncol=2)

## What is cluster 6?
## What are the genes upregulated in cluster 6
vi <- com$cluster == 6
head(vi)
table(vi)

## sanity check, I can use data visualization
df <- data.frame(x = emb[,1],
                 y = emb[,2],
                 col = as.factor(vi))
p <- ggplot(data = df,
            mapping = aes(x = x, y = y)) +
  geom_scattermore(mapping = aes(col = col), 
                   pointsize=3)
p1 <- p + theme_classic()
p1

## is this gene differentially upregualted in c6?
g <- 'Gad1'
t.test(mat[vi, g], mat[!vi, g], 
       alternative='greater')
## p-value low, reject the Ho
## p-value high, can't reject Ho
wilcox.test(mat[vi, g], mat[!vi, g], 
       alternative='greater')
## visualize it for myself
df <- data.frame(x = emb[,1],
                 y = emb[,2],
                 col = mat[, g])
p <- ggplot(data = df,
            mapping = aes(x = x, y = y)) +
  geom_scattermore(mapping = aes(col = col), 
                   pointsize=3)
p1 <- p + theme_classic() + 
  scale_color_gradientn("col", 
                        colours = colorRampPalette(
                          c("darkblue", "white", "red"))(50)
  )
p1

## I'm too lazy to go through every gene manually
## use a loop in R to do that for me
for(i in 1:5) {
  print(i)
  ## push onto a stack
}
x = sapply(1:5, function(i) {
  print(i)
  return(i)
})
## use this for doing all our wilcox test
## save the pvalues
pvs <- sapply(colnames(mat), function(g) {
  x = wilcox.test(mat[vi, g], mat[!vi, g], 
      alternative='two.sided')
  return(x$p.value)
})
## correct for multiple testing
table(p.adjust(pvs) < 0.05)
table(pvs < 0.05)
head(sort(-log10(pvs), decreasing=TRUE))

## if low p-value, then high fold changes?
fcs <- sapply(colnames(mat), function(g) {
  x = mean(mat[vi, g])/mean(mat[!vi, g])
  return(x)
})
pvs["Gad1"]
fcs["Gad1"]
pvs["Htr1b"]
fcs["Htr1b"]

## plot only 100 genes for now
df <- data.frame(name = colnames(mat)[1:100], 
                 pvs = -log10(pvs)[1:100], 
                 fcs = log2(fcs)[1:100])
p <- ggplot(df, mapping = aes(x=fcs, y=pvs)) + 
  geom_point() 
p
## volcano plot
p + geom_label(mapping = aes(label = name)) + theme_classic()
## only label significant things with high fold?






