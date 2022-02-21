# bit.ly/GDV_codex
data <- read.csv('codex_spleen_subset.csv.gz', row.names=1)
head(data)

pos <- data[, 1:2]
area <- data[, 3]
pexp <- data[, 4:ncol(data)]
dim(pexp)
dim(pos)
par(mfrow=c(1,1), mar=rep(2,4))
plot(pos, pch=".")

colnames(pexp)
pp <- pexp[, 'CD8']
names(pp) <- rownames(pexp)
MERINGUE::plotEmbedding(pos, col=pp, cex=0.25)

totalpexp <- rowSums(pexp)
hist(log10(totalpexp+1))
mat <- log10(pexp/totalpexp+1)

plot(area, totalpexp, pch=".")
#plot(log10(area+1), log10(totalpexp+1), pch=".")
#plot(area, log10(totalpexp+1), pch=".")

mat2 <- log10(pexp/area+1)

par(mfrow=c(2,2))
k <- "Vimentin"
hist(mat[, k])
hist(mat2[, k])

## total protein normalized
emb <- Rtsne::Rtsne(mat)
emby <- emb$Y
rownames(emby) <- rownames(mat)

## area normalized
emb2 <- Rtsne::Rtsne(mat2)
emby2 <- emb2$Y
rownames(emby2) <- rownames(mat2)

par(mfrow=c(2,2))
plot(emby, pch=".")
plot(emby2, pch=".")

par(mfrow=c(2,2))
MERINGUE::plotEmbedding(emby, col=area, cex=0.1)
MERINGUE::plotEmbedding(emby2, col=area, cex=0.1)

par(mfrow=c(2,2))
MERINGUE::plotEmbedding(emby, col=totalpexp, cex=0.1)
MERINGUE::plotEmbedding(emby2, col=totalpexp, cex=0.1)

com <- kmeans(mat, centers = 10)$cluster
com2 <- kmeans(mat2, centers = 10)$cluster
par(mfrow=c(2,2))
MERINGUE::plotEmbedding(emby, groups=com, cex=0.1)
MERINGUE::plotEmbedding(emby2, groups=com, cex=0.1)
MERINGUE::plotEmbedding(emby, groups=com2, cex=0.1)
MERINGUE::plotEmbedding(emby2, groups=com2, cex=0.1)

par(mfrow=c(2,2), mar=rep(1,4))
MERINGUE::plotEmbedding(pos, groups=com, cex=0.1)
MERINGUE::plotEmbedding(pos, groups=com2, cex=0.1)
MERINGUE::plotEmbedding(emby, groups=com, cex=0.1)
MERINGUE::plotEmbedding(emby2, groups=com2, cex=0.1)

## prof fan tries things
## like not logging; conclusion: not great; outliers in tSNE
mat <- pexp/totalpexp
mat2 <- pexp/area
## transforming area; conclusion: not great; same issues
mat2 <- pexp/log10(area+1)
## do pca first
mat <- log10(pexp/totalpexp+1)
mat2 <- log10(pexp/area+1)
pcs <- prcomp(mat)
pcs2 <- prcomp(mat2)
par(mfrow=c(2,2))
plot(pcs$sdev, type="l")
plot(pcs2$sdev, type="l")
## scree plots suggest should just use all anyway?
emb <- Rtsne::Rtsne(pcs$x[,1:20])
emby <- emb$Y
rownames(emby) <- rownames(mat)
emb2 <- Rtsne::Rtsne(pcs$x[,1:20])
emby2 <- emb2$Y
rownames(emby2) <- rownames(mat2)
par(mfrow=c(2,2))
plot(emby, pch=".")
plot(emby2, pch=".")
com <- kmeans(pcs$x[,1:15], centers = 10)$cluster
com2 <- kmeans(pcs$x[,1:15], centers = 10)$cluster
par(mfrow=c(2,2))
MERINGUE::plotEmbedding(emby, groups=com, cex=0.1)
MERINGUE::plotEmbedding(emby2, groups=com, cex=0.1)
MERINGUE::plotEmbedding(emby, groups=com2, cex=0.1)
MERINGUE::plotEmbedding(emby2, groups=com2, cex=0.1)
MERINGUE::plotEmbedding(pos, groups=com, cex=0.1)
MERINGUE::plotEmbedding(pos, groups=com2, cex=0.1)
## what if we just use a smaller k
com <- kmeans(mat, centers = 5)$cluster
com2 <- kmeans(mat2, centers = 5)$cluster
par(mfrow=c(2,2))
MERINGUE::plotEmbedding(emby, groups=com, cex=0.1)
MERINGUE::plotEmbedding(emby2, groups=com, cex=0.1)
MERINGUE::plotEmbedding(emby, groups=com2, cex=0.1)
MERINGUE::plotEmbedding(emby2, groups=com2, cex=0.1)
MERINGUE::plotEmbedding(pos, groups=com, cex=0.1)
MERINGUE::plotEmbedding(pos, groups=com2, cex=0.1)

## do cells of similar size group together in space?
par(mfrow=c(1,1))
hist(area)
test = scale(area)[,1]
test[test > 1.5] <- 1.5 ## winsorization
test[test < -1.5] <- -1.5
MERINGUE::plotEmbedding(pos, col=test, cex=0.25)



