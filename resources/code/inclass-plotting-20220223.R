## load data
head(pos)
head(epex)
head(features)

head(pos)
dim(pos)
dim(epex)

range(epex)
table(features$seqnames)
vi <- as.character(features$seqnames) == 'chr19'
table(vi)
subepex <- epex[, vi]
plot(x=1:ncol(subepex), y=colSums(subepex), type="l")
plot(x=1:ncol(subepex), y=subepex[1,], pch=".")
plot(x=1:100, y=subepex[1,1:100], pch=".")
plot(x=1:100, y=subepex[1,1:100], type="l")

par(mfrow=c(1,1), mar=rep(2,4))
MERINGUE::plotEmbedding(pos)
MERINGUE::plotEmbedding(pos, col=rowSums(epex))

