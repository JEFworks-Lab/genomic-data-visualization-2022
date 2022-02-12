data <- read.csv('~/Desktop/class/genomic-data-visualization/data/Visium_Cortex_varnorm.csv.gz')
dim(data)

pos <- data[, c('x', 'y')]
rownames(pos) <- data[,1]

gexp <- data[, 4:ncol(data)]
rownames(gexp) <- data[,1]

## subsaample cells
set.seed(0)
vi <- sample(rownames(pos), 2000)
head(vi)
length(vi)
pos <- pos[vi,]
gexp <- gexp[vi, ]
plot(pos, pch=".")

## subsample genes
set.seed(0)
vi <- sample(colnames(gexp), 1000)
head(vi)
length(vi)
gexp <- gexp[, vi]

## completely shuffle
gexp[1:5,1:5]
set.seed(0)
gexp <- as.matrix(gexp)
x = gexp[,1]
y = sample(gexp[,1], nrow(gexp), replace=FALSE)
length(x)
length(y)
par(mfrow=c(1,1), mar=rep(5,4))
plot(x,y)
hist(x)
hist(y)
rand <- do.call(cbind, lapply(colnames(gexp), function(g) {
  sample(gexp[,g], nrow(gexp), replace=FALSE)
}))
dim(rand)
dim(gexp)
rownames(rand) <- rownames(gexp)
colnames(rand) <- colnames(gexp)
rand[1:5,1:5]
gexp[1:5,1:5]
gexp=rand

## save
save(pos, gexp, file="mystery_data.RData")
