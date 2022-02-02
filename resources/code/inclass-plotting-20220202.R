## from last time
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

# CPM normalize
numgenes <- rowSums(gexp)
normgexp <- gexp/rowSums(gexp)*1e6

# plot 2 genes
library(ggplot2)
c('Gad1', 'Slc32a1') %in% colnames(data)
df <- data.frame(log10(normgexp[, c('Gad1', 'Slc32a1')]+1))
head(df)
ggplot(data = df,
       mapping = aes(x = Gad1, y = Slc32a1)) +
  geom_point()

# What type of data are we visualizing here?
## We are visualizing quantitative data of gene expression for Gad1 and Slc32a1
## log10 transformed counts per million normalized with a pseudocount of 1
## (describing this data visualization to a blind person)

# How are we encoding this information?
## We are encoding this information using the geometric primitives
## of points to represent each Visium spot
## and the visual channels of position on the x axis to encode the 
## gene expression magnitude of Gad1
## position on the y axis to encode the gene expression magnitude of Slc32a1

# what trend in the data are you trying to make more salient with 
## your data visualization?
##"I'm trying to make more salient the fact that most cells express
## both Slc32a1 and Gad1"


# plot 3 or 4 genes?
## I want to create a exploratory data visualization to 
## explore the relationship between 3 genes
library(ggplot2)
c('Gad1', 'Slc32a1', 'Aqp4') %in% colnames(data)
df <- data.frame(log10(normgexp[, c('Gad1', 'Slc32a1', 'Aqp4')]+1))
head(df)
ggplot(data = df,
       mapping = aes(x = Gad1, y = Slc32a1)) +
  geom_point(mapping = aes(size = Aqp4))

# visually channel of size on the geometric point primitive to encode
## information about Aqp4 gene expression magnitude
## 4 genes or 5 genes, we eventually run out of visual channel
## reduce dimensionality -> PCA

colnames(normgexp)
mat <- log10(normgexp+1)
dim(mat)
?prcomp
pcs <- prcomp(mat) ## may take a few minutes for large matrices
save(pcs, file="Visium_pcs.RData") ## going to save for next time
names(pcs)
head(pcs$x)
dim(pcs$x)
plot(pcs$sdev[1:30], type="l")
## loadings
head(sort(pcs$rotation[,1], decreasing=TRUE))
plot(pcs$x[,1:2], pch=".")
# explore and interpret PCs
## I wonder how the PCs relate to spatial information
## How do the PCs relate to gene loading, gene expression value
## How can we make explotatory data visualizations to help us explore

dim(pcs$rotation)
head(pcs$rotation)
pcs$rotation[1:5,1:5]

pcs$rotation['Gad1',1:10]
df <- data.frame(pc1 = pcs$x[,1], 
                 pc2 = pcs$x[,9], 
                 col = log10(normgexp[, 'Gad1']+1))
head(df)
ggplot(data = df,
       mapping = aes(x = pc1, y = pc2)) +
  geom_point(mapping = aes(col=col))

head(pos)
df <- data.frame(x = pos[,1], 
                 y = pos[,2], 
                 col = pcs$x[,1])
ggplot(data = df,
       mapping = aes(x = x, y = y)) +
  geom_point(mapping = aes(col=col))



## Homework tip
## for those having a hard time using prcomp
## can use single value decomposition as faster alternative
## will require installing RSpectra package
## install.packages('RSpectra')
## k is number of pcs to be computed 
## saves time through approximation and not calculating all pcs
svd2 <- RSpectra::svds(A = as.matrix(mat), k=30,
                       opts = list(center = TRUE, 
                                   scale = FALSE, 
                                   maxitr = 2000, 
                                   tol = 1e-10))
pcs2 <- svd2$u %*% diag(svd2$d)
rownames(pcs2) <- rownames(mat)
# plot(pcs$x[,1], pcs2[,1]) ## double check results are same
