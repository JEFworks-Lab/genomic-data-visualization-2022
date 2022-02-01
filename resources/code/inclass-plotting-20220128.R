# read in from github download
# these are pretty big
# MERFISH data will have more cells
# Visium data will have more genes
data <- read.csv('~/Desktop/class/genomic-data-visualization/data/MERFISH_Slice2Replicate2_halfcortex.csv.gz')
dim(data)
head(data)
data[1:5,1:5]
## first col is cell name
## second is x coordinate
## third y coordinate of cell
## rest are genes

pos <- data[, c('x', 'y')]
head(pos)
rownames(pos) <- data[,1]
head(pos)
?plot
plot(pos, pch='.') # change the geometric primitive point

# Read in your dataset
# tinker round and plot as you wish

colnames(data)
library(ggplot2)
c('Gad1', 'Slc32a1') %in% colnames(data)
df <- data.frame(data[, c('Gad1', 'Slc32a1')])
head(df)
ggplot(data = df,
       mapping = aes(x = Gad1, y = Slc32a1)) +
  geom_point()

