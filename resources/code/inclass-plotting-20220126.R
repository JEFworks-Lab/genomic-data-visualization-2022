####################
# Install Rstudio (free)
# coding in R
# (comments in R are demarked with #)
## I like to use 2 ##
####################

# use control+enter to send commands 
# to the command prompt
data(iris)
dim(iris)
head(iris)

class(iris)
iris$Sepal.Length
iris[, 'Sepal.Length'] ## case sensitive
iris[, 1]

x = iris$Sepal.Length
y = iris$Sepal.Width
?plot

## plotting in "base R"
plot(x,y)
head(iris)

## use instead a package ggplot2
## grammar of graphics
## uses geometric primitives
## "aesthetic mappings" describe how data is mapped
## to these visual channels for geometric primitives
#install.packages('ggplot2')
library(ggplot2)
?ggplot
ggplot(data=iris, mapping=aes(x = Sepal.Length,
                              y = Sepal.Width)) +
  geom_point()

ggplot(data=iris, mapping=aes(x = Sepal.Length,
                              y = Sepal.Width)) +
  geom_point(mapping=aes(shape=Species))

ggplot(data=iris, mapping=aes(x = Sepal.Length,
                              y = Sepal.Width)) +
  geom_point(mapping=aes(col=Species))

## if you get ahead, try tinkering around
## come back together as a class and share what we found
## google around ggplot2 cheat sheets
colnames(iris)
ggplot(data=iris, mapping=aes(x = Sepal.Width,
                              y = Petal.Width)) +
  geom_point(mapping=aes(col=Species)) +
  facet_grid(. ~ Species) +
  geom_smooth(method='lm')


