data(iris)

?class
class(iris)

?dim
dim(iris)

?head
head(iris)

iris$Sepal.Length

iris$Sepal.Width

?plot
plot(iris$Sepal.Length, iris$Sepal.Width)

plot(iris$Sepal.Length, iris$Sepal.Width, pch=16)

class(iris$Species)
iris$Species
as.numeric(iris$Species)
c(15, 16, 17)[iris$Species]
c(15, 16, 17)[as.numeric(iris$Species)]

plot(iris$Sepal.Length, iris$Sepal.Width, pch=c(15, 16, 17)[iris$Species])

c('red', 'green', 'blue')[iris$Species]

plot(iris$Sepal.Length, iris$Sepal.Width, col=c('red', 'green', 'blue')[iris$Species])


