data(iris)

install.packages("ggplot2")
library(ggplot2)

# aes==aesthetic mapping describe how variables in the data are mapped to visual properties (aesthetics) of geoms. 
ggplot(data=iris, mapping=aes(x = Sepal.Length, y = Sepal.Width)) + geom_point(mapping=aes(shape=Species))  
ggplot(data=iris, mapping=aes(x = Sepal.Length, y = Sepal.Width)) + geom_point(mapping=aes(color=Species)) 


ggplot(data=iris, mapping=aes(x = Sepal.Length, y = Sepal.Width)) + 
  geom_point(mapping=aes(color=Species)) + 
  facet_grid(. ~ Species) 

ggplot(data=iris, mapping=aes(x = Sepal.Length, y = Sepal.Width)) + 
  geom_point(mapping=aes(color=Species)) + 
  facet_grid(. ~ Species) + 
  geom_smooth(method="lm")

ggplot(data=iris, mapping=aes(x = Sepal.Length, y = Sepal.Width, color=Species)) + 
  geom_point() + 
  facet_grid(. ~ Species) + 
  geom_smooth(method="lm")

ggplot(data=iris, mapping=aes(x = Sepal.Length, y = Sepal.Width, color=Species)) + 
  geom_point() + 
  facet_grid(. ~ Species) + 
  geom_smooth(method="lm") +
  theme_light()

ggplot(data=iris, mapping=aes(x = Sepal.Length, y = Sepal.Width, color=Species)) + 
  geom_point() + 
  facet_grid(. ~ Species) + 
  geom_smooth(method="lm") +
  theme_void()
