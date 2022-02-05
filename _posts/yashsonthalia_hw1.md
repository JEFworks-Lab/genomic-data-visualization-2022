---
layout: post
title:  "Homework 1 Yash Sonthalia"
author: Yash Sonthalia
categories: [ HW1 ]
image: homework/hw1/yashsonthalia_hw1.png
featured: false
---

```{r}
library(tidyverse)
merfish=read.csv('/Users/yashsonthalia/Documents/JHUCourses/Spring2022/genomic-data-visualization/data/MERFISH_Slice2Replicate2_halfcortex.csv.gz',row.names = 1)
merfish=data.frame(merfish)
merfish_data=merfish[,3:485]

cell_counts=rowSums(merfish_data)
merfish_data_norm=log2(merfish_data/cell_counts*1e6+1)
merfish_data_norm$x=merfish$x
merfish_data_norm$y=merfish$y


astrocytes=ggplot(merfish_data_norm,aes(x,y))+geom_point(mapping=aes(color=Aqp4),size=0.4)+
  xlab('X')+ylab('Y')+ggtitle('Spatial locations of astrocytes')
ggsave(filename = 'yashsonthalia_hw1_gdv.png')
```

