---
layout: page
title:  "Lesson 3: Summary Statistics"
---

# Lecture 3

## 3.0 Lesson learning objectives

By the end of this lesson, we should understand what is spatially resolved transcriptomic data, how the data is generated, and how we can begin visualizing and interacting with the data.

## 3.1 Why summarize?

## 3.1.1 Types of summary statistics

- Sum (total)
- Mean (average)
- Spread (variance)
- Dependence (pairwise correlations)

## 3.2 Aesthetic mappings in R

Mapping categorical data to colors:

```r
fac2col <- function (x, 
	level.colors = NULL,
	na.col = "lightgrey") {
	
	x <- as.factor(x)	
	if (is.null(names(x))) {
		names(x) <- seq(length(x))
	}
	
	if (is.null(level.colors)) {
		col <- rainbow(length(levels(x)), s = 1, v = 1)
	}
	else {
		col <- level.colors[1:length(levels(x))]
	}
	names(col) <- levels(x)
	
	y <- col[as.integer(x)]
	names(y) <- names(x)	
	y[is.na(y)] <- na.col
	
	return(y)
}

```

Mapping from quantitative data to colors:

```r
map2col <- function (x, 
	pal = colorRampPalette(c("blue", "grey", "red"))(100), 
	na.col = "lightgrey") {
	
	if (is.null(names(x))) {
		names(x) <- seq(length(x))
	}
	
	original <- x
	x <- na.omit(x)
		
	limits <- range(x)	
	y <- pal[findInterval(x, seq(limits[1], limits[2], length.out = length(pal) + 1), 
		all.inside = TRUE)]
	names(y) <- names(x)
	
	colors <- rep(na.col, length(original))
	names(colors) <- names(original)
	
	colors[names(y)] <- y
	
	return(colors)
}
```

---

# Hands-on component 3

Our in-class hands-on component will analyzing either the MERFISH or Visium dataset to create a data visualization. We will also do an in-class demonstration of the homework submission process. 

---

# Class Lesson Notes 3

Notes from class will be uploaded here.

---

# Homework 3

Make a new data visualization of either the MERFISH or Visium dataset using ggplot in R (do not make the same visualizations that were made in class). Write a description of what you did using vocabulary terms from Lesson 1. Make a pull request to submit your homework (due Thursday).


