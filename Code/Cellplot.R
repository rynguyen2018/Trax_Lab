library(ggplot2)
library(dplyr)
library(animation)
library(gganimate)

setwd("C:/Users/Ryan/Desktop/gitcode/Trax_Lab/Code")
points<- read.csv("Cell_growth_points.csv", header=TRUE)
xs<- points$x
ys<- points$y
frames <- points$frame
p<-ggplot(points, aes(x= xs, y=ys, frame= frames, cumulative = TRUE))+geom_point()
ani.options(interval = 0.06) #animation speed, seconds per frame
gganimate(p)
