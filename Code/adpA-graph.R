library(ggplot2)


setwd("C:/Users/Death Star/Desktop/Trax_Lab/Code")
points<- read.csv("adpA-experimental.csv", header=TRUE)
time <- points$time
ys <- points$IMV.intensity.Mean.Value.

ggplot(points, aes(x=time, y=ys)) +geom_point() + geom_line()
