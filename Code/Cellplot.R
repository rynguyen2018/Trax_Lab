library(ggplot2)
library(dplyr)
library(animation)
library(gganimate)

setwd("C:/Users/Ryan/Desktop/gitcode/Trax_Lab/Code")
points<- read.csv("DLA_points.csv", header=TRUE)
xs<- points$x
ys<- points$y
frames <- points$frame
p<-ggplot(points, aes(x= xs, y=ys, frame= frames, cumulative = TRUE))+geom_point()+theme_bw()+ theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank())
#p<- ggplot(points, aes(x=xs, y=ys)) + geom_point()
#print(p)
ani.options(interval = 0.06) #animation speed, seconds per frame
gganimate(p, "Cell_growth.gif")

