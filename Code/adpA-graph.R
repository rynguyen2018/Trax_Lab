library(ggplot2)

#setwd("C:/Users/Ryan/Desktop/gitcode/Trax_Lab/Code")
setwd("/Users/ryannguyen/Desktop/Trax_code/Code")
#setwd("C:/Users/Death Star/Desktop/Trax_Lab/Code")
points<- read.csv("adpA-experimental_rev2.csv", header=TRUE)
time <- points$Time
ys <- points$Concentration_micromolar

ggplot(points, aes(x=time, y=ys)) +geom_point() + geom_line()+ labs(x="Time(seconds)", y= "[adpA](micromolar)", title= "AdpA Concentration")
