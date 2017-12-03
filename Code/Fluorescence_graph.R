library(deSolve)
library(ggplot2)
library(parallel)
#setwd("/Users/Echo_Base/Desktop/Trax_code/Code")
#setwd("C:/Users/Death Star/Desktop/Trax_Lab/Code")
setwd("/Users/ryannguyen/Desktop/Trax_code/Code")
points<- read.csv("adpA-experimental_rev2.csv", header=TRUE)
time<- points$Time[0:35]
#dyn.load("gene_circuit.dll")
adpAExpression <- points$Average.IMV[0:35]


my_title <- expression(paste("Fluorescence Activity for ", italic(adpA)))
my_y_title <- expression("IMV Fluorescence")
ggplot(data.frame(adpAExpression), aes(x=time, y= adpAExpression))+geom_line()+
  labs(x="Time (seconds)", y= my_y_title)+
  ggtitle(my_title)+
  theme(plot.title = element_text(family = "Trebuchet MS", color="black", face="bold", size=24, hjust=0)) +
  theme(axis.title = element_text(family = "Trebuchet MS", color="black", face="bold", size=18))
