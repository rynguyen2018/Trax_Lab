
library(gplots)
require(VennDiagram)


gene_fold_data= read.csv('C:/Users/Ryan/Desktop/Traxler_Lab/Code/Log2Data.csv', header =T)
# pick the gene fold data in Traxler Lab folder
# make sure it is a csv file! 



#allows for extraction and separation of data by headers within text files 
attach(gene_fold_data)


#Standard Deviation of all of the samples 

#SPB74
SPB74_standev <-sd(SPB74.log.2,na.rm=TRUE)
SPB74_Plus_significant<- subset(id,SPB74.log.2>(mean(SPB74.log.2,na.rm=T)+abs(2*SPB74_standev))) #+2 sigma
SPB74_Neg_significant <- subset(id,SPB74.log.2<(mean(SPB74.log.2,na.rm=T)-abs(2*SPB74_standev))) #-2 sigma

#Virido
virido_standev <-sd(virido.log.2,na.rm=TRUE)
Virido_Plus_significant<- subset(id,virido.log.2>(mean(virido.log.2,na.rm=T)+abs(2*virido_standev))) #+2 sigma
Virido_Neg_significant <- subset(id,virido.log.2<(mean(virido.log.2,na.rm=T)-abs(2*virido_standev))) #-2 sigma

#J1074
J1074_standev <-sd(J1074.log.2,na.rm=TRUE)
J1074_Plus_significant<- subset(id,J1074.log.2>(mean(J1074.log.2,na.rm=T)+abs(2*J1074_standev))) #+2 sigma
J1074_Neg_significant <- subset(id,J1074.log.2<(mean(J1074.log.2,na.rm=T)-abs(2*J1074_standev))) #-2 sigma

#AA4 
AA4_standev <-sd(AA4.log.2,na.rm=TRUE)
AA4_Plus_significant<- subset(id,AA4.log.2>(mean(AA4.log.2,na.rm=T)+abs(2*AA4_standev))) #+2 sigma
AA4_Neg_significant <- subset(id,AA4.log.2<(mean(AA4.log.2,na.rm=T)-abs(2*AA4_standev))) #-2 sigma

standard_dev<- c(SPB74_standev, virido_standev,J1074_standev, AA4_standev)
names(standard_dev)<-c( 
                       "SPB74 Standard Deviation:", 
                       "Virido Standard Deviation:", 
                       "J1074 Standard Deviation: ", 
                       "AA4 Standard Deviation ")
Plus_list= list(SPB74_Plus_significant,Virido_Plus_significant,J1074_Plus_significant,AA4_Plus_significant)
lenB<-length(SPB74_Plus_significant)
lenC<-length(Virido_Plus_significant)
lenD<-length(J1074_Plus_significant)
lenE<- length(AA4_Plus_significant)
n12= length(Reduce(intersect,list(Plus_list[[1]],Plus_list[[2]])))
n13=length(Reduce(intersect,list(Plus_list[[1]],Plus_list[[3]])))
n14=length(Reduce(intersect,list(Plus_list[[1]],Plus_list[[4]])))
#n15=length(Reduce(intersect,list(Plus_list[[1]],Plus_list[[5]])))
n23=length(Reduce(intersect,list(Plus_list[[2]],Plus_list[[3]])))
n24=length(Reduce(intersect,list(Plus_list[[2]],Plus_list[[4]])))
#n25=length(Reduce(intersect,list(Plus_list[[2]],Plus_list[[5]])))
n34=length(Reduce(intersect,list(Plus_list[[3]],Plus_list[[4]])))
#n35=length(Reduce(intersect,list(Plus_list[[3]],Plus_list[[5]])))
#n45=length(Reduce(intersect,list(Plus_list[[4]],Plus_list[[5]])))
n123=length(Reduce(intersect,list(Plus_list[[1]],Plus_list[[2]], Plus_list[[3]])))
n124=length(Reduce(intersect,list(Plus_list[[1]],Plus_list[[2]], Plus_list[[4]])))
#n125=length(Reduce(intersect,list(Plus_list[[1]],Plus_list[[2]], Plus_list[[5]])))
n134= length(Reduce(intersect,list(Plus_list[[1]],Plus_list[[3]], Plus_list[[4]])))
#n135=length(Reduce(intersect,list(Plus_list[[1]],Plus_list[[3]], Plus_list[[5]])))
#n145=length(Reduce(intersect,list(Plus_list[[1]],Plus_list[[4]], Plus_list[[5]])))
n234=length(Reduce(intersect,list(Plus_list[[2]],Plus_list[[3]], Plus_list[[4]])))
#n235=length(Reduce(intersect,list(Plus_list[[2]],Plus_list[[3]], Plus_list[[5]])))
#n245=length(Reduce(intersect,list(Plus_list[[2]],Plus_list[[4]], Plus_list[[5]])))
#n345=length(Reduce(intersect,list(Plus_list[[3]],Plus_list[[4]], Plus_list[[5]])))
n1234=length(Reduce(intersect,list(Plus_list[[1]],Plus_list[[2]], Plus_list[[3]], Plus_list[[4]])))
#n1235=length(Reduce(intersect,list(Plus_list[[1]],Plus_list[[2]], Plus_list[[3]], Plus_list[[5]])))
#n1245=length(Reduce(intersect,list(Plus_list[[1]],Plus_list[[2]], Plus_list[[4]], Plus_list[[5]])))
#n1345= length(Reduce(intersect,list(Plus_list[[1]],Plus_list[[3]], Plus_list[[4]], Plus_list[[5]])))
#n2345=length(Reduce(intersect,list(Plus_list[[2]],Plus_list[[3]], Plus_list[[4]], Plus_list[[5]])))
#n12345=length(Reduce(intersect,list(Plus_list[[1]],Plus_list[[2]], Plus_list[[3]], Plus_list[[4]], Plus_list[[5]])))
money_shot_positive<-(Reduce(intersect,list(Plus_list[[1]],Plus_list[[2]], Plus_list[[3]], Plus_list[[4]])))

vplot<-draw.quad.venn(area1=lenB, 
                      area2=lenC, 
                      area3=lenD, 
                      area4=lenE, 
                      n12=n12, 
                      n13=n13, 
                      n14=n14, 
                      n23=n23, 
                      n24=n24,
                      n34=n34, 
                      n123=n123, 
                      n124=n124, 
                      n134=n134, 
                      n234=n234, 
                      n1234=n1234, 
                      category = c("SPB74", "Virido", "J1074", "AA4"), 
                      fill = c("orange", "red", "green", "blue"),
                      lty = "dashed",
                      cex = 2,
                      cat.cex = 2,
                      cat.col = c("orange", "red", "green", "blue"))



