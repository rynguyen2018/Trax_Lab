
#library(gplots)
require(VennDiagram)


gene_fold_data= read.csv('C:/Users/Ryan/Desktop/Traxler_Lab/Code/Log2Data.csv', header =T)
# pick the gene fold data in Traxler Lab folder
# make sure it is a csv file! 



 #allows for extraction and separation of data by headers within text files 
attach(gene_fold_data)


#Standard Deviation of all of the samples 
#E14
E14_standev <-sd(E14.log.2,na.rm=TRUE)
E14_Plus_significant<- subset(id,E14.log.2>(mean(E14.log.2,na.rm=T)+abs(2*E14_standev))) #+2 sigma
E14_Neg_significant <- subset(id,E14.log.2<(mean(E14.log.2,na.rm=T)-abs(1.8*E14_standev))) #-2 sigma

#SPB74
SPB74_standev <-sd(SPB74.log.2,na.rm=TRUE)
SPB74_Plus_significant<- subset(id,SPB74.log.2>(mean(SPB74.log.2,na.rm=T)+abs(2*SPB74_standev))) #+2 sigma
SPB74_Neg_significant <- subset(id,SPB74.log.2<(mean(SPB74.log.2,na.rm=T)-abs(1.8*SPB74_standev))) #-2 sigma

#Virido
virido_standev <-sd(virido.log.2,na.rm=TRUE)
Virido_Plus_significant<- subset(id,virido.log.2>(mean(virido.log.2,na.rm=T)+abs(2*virido_standev))) #+2 sigma
Virido_Neg_significant <- subset(id,virido.log.2<(mean(virido.log.2,na.rm=T)-abs(1.8*virido_standev))) #-2 sigma

#J1074
J1074_standev <-sd(J1074.log.2,na.rm=TRUE)
J1074_Plus_significant<- subset(id,J1074.log.2>(mean(J1074.log.2,na.rm=T)+abs(2*J1074_standev))) #+2 sigma
J1074_Neg_significant <- subset(id,J1074.log.2<(mean(J1074.log.2,na.rm=T)-abs(1.8*J1074_standev))) #-2 sigma

#AA4 
AA4_standev <-sd(AA4.log.2,na.rm=TRUE)
AA4_Plus_significant<- subset(id,AA4.log.2>(mean(AA4.log.2,na.rm=T)+abs(2*AA4_standev))) #+2 sigma
AA4_Neg_significant <- subset(id,AA4.log.2<(mean(AA4.log.2,na.rm=T)-abs(1.8*AA4_standev))) #-2 sigma

standard_dev<- c(E14_standev, SPB74_standev, virido_standev,J1074_standev, AA4_standev)
names(standard_dev)<-c("E14 Standard Deviation:", 
                       "SPB74 Standard Deviation:", 
                       "Virido Standard Deviation:", 
                       "J1074 Standard Deviation: ", 
                       "AA4 Standard Deviation ")
#print(standard_dev)
Neg_list= list(SPB74_Neg_significant,Virido_Neg_significant,J1074_Neg_significant,AA4_Neg_significant)
lenA<-length(E14_Neg_significant)
lenB<-length(SPB74_Neg_significant)
lenC<-length(Virido_Neg_significant)
lenD<-length(J1074_Neg_significant)
lenE<- length(AA4_Neg_significant)
n12= length(Reduce(intersect,list(Neg_list[[1]],Neg_list[[2]])))
n13=length(Reduce(intersect,list(Neg_list[[1]],Neg_list[[3]])))
n14=length(Reduce(intersect,list(Neg_list[[1]],Neg_list[[4]])))
#n15=length(Reduce(intersect,list(Neg_list[[1]],Neg_list[[5]])))
n23=length(Reduce(intersect,list(Neg_list[[2]],Neg_list[[3]])))
n24=length(Reduce(intersect,list(Neg_list[[2]],Neg_list[[4]])))
#n25=length(Reduce(intersect,list(Neg_list[[2]],Neg_list[[5]])))
n34=length(Reduce(intersect,list(Neg_list[[3]],Neg_list[[4]])))
#n35=length(Reduce(intersect,list(Neg_list[[3]],Neg_list[[5]])))
#n45=length(Reduce(intersect,list(Neg_list[[4]],Neg_list[[5]])))
n123=length(Reduce(intersect,list(Neg_list[[1]],Neg_list[[2]], Neg_list[[3]])))
n124=length(Reduce(intersect,list(Neg_list[[1]],Neg_list[[2]], Neg_list[[4]])))
#n125=length(Reduce(intersect,list(Neg_list[[1]],Neg_list[[2]], Neg_list[[5]])))
n134= length(Reduce(intersect,list(Neg_list[[1]],Neg_list[[3]], Neg_list[[4]])))
#n135=length(Reduce(intersect,list(Neg_list[[1]],Neg_list[[3]], Neg_list[[5]])))
#n145=length(Reduce(intersect,list(Neg_list[[1]],Neg_list[[4]], Neg_list[[5]])))
n234=length(Reduce(intersect,list(Neg_list[[2]],Neg_list[[3]], Neg_list[[4]])))
#n235=length(Reduce(intersect,list(Neg_list[[2]],Neg_list[[3]], Neg_list[[5]])))
#n245=length(Reduce(intersect,list(Neg_list[[2]],Neg_list[[4]], Neg_list[[5]])))
#n345=length(Reduce(intersect,list(Neg_list[[3]],Neg_list[[4]], Neg_list[[5]])))
n1234=length(Reduce(intersect,list(Neg_list[[1]],Neg_list[[2]], Neg_list[[3]], Neg_list[[4]])))
#n1235=length(Reduce(intersect,list(Neg_list[[1]],Neg_list[[2]], Neg_list[[3]], Neg_list[[5]])))
#n1245=length(Reduce(intersect,list(Neg_list[[1]],Neg_list[[2]], Neg_list[[4]], Neg_list[[5]])))
#n1345= length(Reduce(intersect,list(Neg_list[[1]],Neg_list[[3]], Neg_list[[4]], Neg_list[[5]])))
#n2345=length(Reduce(intersect,list(Neg_list[[2]],Neg_list[[3]], Neg_list[[4]], Neg_list[[5]])))
#n12345=length(Reduce(intersect,list(Neg_list[[1]],Neg_list[[2]], Neg_list[[3]], Neg_list[[4]], Neg_list[[5]])))


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



# vplot<-draw.quintuple.venn(area1=lenA, 
#                            area2=lenB, 
#                            area3=lenC, 
#                            area4=lenD, 
#                            area5=lenE, 
#                            n12=n12, 
#                            n13=n13, 
#                            n14=n14, 
#                            n15=n15,
#                            n23=n23, 
#                            n24=n24, 
#                            n25=n25, 
#                            n34=n34, 
#                            n35=n35, 
#                            n45=n45, 
#                            n123=n123, 
#                            n124=n124, 
#                            n125=n125, 
#                            n134=n134,
#                            n135=n135, 
#                            n145=n145, 
#                            n234=n234, 
#                            n235=n235, 
#                            n245=n245, 
#                            n345=n345, 
#                            n1234=n1234, 
#                            n1235=n1235,
#                            n1245=n1245, 
#                            n1345=n1345, 
#                            n2345=n2345, 
#                            n12345=n12345,
#                            category = c("A", "B", "C", "D", "E"),
#                            fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
#                            cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
#                            cat.cex = 2,
#                            margin = 0.05,
#                            cex = c(1.5, 1.5, 1.5, 1.5, 1.5, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8,
#                                    1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 1, 1, 1, 1, 1.5),
#                            ind = TRUE)


#v.table<- venn(Neg_list,show.plot= FALSE)
#print(v.table)


