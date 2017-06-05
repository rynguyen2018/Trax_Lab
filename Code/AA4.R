
gene_fold_data= read.csv('C:/Users/Ryan/Desktop/Traxler_Lab/Code/Log2Data.csv', header =T)
# pick the gene fold data in Traxler Lab folder
# make sure it is a csv file! 



 #allows for extraction and separation of data by headers within text files 
attach(gene_fold_data)

#sort DNA based on length
sorted_DNA<- id  
AA4cluster_data<- AA4.log.2


#exclude DNA with cluster= N/A 
x<- sorted_DNA
y<-AA4cluster_data
#better.y<-as.numeric(as.character(y))

plot(x,y, main= "AA4 Gene Fold Data", ylab = "Log_2 Gene Fold", xlab="Strep Gene")
abline(h=mean(gene_fold_data$AA4.log.2, na.rm=T),col="red") 
xcoord<- identify(x)
points<-c(xcoord,y[xcoord])
print(points)
