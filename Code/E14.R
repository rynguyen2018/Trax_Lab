
gene_fold_data= read.csv('C:/Users/Ryan/Desktop/Traxler_Lab/Log2Data.csv', header =T)
# pick the gene fold data in Traxler Lab folder
# make sure it is a csv file! 



 #allows for extraction and separation of data by headers within text files 
attach(gene_fold_data)

#sort DNA based on length
sorted_DNA<- id  
E14cluster_data<- E14.log.2


#exclude DNA with cluster= N/A 
x<- sorted_DNA
y<-E14cluster_data
#better.y<-as.numeric(as.character(y))

plot(x,y, main= "E14 Gene Fold Expression as a Function of DNA Position", ylab = "Log_2 Gene Fold", xlab="DNA position(bp)")
points<- identify(x,y)
print(points)
print(y[points])
 

