
gene_fold_data= read.csv('C:/Users/Ryan/Desktop/Traxler_Lab/all_4_fold_venn_JGI_annot.csv', header =T)
# pick the gene fold data in Traxler Lab folder
# make sure it is a csv file! 



 #allows for extraction and separation of data by headers within text files 
attach(gene_fold_data)

#sort DNA based on length
sorted_DNA<- sort(DNA.Sequence.Length..bp.)
cluster_data<- venn.cat.1


#exlude DNA with cluster= N/A 
x<- subset(sorted_DNA, venn.cat != '#N/A') 
y<- subset(cluster_data, venn.cat!='#N/A')
better.y<-as.numeric(as.character(y))

plot(x,better.y, main= "Gene Fold Expression as a Function of DNA Position", ylab = "Log_2 Gene Fold", xlab="DNA position(bp)")

locator(n=200, type='p')
 

