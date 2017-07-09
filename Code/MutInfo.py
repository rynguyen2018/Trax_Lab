from __future__ import division
from collections import Counter 
import sys 
import csv
from math import log
import heapq
from operator import itemgetter
import argparse
import matplotlib.pyplot as plt
from collections import Counter

"""
Task 1) Read in RNA seq data in csv format  
"""


#get single columns 
def column(matrix, i): 
	return [row[i] for row in matrix]

"""
Task 2) Calculate probability of finding particular regulation in a column 

Note: calculating the Shannon entropy associated with a column should possibly be generalized.
In other words, there should be a method created that lumps the probablity of finding a 2 or a 3 
into one probability as to generalize up regulation. The same should be true for down regulatory genes.  
"""
#calculates probability of particular column
def findProb(col):
	col_length= len(col)
	prob={}
	# counts instances of log 2 regulation and calculates regulatory instance probability  
	prob_temp= Counter(col)
	for key, value in prob_temp.iteritems(): 
		prob[key]= value/col_length
	return prob

#paired column calculation:
#get paired columns 
def pairColumnProb(matrix):
	pairColProb={}
	seqList={}
	for x in range(0, len(matrix[0])): 
		col1= [str(col[x]) for col in matrix] 
		for y in range(x+1, len(matrix[0])): 
			col2= [str(col2[y]) for col2 in matrix]
			seq= map("".join,zip(col1,col2)) 
			col_cat= str(x)+","+str(y)
			seqList[col_cat]= seq
	for position, reg_pair in seqList.iteritems(): 
		pairColProb[position]= findProb(reg_pair)		
	return pairColProb

""" 
Task 3) Calculate the Shannon entropy of each column 
"""
def findEntropy(prob_col):
	entropy_value=0 
	for key, value in prob_col.iteritems(): 
		if value == "null": 
			return -9999999999
		else: 
			entropy_value+= (-1*prob_col[key]*log(prob_col[key],2))	

	return round(entropy_value,6)


"""
Task 4) Calculate the mutual information 
"""

def mutInfo(matrix, entropy_list, ij_entropy):
	mut_info={}
	for i in range(0, len(matrix[0])): 
		h_i= entropy_list[i]  
		for j in range(i+1, len(matrix[0])): 
			ij= str(i)+","+str(j)	
			h_j= entropy_list[j]
			h_ij= ij_entropy[ij]
			mut_info[ij]= h_i+ h_j - h_ij
	return mut_info


def main(): 	
	parser = argparse.ArgumentParser(description='Process some RNAseq data.')
	parser.add_argument('-i','--input', help= "csv file of RNA seq data")
	parser.add_argument('-max', help= "max gene with mutual information of interest")
	parser.add_argument('-min', help= 'min gene with mutual information of interest')
	parser.add_argument('-a','--annotate', help= "compares high scorers of mutual information with annotated gene base")

	args = parser.parse_args()
	file= args.input
	gene_interest_max= int(args.max) 
	gene_interest_min= int(args.min)
	annotation_file= args.annotate 


	with open(file) as csvfile: 
		reader=csv.reader(csvfile, delimiter=",")
		header= next(csvfile)
		temp_data= []
		for row in reader: 
			row= ["0" if(value=="") else value for value in row]
			temp_data.append(row)
	#Turn data from string to numbers 
	data= []
	for row in temp_data:
		if (row[x]== "null" for x in row):
			data.append(row)
		else: 
			row=[round(float(i)*2)/2 for i in row]
			data.append(row)
	temp_data=None

	#transpose data 
	data= map(list, zip(*data))

	#find entropy
	entropy={}
	for col in range(0, len(data[0])):
		entropy[col]= findEntropy(findProb(column(data,col)))
	print "Entropy values found!"
	# vals= []
	# for key, value in entropy.iteritems(): 
	# 	vals.append((key,value))
	
	#vals= sorted(vals,  key=itemgetter(1))

	#find pairwise column prob
	ij_prob= pairColumnProb(data)

	ij_entropy={}
	for pair, value in ij_prob.iteritems(): 
		ij_entropy[pair]= findEntropy(value)

	#get mutual information	
	mutual_information= mutInfo(data, entropy, ij_entropy)
	# print "Mutual Information: "
	# for key,value in heapq.nlargest(200, mutual_information.items(), key=itemgetter(1)):
	# 	print key
	#filter mutual information based on RNA seq up/down regulation 
		#find the conservons!!!! 
	mut_info_vals=[]
	for key, value in mutual_information.iteritems(): 
		#print type(int(key.split(",")[0]))
		if (int(key.split(",")[0]) >= gene_interest_min and int(key.split(",")[0])<= gene_interest_max): 
			mut_info_vals.append((int(key.split(",")[1]), value))
		if (int(key.split(",")[1])>= gene_interest_min and int(key.split(",")[1])<= gene_interest_max): 
			mut_info_vals.append((int(key.split(",")[0]), value))
			
	print "Mutual Information found! "
	mut_info_vals= sorted(mut_info_vals,  key=lambda x: x[1])
	#print mut_info_vals
	top_scores=mut_info_vals[-100:]
	top_genes= [int(i[0]) for i in top_scores]
	corresponding_genes_of_interest_list= [int(i[1]) for i in top_scores]
	#print len(top_genes)
	#raw_input("Press Enter to continue")
	# plt.hist(top_genes, bins=len(top_genes)*2)
	# plt.title("Mutual Information Histogram for Genes of Interest")
	# plt.xlabel("Value")
	# plt.ylabel("Frequency")
	# plt.show()

	total_gene_list= []
	annotation_list=[]
	#annotate that ish bruhhhhhh
	with open(annotation_file) as csvfile: 
		reader=csv.reader(csvfile, delimiter= ",")
		header= next(csvfile)
		for row in reader:
			total_gene_list.append(row[7])
			annotation_list.append(row[11])
	for gene in range(0,len(total_gene_list)):
		temp_gene= total_gene_list[gene].strip("SCO")
		temp_gene= temp_gene.lstrip("0") 
		total_gene_list[gene]=int(temp_gene) 
	result_annotation= {}
	f = open(str(gene_interest_min) + "," + str(gene_interest_max) +"-" 'mut_annotation.txt', 'w')
	annotation_list= (zip(total_gene_list, annotation_list))
	#print annotation_list
	for i in range(0,len(annotation_list)):
		if annotation_list[i][0] in top_genes:
			f.write(str(annotation_list[i][0]) + "\t" + str(annotation_list[i][1]) +"\t"+ corresponding_genes_of_interest_list[i] + "\t"+top_scores[i] + "\n")
	#print result_annotation
	f.close()

	print "annotation done! Now printing statistics"
	plt.hist(top_genes, bins=len(top_genes)*2)
	plt.title("Mutual Information Histogram for Genes of Interest")
	plt.xlabel("Value")
	plt.ylabel("Frequency")
	plt.show()


	#top_scores= sorted(top_scores, key=lambda x: x[0])
	#d = {x:top_scores.count(x) for x in top_scores}
	#od =sorted(d.items(), key=lambda x: x[1])
	#print Counter(elem[0] for elem in top_scores)

 #   	big_regulation_list=[]
	# for col in range(0, len(data[0])):
	# 	temp_col=column(data,col)
	# 	average=abs(sum(list(map(float, temp_col)))/len(temp_col))
	# 	if average> 1.5: 
	# 		big_regulation_list.append(col)
	# new_mut_info= {k:v for k,v in mutual_information.iteritems() if int(k.split(",")[0]) in big_regulation_list or int(k.split(",")[1])}
			
	# for key,value in heapq.nlargest(200, new_mut_info.items(), key=itemgetter(1)):
	# 	print key

	#find the conservons!!!! 
	# for key, value in new_mut_info.iteritems(): 
	# 	if (int(key.split(",")[0])<= 1166) or (int(key.split(",")[1])<= 1166): 
	# 		print key 
if __name__ == '__main__':
	main()
