from __future__ import division
from collections import Counter 
import sys 
import csv
from math import log
import heapq
from operator import itemgetter

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
			col2= [str(col2[x]) for col2 in matrix]
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
			ij= str(x)+","+str(y)	
			h_j= entropy_list[j]
			h_ij= ij_entropy[ij]
			mut_info[ij]= h_i+ h_j - h_ij
	return mut_info


def main(): 	
	file= sys.argv[1]

	with open(file) as csvfile: 
		reader=csv.reader(csvfile, delimiter= ",")
		header= next(csvfile)
		temp_data= []
		for row in reader: 
			row= ["0" if(value=="") else value for value in row]
			temp_data.append(row)

	#Turn data from string to numbers 
	data= []
	for row in temp_data: 
		row=[round(float(i)*2)/2 for i in row]
		data.append(row)
	temp_data=None

	#transpose data 
	data= map(list, zip(*data))


	entropy={}
	for col in range(0, len(data[0])):
		entropy[col]= findEntropy(findProb(column(data,col)))
	#print entropy
	#for val in entropy:
	#	print val 
	#	raw_input("Bing")
	vals= []
	for key, value in entropy.iteritems(): 
		vals.append((key,value))
	
	vals= sorted(vals,  key=itemgetter(1))

	#for key,value in heapq.nsmallest(10, entropy.items(), key=itemgetter(1)):
	#	print key,value


	ij_prob= pairColumnProb(data)
	ij_entropy={}
	for pair, value in ij_prob.iteritems(): 
		ij_entropy[pair]= findEntropy[value]

	mutual_information= mutInfo(data, entropy_list, ij_entropy)

	for key,value in heapq.nlargest(50, mutual_information.items(), key=itemgetter(1)):
		print key
if __name__ == '__main__':
	main()