from __future__ import division
from collections import Counter 
import sys 
import csv
"""
Task 1) Read in RNA seq data in csv format  
"""
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
	row=[round(float(i)) for i in row]
	data.append(row)
temp_data=None

#transpose data 
data= map(list, zip(*data))


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
	return col_length
	


