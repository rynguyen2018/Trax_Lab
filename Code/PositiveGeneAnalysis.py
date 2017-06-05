import csv 

my_file=open('Log2Data.csv')
my_data= csv.reader(my_file)
next(my_data,None)


E14_data= [(row[2]) for row in my_data]


for row in E14_data:
	new = [float(x) for x in row if x != '']


E14logdata=[] 

my_file.close()