import re

file_num=0 
with open('ALL_GNPS.mgf') as f: 
	for line in f: 
		if "BEGIN IONS" in line: 
			file_num+=1 
			title= "GNPS_Output_" + str(file_num) +".txt" 
			f= open(title,"w")
		if "END IONS" not in line: 
			f.write(line)
		
f.close()
							 
#(file_num)