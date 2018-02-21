import argparse
import json

"""
	This script is designed to take the output files of antismash_sim.py and output them into a json file to be processed
	for D3 visualization 
"""

def storeFastaToArray(x):
	x_array = []
	with open(x) as file: 
		for line in file: 
			if ">" in line:
				x_array.append(line[1:])
	return x_array

def storeTxtToArray(y):
	
	f = open(y)
	message = f.read()
	#print(message)
	output_genes = message.split(",")
	for i in range(0, len(output_genes)): 
		output_genes[i] = output_genes[i].strip("'")
		output_genes[i] = output_genes[i].strip("[")
		output_genes[i] = output_genes[i].strip("]")
		output_genes[i] = output_genes[i].strip(",")
		output_genes[i] = output_genes[i].strip('"')
		output_genes[i] = output_genes[i].strip()
		output_genes[i] = output_genes[i].strip("'")
	f.close()

	return output_genes 



def main(): 
	parser = argparse.ArgumentParser(description='JSON maker')
	parser.add_argument('-x', '--input', help= "input fastaa ")
	parser.add_argument('-y', '--output', help= "output genes")

	args = parser.parse_args()
	x = args.input
	y = args.output 

	x_array = storeFastaToArray(x)
	y_array = storeTxtToArray(y)




	output = "var data =[" 

	dict_array = [] 
	for j in range(0, len(x_array)): 
		temp_dict_forward = {}
		temp_dict_reverse = {}
		temp_array_forward = []
		temp_array_reverse = []
		temp_dict_forward["name"] = "flare.c3ba." + str(x_array[j]).strip("\n")
		temp_array_reverse.append(temp_dict_forward["name"])
		temp_dict_reverse["name"] = "flare.viola." + str(y_array[j]).strip("\n") 
		temp_array_forward.append(temp_dict_reverse["name"])
		
		temp_dict_forward["size"] = 300
		temp_dict_reverse["size"] = 3000
		
		temp_dict_forward["imports"] = temp_array_forward
		temp_dict_reverse["imports"] = temp_array_reverse

		dict_array.append(temp_dict_forward)
		dict_array.append(temp_dict_reverse)
	with open("output.json", "w") as outfile: 
		json.dump(dict_array,outfile)
	#print(json.dumps(dict_array, ensure_ascii=False))


	# for j in range(0, len(x_array)): 
	# 	input_node = "{'name': genes.c3ba." + str(x_array[j]) + "',"
	# 	output_node = "'imports': ['genes.viola." + str(y_array[j]) + "']}"  
	# 	if j != len(x_array) -1: 
	# 		ender = ","
	# 	else: 
	# 		ender = ";"
	# 	output += input_node + output_node + ender

	# d = json.loads(output[0])
	#print(d['names'])
if __name__ == '__main__':
	main()