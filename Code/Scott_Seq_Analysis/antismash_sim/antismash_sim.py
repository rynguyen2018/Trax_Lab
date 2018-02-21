import sys
import os 
import argparse
#from Bio import pairwise2 as pw2
import csv
#from Bio import SeqIO

""" 
	1) Reads in GenBank file of "surugamide" clusters from C3BA antismash 
	2) Blast C3BA genes and take first hit protein sequence from each protein sequence in cluster 
	3) Compare hits to each protein in each cluster in viola antismash result
	4) When >90% similarity, store gene into array as well as similarity score 
	5) Output Top 3 candidates 
	6) Compare candidate similarity to protein sequence
	7) Store normalized similarity into some sort of data structure that allows for pretty data viz haha 

"""



def renameDuplicates(fasta_file):
	##renames dulicate seqIDs in fasta file because antismash couldn't figure out how to do that themselves... selfish pricks
	from shutil import copyfile
	from os import remove

	temp_file = fasta_file +"_new"
	name_dict = {}
	with open(fasta_file) as old_file, open(temp_file, "w") as new_file:
		for line in old_file:
			if ">" in line:
				if line in name_dict.keys():
					name_dict[line] += 1
					line = line.rstrip() + "." + str(name_dict[line]) + "\n"
				else:
					name_dict[line] = 1
			new_file.write(line)
	copyfile(temp_file,fasta_file)
	remove(temp_file)


def GenbankToFasta(gen_file, out_file):
	print("Now Converting Genbank to Fasta")
	CBA_fasta = open(out_file, "w") 
	prot_seq = ""
	prot_tag = ""
	with open(gen_file) as file: 
		for line in file: 
			if "/old_locus_tag" in line:
				next(file)
			if "/locus_tag" in line: 
				prot_tag = line
				prot_tag.strip()
				prot_tag = prot_tag.replace('"', "")
				prot_tag = prot_tag.split("=")[1]
				prot_tag = ">" + prot_tag
				CBA_fasta.write(prot_tag)
				
			if "/translation" in line:
				prot_seq += line.strip()
				for line in file:
					if "CDS" not in line and "aSDomain" not in line and "ORIGIN" not in line and "/" not in line and "gene" not in line and "rRNA" not in line and "tRNA" not in line:
						prot_seq += line.strip()
					elif "ORIGIN" in line: 
						prot_seq.strip()
						prot_seq = prot_seq.replace('"',"")
						prot_seq = prot_seq.split("=")[1]
						CBA_fasta.write(prot_seq)
						break
					else: #"CDS" in line:
						prot_seq.strip()
						prot_seq = prot_seq.replace('"',"")
						prot_seq = prot_seq.split("=")[1]
						CBA_fasta.write(prot_seq + "\n")
						break

	CBA_fasta.close()
def BlastP(fasta_file, output_file= "blast.csv", ref_db_file = "viola_database.fasta", output_fasta_file = "viola_hit.fasta"):
	print("Making Blast database")
	
	string = "blastp -query " + fasta_file + " -out " + output_file + " -db "+ ref_db_file+ " -max_target_seqs 1 -outfmt 10"
	os.system(string)
	#print(string)
	#sys.exit()
	#seq_fasta = open("viola_hit.fasta", "w")
	outfile = open(output_fasta_file, "w") 

	if output_fasta_file != "BOOTYBOOTY":
		with open(output_file) as query_file: 
			reader = csv.reader(query_file)
			for row in reader: 
				seq_id = row[1]
				ref_file = open(ref_db_file) 
				for line in ref_file:
					if seq_id in line:
						outfile.write(line)
						seq = ""
						for line in ref_file:
							if line != '\n':
				 				seq += line
							else:
								outfile.write(seq)
								break
				ref_file.close()
		outfile.close()
	print("BlastP finished")

# def compareSequences(query_seq, compare_seq): 
# 	print("Now comparing sequences")
# 	best_score_list = []
# 	best_name_list = []
# 	for record in SeqIO.parse(query_seq, "fasta"): 
# 		best_score = -99999999
# 		best_name = ""
# 		first_seq = record.seq
# 		#print(record.id)
# 		for record2 in SeqIO.parse(compare_seq, "fasta"):
# 			second_seq = record2.seq
# 			global_align = pw2.align.localxx(first_seq, second_seq, 2, -1, -2, score_only = True)
# 			#print(len(first_seq), len(second_seq))
# 			seq_length = min(len(first_seq), len(second_seq))
# 			if seq_length !=0 :
# 				if type(global_align) == list: 
# 					global_align = 0
# 				temp_score = global_align/seq_length
# 				if temp_score > best_score: 
# 					best_score = temp_score
# 					best_name = record2.id
# 				if temp_score >= 0.90:
# 					break
# 		best_score_list.append(best_score)
# 		best_name_list.append(best_name)
# 	return best_score_list, best_name_list

def main(): 
	parser = argparse.ArgumentParser(description='Process some protein data from antismash. swag')
	parser.add_argument('-i','--input', help = "genbank file with AntiSmash genes from organism of interest")
	parser.add_argument('-o','--output', help = "name of fasta file to make")

	parser.add_argument('--smash_file', help = "antismash genbank file to compare with")
	parser.add_argument('--smash_output', help = "name of comparator antismash fastaa file to make")

	args = parser.parse_args()
	gen_file = args.input
	output_file = args.output
	compare_input_file = args.smash_file
	compare_output_file = args.smash_output

	print("GO")
	#print(output_file)
	if os.path.isfile(output_file) == False:
		GenbankToFasta(gen_file, output_file)
		BlastP(output_file)

	#GenbankToFasta(compare_input_file, compare_output_file)
	

	if os.path.isfile(compare_output_file) == False:
		GenbankToFasta(compare_input_file, compare_output_file)
		renameDuplicates(compare_output_file)
		os.system("makeblastdb -in " + compare_output_file + " -parse_seqids -dbtype prot")
		BlastP("viola_hit.fasta", "final_blast.csv", ref_db_file = "viola_smash.fasta", output_fasta_file = "BOOTYBOOTY")
	# best_score_list, best_name_list = compareSequences(output_file,compare_output_file)

	# print(best_name_list)

if __name__ == '__main__':
	main()


