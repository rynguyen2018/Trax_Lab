import antismash_sim 
from Bio import pairwise2 as pw2
import csv
from Bio import SeqIO

gbk_filename = "viola_hit.fasta"
faa_filename = "viola_smash.fasta"



best_score_list, best_name_list = antismash_sim.compareSequences(gbk_filename, faa_filename)

print(best_name_list)

#with open("viola_smash_hit.csv", 'wb') as myfile:
#	wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
#	wr.writerow(best_name_list)
#	wr.writerow(best_score_list)
#for seq_record in SeqIO.parse(input_handle, "genbank") :
    #print(seq_record.id)
    #print("Dealing with GenBank record %s" % seq_record.id)
    #output_handle.write(">%s %s\n%s\n" % (
    #      seq_record.id,
    #      seq_record.description,
    #      seq_record.seq))

#output_handle.close()
#input_handle.close()