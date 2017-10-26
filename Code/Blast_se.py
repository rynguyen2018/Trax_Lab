from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import numpy as np
###nucleotide align
num_genes=0
seq_array= []
with open("test_seq.fa") as gen_seq:
    for line in gen_seq:
        if ">" in line:
            num_genes+=1 
        else: 
            seq_array.append(line)
num_sur_genes=0
sur_seq_array=[]
with open("sur_sequence.fa") as sur_seq:
    for line in sur_seq:
        if ">" in line:
            num_sur_genes+=1 
        else:
            sur_seq_array.append(line)

print("Sequences ready for analysis")
sur_gene_number= 1
sur_score_dict={}
for sur_gene in sur_seq_array:
    sur_gene_name= "Sur Gene " + str(sur_gene_number)
    sur_gene_alignment_score=[]
    print(sur_gene_alignment_score)
    seq_gene_nums=[]
    for i in range(0, len(seq_array)):
        #print(sur_gene)
        alignment_score1 = pairwise2.align.localxx(sur_gene, seq_array[i], score_only=True)
        alignment_score2 = pairwise2.align.localxx(seq_array[i], sur_gene,score_only=True)
        alignment_score= max(alignment_score1, alignment_score2)
        sur_gene_alignment_score.append(alignment_score)
        seq_gene_nums.append(str(i))
        #input("...")
    print(sur_gene_alignment_score)
    print("Now finding best matches for ", sur_gene_name, "...")
    #print(list(zip(sur_gene_alignment_score,seq_gene_nums))[45])
    print(sorted(zip(sur_gene_alignment_score, seq_gene_nums), reverse=True))
    print("Onto next gene ....s")
    sur_gene_alignment_score=None
    input("...")
    sur_gene_number+=1


