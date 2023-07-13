import argparse
import json
from Bio import SeqIO, Phylo
import pandas as pd
import matplotlib.pyplot as plt
from collections import OrderedDict
import numpy as np


def cumulative(dictionary):
    """
    Returns cumulative values for each element of an ordered dictionary
    Input:  
        dictionary with cumulative values for each ordered key 
    """
    od = OrderedDict(sorted(dictionary.items()))
    x = list(od.values())
    res = np.cumsum(x)
    cumul= dict()
    for i, j in zip(od.keys(), res): 
        cumul[i] = j
    return(cumul)


def finding_duplications(nwktree, reconstructeddupl, which):

    """returns earliest sequence of duplication 1 duplication 2 or preduplication based on phylogenetic nwk tree and reconstructed dupl sequences"""

    tree_file  = Phylo.read(nwktree, "newick")
    tree_file.root_at_midpoint()
    tree_file.find_clades()
    duplication_file = SeqIO.parse(reconstructeddupl, "fasta")
    seq_dict = dict()
    for record in duplication_file: seq_dict[record.id] = record.seq

    for branch in tree_file.get_nonterminals(order='preorder'):
        if pd.isna(branch.name) == False:
            if '-'*int(72) not in seq_dict[branch.name]:
                first_dupl_sequence = seq_dict[branch.name]
                break  
    for branch in tree_file.get_nonterminals(order='preorder'):
        if pd.isna(branch.name) == False:
                predupl_sequence = seq_dict[branch.name][72:][1:-2]
                break  
    post_dupl_1 = first_dupl_sequence[:72][1:-2]
    post_dupl_2 = first_dupl_sequence[72:][1:-2]

    if which == "1":return(post_dupl_1)
    if which == "2": return(post_dupl_2)
    if which == "pre": return(predupl_sequence)



def dictionary_of_mutations(duplication, mutation_matrix, type_):
    sum_of_rows = mutation_matrix.sum(axis = 1)
    translations = {'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'], 'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'], 'C': ['TGT', 'TGC'], 'W': ['TGG'], 'E': ['GAA', 'GAG'], 'D': ['GAT', 'GAC'], 'P': ['CCT', 'CCC', 'CCA', 'CCG'], 'V': ['GTT', 'GTC', 'GTA', 'GTG'], 'N': ['AAT', 'AAC'], 'M': ['ATG'], 'K': ['AAA', 'AAG'], 'Y': ['TAT', 'TAC'], 'I': ['ATT', 'ATC', 'ATA'], 'Q': ['CAA', 'CAG'], 'F': ['TTT', 'TTC'], 'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'], 'T': ['ACT', 'ACC', 'ACA', 'ACG'], '*': ['TAA', 'TAG', 'TGA'], 'A': ['GCT', 'GCC', 'GCA', 'GCG'], 'G': ['GGT', 'GGC', 'GGA', 'GGG'], 'H': ['CAT', 'CAC']}
    each_position = dict()
    if type_ == "point_mut":
        for i, char in enumerate(duplication):
            if i%3 == 0:
                codon = duplication[i:i+3]
                for key, entry in translations.items():
                    if codon in entry:
                            if len(entry) == 4: each_position[i+2] = sum_of_rows[codon[-1]]
                            elif len(entry) == 2 or len(entry) == 3:
                                sum_of_muts = 0
                                for ele in entry: sum_of_muts += mutation_matrix.at[codon[-1], ele[-1]]
                                each_position[i+2] = sum_of_muts
                            elif len(entry) == 6:
                                each_position[i+2] = sum_of_rows[codon[-1]]
                                for ele in entry:
                                    if ele[0] != codon[0]:
                                        sum_ = float(mutation_matrix.at[codon[0], ele[0]])
                                        each_position[i] = sum_
    else:
        for i, char in enumerate(duplication):
            if i%3 == 0 and i !=0:
                codon_plus_prev = str(duplication[i-1:i+3])
                codon = str(duplication[i:i+3])
                for key, entry in translations.items():
                    if codon in entry:
                            if len(entry) == 4: each_position[i+2] = sum_of_rows[codon[1:]]
                            if len(entry) == 2 or len(entry) == 3:
                                sum_of_muts = 0
                                for ele in entry: sum_of_muts += mutation_matrix.at[codon[1:], ele[-1]]
                                each_position[i+2] = sum_of_muts
                            if len(entry) == 6:
                                each_position[i+2] = sum_of_rows[codon[1:]]
                                for ele in entry:
                                    if ele[0] != codon[0]:
                                        sum_ = float(mutation_matrix.at[codon_plus_prev[:2], ele[0]])
                                        each_position[i] = sum_
    cumulative_= cumulative(each_position)
    return(cumulative_)


if __name__=="__main__":
    parser = argparse.ArgumentParser(
        description="Construct synynymous mutation matrix",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument('--matrix', required=True, type=str, help="csv file with scaled gene mutations")
    parser.add_argument('--tree', required=True, type=str, help="Tree nwk")
    parser.add_argument('--ref', type=str, help="genbank reference file")
    parser.add_argument('--output', type=str, help="graph png file")
    parser.add_argument('--type', required=True, help="context of the mutation")
    parser.add_argument('--duplicationseq', required=True, help="fasta file with reconstructed sequences")

    args = parser.parse_args()

    mut_matrix = pd.read_csv(args.matrix)
    mut_matrix = mut_matrix.set_index('Unnamed: 0')

    post_1 = finding_duplications(args.tree, args.duplicationseq, "1")    
    post_2 = finding_duplications(args.tree, args.duplicationseq, "2")   
    pre = finding_duplications(args.tree, args.duplicationseq, "pre")    

    cumulative_1 = dictionary_of_mutations(post_1, mut_matrix, args.type)
    cumulative_2 = dictionary_of_mutations(post_2, mut_matrix, args.type)
    cumulative_pre = dictionary_of_mutations(pre, mut_matrix, args.type)

    plt.step(cumulative_1.keys(), cumulative_1.values(), label= f'1st copy postduplication', where='post')
    plt.step(cumulative_2.keys(), cumulative_2.values(), label=f'2nd copy postduplication', where='post' )
    plt.step(cumulative_pre.keys(), cumulative_pre.values(), label= f'preduplication', where='post' )
    plt.legend(loc='lower right')
    plt.xlabel("gene location")
    plt.ylabel("cumulative sum of mutations")
    plt.savefig(args.output)
