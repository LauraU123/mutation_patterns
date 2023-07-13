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




sum_of_rows = scaled_and_normalized.sum(axis = 1)


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


translations = {'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'], 'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'], 'C': ['TGT', 'TGC'], 'W': ['TGG'], 'E': ['GAA', 'GAG'], 'D': ['GAT', 'GAC'], 'P': ['CCT', 'CCC', 'CCA', 'CCG'], 'V': ['GTT', 'GTC', 'GTA', 'GTG'], 'N': ['AAT', 'AAC'], 'M': ['ATG'], 'K': ['AAA', 'AAG'], 'Y': ['TAT', 'TAC'], 'I': ['ATT', 'ATC', 'ATA'], 'Q': ['CAA', 'CAG'], 'F': ['TTT', 'TTC'], 'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'], 'T': ['ACT', 'ACC', 'ACA', 'ACG'], '*': ['TAA', 'TAG', 'TGA'], 'A': ['GCT', 'GCC', 'GCA', 'GCG'], 'G': ['GGT', 'GGC', 'GGA', 'GGG'], 'H': ['CAT', 'CAC']}


def dictionary_of_mutations(duplication, mutation_matrix):
    translations = {'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'], 'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'], 'C': ['TGT', 'TGC'], 'W': ['TGG'], 'E': ['GAA', 'GAG'], 'D': ['GAT', 'GAC'], 'P': ['CCT', 'CCC', 'CCA', 'CCG'], 'V': ['GTT', 'GTC', 'GTA', 'GTG'], 'N': ['AAT', 'AAC'], 'M': ['ATG'], 'K': ['AAA', 'AAG'], 'Y': ['TAT', 'TAC'], 'I': ['ATT', 'ATC', 'ATA'], 'Q': ['CAA', 'CAG'], 'F': ['TTT', 'TTC'], 'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'], 'T': ['ACT', 'ACC', 'ACA', 'ACG'], '*': ['TAA', 'TAG', 'TGA'], 'A': ['GCT', 'GCC', 'GCA', 'GCG'], 'G': ['GGT', 'GGC', 'GGA', 'GGG'], 'H': ['CAT', 'CAC']}
    each_position = dict()
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
    cumulative_= cumulative(each_position)
    return(cumulative_)


plt.step(cumulative_1.keys(), cumulative_1.values(), label= f'1st copy postduplication', where='post')
plt.step(cumulative_2.keys(), cumulative_2.values(), label=f'2nd copy postduplication', where='post' )
plt.step(cumulative_pre.keys(), cumulative_pre.values(), label= f'preduplication', where='post' )
plt.legend(loc='lower right')
plt.suptitle('Expected Synonymous Mutations in RSV-A')
plt.xlabel("gene location")
plt.ylabel("cumulative sum of mutations")
plt.savefig("expected_synonymous.png")


