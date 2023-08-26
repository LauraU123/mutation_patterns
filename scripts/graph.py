import argparse
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

def finding_duplications(nwktree, reconstructeddupl, subtype, which):

    """returns earliest sequence of duplication 1 duplication 2 or preduplication based on phylogenetic nwk tree and reconstructed dupl sequences"""

    tree_file  = Phylo.read(nwktree, "newick")
    tree_file.root_at_midpoint()
    tree_file.find_clades()
    duplication_file = SeqIO.parse(reconstructeddupl, "fasta")
    seq_dict = dict()
    for record in duplication_file: seq_dict[record.id] = record.seq

    for branch in tree_file.get_nonterminals(order='preorder'):
        if pd.isna(branch.name) == False:
            if '-' not in seq_dict[branch.name]:
                first_dupl_sequence = seq_dict[branch.name]
                #print(first_dupl_sequence, branch.name)
                break  
    for branch in tree_file.get_nonterminals(order='preorder'):
        if pd.isna(branch.name) == False:
                predupl_sequence = seq_dict[branch.name][int(subtype):]
                #print(predupl_sequence)
                break  
    post_dupl_1 = first_dupl_sequence[:int(subtype)]
    post_dupl_2 = first_dupl_sequence[int(subtype):]

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
            if i%3 == 0 and i !=0 :
                codon_plus_prev = str(duplication[i-1:i+3])
                codon_plus_after = str(duplication[i:i+4])
                codon_plus_both = str(duplication[i-1:i+4])
                codon = str(duplication[i:i+3])
                if type_ == "one_before":
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
                                            
                elif type_ == "one_after" and len(codon_plus_after)==4:
                    for key, entry in translations.items():
                        if codon in entry:
                                if len(entry) == 4:  
                                    each_position[i+2] = sum_of_rows[codon_plus_after[2:]]
                                if len(entry) == 2 or len(entry) == 3:
                                    sum_of_muts = 0
                                    for ele in entry: sum_of_muts += mutation_matrix.at[codon_plus_after[2:], ele[-1]]
                                    each_position[i+2] = sum_of_muts
                                if len(entry) == 6:
                                    each_position[i+2] = sum_of_rows[codon_plus_after[2:]]
                                    for ele in entry:
                                        if ele[0] != codon[0]:
                                            sum_ = float(mutation_matrix.at[codon_plus_after[:-2], ele[0]])
                                            each_position[i] = sum_
                                            
                elif type_ == "before_after" and len(codon_plus_both)==5:
                    for key, entry in translations.items():
                        if codon in entry:
                                if len(entry) == 4: 
                                    each_position[i+2] = sum_of_rows[codon_plus_both[2:]]
                                if len(entry) == 2 or len(entry) == 3:
                                    sum_of_muts = 0
                                    for ele in entry: sum_of_muts += mutation_matrix.at[codon_plus_both[2:], ele[-1]]
                                    each_position[i+2] = sum_of_muts
                                if len(entry) == 6:
                                    each_position[i+2] = sum_of_rows[codon_plus_both[2:]]
                                    for ele in entry:
                                        if ele[0] != codon[0]:
                                            sum_ = float(mutation_matrix.at[codon_plus_both[:-2], ele[0]])
                                            each_position[i] = sum_

    cumulative_= cumulative(each_position)
    return(cumulative_)

    
def nonsynonymous(duplication, mutation_matrix, type_):
    translations = {'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'], 'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'], 'C': ['TGT', 'TGC'], 'W': ['TGG'], 'E': ['GAA', 'GAG'], 'D': ['GAT', 'GAC'], 'P': ['CCT', 'CCC', 'CCA', 'CCG'], 'V': ['GTT', 'GTC', 'GTA', 'GTG'], 'N': ['AAT', 'AAC'], 'M': ['ATG'], 'K': ['AAA', 'AAG'], 'Y': ['TAT', 'TAC'], 'I': ['ATT', 'ATC', 'ATA'], 'Q': ['CAA', 'CAG'], 'F': ['TTT', 'TTC'], 'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'], 'T': ['ACT', 'ACC', 'ACA', 'ACG'], '*': ['TAA', 'TAG', 'TGA'], 'A': ['GCT', 'GCC', 'GCA', 'GCG'], 'G': ['GGT', 'GGC', 'GGA', 'GGG'], 'H': ['CAT', 'CAC']}
    each_position = dict()
    #print("nonsyn")
    for i, char in enumerate(duplication):
        if i%3 == 0:
            codon = duplication[i:i+3]
            codon_with_prev = duplication[i-1:i+3]
            codon_with_after = duplication[i:i+4]
            codon_before_after = duplication[i-1:i+4]
            sum_of_muts_1, sum_of_muts_2, sum_of_muts_3 = (0 for j in range(3))
            new_codons_1, new_codons_2, new_codons_3 = ([] for lst in range(3))
            for pos, char in enumerate(codon):
                sum_of_muts_1, sum_of_muts_2, sum_of_muts_3 = (0 for j in range(3))
                for mut in ["A", "C", "G", "T"]:
                    if pos == 0 and codon[0] != mut: new_codons_1.append(f"{mut}{codon[1]}{codon[2]}")
                    if pos == 1 and codon[1] != mut: new_codons_2.append(f"{codon[0]}{mut}{codon[2]}")
                    if pos == 2 and codon[2] != mut: new_codons_3.append(f"{codon[0]}{codon[1]}{mut}")
                for key, entry in translations.items():
                    if codon in entry:
                        codons_diff_1 = set(new_codons_1).difference(set(entry))
                        codons_diff_2 = set(new_codons_2).difference(set(entry))
                        codons_diff_3 = set(new_codons_3).difference(set(entry))
                if type_ == "point_mut":
                    for ele in codons_diff_1:  sum_of_muts_1 += mutation_matrix.at[codon[0], ele[0]]
                    for ele in codons_diff_2: sum_of_muts_2 += mutation_matrix.at[codon[1], ele[1]]
                    for ele in codons_diff_3: sum_of_muts_3 += mutation_matrix.at[codon[2], ele[2]]
                    each_position[i] = sum_of_muts_1
                    each_position[i+1] = sum_of_muts_2
                    each_position[i+2] = sum_of_muts_3
                elif type_ == "one_before" and len(codon_with_prev) == 4:        
                    for ele in codons_diff_1: sum_of_muts_1 += mutation_matrix.at[f"{codon_with_prev[0]}{codon[0]}", ele[0]]
                    for ele in codons_diff_2:sum_of_muts_2 += mutation_matrix.at[f"{codon[0]}{codon[1]}", ele[1]]
                    for ele in codons_diff_3: sum_of_muts_3 += mutation_matrix.at[f"{codon[1]}{codon[2]}", ele[2]]
                    each_position[i] = sum_of_muts_1
                    each_position[i+1] = sum_of_muts_2
                    each_position[i+2] = sum_of_muts_3
                elif type_ == "one_after" and len(codon_with_after)==4:        
                    for ele in codons_diff_1: sum_of_muts_1 += mutation_matrix.at[f"{codon[0]}{codon[1]}", ele[0]]
                    for ele in codons_diff_2: sum_of_muts_2 += mutation_matrix.at[f"{codon[1]}{codon[2]}", ele[1]]
                    for ele in codons_diff_3: sum_of_muts_3 += mutation_matrix.at[f"{codon[2]}{codon_with_after[-1]}", ele[2]]
                    each_position[i] = sum_of_muts_1
                    each_position[i+1] = sum_of_muts_2
                    each_position[i+2] = sum_of_muts_3
                elif type_ == "before_after" and len(codon_before_after)==5:        
                    for ele in codons_diff_1: sum_of_muts_1 += mutation_matrix.at[f"{codon_before_after[0]}{codon[0]}{codon[1]}", ele[0]]
                    for ele in codons_diff_2: sum_of_muts_2 += mutation_matrix.at[f"{codon}", ele[1]]
                    for ele in codons_diff_3: sum_of_muts_3 += mutation_matrix.at[f"{codon[1]}{codon[2]}{codon_before_after[-1]}", ele[2]]
                    each_position[i] = sum_of_muts_1
                    each_position[i+1] = sum_of_muts_2
                    each_position[i+2] = sum_of_muts_3
    cumulative_= cumulative(each_position)
    return(cumulative_)

if __name__=="__main__":
    parser = argparse.ArgumentParser(
        description="Construct synynymous mutation matrix",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument('--matrix', required=True, type=str, help="csv file with scaled gene mutations")
    parser.add_argument('--tree', required=True, type=str, help="Tree nwk")
    parser.add_argument('--ref', type=str, help="genbank reference file")
    parser.add_argument('--syn', type=str, help="graph png file")
    parser.add_argument('--nonsyn', type=str, help="graph png file")
    parser.add_argument('--type', required=True, help="context of the mutation")
    parser.add_argument("--subtype", required=True, help="A or B - string")
    parser.add_argument('--duplicationseq', required=True, help="fasta file with reconstructed sequences")

    args = parser.parse_args()
    mut_matrix = pd.read_csv(args.matrix)
    mut_matrix = mut_matrix.set_index('Unnamed: 0')

    post_1 = finding_duplications(args.tree, args.duplicationseq, args.subtype,  "1")    
    post_2 = finding_duplications(args.tree, args.duplicationseq, args.subtype, "2")   
    pre = finding_duplications(args.tree, args.duplicationseq, args.subtype, "pre")    

    cumulative_1 = dictionary_of_mutations(post_1, mut_matrix, args.type)
    cumulative_2 = dictionary_of_mutations(post_2, mut_matrix, args.type)
    cumulative_pre = dictionary_of_mutations(pre, mut_matrix, args.type)

    cumulative_1_nonsyn = nonsynonymous(post_1, mut_matrix, args.type)
    cumulative_2_nonsyn = nonsynonymous(post_2, mut_matrix, args.type)
    cumulative_pre_nonsyn = nonsynonymous(pre, mut_matrix, args.type)

    plt.figure()
    plt.step(cumulative_1_nonsyn.keys(), cumulative_1_nonsyn.values(), label= f'1st copy postduplication', where='post')
    plt.step(cumulative_2_nonsyn.keys(), cumulative_2_nonsyn.values(), label=f'2nd copy postduplication', where='post' )
    plt.step(cumulative_pre_nonsyn.keys(), cumulative_pre_nonsyn.values(), label= f'preduplication', where='post' )

    plt.legend(loc='lower right')
    plt.xlabel("gene location")
    plt.ylabel("cumulative sum of mutations")
    plt.savefig(args.nonsyn)

    cumulative_1 = dictionary_of_mutations(post_1, mut_matrix, args.type)
    cumulative_2 = dictionary_of_mutations(post_2, mut_matrix, args.type)
    cumulative_pre = dictionary_of_mutations(pre, mut_matrix, args.type)

    plt.figure()
    plt.step(cumulative_1.keys(), cumulative_1.values(), label= f'1st copy postduplication', where='post')
    plt.step(cumulative_2.keys(), cumulative_2.values(), label=f'2nd copy postduplication', where='post' )
    plt.step(cumulative_pre.keys(), cumulative_pre.values(), label= f'preduplication', where='post' )

    plt.legend(loc='lower right')
    plt.xlabel("gene location")
    plt.ylabel("cumulative sum of mutations")
    plt.savefig(args.syn)