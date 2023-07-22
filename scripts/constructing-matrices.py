import argparse
import json
from Bio import SeqIO
import pandas as pd
from collections import Counter, defaultdict


#Translation Matrix 
translations = {'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'], 'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'], 'C': ['TGT', 'TGC'], 'W': ['TGG'], 'E': ['GAA', 'GAG'], 'D': ['GAT', 'GAC'], 'P': ['CCT', 'CCC', 'CCA', 'CCG'], 'V': ['GTT', 'GTC', 'GTA', 'GTG'], 'N': ['AAT', 'AAC'], 'M': ['ATG'], 'K': ['AAA', 'AAG'], 'Y': ['TAT', 'TAC'], 'I': ['ATT', 'ATC', 'ATA'], 'Q': ['CAA', 'CAG'], 'F': ['TTT', 'TTC'], 'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'], 'T': ['ACT', 'ACC', 'ACA', 'ACG'], '*': ['TAA', 'TAG', 'TGA'], 'A': ['GCT', 'GCC', 'GCA', 'GCG'], 'G': ['GGT', 'GGC', 'GGA', 'GGG'], 'H': ['CAT', 'CAC']}


def CDS_finder_gbk(reference):
    #this function finds CDS location 
    cds_ = dict()
    for feature in reference.features:
        if feature.type == 'CDS':  cds_[feature.qualifiers['gene'][0]] = (list(feature.location))
    return(cds_)
"""
def Synonymous_Mutations(reffile, node, dictionary_=None, new_=None):
    #Finds Synonymous mutations in CDS regions. Input:nested json
    ref_file = SeqIO.read(reffile, "genbank")
    gene_cds = CDS_finder(ref_file)
    if new_ is None: new_ = []
    if dictionary_ is None: dictionary_ = dict()
    if 'mutations' in node['branch_attrs']:
        aa_mutations, new_, in_it = ([] for i in range(3))
        if 'nuc' in node['branch_attrs']['mutations']:
            for gene, loc in gene_cds.items():
                if gene in node['branch_attrs']['mutations']:
                    for mut in node['branch_attrs']['mutations'][gene]:
                        aa_mutations.append(int(mut[1:-1])*3+loc[0]) # converting the amino acid location to nucleotide
            for mut in node['branch_attrs']['mutations']['nuc']:
                if '-' not in mut and '*' not in mut and 'N' not in mut and "R" not in mut and "Y" not in mut and "M" not in mut and "D" not in mut:
                    print(mut, aa_mutations)
                    #make sure the codon is correct. 
                    if int(mut[1:-1]) not in aa_mutations and int(mut[1:-1])+2 not in aa_mutations and int(mut[1:-1])+1 not in aa_mutations:  #if the mutation is not in the same codon as a aa mutation
                        new_.append(mut)
                    else: in_it.append(mut[1:-1])
    if 'name' in node: dictionary_[node['name']] = new_
    if 'children' in node:
        for child in node['children']: Synonymous_Mutations(reffile, child, dictionary_, new_=None)
    return(dictionary_)
"""

def CDS_finder(jsonfile):
    CDS_locations = dict()
    for gene, data in jsonfile['meta']['genome_annotations'].items(): CDS_locations[gene] = [i for i in range(data['start'], data['end']+1)]
    CDS_locations.pop('nuc')
    return(CDS_locations)


def all_loc_CDS(jsonfile):
    gene_cds = CDS_finder(jsonfile)
    all_loc_CDS_ = []
    for gene, locations in gene_cds.items(): all_loc_CDS_.extend(locations)
    return(all_loc_CDS_)
    
    

def Synonymous_Mutations(f, node, dictionary_=None, new_=None):
    """ Finds Synonymous mutations in CDS regions. Input:nested json"""
    gene_cds = CDS_finder(f)
    all_cds = all_loc_CDS(f)
    if new_ is None: new_ = []
    if dictionary_ is None: dictionary_ = dict()
    if 'mutations' in node['branch_attrs']:
        aa_mutations, new_ = ([] for i in range(2))
        if 'nuc' in node['branch_attrs']['mutations']:
            for gene, loc in gene_cds.items():
                if gene in node['branch_attrs']['mutations']:
                    for mut in node['branch_attrs']['mutations'][gene]:
                        #each possible codon location
                        aa_mutations.append(int(mut[1:-1])*3+ loc[0]-1) 
                        aa_mutations.append(int(mut[1:-1])*3+ loc[0]-2) 
                        aa_mutations.append(int(mut[1:-1])*3+ loc[0]-3) 
            for nucl in node['branch_attrs']['mutations']['nuc']:
                if int(nucl[1:-1]) not in aa_mutations: 
                    if '-' not in nucl and int(nucl[1:-1]) in all_cds: new_.append(nucl)
    if 'name' in node: dictionary_[node['name']] = new_
    if 'children' in node:
        for child in node['children']: Synonymous_Mutations(f, child, dictionary_, new_=None)
    return(dictionary_)


def mutations_matrix_unscaled(synonymous, reconstructed, which, type_):
    """Constructing matrix from synonymous mutations"""
    mut_by_branch_CDS, all_dinucleotides = (defaultdict(list) for i in range(2))
    all_muts = []

    for branch, muts in synonymous.items():
        if branch != []:
            for mut in muts:
                if mut[0] and mut[-1] in ["A", "T", "C", "G"]: 
                    all_muts.append(f'{mut[0]}{mut[-1]}')
                    mut_by_branch_CDS[branch].append(mut)

    if type_ != "point_mut":
        aligned_for_tree = SeqIO.parse(reconstructed, "fasta")
        for entry in aligned_for_tree:
            for i in mut_by_branch_CDS[entry.id]:
                location_of_interest = int(i[1:-1])
                if type_ == "one_before":
                    if entry.seq[location_of_interest-2] != '-':
                        all_dinucleotides[f'{entry.seq[location_of_interest-2]}{i[0]}'].append(i[-1])
                elif type_ == "one_after":
                    if entry.seq[location_of_interest+1] != '-': all_dinucleotides[f'{i[0]}{entry.seq[location_of_interest+1]}'].append(i[-1])
                elif type_ == "before_after":
                    print(i, entry.id)
                    print(entry.seq[location_of_interest-1])
                    #print(entry.seq[location_of_interest+1])

                    
                    if entry.seq[location_of_interest+1] != '-' and entry.seq[location_of_interest-2] != '-': all_dinucleotides[f'{entry.seq[location_of_interest-2]}{i[0]}{entry.seq[location_of_interest+1]}'].append(i[-1])

        with_counters = dict()
        for type, mut in all_dinucleotides.items():
            with_counters[type] = Counter(mut)


    all_muts_counter = Counter(all_muts)
    if type_ == "point_mut": df = pd.DataFrame(index=['A', 'C', 'G', 'T'], columns=['A', 'C', 'G', 'T'])

    elif type_ == "one_before": df = pd.DataFrame(index=['AA', 'CA', 'GA', 'TA', 'AC', 'CC', 'GC', 'TC', 'AG', 'CG', 'GG', 'TG' ,'AT', 'CT', 'GT', 'TT'], columns=['A', 'C', 'G', 'T'])   
    elif type_ == "one_after": df = pd.DataFrame(index=['AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT', 'GA', 'GC', 'GG', 'GT', 'TA', 'TC', 'TG', 'TT'], columns=['A', 'C', 'G', 'T'])   
    elif type_ == "before_after": df = pd.DataFrame(index=['AAA', 'AAC', 'AAG', 'AAT', 'CAA', 'GAA', 'TAA', 'CAC', 'CAG', 'CAT', 'GAC', 'GAG', 'GAT', 'TAC', 'TAG', 'TAT' ,
                                                           'ACA', 'ACC', 'ACG', 'ACT', 'CCA', 'GCA', 'TCA', 'CCC', 'CCG', 'CCT', 'GCC', 'GCG', 'GCT', 'TCC', 'TCG', 'TCT' ,
                                                           'AGA', 'AGC', 'AGG', 'AGT', 'CGA', 'GGA', 'TGA', 'CGC', 'CGG', 'CGT', 'GGC', 'GGG', 'GGT', 'TGC', 'TGG', 'TGT' ,
                                                           'ATA', 'ATC', 'ATG', 'ATT', 'CTA', 'GTA', 'TTA', 'CTC', 'CTG', 'CTT', 'GTC', 'GTG', 'GTT', 'TTC', 'TTG', 'TTT'] , columns=['A', 'C', 'G', 'T'])   
    
    if type_ == "point_mut":
        for mutation, nr in all_muts_counter.items(): df.at[mutation[0], mutation[-1]] = int(nr)
        return(df)
        
    else:
        for mutation, count in with_counters.items():
            if 'N' not in mutation:
                for type, c in count.items():df.at[mutation, type] = c
        return(df)


def possible_syn_mut_locations(reference):
    """finding ratio of synonymous to nonsynonymous mutation locations in the genome"""
    ref_file = SeqIO.read(reference, "genbank")
    gene_cds = CDS_finder_gbk(ref_file)
    sequence_ref_cds = dict()
    whole_seq_CDS = ""
    for gene, cds in gene_cds.items(): 
        sequence_ref_cds[gene] = ref_file.seq[cds[0]:cds[-1]]
        whole_seq_CDS = whole_seq_CDS+ref_file.seq[cds[0]:cds[-1]]

    synonymous_possibilities, nonsynonymous_possibilities = (0 for i in range(2))
    for gene, sequence in sequence_ref_cds.items():
        for i, letter in enumerate(sequence):
            if i%3 == 0:
                codon = sequence[i: i+3]
                for key, entry in translations.items():
                    if codon in entry:
                        synonymous = len(entry)
                        synonymous_possibilities += synonymous-1 
                        nonsynonymous = 9 - len(entry)-1
                        nonsynonymous_possibilities += nonsynonymous
    syn_ratio = synonymous_possibilities/(nonsynonymous_possibilities+synonymous_possibilities)
    return(syn_ratio)



#now have to divide by the nr of locations where the mutation can occur


def scaled_by_ratio(dataframe, ratio):
    """dataframe scaled by a given ratio"""
    scaled = dataframe.divide(ratio)
    scaled = scaled.fillna(0)
    return(scaled)

def count_of_nucleotides_in_syn_positions(reference, type_):
    ref_file = SeqIO.read(reference, "genbank")
    gene_cds = CDS_finder_gbk(ref_file)
    sequence_ref_cds = dict()
    for gene, cds in gene_cds.items(): sequence_ref_cds[gene] = ref_file.seq[cds[0]:cds[-1]]
    list_all= []
    for gene, sequence in sequence_ref_cds.items():
        for i, letter in enumerate(sequence):
            if i%3 == 0:
                codon = sequence[i: i+3]
                for key, entry in translations.items():
                    if codon in entry:
                        if type_ == "point_mut":
                            if 1< len(entry) <= 4: list_all.append(codon[-1])
                            elif len(entry) > 4:
                                list_all.append(codon[0])
                                list_all.append(codon[0])
                            elif len(entry) == 1: continue
                        elif type_ == "one_before":
                            if 1< len(entry) <= 4: list_all.append(str(codon[1:]))
                            elif len(entry) > 4:
                                list_all.append(str(sequence[i-1:i+1]))
                                list_all.append(str(codon[1:]))
                            elif len(entry) == 1: continue
                        elif type_ == "one_after":
                            if 1< len(entry) <= 4: list_all.append(str(sequence[i+2:i+4]))
                            elif len(entry) > 4:
                                list_all.append(str(codon[:-1]))
                                list_all.append(str(sequence[i+2:i+4]))
                            elif len(entry) == 1: continue
                        elif type_ == "before_after":
                            if 1< len(entry) <= 4: list_all.append(str(sequence[i+1:i+4]))
                            elif len(entry) > 4:
                                list_all.append(str(sequence[i-1:i+2]))
                                list_all.append(str(sequence[i+1:i+4]))
                            elif len(entry) == 1: continue

    counter = Counter(list_all)
    return(counter)



def scaled_by_nucleotides_and_normalized(dataframe, nucl_dict, type_):
    """scaling dataframe by nucleotides at synonymous positions in reference"""
    total = 0
    if type_ == "point_mut": nuc = ["A", "C", "G", "T"]
    if type_ == "one_before": nuc = ['AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT', 'GA', 'GC', 'GG', 'GT', 'TA', 'TC', 'TG', 'TT']
    if type_ == "one_after": nuc = ['AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT', 'GA', 'GC', 'GG', 'GT', 'TA', 'TC', 'TG', 'TT']
    if type_ == "before_after": nuc = ['AAA', 'AAC', 'AAG', 'AAT', 'CAA', 'GAA', 'TAA', 'CAC', 'CAG', 'CAT', 'GAC', 'GAG', 'GAT', 'TAC', 'TAG', 'TAT' ,
                                                           'ACA', 'ACC', 'ACG', 'ACT', 'CCA', 'GCA', 'TCA', 'CCC', 'CCG', 'CCT', 'GCC', 'GCG', 'GCT', 'TCC', 'TCG', 'TCT' ,
                                                           'AGA', 'AGC', 'AGG', 'AGT', 'CGA', 'GGA', 'TGA', 'CGC', 'CGG', 'CGT', 'GGC', 'GGG', 'GGT', 'TGC', 'TGG', 'TGT' ,
                                                           'ATA', 'ATC', 'ATG', 'ATT', 'CTA', 'GTA', 'TTA', 'CTC', 'CTG', 'CTT', 'GTC', 'GTG', 'GTT', 'TTC', 'TTG', 'TTT']

    for i in nucl_dict.values(): total += i

    df_ratios = pd.DataFrame.from_dict(nucl_dict, orient='index').astype(int).T
    df_ratios = df_ratios.divide(total)
    for n in nuc: dataframe.loc[[n]] = dataframe.loc[[n]].div(float(df_ratios[n]))

    #normalizing step
    sum_ = dataframe.to_numpy().sum()
    scaled_and_normalized = dataframe.divide(sum_)
    return(scaled_and_normalized)


if __name__=="__main__":
    parser = argparse.ArgumentParser(
        description="Construct synynymous mutation matrix",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument('--ref', required=True, type=str, help="reference file, genbank format")
    parser.add_argument('--tree', required=True, type=str, help="Tree json annotated with amino acid and nucleotide mutations")
    parser.add_argument('--output', type=str, help="output CSV file")
    parser.add_argument('--rsvsubtype', type=str, help="a or b")
    parser.add_argument("--reconstructedseq", type=str, required=True)
    parser.add_argument('--type', type=str, help="type of context of mutation")
    args = parser.parse_args()

    with open (args.tree) as file_:
        f = json.load(file_)  

    CDS_dict = CDS_finder(f)

    synonymous = Synonymous_Mutations(f, f['tree'])

    mutation_matrix = mutations_matrix_unscaled(synonymous, args.reconstructedseq, args.rsvsubtype, args.type)

    syn_ratio = possible_syn_mut_locations(args.ref)    

    scaled_by_ratio_ = scaled_by_ratio(mutation_matrix, syn_ratio)

    syn_mut_count_reference = count_of_nucleotides_in_syn_positions(args.ref, args.type)

    scaled_normalized = scaled_by_nucleotides_and_normalized(scaled_by_ratio_, syn_mut_count_reference, args.type)

    scaled_normalized.to_csv(args.output)


