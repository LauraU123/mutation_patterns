import argparse
import pandas as pd
import matplotlib.pyplot as plt
from collections import OrderedDict

if __name__=="__main__":
    parser = argparse.ArgumentParser(
        description="Construct synynymous mutation spetra",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument('--matrix', required=True, type=str, help="matrix type")
    parser.add_argument('--output', type=str, help="spectra png file")
    parser.add_argument('--rsvsubtype', type=str, help="a or b")
    args = parser.parse_args()

    lst = ["A", "C", "G", "T"]
    lst_double = ['AA', 'CA', 'GA', 'TA', 'AC', 'CC', 'GC', 'TC', 'AG', 'CG', 'GG', 'TG' ,'AT', 'CT', 'GT', 'TT']
    lst_triple =['AAA', 'AAC', 'AAG', 'AAT', 'CAA', 'GAA', 'TAA', 'CAC', 'CAG', 'CAT', 'GAC', 'GAG', 'GAT', 'TAC', 'TAG', 'TAT' ,
                                                           'ACA', 'ACC', 'ACG', 'ACT', 'CCA', 'GCA', 'TCA', 'CCC', 'CCG', 'CCT', 'GCC', 'GCG', 'GCT', 'TCC', 'TCG', 'TCT' ,
                                                           'AGA', 'AGC', 'AGG', 'AGT', 'CGA', 'GGA', 'TGA', 'CGC', 'CGG', 'CGT', 'GGC', 'GGG', 'GGT', 'TGC', 'TGG', 'TGT' ,
                                                           'ATA', 'ATC', 'ATG', 'ATT', 'CTA', 'GTA', 'TTA', 'CTC', 'CTG', 'CTT', 'GTC', 'GTG', 'GTT', 'TTC', 'TTG', 'TTT'] 


    mut_matrix = pd.read_csv(f"results/{args.matrix}_{args.rsvsubtype}_matrix.csv")
    mut_matrix = mut_matrix.set_index('Unnamed: 0')

    if args.matrix == "point_mut":
        dct = dict()
        for nuc in lst:
            for nucl in lst: 
                if nuc != nucl: dct[f'{nuc}->{nucl}'] = mut_matrix.at[nuc, nucl]
        ordered = OrderedDict(dct)
        plt.plot(dct.keys(), dct.values(), 'o', markersize=10)
        plt.xticks(rotation=90)
        plt.xticks(fontsize= 20, weight = 'bold')
        plt.yticks(fontsize=20)
        plt.tight_layout()
        plt.savefig(args.output)
    
    elif args.matrix == "one_before":
        dct = dict()
        for i in lst_double:
            for j in lst: 
                if i[-1] != j: dct[f'{i}->{j}'] = mut_matrix.at[i, j]
        
        plt.figure(figsize=(30,15))
        plt.plot(dct.keys(), dct.values(), 'o', markersize=20)
        plt.xticks(rotation=90)
        plt.xticks(fontsize=40, weight = 'bold')
        plt.yticks(fontsize=40)
        plt.tight_layout()
        plt.savefig(args.output)

    elif args.matrix == "one_after":
        dct = dict()
        for i in lst_double:
            for j in lst: 
                if i[0] != j: dct[f'{i}->{j}'] = mut_matrix.at[i, j]
        ordered = OrderedDict(dct)
        plt.figure(figsize=(30,15))
        plt.plot(ordered.keys(), ordered.values(), 'o', markersize=20)
        plt.xticks(rotation=90)
        plt.xticks(fontsize= 40, weight = 'bold')
        plt.yticks(fontsize=40)
        plt.tight_layout()
        plt.savefig(args.output)
    
    elif args.matrix == "before_after":
        dct = dict()
        for i in lst_triple:
            for j in lst: 
                if i[1] != j: dct[f'{i}->{j}'] = mut_matrix.at[i, j]
        plt.figure(figsize=(30,15))
        plt.plot(dct.keys(), dct.values(), 'o', markersize=20)
        plt.xticks(rotation=90)
        plt.xticks(fontsize= 40, weight = 'bold')
        plt.yticks(fontsize=40)
        plt.tight_layout()
        plt.savefig(args.output)

         #A, C, G, T       
        for base in lst:
            dct_for_subtype = dict()
            for entry, code in dct.items():
                #for each entry in the dict, if the middle value is the base in question
                    if entry[1] == base:dct_for_subtype[entry] = code
            for base_ in lst:
                dct_for_base = dict()
                for entry_, code_ in dct_for_subtype.items():
                    if entry_[-1] == base_: 
                        print(entry_, entry_[-1])
                        dct_for_base[entry_] = code_

                plt.figure(figsize=(20,14))
                plt.plot(dct_for_base.keys(), dct_for_base.values(), 'o', markersize=20)
                plt.xticks(rotation=90)
                plt.xticks(fontsize= 40, weight = 'bold')
                plt.yticks(fontsize=40)
                plt.tight_layout()
                plt.savefig(args.rsvsubtype + "_before_after" + f"{base}_to_{base_}.png")
