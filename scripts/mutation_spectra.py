import argparse
import pandas as pd
import matplotlib.pyplot as plt



if __name__=="__main__":
    parser = argparse.ArgumentParser(
        description="Construct synynymous mutation spetra",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument('--matrix', required=True, type=str, help="matrix type")
    parser.add_argument('--output', type=str, help="spectra png file")
    args = parser.parse_args()


    lst = ["A", "C", "G", "T"]
    lst_double = ['AA', 'CA', 'GA', 'TA', 'AC', 'CC', 'GC', 'TC', 'AG', 'CG', 'GG', 'TG' ,'AT', 'CT', 'GT', 'TT']
    lst_triple =['AAA', 'AAC', 'AAG', 'AAT', 'CAA', 'GAA', 'TAA', 'CAC', 'CAG', 'CAT', 'GAC', 'GAG', 'GAT', 'TAC', 'TAG', 'TAT' ,
                                                           'ACA', 'ACC', 'ACG', 'ACT', 'CCA', 'GCA', 'TCA', 'CCC', 'CCG', 'CCT', 'GCC', 'GCG', 'GCT', 'TCC', 'TCG', 'TCT' ,
                                                           'AGA', 'AGC', 'AGG', 'AGT', 'CGA', 'GGA', 'TGA', 'CGC', 'CGG', 'CGT', 'GGC', 'GGG', 'GGT', 'TGC', 'TGG', 'TGT' ,
                                                           'ATA', 'ATC', 'ATG', 'ATT', 'CTA', 'GTA', 'TTA', 'CTC', 'CTG', 'CTT', 'GTC', 'GTG', 'GTT', 'TTC', 'TTG', 'TTT'] 


    mut_matrix = pd.read_csv(f"results/{args.matrix}_matrix.csv")
    mut_matrix = mut_matrix.set_index('Unnamed: 0')


    if args.matrix == "point_mut":
        dct = dict()
        for nuc in lst:
            for nucl in lst: 
                if nuc != nucl: dct[f'{nuc}->{nucl}'] = mut_matrix.at[nuc, nucl]
        plt.plot(dct.keys(), dct.values(), 'o')
        plt.xticks(rotation=90)
        plt.savefig(args.output)
    
    elif args.matrix == "one_before":
        dct = dict()
        for i in lst_double:
            for j in lst: 
                if i[-1] != j: dct[f'{i}->{j}'] = mut_matrix.at[i, j]
        plt.figure(figsize=(20,5))
        plt.plot(dct.keys(), dct.values(), 'o')
        plt.xticks(rotation=90)
        plt.savefig(args.output)

    elif args.matrix == "one_after":
        dct = dict()
        for i in lst_double:
            for j in lst: 
                if i[0] != j: dct[f'{i}->{j}'] = mut_matrix.at[i, j]
        plt.figure(figsize=(20,5))
        plt.plot(dct.keys(), dct.values(), 'o')
        plt.xticks(rotation=90)
        plt.savefig(args.output)
    
    elif args.matrix == "before_after":
        dct = dict()
        for i in lst_triple:
            for j in lst: 
                if i[1] != j: dct[f'{i}->{j}'] = mut_matrix.at[i, j]
        plt.figure(figsize=(40,7))
        plt.plot(dct.keys(), dct.values(), 'o')
        plt.xticks(rotation=90)
        plt.savefig(args.output)

        
    