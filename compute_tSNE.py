#!/home/ferrari/anaconda3/bin/python

from sklearn.manifold import TSNE
import numpy as np
import pandas as pd 
import argparse

def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("-t","--chromatin_table", help="table output of 'from_segmentation_to_final_table.py'", dest="final_table", required=True)
    parser.add_argument("-p","--perplexity", type=int, help="perplexity of tSNE", default=30)
    parser.add_argument("--seed",type=int, help="seed for tSNE", default=2)
    
    return parser

def main():
    parser = parse_args()
    args = parser.parse_args()
    print(args)

    data_df = pd.read_csv(args.final_table, sep="\t", index_col=0)
    #print(data_matrix.head())
    data_matrix = data_df.as_matrix()
    #print(type(data_matrix))
    #print(data_matrix.shape)

    #run tSNE
    X_embedded = TSNE(n_components=2, random_state=args.seed, perplexity=args.perplexity).fit_transform(data_matrix)

    #save results to file
    to_be_printed=pd.DataFrame(X_embedded, data_df.index, columns=["Dim1","Dim2"])
    to_be_printed.to_csv("./tSNE_COORDINATES_perp{}-seed{}.tsv".format(args.perplexity, args.seed), sep="\t")

if __name__ == "__main__":
    main()

