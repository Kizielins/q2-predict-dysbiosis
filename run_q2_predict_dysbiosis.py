from q2_predict_dysbiosis import calculate_index 
import pandas as pd


taxa = pd.read_csv("test_files/species.txt", sep="\t", index_col=0)
paths_strat = pd.read_csv("test_files/pathways_stratified.txt", sep="\t", index_col=0)
paths_unstrat = pd.read_csv("test_files/pathways_unstratified.txt", sep="\t", index_col=0)

mode = "full" # if you want to run without taxonomic features, change this to "no_taxonomy"

print(calculate_index(taxa, paths_strat, paths_unstrat, mode))