from q2_predict_dysbiosis import calculate_index 
import pandas as pd


taxa = pd.read_csv("test_files/species.txt", sep="\t", index_col=0)
paths_strat = pd.read_csv("test_files/pathways_stratified.txt", sep="\t", index_col=0)
paths_unstrat = pd.read_csv("test_files/pathways_unstratified.txt", sep="\t", index_col=0)


print(calculate_index(taxa, paths_strat, paths_unstrat))