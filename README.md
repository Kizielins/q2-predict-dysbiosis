# q2-predict-dysbiosis

This repository contains our function-based gut microbiome health index as well as steps to reproduce analyses contained within our manuscript, linked at the bottom of the page.
The model can be ran in two modes:
(a) full model 
(b) the model without additional taxonomic features (better at predicting diseases that do not originate in the gut)

For more details please see our manuscript.

## Input prep:
Sample inputs can be found in the "test_files" folder.

Taxonomy table: standard QIIME 2 *qza feature table, collapsed to species level, with removed "s__" and underscores instead of spaces (ie "Escherichia_coli")

Stratified pathways table: standard QIIME 2 *qza feature table, produced by HUMAnNN, collapsed to species level, with underscores instead of spaces (ie ANAEROFRUCAT-PWY:_homolactic_fermentation|g__Citrobacter.s__Citrobacter_freundii)

Unstratified pathways table: standard QIIME 2 *qza feature table, produced by HUMAnNN, with underscores instead of spaces (ie AEROBACTINSYN-PWY:_aerobactin_biosynthesis)

Metadata: standard QIIME 2 metadata format, with "id" and columns representing sample IDs and labelling.
The values in all tables should be expressed as relative abundance.

## Running q2-predict-dysbiosis

The easiest way to run the script is to clone the repository, modify the "run_q2_predict_dysbiosis.py" accordingly and run! :)

## Original publication / citation
If you want to learn more about this method, or to cite it, please refer to our article: https://www.biorxiv.org/content/10.1101/2023.12.04.569909v4

