# q2-predict-dysbiosis

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

A full script to reproduce all figures in the article will be available shortly.

## Acknowledgements
We would like to acknowledge the Authors of the q2-health-index plugin, whose scripts formed the foundation of our work.
