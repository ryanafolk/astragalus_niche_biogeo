# astragalus_niche_biogeo
Analyses for the paper entitled, "Anatomy of a mega-radiation: Biogeography and niche evolution in *Astragalus*". Folders are as follows in subsections.

##all_trees_undated


## BAMM
Scripts, input, and output for BAMM. `Tree_check.R` is basic input tree and prior checking (output `myPriors.txt`). `config.txt` gives settings for the BAMM run. `bammtools_diversification.R` is postprocessing of BAMM results in BAMMtools. `sample_fractions.txt` is the file that encodes missing taxon proportions. `chromosomes_raw.txt` gives input chromosome data, recoded in file `chromosomes_cleaned.txt` per the Methods section. Otherwise all `.csv` files are input data for trait-associated diversification and all `.svg` and `.pdf` files in folder `figs`are figures from the output.


## biogeography_manuallycorrected_trimmedtree


## categorical_reconstructions


## chromosome_reconstruction


## continuous_ancestral_reconstructions


## dated_trees


## DR_rates
`tip_DR.csv` gives calculated DR rates; `figs` contains figures from outputs comparing DR to other data.


## environmental_data
CSVs representing environmental means (soil types in a separate directory).

## nquire


## occurrences
CSVs representing raw occurrences.


## soil_types
CSVs representing soil types, given both as full names and numerical coding (see `.xlsx` for number-name correspondence).


## tree_quality_control


## Tree_with_labels
