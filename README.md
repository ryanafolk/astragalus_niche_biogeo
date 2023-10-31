# astragalus_niche_biogeo
Analyses for the paper entitled, "Anatomy of a mega-radiation: Biogeography and niche evolution in *Astragalus*". Folders are as follows in subsections.

## all_trees_undated
Trees before dating, branches in coalescent units. `astragalus.hybpiper.astral.scored.sectionnames.afterprunning.tre` has uncorrected names (conformant to the NitFix taxonomy prior to miscellaneous changes) and `astragalus.hybpiper.astral.scored.sectionnames.afterprunning.taxonomycorrection.tre` has corrected names. `figs` contains PDFs of plotted trees. `bs10` indicates the result from an additional filtering step on the gene trees.

## BAMM
Scripts, input, and output for BAMM. `Tree_check.R` is basic input tree and prior checking (output `myPriors.txt`). `config.txt` gives settings for the BAMM run. `bammtools_diversification.R` is postprocessing of BAMM results in BAMMtools. `sample_fractions.txt` is the file that encodes missing taxon proportions. `chromosomes_raw.txt` gives input chromosome data, recoded in file `chromosomes_cleaned.txt` per the Methods section. Otherwise all `.csv` files are input data for trait-associated diversification and all `.svg` and `.pdf` files in folder `figs`are figures from the output.

## biogeography_manuallycorrected_trimmedtree
BioGeoBEARS directory. `biogeo_astragalus.R` is the analysis script, `BIOGEOCODING.formatted.biogeobears.treematched.curated.newtrim23.txt` contains the geographic coding, and `astragalus.hybpiper.norogues.originalnames.nosections.optimize.calibrated.no.root.constraint.noplatenames.nodups.localitymatched.newtrim23.tre` is the precise tree used. `figs` contains plotted results, `AIC` contains tables of AIC results, and `Rdata_objects` contains object saves for steps that take a long time to compute.

## categorical_reconstructions
`new.astragalus.reconstructions.mar23.R` is the analysis script. `astragalus.hybpiper.norogues.originalnames.nosections.optimize.calibrated.no.root.constraint.noplatenames.nodups.forsoilrec.newtrim23.tre` is the precise tree used for analysis (the other tree in this directory shows the tree prior to name cleaning). `astragalus.aic.table.csv` contains AIC output. `figs` contains plotted figures.

## chromosome_reconstruction
`chromevol_astragalus.Rev` is the analysis script (RevBayes), `ancstatesrecon.R` is the plotting script (RevGadgets), `chrom.ordered.tsv` contains the chromosome data, `tree.for.chromevol.nex` is the tree used for the analysis, `ChromEvol_simple_anc_states_50K.log` is the log file, `ChromEvol_simple_final_50K.tree` is the output (NEXUS), and `figs` contains plotted output.

## continuous_ancestral_reconstructions
`niche_conservatism_ancestral_recon.r` is the analysis script (also includes niche conservatism tests), `traitDependent_functions.R` is a module from [https://github.com/macroevolution/fisse], *p*-values and other outputs are summarized in `niche_conserve.xlsx` and `MANOVA_results.xlsx`, and `figs` contains plotted figures.

## dated_trees
Four trees after dating in TreePL (`.tre` extension). Naming convention: "localitymatched" is a version of the tree with taxa excluded that have no occurrence records. "newtrim23" indicates a version of the tree with taxa dropped per the main text and summarized in `trimlist.Mar23.txt`. The full list of dropped taxa is given in `complete.trimlist.csv`. `astragalus_treepl_config.txt` is the analysis file for TreePL. 

## DR_rates
`tip_DR.csv` gives calculated DR rates; `figs` contains figures from outputs comparing DR to other data.


## environmental_data
CSVs representing environmental means (soil types in a separate directory).

## nquire
`ploidy.call.csv` contains calls of ploidy level by nQuire; `ploidy.call.original.csv` is the same but including NitFix plate and well numbers. `figs` contains plotted analysis results.

## occurrences
CSVs representing raw occurrences.

## soil_types
CSVs representing soil types, given both as full names and numerical coding (see `.xlsx` for number-name correspondence).

## tree_quality_control
Miscellaneous QC. `HybPiper_test_stats.xlsx` contains standard HybPiper quality control output. `astragalus_treeQC.R` is an analysis script for plotting QC data on a tree; `astragalus.wBranchLenghts.outliersremoved.renamed.tre` is the plotting tree; `all_gene_trees.bs10.output.quartet.tre` has quartet scores plotted on the branches; `qc_astragalus_rename.csv` is slightly reformatted QC data for plotting. `figs` contains output figures.
