# multispecies_rna-seq

###This repository includes data and code to accompany:

###"The transcriptional response to low temperature is weakly conserved across the Enterobacteriaceae" by Hoang and Stoebel

###Now published in mSystems

###[https://doi.org/10.1128/msystems.00785-24]

## Outputs

The directory `outputs` includes:
- all of the figures used in the paper
- all of the outputs of from DEseq (in `/sig_tables`)
    - files of the form `{species name genotype}_significance_of_temp` compare 37°C vs 15°C samples of that strain
    - files of the form `{species name temperature}_significance_of_genotype` compare wt vs ∆rpoS samples of that species at that temperature.
- all of the outputs of GO enrichment analysis (in `/GO`) 
    - files of the form `enriched_GO_terms_in_{strain}` are the results of the GO enrichment of the DE genes between 37°C and 15°C in that strain.
    - files of the form `{strain}_GO_terms_for_topGO` are intermediate files used by data processing but may not be useful for others
- a pdf of expression of the genes for each GO term that was enriched in at least one species (`all_enriched_plots.pdf`)
- The growth rate parameters estimated for each strain, which appear in Table 1 in the paper.

## Data analysis files

### Four files do the analysis of RNA-seq data that are the majority of the data. 

- `basic_analysis.R` does the fundamental analysis of the RNA-seq data using the software DEseq. This output of the file is saved to `basic_analysis.RData` so that it doesn't need to be run every time. 

- `wt differentiall expression.Rmd` takes care of the analysis of RNA-seq data on the six wild-type species, including creating the upset plots and PCA plots in the paper.
- `GO enrichment.Rmd` does all of the GO enrichment analysis on the *wild-type* strains
- `wt vs rpoS data.Rmd` does all of the analysis on the ∆rpoS strains, including both the comparison of significant genes in the two species as well as the comparison of GO enrichment.

### Two files for other data analysis

- `analysis of all western blots.Rmd` analyses the RpoS western blots, creates figure S1
- `growth_rates.Rmd` analyses all growth rate data and produces the data for Table 1

### One file justifies analysis choices

- `quality control.Rmd` provides justification for dropping four samples

## Data

- RNA-seq count files are in `RNAseq_data`. (Raw data are in GEO and linked in the paper)
- Growth rate data are in `growth_rates`
- Western blot data are in `western_blots`
- The phylogeny inferred by ASTRAL is in `trees`
- The genome analysis output from xenoGI is in `xenoGI_outputs`
- The GO mapping from eggNOG-mapper is in `eggNOG-mapper`
- The data used to compare these RNA-seq data to those from Adams et al. 2023 is in `adamseteal_combined_deseq_analysis.csv`

## Package managagement

This project used `renv` to manage packages, so the directory `renv` and the file `renv.lock` will set those up.


