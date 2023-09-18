# HGT detection pipeline v.1

# Usage

## Add your own (sub)-taxonomy

Insert tutorial

### Installation

#### Dependencies

* snakemake - conda
* hgtector2 - pip
* taxonkit - conda
* blast suite - conda
* mafft - conda
* proteinortho - conda
* famsa - conda
* seqkit - conda
* trimal - conda
* fasttree - conda
* FastRoot - pip
* tidyverse R package - CRAN
* phytools -CRAN
* ggtree - Bioconductor
* ggtreeExtra - Bioconductor
* Biostrings - Bioconductor
* ComplexHeatmap - Bioconductor

You can use conda to create this environment, ideally use mamba:

```
conda create -c conda-forge -c bioconda -n hgt snakemake
conda activate hgt

# install hgtector
conda install -c conda-forge pyyaml pandas matplotlib scikit-learn bioconda::diamond
pip install git+https://github.com/qiyunlab/HGTector.git

# install other dependencies
conda install -c bioconda mafft trimal fasttree taxonkit proteinortho famsa seqkit

# install R if you need
conda install -c conda-forge r-base
conda install -c conda-forge -c bioconda r-tidyverse r-optparse r-magic r-phytools bioconductor-ggtree bioconductor-ggtreeextra bioconductor-complexheatmap 
bioconductor-biostrings

# install fastroot
python3 -m pip install FastRoot
fr=$(which FastRoot.py)
dos2unix $fr
```


# TODOs

* mask proteins and then do dataframe with protein length, %masked proportion to eventually filter downstream results
* add link to deps
* what to do if less than 2 seqs in plot_Tree?
* more elegant way to get sequences to do the tree 
* add QC plots
* what to do with foldseek?
* matreex plots?
* mark HGT tree by topological rule
* get df of tree with all different properties so user can filter
* when deleting hits from the same species you may remove paralogs that are not identified by proteinortho but still may be close in the tree
* annotate domains and protein functions
* maybe one strategy would be to compute big trees in a fast way then reduce taxonomic redundancy and compute a more accurate tree
* remove fastroot dependency

* option to define "close" outgroup and compute other indexes
* identify potential donors from blast results (recurring weird taxons)
* one day 3 modalities
    * --interdomain: easy mode for proka2euka lgt
    * --donors: user-defined putative donors
    * --reconcile: with small curated database and species tree run reconciliation
