# Broad-scale climatic gradients drive multiple facets of scorpion beta diversity in northeastern Brazil

This repository accompanies the article Broad-scale climatic gradients drive multiple facets of scorpion beta diversity in northeastern Brazil, published on *Journal of Biogeography*.

## File Description

### R Scripts

- `scripts/data-analysis.R` contains the R code to reproduce all results in the manuscript.
- `scripts/fit-beta-div.R` custom R function that takes a community matrix, site-level environmental variables, spatial descriptors (MEMs), and a tree (phylogenetic or functional) and calculates 1) beta-diversity metrics, including relative contributions of species replacement and richness differences to total beta-diversity, variation partitioning of beta-diversity matrices (e.g., replacement, richness difference) as well MEMs retained by forward selection. Use this function without a tree if you want to calculate taxonomic beta-diversity.

### Data Files

- `data/scorp-comm.xlsx` community matrix with species abundances per site. Also contains a few info about the sites, such as forest type (Caatinga or Atlantic Forest), sampling year, longitude and latitude.
- `data/scorp-traits.csv` species-level trait data (species means), including log<sub>10</sub> carapace length, chela measurements (length, width and height) relative to carapace length, and functional group.
- `data/site-envs.csv` site-level environmental variables used in the study.
- `data/scorp-pruned-tree.tre` time-calibrated molecular phylogeny from Foerster ([2026](https://doi.org/10.1093/evolut/qpag097)), pruned to the sampled species.
- `data/results.xlsx` beta-diversity and proportion of beta-diversity components (used for plotting only).
- `shp/` shapefiles used for plotting.

### Data usage notice

The code in this repository is released under the MIT License. However, the data files are **not** covered by the MIT License.

The data are provided for transparency and reproducibility purposes only. Users wishing to reuse, redistribute, or publish analyses based on these data must first obtain permission from the author.

Please contact: stenioit\@gmail.com
