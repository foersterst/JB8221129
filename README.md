# Broad-scale climatic gradients drive multiple facets of scorpion beta diversity in northeastern Brazil

Please make sure to cite the following paper if you use any data from this repository:

Foerster, S. Í. A., & Lira, A. F. A. (2026). Broad-Scale Climatic Gradients Drive Multiple Facets of Scorpion Beta Diversity in Northeastern Brazil. *Journal of Biogeography*, 53(6), e70271. <https://doi.org/10.1111/jbi.70271>

## File Description

### R Scripts

- `scripts/data-analysis.R` contains the R code to reproduce all results in the manuscript.
- `scripts/fit-beta-div.R` custom R function that takes a community matrix, site-level environmental variables, spatial descriptors (MEMs), and an optional phylogenetic or functional tree to calculate beta-diversity metrics. The function quantifies total beta-diversity and its relative components, including species replacement and richness differences, performs variation partitioning on beta-diversity matrices (e.g., replacement and richness difference), and identifies spatial descriptors (MEMs) retained through forward selection. When no tree is provided, the function calculates taxonomic beta-diversity.

### Data Files

- `data/scorp-comm.xlsx` community matrix with species abundances per site. Also contains a few info about the sites, such as forest type (Caatinga or Atlantic Forest), sampling year, longitude and latitude.
- `data/scorp-traits.csv` species-level trait data (species means), including log<sub>10</sub> carapace length, chela measurements (length, width and height) relative to carapace length, and functional group.
- `data/site-envs.csv` site-level environmental variables used in the study.
- `data/scorp-pruned-tree.tre` time-calibrated molecular phylogeny from Foerster ([2026](https://doi.org/10.1093/evolut/qpag097)), pruned to the sampled species.
- `data/results.xlsx` beta-diversity and proportion of beta-diversity components (used for plotting only).
- `shp/` shapefiles used for plotting.

### Contact

I welcome inquiries and opportunities for future collaboration. Contact details and information about my research are available at [https://foersterst.github.io](https://foersterst.github.io/)
