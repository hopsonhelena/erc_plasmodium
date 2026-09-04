# ERC Plasmodium

This [R Shiny application](https://lefflerlab.chpc.utah.edu/erc_plasmodium/) accompanies [Hopson et al. (2026)](https://doi.org/10.1093/gbe/evag203) and enables exploration of evolutionary rate covariation (ERC) across the genomes of 22 *Plasmodium* species: 

- **Group of genes:** inspect pairwise ERC values, visualize a clustered
  heatmap, and test whether the group's mean ERC is greater than that of
  randomly sampled gene sets of the same size.
- **Single gene:** retrieve genes with the highest ERC values for a selected
  gene.
- **1:1 orthologs:** find members of an inferred ortholog group across selected
  *Plasmodium* species.

## Data

The source data are avaliable on Dryad:
[doi:10.5061/dryad.cfxpnvxmw](https://doi.org/10.5061/dryad.cfxpnvxmw).

See [`data/README.md`](data/README.md) for expected inputs.

## Installation

The project uses
[`renv`](https://rstudio.github.io/renv/) to record its R package environment.
From the project root, install `renv` (if needed) and restore the recorded package environment:

```r
install.packages("renv")
renv::restore()
```

The application uses these R packages: `shiny`, `DT`, `bslib`,
`ggplot2`, `writexl`, `dplyr`, `reshape2`, `readxl`, `colourpicker`,
`shinycssloaders`, and `openxlsx`. Their dependencies are captured in
`renv.lock`.

## Running the application
With files in `data/`, run: 

```r
shiny::runApp()
```

## Repository structure

```text
app.R                          Shiny user interface, server, and data loading
modules/gene_input_module.R    Gene-list input
modules/result_matrix_module.R Pairwise matrix and permutation-test output
modules/heatmap_module.R       Clustered ERC heatmap
modules/top_gene_clusters.R    Ranked single-gene results
modules/select_species.R       Ortholog lookup
modules/utils.R                Shared ERC functions
data/README.md                 Data documentation
```

## Citation
The underlying R functions were adapted from the
[nclark-lab/erc](https://github.com/nclark-lab/erc) project.

If you use this application please cite:

> Hopson HD, Omelianczyk RI, Ramirez A, Little JH, Clark N, Sigala PA,
> Leffler EM. Evolutionary rate covariation across malaria parasite species
> enables inference of protein interactions. *Genome Biology and Evolution*.
> 2026. <https://doi.org/10.1093/gbe/evag203>

## Contact

This application was created by the Leffler Lab at the University of Utah.
For more information, visit <https://lefflerlab.org/>.
