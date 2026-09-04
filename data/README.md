# Application data

The full ERC data are excluded from Git because the directory is
large and contains tens of thousands of generated files.

The application expects:

| Path | Purpose |
| --- | --- |
| `cluster_modified_to_cluster.rds` | ERC matrix indexed by ortholog cluster |
| `gene_to_cluster_long.rds` | Mapping among genes, species, and ortholog clusters |
| `per_gene_rds/<gene identifier>.rds` | Precomputed ranked ERC results |

## Data availability

The source dataset is archived on Dryad:
[doi:10.5061/dryad.cfxpnvxmw](https://doi.org/10.5061/dryad.cfxpnvxmw).

The Dryad files used to prepare the application data are:

| Dryad file | Application use |
| --- | --- |
| `Plasmodium_ftERC.txt` | Source for `cluster_modified_to_cluster.rds` |
| `Plasmodium_OGs_singlecopy.xlsx` | Source for `gene_to_cluster_long.rds` |
| `Species_name_key.xlsx` | Species names and abbreviations |