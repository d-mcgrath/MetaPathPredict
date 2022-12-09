# MetaPredict

The MetaPredict R package utilizes machine learning models (stacked neural network/XGBoost ensemble and standalone neural network architectures) to predict the presence or absence of KEGG metabolic modules in bacterial genomes recovered from environmental sequencing efforts.

## Installation

``` r
# GitHub install of recipeselectors package separately
devtools::install_github("stevenpawley/recipeselectors")

# GitHub install of MetaPredict
devtools::install_github("d-mcgrath/MetaPredict/MetaPredict")
```

## Functions

The following functions can be implemented to run MetaPredict:

- `read_data` parses one or more input KEGG Ortholog genome (MAG, SAG, isolate) gene annotation datasets (currently only bacterial genome data is supported)

- `metapredict` reconstructs KEGG modules within the input annotation dataset and predicts the presence or absence of incomplete or missing KEGG modules.

## Basic usage

``` r
input <- read_data(
  metadata_df = data.frame(
  filepath = 'path/to/genome/annotation/file',
  genome_name = 'test_genome'))

results <- metapredict(input)

```
