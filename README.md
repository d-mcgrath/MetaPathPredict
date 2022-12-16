# MetaPredict

The MetaPredict R package utilizes machine learning models (stacked neural network/XGBoost ensemble and standalone neural network architectures) to predict the presence or absence of KEGG metabolic modules in bacterial genomes recovered from environmental sequencing efforts.

## Installation

``` r
# GitHub install of recipeselectors package separately
devtools::install_github("stevenpawley/recipeselectors")

# GitHub install of MetaPredict
devtools::install_github("d-mcgrath/MetaPredict/MetaPredict")
```

#### Download MetaPredict's SQL database
The SQL database containing MetaPredict's models can be downloaded [here](https://zenodo.org/record/7419289).


## Functions

The following functions can be implemented to run MetaPredict:

- `read_data` parses one or more input KEGG Ortholog genome (MAG, SAG, isolate) gene annotation datasets (currently only bacterial genome data is supported). This function currently accepts as input the output files from KofamScan and DRAM gene annotation command line tools. Run either of these tools first and then use their output .tsv files as input for MetaPredict. A sample gene annotation output file from KofamScan is included in the repository, and is called `genome_annotations.tsv`. To test MetaPredict or see a sample input, download this file and use it as input.

- `metapredict` reconstructs KEGG modules within the input annotation dataset and predicts the presence or absence of incomplete or missing KEGG modules. Be sure to include the path to the downloaded SQL database with the `db_path` argument.

## Basic usage

``` r
input <- read_data(
  metadata_df = data.frame(
  filepath = 'path/to/genome/annotation/file',
  genome_name = 'test_genome'))

results <- metapredict(input, db_path = "/path/to/sql/database")

```
