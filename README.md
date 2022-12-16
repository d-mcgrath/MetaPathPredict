# MetaPredict

The MetaPredict R package utilizes machine learning models (stacked neural network/XGBoost ensemble and standalone neural network architectures) to predict the presence or absence of KEGG metabolic modules in bacterial genomes recovered from environmental sequencing efforts.

## Installation

``` r
# Run the following commands in an R session:

# GitHub install of recipeselectors package separately
devtools::install_github("stevenpawley/recipeselectors")

# GitHub install of MetaPredict
devtools::install_github("d-mcgrath/MetaPredict/MetaPredict")
```

#### Download MetaPredict's SQL database
MetaPredict requires an SQL database that contains its machine learning models in order to make predictions. This SQL database can be downloaded [here](https://zenodo.org/record/7419289).


## Functions

The following functions can be implemented to run MetaPredict:

- `read_data` parses one or more input KEGG Ortholog genome (MAG, SAG, isolate) gene annotation datasets (currently only bacterial genome data is supported). This function currently accepts as input the output files from KofamScan and DRAM gene annotation command line tools. Run either of these tools first and then use their output .tsv files as input for MetaPredict. A sample gene annotation output file from KofamScan is included in the repository, and is called `genome_annotations.tsv`. To test MetaPredict or see a sample input, download this file and use it as input.

- `metapredict` reconstructs KEGG modules within the input annotation dataset and predicts the presence or absence of incomplete or missing KEGG modules. Be sure to include the path to the downloaded SQL database with the `db_path` argument. To specify a specific module or modules, include the module identifier or identifiers as a character vector for the argument `module_vector`. 

- To view which KEGG modules MetaPredict has models for, run the following: `available_modules()`.

## Basic usage

``` r
input <- read_data(
  metadata_df = data.frame(
  filepath = 'path/to/genome/annotation/file',
  genome_name = 'test_genome'))

# make presence/absence predictions for KEGG modules M00001 and M00005
results <- metapredict(
  input,
  module_vector = c('M00001', 'M00005'),
  db_path = "/path/to/sql/database")

```
