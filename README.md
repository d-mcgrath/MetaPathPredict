# MetaPathPredict

The MetaPathPredict Python module utilizes PyTorch deep learning models to predict the presence or absence of KEGG metabolic modules in bacterial genomes recovered from environmental sequencing efforts.

<br>

## Installation

```
# Run the following commands from a terminal window:
pip install MetaPathPredict [not available yet]

# GitHub install of MetaPathPredict
[insert command here]
```

## Functions

The following functions can be implemented to run MetaPathPredict:

- `MetaPathPredict` parses one or more input KEGG Ortholog genome (MAG, SAG, isolate) gene annotation datasets (currently only bacterial genome data is supported). This function currently accepts as input the output files from KofamScan and DRAM gene annotation command line tools. Run either of these tools first and then use their output .tsv files as input for MetaPathPredict. A sample gene annotation output file from KofamScan is included in the repository, and is called `genome_annotations.tsv`. To test MetaPathPredict or see a sample input, download this file and use it as input.

- This function reconstructs KEGG modules within the input annotation dataset and predicts the presence or absence of incomplete or missing KEGG modules. To specify a specific module or modules, include the module identifier or identifiers as a comma-separated list to the argument `--modules`. 

- To view which KEGG modules MetaPathPredict has models for, run the following: `available_modules()`.

<br>

## Basic usage

```
# predict method for making KEGG module presence/absence predictions on input gene annotations
usage: MetaPathPredict [-h] --input INPUT --output OUTPUT --model-file MODEL_IN

options:
  -h, --help            show this help message and exit
  --input INPUT, -i INPUT
                        input file name [required]
  --annotation-format, -a ANNOTATION_FORMAT
                        annotation format [kofamscan, dram, koala; default: kofamscan]
  --output OUTPUT, -o OUTPUT
                        output file name [required; no default folder created]
  --model-file MODEL_IN, -m MODEL_IN
                        input model file names [default: models/ directory of the MetaPathPredict folder]
```

<br>

## Example usage with sample datasets

```
# with one KofamScan gene annotation dataset
MetaPathPredict -i kofamscan_annotation.tsv -a kofamscan -o /results/predictions.tsv

# with multiple KofamScan datasets
MetaPathPredict -i annotations/*.tsv -a kofamscan -o /results/predictions.tsv

# with one DRAM gene annotation dataset
MetaPathPredict -i dram_annotation.tsv -a dram -o /results/predictions.tsv

# with multiple DRAM datasets
MetaPathPredict -i annotations/*.tsv -a dram -o /results/predictions.tsv
```

<br>

## Developer usage

```
usage: MetaPathTrain [-h] --train-targets TRAIN_TARGETS --train-features TRAIN_FEATURES [--num-epochs NUM_EPOCHS] --model-out
                            MODEL_OUT [--use-gpu] [--num-cores NUM_CORES] [--num-hidden-layers NUM_HIDDEN_LAYERS]
                            [--hidden-nodes-per-layer HIDDEN_NODES_PER_LAYER] [--num-features NUM_FEATURES]

options:
  -h, --help            show this help message and exit
  --train-targets TRAIN_TARGETS
                        training targets file
  --train-features TRAIN_FEATURES
                        training features
  --num-epochs NUM_EPOCHS
                        number of epochs
  --model-out MODEL_OUT, -m MODEL_OUT
                        model file name output
  --use-gpu             use GPU if available
  --num-cores NUM_CORES
                        Number of cores for parallel processing

Neural Net parameters:
  --num-hidden-layers NUM_HIDDEN_LAYERS
                        number of hidden layers
  --hidden-nodes-per-layer HIDDEN_NODES_PER_LAYER
                        number of nodes in each hidden layer
  --num-features NUM_FEATURES
                        number of features to retain from training data
```
