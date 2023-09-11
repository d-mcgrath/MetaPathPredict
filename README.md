# MetaPathPredict

The MetaPathPredict Python module utilizes PyTorch deep learning models to predict the presence or absence of KEGG metabolic modules in bacterial genomes recovered from environmental sequencing efforts.

## Installation

```
# Run the following commands from a terminal window:
pip install MetaPathPredict [not available yet]

# GitHub install of MetaPathPredict
[insert command here]
```

## Functions

The following functions can be implemented to run MetaPathPredict:

- `MetaPathPredict` parses one or more input KEGG Ortholog gene annotation datasets (currently only bacterial genome data is supported). This function currently accepts as input the output files from KofamScan, DRAM, BlastKoala, and GhostKoala gene annotation tools. Run any of these tools first and then use one or more of their output gene annotation files as input to MetaPathPredict. A sample of each annotation output file form is included in this repository. To test MetaPathPredict or see a sample input file for use with this tool, download the repository, install MetaPathPredict, and run predictions on one or more of the sample gene annotation files.

- The `MetaPathPredict` function reconstructs KEGG modules within the input annotation dataset and predicts the presence or absence of incomplete or missing KEGG modules. To specify a specific KEGG module or modules to reconstruct and predict, include the module identifier (e.g., M00001) or identifiers as a comma-separated list to the argument `--kegg-modules`. 

- To view which KEGG modules MetaPathPredict can reconstruct and make predictions for, run the following command: `AvailableModules`.

<br>

## Basic usage

```
# predict method for making KEGG module presence/absence predictions on input gene annotations
usage: MetaPathPredictPredict [-h] --input INPUT [INPUT ...] --annotation-format
                              ANNOTATION_FORMAT --output OUTPUT --model-files MODEL_IN
                              [MODEL_IN ...] --scaler-files SCALER_IN [SCALER_IN ...]
                              [--kegg-modules KEGG_MODULES [KEGG_MODULES ...]]

options:
  -h, --help            show this help message and exit
  --input INPUT [INPUT ...], -i INPUT [INPUT ...]
                        input file path(s) [required]
  --annotation-format ANNOTATION_FORMAT, -a ANNOTATION_FORMAT
                        annotation format [kofamscan, dram, koala; default: kofamscan]
  --output OUTPUT, -o OUTPUT
                        output file path and name [required; no default output
                        file name or folder created]
  --model-files MODEL_IN [MODEL_IN ...], -m MODEL_IN [MODEL_IN ...]
                        input model files location [default: MetaPathPredict folder]
  --scaler-files SCALER_IN [SCALER_IN ...], -s SCALER_IN [SCALER_IN ...]
                        input scaler files location [default: MetaPathPredict folder]
  --kegg-modules KEGG_MODULES [KEGG_MODULES ...], -k KEGG_MODULES [KEGG_MODULES ...]
                        KEGG modules to predict [default: MetaPathPredict KEGG
                        modules]
```

<br>

## Examples with sample datasets

```
# One KofamScan gene annotation dataset
MetaPathPredict -i kofamscan_annotations_1.tsv -a kofamscan -o /results/predictions.tsv

# Three KofamScan gene annotation datasets, with predictions for modules M00001 and M00003
MetaPathPredict \
-i kofamscan_annotations_1.tsv kofamscan_annotations_2.tsv kofamscan_annotations_3.tsv \
-a kofamscan \
-o /results/predictions.tsv

# Multiple KofamScan datasets in a directory
MetaPathPredict -i annotations/*.tsv -a kofamscan -o /results/predictions.tsv

# One DRAM gene annotation dataset
MetaPathPredict -i dram_annotation.tsv -a dram -o /results/predictions.tsv

# Multiple DRAM datasets in a directory
MetaPathPredict -i annotations/*.tsv -a dram -o /results/predictions.tsv
```

<br>

## Developer usage

```
usage: MetaPathPredictTrain [-h] --train-targets TRAIN_TARGETS --train-features
                            TRAIN_FEATURES [--num-epochs NUM_EPOCHS] --model-out
                            MODEL_OUT [--use-gpu] [--num-cores NUM_CORES]
                            [--num-hidden-layers NUM_HIDDEN_LAYERS]
                            [--hidden-nodes-per-layer HIDDEN_NODES_PER_LAYER]
                            [--num-features NUM_FEATURES] [--threshold THRESHOLD]

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
  --threshold THRESHOLD
                        threshold for SelectKBest feature selection
```
