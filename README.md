# MetaPathPredict

The MetaPathPredict Python module utilizes deep learning models to predict the presence or absence of KEGG metabolic modules in bacterial genomes recovered from environmental sequencing efforts.

## Installation

```
# Run the following commands from a terminal window:
pip install MetaPathPredict [not available yet]

# GitHub install of MetaPathPredict
[insert command here]
```

## Functions

The following functions can be implemented to run MetaPathPredict:

- `MetaPathPredictPredict` parses one or more input KEGG Ortholog genome (MAG, SAG, isolate) gene annotation datasets (currently only bacterial genome data is supported). This function currently accepts as input the output files from KofamScan and DRAM gene annotation command line tools. Run either of these tools first and then use their output .tsv files as input for MetaPathPredict. A sample gene annotation output file from KofamScan is included in the repository, and is called `genome_annotations.tsv`. To test MetaPathPredict or see a sample input, download this file and use it as input.

- This function reconstructs KEGG modules within the input annotation dataset and predicts the presence or absence of incomplete or missing KEGG modules. To specify a specific module or modules, include the module identifier or identifiers as a comma-separated list to the argument `--modules`. 

- To view which KEGG modules MetaPathPredict has models for, run the following: `available_modules()`.

## Basic usage

```
# predict method for making KEGG module presence/absence predictions on input gene annotations
usage: MetaPathPredictPredict [-h] --input INPUT --output OUTPUT --model-file MODEL_IN

options:
  -h, --help            show this help message and exit
  --input INPUT, -i INPUT
                        input file
  --output OUTPUT, -o OUTPUT
                        output file
  --model-file MODEL_IN, -m MODEL_IN
                        input model file name




# for development or advanced users: training method for training deep learning model to make KEGG module presence/absence predictions

usage: MetaPathPredictTrain [-h] --train-targets TRAIN_TARGETS --train-features TRAIN_FEATURES [--num-epochs NUM_EPOCHS] --model-out MODEL_OUT [--use-gpu]
                            [--num-hidden-layers NUM_HIDDEN_LAYERS] [--hidden-nodes-per-layer HIDDEN_NODES_PER_LAYER]

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

Neural Net parameters:
  --num-hidden-layers NUM_HIDDEN_LAYERS
                        number of hidden layers
  --hidden-nodes-per-layer HIDDEN_NODES_PER_LAYER
                        number of nodes in each hidden layer
```
