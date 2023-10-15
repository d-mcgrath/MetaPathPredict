# MetaPathPredict

The MetaPathPredict Python module utilizes deep learning models to predict the presence or absence of KEGG metabolic modules in bacterial genomes recovered from environmental sequencing efforts.

## Installation

To run MetaPathPredict, download this repository and install it as a Python module (see download and installation instructions below):


### GitHub install:

NOTE: [Conda](https://docs.conda.io/en/latest/) is required for this installation.

1. Open a Terminal/Command Prompt window and run the following command to download the
GitHub repository to the desired location (note: change your current working directory first
to the desired download location, e.g., `~/Downloads` on MacOS):
`git clone https://github.com/d-mcgrath/MetaPathPredict.git`

    1. NOTE: You can also download the repository zip file from GitHub

2. In a Terminal/Command Prompt window, run the following commands from the parent directory the MetaPathPredict repository was cloned to:
```bash
conda create -n MetaPathPredict python=3.10.6 scikit-learn=1.3.0 tensorflow=2.10.0 numpy=1.23.4 pandas=1.5.2 
keras=2.10.0
```
NOTE: You will be prompted (y/n) to confirm creating this conda environment. Now activate it:

```bash
conda activate MetaPathPredict
```

3. pip install pyxet:
```bash
pip install pyxet
```

4. Install Pytorch. To do this, you must pick which [version](https://pytorch.org/get-started/locally/) to install from conda.

To install the Linux CPU version, use the following:

```bash
conda install pytorch torchvision torchaudio cpuonly -c pytorch
```

To install the Mac default version, use the following:

```bash
conda install pytorch::pytorch torchvision torchaudio -c pytorch
```

To install the Windows CPU version, use the following:

```bash
conda install pytorch torchvision torchaudio cpuonly -c pytorch
```

5. Once complete, pip install MetaPathPredict:
```bash
pip install MetaPathPredict/package
```

6. Download MetaPathPredict's models by running the following command:
```bash
DownloadModels
```

Note: MetaPathPredict is now installed in the `MetaPathPredict` conda environment. Activate the conda environment prior to any use of MetaPathPredict.

### XetHub install:
Follow the instructions from the MetaPathPredict [XetHub](https://xethub.com/dgellermcgrath/MetaPathPredict) repository.

### pip install:
[not available yet]

<br>

## Functions

The following functions can be implemented to run MetaPathPredict on the command line:

- `MetaPathPredict` parses one or more input KEGG Ortholog gene annotation datasets (currently only bacterial genome data is supported) and predicts the presence or absence of [KEGG Modules](https://www.genome.jp/kegg/module.html). This function takes as input the .tsv output files from the [KofamScan](https://github.com/takaram/kofam_scan) and [DRAM](https://github.com/WrightonLabCSU/DRAM) gene annotation tools as well as the KEGG KOALA online annotation platforms [blastKOALA](https://www.kegg.jp/blastkoala/), [ghostKOALA](https://www.kegg.jp/ghostkoala/), and [kofamKOALA](https://www.genome.jp/tools/kofamkoala/). Run any of these tools first and then use one or more of their output .tsv files as input to MetaPathPredict.
    - A single file or multiple space-separated files can be specified to the `--input` parameter, or use a wildcard (e.g., /results/*.tsv). Include full or relative paths to the input file(s). A sample of each annotation file format that MetaPathPredict can process is included in this repository in the [annotatation_examples](annotatation_examples) folder. The sample annotation files in [annotatation_examples](annotatation_examples) can optionally be used as input to test the installation.
    - The format of the gene annotation files (kofamscan, kofamkoala, dram, or koala) that is used as input must be specified with the `--annotation-format` parameter. Currently, only one input type can be specified at a time.
    - The full or relative path to the desired destination for MetaPathPredict's output .tsv file must be specified, as well as a name for the file. The output file path and name can be specified using the `--output` parameter. By default, MetaPathPredict does not create any default output directory nor does the output have a default file name.
    - To specify a specific KEGG module or modules to reconstruct and predict, include the module identifier (e.g., M00001) or identifiers as a space-separated list to the argument `--kegg-modules`.

- To view which KEGG modules MetaPathPredict can reconstruct and make predictions for, run the following on the command line: `MetaPathModules`.

<br>

## Basic usage

```
# predict method for making KEGG module presence/absence predictions on input gene annotations

usage: MetaPathPredict [-h] --input INPUT [INPUT ...] --annotation-format ANNOTATION_FORMAT
                       [--kegg-modules KEGG_MODULES [KEGG_MODULES ...]] --output OUTPUT

options:
  -h, --help            show this help message and exit
  --input INPUT [INPUT ...], -i INPUT [INPUT ...]
                        input file path(s) and name(s) [required]
  --annotation-format ANNOTATION_FORMAT, -a ANNOTATION_FORMAT
                        annotation format (kofamscan, kofamkoala, dram, or koala) [default:
                        kofamscan]
  --kegg-modules KEGG_MODULES [KEGG_MODULES ...], -k KEGG_MODULES [KEGG_MODULES ...]
                        KEGG modules to predict [default: MetaPathPredict KEGG modules]
  --output OUTPUT, -o OUTPUT
                        output file path and name [required]
```

<br>

## Examples with sample datasets

```
# One KofamScan gene annotation dataset
MetaPathPredict -i /path/to/kofamscan_annotations_1.tsv -a kofamscan -o /results/predictions.tsv

# Three KofamScan gene annotation datasets, with predictions for modules M00001 and M00003
MetaPathPredict \
-i kofamscan_annotations_1.tsv kofamscan_annotations_2.tsv kofamscan_annotations_3.tsv \
-a kofamscan \
-k M00001 M00003 \
-o /results/predictions.tsv

# Multiple KofamScan datasets in a directory
MetaPathPredict -i annotations/*.tsv -a kofamscan -o /results/predictions.tsv

# One DRAM gene annotation dataset
MetaPathPredict -i dram_annotation.tsv -a dram -o /results/predictions.tsv

# Multiple DRAM datasets in a directory
MetaPathPredict -i annotations/*.tsv -a dram -o /results/predictions.tsv
```

<br>

## Understanding the output

The output of running `MetaPathPredict` is a table. The first column, `file`, displays the full file name of each input gene annotation file. The remaining columns give the class predictions (module present = 1; module absent = 0) of KEGG modules. Each KEGG module occupies a single column in the table and is labelled by its module identifier. See a sample output below of four KEGG module predictions for three input annotation files:

| file                                 | M00001 | M00002 | M00003 | M00004 |
|--------------------------------------|--------|--------|--------|--------|
| /path/to/kofamscan_annotations_1.tsv | 1      | 1      | 0      | 1      |
| /path/to/kofamscan_annotations_2.tsv | 0      | 1      | 0      | 0      |
| /path/to/kofamscan_annotations_3.tsv | 1      | 0      | 0      | 0      |

<br>

## Developer usage

```
# training method for MetaPathPredict's internal models

usage: MetaPathTrain [-h] --train-targets TRAIN_TARGETS --train-features TRAIN_FEATURES
                     [--num-epochs NUM_EPOCHS] --model-out MODEL_OUT [--use-gpu]
                     [--num-cores NUM_CORES] [--num-hidden-layers NUM_HIDDEN_LAYERS]
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
