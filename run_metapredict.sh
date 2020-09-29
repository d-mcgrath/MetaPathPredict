#!/bin/bash

#
CLEAR='\033[0m'
RED='\033[0;31m'

function usage() {
  if [ -n "$1" ]; then
    echo -e "${RED}👉 $1${CLEAR}\n";
  fi
    echo "Usage: $0 [-w -working-dir] [-o output-dir] [-g genus-file] [-c cores] [-f input-files]"
    echo "  -w, --working-dir  : The directory containing input files [default current directory]"
    echo "  -o, --output-dir   : The directory for output files. Will be created in working directory if it does not exist [default metapredict_output]"
    echo "  -g, --genus-file   : Full path to a CSV file containing the genus of each input fasta file, in the same order as fasta files in the working directory"
    echo "  -c, --cores        : Number of CPU cores to use [default 1]"
    echo "  -f, --files        : One or more input genome fasta files"
    echo ""
    echo "Example with long flags: $0 --working-dir /path/to/input_files --output-dir /path/to/output_dir --genus-file /full/path/to/genus_file.csv --cores 1 --files *.fa"
    echo ""
    echo "Example with short flags: $0 -w /path/to/input_files -o /path/to/output_dir -g /full/path/to/genus_file.csv -c 1 -f *.fa"
    exit 1
}

#parse parameters
while [[ "$#" > 0 ]]; do case $1 in
  -w|--working-dir) WORKING_DIR="$2"; shift;shift;;
  -o|--output-dir) OUTPUT_DIR="$2"; shift;shift;;
  -g|--genus-file) GENUS_FILE="$2"; shift;shift;;
  -c|--cores) CORES="$2"; shift;shift;;
  -f|--files) FILES="$2"; shift;shift;;
  *) usage "Unknown parameter passed $1"; shift;shift;;
esac; done

#verify parameters
if [ -z "${WORKING_DIR}" ]; then WORKING_DIR=$(pwd); fi;
if [ -z "${OUTPUT_DIR}" ]; then OUTPUT_DIR=$(echo metapredict_output); fi;
if [ -z "${GENUS_FILE}" ]; then usage "Error: no genus file detected. A genus CSV file is required. Did you provide the full path to the genus file?"; fi;
if [ -z "${CORES}" ]; then CORES=$(echo 1); fi;
if [ -z "${FILES}" ]; then usage "Error: no input files detected. One or more input genome fasta files is required"; fi;

#load conda environment containing Prodigal, Kofamscan, and required R packages
module load anaconda/5.1
source activate metapredict-r-env
cd "${WORKING_DIR}"

#
if [ -d "${OUTPUT_DIR}" ]; then cd "${OUTPUT_DIR}"
else mkdir "${OUTPUT_DIR}"; cd "${OUTPUT_DIR}"; fi

#run prodigal using parallel processing
find "${WORKING_DIR}" -maxdepth 1 -name \""${FILES}"\" -printf "%f\n" | parallel --jobs "${CORES}" "echo Running Prodigal on {}; prodigal -i "${WORKING_DIR}"/{} -a {.}-genes.fa; printf "\nDone with %s\n\n" {}"

#run kofamscan
for FILE in $(ls *-genes.fa)
do echo Running Kofamscan on "${FILE}"...
PREFIX=$(echo "${FILE}" | sed -r 's/(.*)-genes.fa/\1/')
exec_annotation -o "${PREFIX}"-ko.tsv --cpu "${CORES}" -f detail-tsv -p ../profiles/ -k ../ko_list "${FILE}"
#cat "${PREFIX}"-ko.tsv | awk -F"\t" '{print $3,$6}' | sed '/^--/d' > "${PREFIX}"-ko-only.txt
#echo Parsing HMM hits and E-values into file "${PREFIX}"-ko-only.txt...
printf "Done with %s\n\n" "${FILE}"
done

#run metapredict
echo Starting MetaPredict on files: "$(ls *-ko.tsv | paste -s -d ,)"...
Rscript ../full_metapredict.R -p ./ -g ../"${GENUS_FILE}" -c "${CORES}"

#all done
