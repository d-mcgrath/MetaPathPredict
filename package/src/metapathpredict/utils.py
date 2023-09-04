import csv
import pandas as pd
import re


def read_kofamscan_detailed_tsv(files):
  """Reads in multiple .tsv files, each with columns: 0: "surpassed_threshold", 
  1: 'gene_name', 2: "k_number", 3: "adaptive_threshold", 4: "score",
  5: "evalue", 6: "definition". Keeps only rows where "surpassed_threshold" is 
  equal to "*". When there are duplicate values in "gene name", keeps the 
  row containing the highest value in the "score" column. If column "gene name" 
  contains multiple rows with the same maximum value, calculates the 
  score-to-adaptive-threshold ratio, and picks the annotation with the highest
  ratio.

  Args:
    files: A list of .tsv file names.

  Returns:
    A list of lists, where each inner list is the annotation data from one file.
  """
  kofam_list = []
  
  for file in files:
    lines = []
    with open(tsv, "rb") as f:
      for row in f:
        if row.decode().split("\t")[0] == "*":
          lines.append(row.decode().split("\t"))
    
    data = pd.DataFrame(lines)
    data.rename(columns={0: "surpassed_threshold", 1: 'gene_name',
    2: "k_number", 3: "adaptive_threshold", 4: "score",
    5: "evalue", 6: "definition"}, inplace=True)
  
    data[["adaptive_threshold", "score", "evalue"]] = data[["adaptive_threshold", "score", "evalue"]].apply(pd.to_numeric, axis=1)
    data = data.groupby("gene_name").apply(lambda group: group.loc[group["score"] == group["score"].max()]).reset_index(level=0, drop=True)
  
    data["group_size"] = data.groupby(["gene_name"]).transform("size")
  
    if data["group_size"].max() > 1:
        n_genes = (data[['gene_name', 'group_size']].drop_duplicates()['group_size'] > 1).sum()
        print(f"""{n_genes} gene(s) contained multiple annotations that surpassed the adaptive threshold. 
        Picking the annotation with the highest score-to-adaptive_threshold ratio for these genes.""")
    
        data["ratio"] = data["score"] / data["adaptive_threshold"]
        data = data.groupby("gene_name").apply(lambda group: group.loc[group["ratio"] == group["ratio"].max()]).reset_index(level=0, drop=True)
  
        data = data.drop(["ratio"], axis=1)
  
    data = data[["gene_identifier", "k_number", "definition"]]
    
    kofam_list.append(data)
    
  return data





def read_dram_annotation_tsv(files):
  """Reads in multiple DRAM annotation.tsv files, keeping the "gene_identifer" 
  as column 0, "k_number"" as column 1, and "definition" as column 2. Keeps 
  only rows where a gene had a KEGG Ortholog annotation.

  Args:
    files: A list of .tsv file names.

  Returns:
    A list of lists, where each inner list is the annotation data from one file.
  """
  pattern = "K[0-9]{5}"
  dram_list = []
  
  for file in files:
    lines = []
    with open(file, "rb") as f:
      for row in f:
        if re.match(pattern, row.decode().split("\t")[8]):
          lines.append(row.decode().split("\t"))
    
    data = pd.DataFrame(lines)[[0,8,9]]
    data.rename(columns={0: "gene_identifer", 8: 'k_number',
    9: "definition"}, inplace=True)
    
    kofam_list.append(data)
    
  return dram_list




def read_koala_tsv(files):
  """Reads in multiple blastKoala or ghostKoala .tsv files, keeping the 
  "gene_identifer" as column 0, "k_number"" as column 1, and "definition" as 
  column 2. Keeps only rows where a gene had a KEGG Ortholog annotation.

  Args:
    files: A list of .tsv file names.

  Returns:
    A list of lists, where each inner list is the annotation data from one file.
  """
  pattern = "K[0-9]{5}"
  dram_list = []
  
  for file in files:
    lines = []
    with open(file, "rb") as f:
      for row in f:
        if re.match(pattern, row.decode().split("\t")[0]):
          lines.append(row.decode().split("\t"))
    
    data = pd.DataFrame(lines)[[0,1,2]]
    data.rename(columns={0: "gene_identifer", 1: 'k_number',
    2: "definition"}, inplace=True)
    
    kofam_list.append(data)
    
  return dram_list








def create_feature_matrix(annotation_list):
  """Converts as list of annotations into a Pandas feature DataFrame.

  Args:
    x: A list of DataFrames.

  Returns:
    A Pandas DataFrame.
  """
  data = pd.DataFrame()
  for df in x:
    df["count"] = 1
    data = data.append(df)

  data = data.groupby(["file_name", "k_number"]).agg(count=("count", "sum")).unstack()
  data.fillna(0, inplace=True)

  return data



