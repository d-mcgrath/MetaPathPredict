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
    with open(file, "rb") as f:
      for row in f:
        if row.decode().split("\t")[0] == "*":
          lines.append(row.decode().split("\t"))
    
    data = pd.DataFrame(lines)
    data.rename(columns={0: "surpassed_threshold", 1: 'gene_identifier',
    2: "k_number", 3: "adaptive_threshold", 4: "score",
    5: "evalue", 6: "definition"}, inplace=True)
  
    data[["adaptive_threshold", "score", "evalue"]] = data[["adaptive_threshold", "score", "evalue"]].apply(pd.to_numeric, axis = 1)
    data = data.groupby("gene_identifier").apply(lambda group: group.loc[group["score"] == group["score"].max()]).reset_index(level = 0, drop = True)
  
    data["group_size"] = data.groupby(["gene_identifier"]).transform("size")
  
    if data["group_size"].max() > 1:
        n_genes = (data[['gene_identifier', 'group_size']].drop_duplicates()['group_size'] > 1).sum()
        print(f"""{n_genes} gene(s) contained multiple annotations that surpassed the adaptive threshold. 
        Picking the annotation with the highest score-to-adaptive_threshold ratio for these genes.""")
    
        data["ratio"] = data["score"] / data["adaptive_threshold"]
        data = data.groupby("gene_identifier").apply(lambda group: group.loc[group["ratio"] == group["ratio"].max()]).reset_index(level = 0, drop = True)
  
        data = data.drop(["ratio"], axis = 1)
  
    data["file_name"] = file
    data = data[["file_name", "gene_identifier", "k_number", "definition"]]
    
    kofam_list.append(data)
    
  return kofam_list





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
    data["file_name"] = file
    data = data[["file_name", "gene_identifier", "k_number", "definition"]]

    
    dram_list.append(data)
    
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
  koala_list = []
  
  for file in files:
    lines = []
    with open(file, "rb") as f:
      for row in f:
        if re.match(pattern, row.decode().split("\t")[0]):
          lines.append(row.decode().split("\t"))
    
    data = pd.DataFrame(lines)[[0,1,2]]
    data.rename(columns={0: "gene_identifer", 1: 'k_number',
    2: "definition"}, inplace=True)
    data["file_name"] = file
    data = data[["file_name", "gene_identifier", "k_number", "definition"]]

    
    koala_list.append(data)
    
  return koala_list





def create_feature_df(annotation_list):
  """Converts as list of annotations into a Pandas feature DataFrame.

  Args:
    annotation_list: A list of DataFrames.

  Returns:
    A Pandas DataFrame.
  """
  data = pd.DataFrame()
  for df in annotation_list:
    df["count"] = 1
    data = pd.concat([data, df], axis = 0)
  
  data = data.groupby(["file_name", "k_number"]).agg(count=("count", "sum")).reset_index().pivot_table(
  index = "file_name",
  columns = "k_number",
  values = "count",
  aggfunc = "first")
  
  data = data.replace(np.NaN, 0)
  data = data.where(data <= 1, 1)

  return data





def check_feature_columns(required_columns, feature_df):
  """Checks that all required columns are present for both of MetaPathPredict's models.

  Args:
    required_columns: A list of all required column names.
    feature_df: A DataFrame containing all predictor columns for both of MetaPathPredict's models.

  Returns:
    A Pandas DataFrame.
  """

  for column in required_columns
    if column not in feature_df.columns
      feature_df[column] = 0
  
  for column in feature_df.columns
    if column not in required_columns
      feature_df.drop([column], axis = 1)

  return feature_df





def select_model_features(required_columns, feature_df):
  """Selects all required columns for the specified MetaPathPredict model (either model #1 or model #2).

  Args:
    required_columns: A list of all required column names.
    feature_df: A DataFrame containing all predictor columns for both of MetaPathPredict's models.

  Returns:
    A Pandas DataFrame.
  """

  feature_df = feature_df[required_columns]

  return feature_df






## testing the functions

# fs = ["/Users/davidgeller-mcgrath/Downloads/MAG_172-ko2.tsv", 
# "/Users/davidgeller-mcgrath/Downloads/MAG_130-ko.tsv", 
# "/Users/davidgeller-mcgrath/Downloads/MAG_172-kofamscan.tsv"]
# x = read_kofamscan_detailed_tsv(fs)


# g = create_feature_matrix(x)
# print(g. iloc[:,:10])
# g.to_csv('/Users/davidgeller-mcgrath/Downloads/sample.tsv', sep="\t", index = False)







