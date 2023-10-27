import csv
import re
import gzip
import numpy as np
import pandas as pd


class InputData:
  
    """Data parsing functions of input data"""


    def __init__(self, files, annotations = []):
      self.files = files
      self.annotations = annotations

    def read_kofamscan_detailed_tsv(self):
      """Reads in multiple .tsv files, each with columns: 0: "surpassed_threshold", 
      1: 'gene_name', 2: "k_number", 3: "adaptive_threshold", 4: "score",
      5: "evalue", 6: "definition". Keeps only rows where "surpassed_threshold" is 
      equal to "*". When there are duplicate values in "gene name", keeps the 
      row containing the highest value in the "score" column. If column "gene name" 
      contains multiple rows with the same maximum value, calculates the 
      score-to-adaptive-threshold ratio, and picks the annotation with the highest
      ratio.
    
      Returns:
        A list of lists, where each inner list is the annotation data from one file.
      """

      if type(self.files) is str:
        self.files = [self.files]
      
      for file in self.files:
        lines = []
        
        if file.endswith(".gz"):
          with gzip.open(file, "rb") as f:
            for row in f:
              if row.decode().split("\t")[0] == "*":
                lines.append(row.decode().split("\t"))
        else:
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
        
        self.annotations.append(data)
        
        
        
        
    def read_kofamkoala(self):
      """Reads in multiple .tsv files, each with columns: 0: "gene_identifier", 
      1: 'k_number', 2: "adaptive_threshold", 3: "score", 4: "evalue", 
      5: "definition", 6: "definition_2". Keeps only rows where 
      "surpassed_threshold" is equal to "*". When there are duplicate values in 
      "gene name", keeps the row containing the highest value in the "score" 
      column. If column "gene name" contains multiple rows with the same maximum 
      value, calculates the score-to-adaptive-threshold ratio, and picks the 
      annotation with the highest ratio.
    
      Returns:
        A list of lists, where each inner list is the annotation data from one file.
      """

      if type(self.files) is str:
        self.files = [self.files]
      
      for file in self.files:
        lines = []
        
        if file.endswith(".gz"):
          with gzip.open(file, "rb") as f:
             for row in f:
              if row.decode().split("\t")[0] == "gene":
                continue
              elif row.decode().split("\t")[3] == "-":
                continue
              elif row.decode().split("\t")[2] == "-":
                if float(row.decode().split("\t")[4]) <= 1e-50:
                  lines.append(row.decode().split("\t"))
                else:
                  continue
              else:
                if float(row.decode().split("\t")[3]) > float(row.decode().split("\t")[2]):
                  lines.append(row.decode().split("\t"))
        else:
          with open(file, "rb") as f:
            for row in f:
              if row.decode().split("\t")[0] == "gene":
                continue
              elif row.decode().split("\t")[3] == "-":
                continue
              elif row.decode().split("\t")[2] == "-":
                if float(row.decode().split("\t")[4]) <= 1e-50:
                  lines.append(row.decode().split("\t"))
                else:
                  continue
              else:
                if float(row.decode().split("\t")[3]) > float(row.decode().split("\t")[2]):
                  lines.append(row.decode().split("\t"))
        
        data = pd.DataFrame(lines)
        data.rename(columns={0: "gene_identifier", 1: 'k_number',
        2: "adaptive_threshold", 3: "score", 4: "evalue",
        5: "definition", 6: "definition_2"}, inplace=True)
        
        data.loc[data["adaptive_threshold"] == "-", "adaptive_threshold"] = 1
              
        data[["adaptive_threshold", "score", "evalue"]] = data[["adaptive_threshold", "score", "evalue"]].apply(pd.to_numeric, axis = 1)
        data = data.groupby("gene_identifier").apply(lambda group: group.loc[group["score"] == group["score"].max()]).reset_index(level = 0, drop = True)
              
        data["group_size"] = data.groupby(["gene_identifier"]).transform("size")
              
        if data["group_size"].max() > 1:
          n_genes = (data[['gene_identifier', 'group_size']].drop_duplicates()['group_size'] > 1).sum()
          print(f"""{n_genes} gene(s) contained multiple annotations that surpassed the adaptive threshold. 
          Picking the annotation with the highest score-to-adaptive_threshold ratio for these genes.""")
                
        data["ratio"] = data["score"] / data["adaptive_threshold"]
        data = data.groupby("gene_identifier", group_keys = False).apply(lambda group: group.loc[group["ratio"] == group["ratio"].max()]).reset_index(level = 0, drop = True)
              
        data = data.drop(["ratio"], axis = 1)
      
        data["file_name"] = file
        data = data[["file_name", "gene_identifier", "k_number", "definition"]]
        
        self.annotations.append(data)

    
    
    def read_dram_annotation_tsv(self):
      """Reads in multiple DRAM annotation.tsv files, keeping the "gene_identifier" 
      as column 0, "k_number"" as column 1, and "definition" as column 2. Keeps 
      only rows where a gene had a KEGG Ortholog annotation.
    
      Returns:
        A list of lists, where each inner list is the annotation data from one file.
      """
      
      pattern = "K[0-9]{5}"
      
      if type(self.files) is str:
        self.files = [self.files]
      
      for file in self.files:
        lines = []
        if file.endswith(".gz"):
          with gzip.open(file, "rb") as f:
            for row in f:
              if re.match(pattern, row.decode().split("\t")[8]):
                lines.append(row.decode().split("\t"))
        else:
          with open(file, "rb") as f:
            for row in f:
              if re.match(pattern, row.decode().split("\t")[8]):
                lines.append(row.decode().split("\t"))
        
        data = pd.DataFrame(lines)[[0,8,9]]
        data.rename(columns={0: "gene_identifier", 8: 'k_number',
        9: "definition"}, inplace=True)
        data["file_name"] = file
        data = data[["file_name", "gene_identifier", "k_number", "definition"]]
    
        
        self.annotations.append(data)
        


    def read_koala_tsv(self):
      """Reads in multiple blastKoala or ghostKoala .tsv files, keeping the 
      "gene_identifier" as column 0, "k_number"" as column 1, and "definition" as 
      column 2. Keeps only rows where a gene had a KEGG Ortholog annotation.
    
      Returns:
        A list of lists, where each inner list is the annotation data from one file.
      """
      
      pattern = "K[0-9]{5}"

      if type(self.files) is str:
        self.files = [self.files]
      
      for file in self.files:
        lines = []
        if file.endswith(".gz"):
          with gzip.open(file, "rb") as f:
            for row in f:
              if re.match(pattern, row.decode().split("\t")[1]):
                lines.append(row.decode().split("\t"))
        else:
          with open(file, "rb") as f:
            for row in f:
              if re.match(pattern, row.decode().split("\t")[1]):
                lines.append(row.decode().split("\t"))
        
        data = pd.DataFrame(lines)[[0,1,2]]
        data.rename(columns={0: "gene_identifier", 1: 'k_number',
        2: "definition"}, inplace=True)
        data["file_name"] = file
        data = data[["file_name", "gene_identifier", "k_number", "definition"]]
    
        self.annotations.append(data)
        


class AnnotationList:
  
    """Data formatting functions to feed formatted data to the MetaPathPredict function"""

  
    def __init__(self, requiredColumnsAll, requiredColumnsModel0, requiredColumnsModel1, annotations, feature_df = pd.DataFrame()): 
      self.requiredColumnsAll = requiredColumnsAll # all required columns for model #1 and model #2
      self.requiredColumnsModel0 = requiredColumnsModel0 # list of all required columns for model #1
      self.requiredColumnsModel1 = requiredColumnsModel1 # list of all required columns for model #2
      self.annotations = annotations
      self.feature_df = feature_df



    def create_feature_df(self):
      """Converts as list of annotations into a Pandas feature DataFrame.
    
      Returns:
        A Pandas DataFrame.
      """

      for df in self.annotations:
        df["count"] = 1
        self.feature_df = pd.concat([self.feature_df, df], axis = 0)
      
      self.feature_df = self.feature_df.groupby(["file_name", "k_number"]).agg(count=("count", "sum")).reset_index().pivot_table(
      index = "file_name",
      columns = "k_number",
      values = "count",
      aggfunc = "first")
      
      self.feature_df = self.feature_df.replace(np.NaN, 0)
      self.feature_df = self.feature_df.where(self.feature_df <= 1, 1)
    


    def check_feature_columns(self):
      """Checks that all required columns are present for both of MetaPathPredict's models.
    
      Returns:
        A Pandas DataFrame.
      """
      
      cols_to_add = [col for col in self.requiredColumnsAll if col not in self.feature_df.columns]
      #self.feature_df.loc[:, cols_to_add] = 0
      col_dict = dict.fromkeys(cols_to_add, 0)
      temp_df = pd.DataFrame(col_dict, index = self.feature_df.index)
      self.feature_df = pd.concat([self.feature_df, temp_df], axis = 1)

      cols_to_drop = [col for col in self.feature_df.columns if col not in self.requiredColumnsAll]
      self.feature_df.drop(cols_to_drop, axis = 1, inplace = True)
      
      self.feature_df = self.feature_df.reindex(self.requiredColumnsAll, axis = 1)
      
      self.feature_df = [self.feature_df, self.feature_df]
      


    # def select_model_features(self):
    #   """Selects all required columns for the specified MetaPathPredict model (both model #1 and model #2).
    # 
    #   Returns:
    #     A Pandas DataFrame.
    #   """
    # 
    #   self.feature_df[0] = self.feature_df[0][self.requiredColumnsModel0]
    #   self.feature_df[0] = self.feature_df[0].reindex(self.requiredColumnsModel0, axis = 1)
    # 
    #   self.feature_df[1] = self.feature_df[1][self.requiredColumnsModel1]
    #   self.feature_df[1] = self.feature_df[1].reindex(self.requiredColumnsModel1, axis = 1)



    # def transform_model_features(self, scaler_0, scaler_1):
    #   """Transforms all required columns for the specified MetaPathPredict model (both model #1 and model #2).
    # 
    #   Returns:
    #     A Pandas DataFrame.
    #   """
    # 
    #   scaled_features_0 = scaler_0.transform(self.feature_df[0])
    #   self.feature_df[0] = pd.DataFrame(scaled_features_0, index = self.feature_df[0].index, columns = self.feature_df[0].columns)
    # 
    #   scaled_features_1 = scaler_1.transform(self.feature_df[1])
    #   self.feature_df[1] = pd.DataFrame(scaled_features_1, index = self.feature_df[1].index, columns = self.feature_df[1].columns)
