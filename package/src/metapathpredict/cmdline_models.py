"""
Command Line Interface for MetaPathPredict Tools:
====================================

.. currentmodule:: metapathpredict

class methods:
   MetaPathPredict methods
"""

import logging
import argparse
import datetime
import pickle
import os
import sys
import re
import math

import sklearn
import numpy as np
import pandas as pd
from typing import Iterable, List, Dict, Set, Optional, Sequence
from itertools import chain

from torchvision import transforms
import torch.optim as optim
from torch.utils.data import Dataset, DataLoader, TensorDataset
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.metrics import classification_report
import torch
import torch.nn as nn

from metapathpredict.utils import InputData
from metapathpredict.utils import AnnotationList



# CUDA for PyTorch
use_cuda = torch.cuda.is_available()
device = torch.device("cuda:0" if use_cuda else "cpu")
# device = "cpu"

torch.backends.cudnn.benchmark = True

# Parameters
params = {"batch_size": 64, "shuffle": True, "num_workers": 6}
# Configure the logging system
logging.basicConfig(
    filename="metapathpredict.log",
    level=logging.DEBUG,
    format="%(asctime)s - %(levelname)s - %(message)s",
)



class CustomDataset(Dataset):
    def __init__(self, data, targets, transform=None):
        print("type", type(data), data.shape)
        self.data = torch.tensor(data, dtype=torch.float32)
        self.targets = torch.tensor(targets, dtype=torch.float32)
        self.transform = transform

    def __len__(self):
        return len(self.data)

    def __getitem__(self, idx):
        features, target = self.data[idx], self.targets[idx]

        if self.transform:
            sample = self.transform(sample)

        return features, target



class CustomModel(nn.Module):
    def __init__(self, num_hidden_nodes_per_layer=1024, num_hidden_layers=5):
        super(CustomModel, self).__init__()
        NUM_HIDDEN_NODES = num_hidden_nodes_per_layer
        self.NUM_HIDDEN_LAYERS = num_hidden_layers

        self.fc1 = nn.Linear(2000, NUM_HIDDEN_NODES)
        self.relu = nn.ReLU()
        self.dropout = nn.Dropout(0.1)

        # array of hidden layers
        self.fcs = [
            nn.Linear(NUM_HIDDEN_NODES, NUM_HIDDEN_NODES)
            for i in range(num_hidden_layers)
        ]

        self.output_layer = nn.Linear(NUM_HIDDEN_NODES, 94)
        self.sigmoid = nn.Sigmoid()

    def forward(self, x):
        x = self.fc1(x)
        x = self.relu(x)
        x = self.dropout(x)

        for i in range(self.NUM_HIDDEN_LAYERS - 1):
            x = self.fcs[i](x)
            x = self.relu(x)
            x = self.dropout(x)

        x = self.fcs[self.NUM_HIDDEN_LAYERS - 1](x)
        x = self.relu(x)

        x = self.output_layer(x)
        x = self.sigmoid(x)
        return x



class Models:

    """Platform-agnostic command line functions available in MetaPathPredict tools."""

    @classmethod
    def train(cls, args: Iterable[str] = None) -> int:
        """Train a model from the input data .

        Writes out a DNN model in the keras forma

        Parameters
        ----------
        args : Iterable[str], optional
            value of None, when passed to `parser.parse_args` causes the parser to
            read `sys.argv`

        Returns
        -------
        return_call : 0
            return call if the program completes successfully

        """
        parser = argparse.ArgumentParser()

        parser.add_argument(
            "--train-targets",
            dest="train_targets",
            required=True,
            help="training targets file",
        )
        parser.add_argument(
            "--train-features",
            dest="train_features",
            required=True,
            help="training features",
        )
        parser.add_argument(
            "--num-epochs",
            dest="num_epochs",
            required=False,
            default=100,
            type=int,
            help="number of epochs",
        )
        parser.add_argument(
            "--model-out",
            "-m",
            dest="model_out",
            required=True,
            help="model file name output",
        )
        parser.add_argument(
            "--use-gpu",
            dest="use_gpu",
            required=False,
            action="store_true",
            help="use GPU if available",
        )
        parser.add_argument(
            "--num-cores",
            dest="num_cores",
            required=False,
            default=10,
            type=int,
            help="Number of cores for parallel processing",
        )
        neural_net_params = parser.add_argument_group("Neural Net parameters")
        neural_net_params.add_argument(
            "--num-hidden-layers",
            default=5,
            required=False,
            type=int,
            help="number of hidden layers",
        )
        neural_net_params.add_argument(
            "--hidden-nodes-per-layer",
            type=int,
            required=False,
            default=1024,
            help="number of nodes in each hidden layer",
        )
        neural_net_params.add_argument(
            "--num-features",
            dest="num_features",
            default=2000,
            required=False,
            type=int,
            help="number of features to retain from training data",
        )
        neural_net_params.add_argument(
            "--threshold",
            dest="threshold",
            default=6432,
            required=False,
            type=float,
            help="threshold for SelectKBest feature selection",
        )


        args = parser.parse_args()

        # CUDA for PyTorch
        device = "cpu"
        if args.use_gpu:
            use_cuda = torch.cuda.is_available()
            device = torch.device("cuda:0" if use_cuda else "cpu")

        logging.info(f"Using device: {device}")

        # read in features
        features = pd.read_table(args.train_features, compression="gzip")
        logging.info(f"reading input features of shape: {features.shape[0]} x {features.shape[1]}")

        # read in labels
        targets = pd.read_table(args.train_targets, compression="gzip")
        logging.info(f"reading input labels of shape: {targets.shape[0]} x {targets.shape[1]}")

        # split the data into training and test sets
        test_size = 0.25
        x, x_test, y, y_test = train_test_split(
            features,
            targets,
            stratify=targets,
            shuffle=True,
            test_size= test_size,
            random_state=111,
        )
        logging.info(f"creating test size of: {test_size}%")

        # Split the remaining data to train and validation
        x_train, x_val, y_train, y_val = train_test_split(
            x, y, stratify=y, test_size=0.2, shuffle=True, random_state=111
        )

        print("features size", features.shape)
        print("targets size", targets.shape)

        print("x_test", x_test.shape, " y_test ", y_test.shape)
        print("x", x.shape, " y ", y.shape)

        print("x_train", x_train.shape, " y_train ", y_train.shape)
        print("x_val", x_val.shape, " y_val ", y_val.shape)
        print("x_test", x_test.shape, " y_test ", y_test.shape)
        
        
        
        # Initialize the StandardScaler
        scaler = StandardScaler()
        
        # Fit the scaler to training data and transform it
        # and then transform val and test data w/ the fitted scaler object 
        # (std. dev., variance, etc. are based on training data columns)
        scaled_features = scaler.fit_transform(x_train)
        x_train = pd.DataFrame(scaled_features, index = x_train.index, columns = x_train.columns)
        x_val = pd.DataFrame(scaler.transform(x_val), index = x_val.index, columns = x_val.columns) 
        x_test = pd.DataFrame(scaler.transform(x_test), index = x_test.index, columns = x_test.columns) 
        logging.info(f"normalizing the training input features")
        
        

        # feature selection based only on the training data
        # Select features according to the k highest F-values
        # from running ANOVA on y_train and x_train
        selected_features = []
        for label in y_train:
            selector = SelectKBest(f_classif, k = 'all')
            selector.fit(x_train, y_train[label])
            selected_features.append(list(selector.scores_))

        # select threshold that retains 2000 features
        threshold = args.threshold

        # # MeanCS
        logging.info(f"total number of features in input: {x_train.shape[1]}")
        selected_features2 = np.mean(selected_features, axis = 0) > threshold
        logging.info(f"number of features selected for training: {sum(selected_features2)}")

        # create new training, validation, and test datasets retaining only the 2000 top features
        # determined from the training data
        x_train2 = x_train.loc[:, selected_features2]
        x_val2 = x_val.loc[:, selected_features2]
        x_test2 = x_test.loc[:, selected_features2]
        features_used = x_train2.columns.values
        labels_used = y_val.columns.values

        logging.info(f"Using features : {str(features_used)}")
        logging.info(f"Using labels : {str(labels_used)}")

        # Initialize the StandardScaler
        #scaler = StandardScaler()

        # Fit the scaler to your data and transform it
        #x_train2 = scaler.fit_transform(x_train2)
        #x_val2 = scaler.fit_transform(x_val2)
        #logging.info(f"normalizing the training input features")

        y_train = np.asarray(y_train.values)
        y_val = np.asarray(y_val.values)

        print()
        print("x_train2", x_train2.shape)
        print("x_val2", x_val2.shape)
        print("x_test2", x_test2.shape)

        # outline the neural network architecture - multilable classifier
        # 1 input layer, 5 hidden layers, 1 output layer
        # inclue dropout for all hidden layers
        model = CustomModel(
            num_hidden_nodes_per_layer=args.hidden_nodes_per_layer,
            num_hidden_layers=args.num_hidden_layers,
        ).to(device)

        # Define loss function and optimizer
        criterion = nn.BCELoss()
        optimizer = optim.Adam(model.parameters(), lr=0.001)
        logging.info(f"optimizer Adam with learning rate: 0.001")

        # Define early stopping
        early_stopping = torch.optim.lr_scheduler.ReduceLROnPlateau(
            optimizer, "min", patience=10
        )

        # Create an empty transform
        no_transform = transforms.Compose([])

        # dataset DataLoader
        x_train2 = np.asarray(x_train2)
        x_val2 = np.asarray(x_val2)
        print("xtrain2", x_train2.shape, y_train.shape)

        logging.info(f"loading training dataset into dataloader")
        dataset = CustomDataset(data=x_train2, targets=y_train, transform=None)

        batch_size = 10000
        train_data_loader = DataLoader(
            dataset, batch_size=batch_size, num_workers=args.num_cores, shuffle=True
        )

        logging.info(f"loading testing dataset into dataloader")
        val_dataset = CustomDataset(data=x_val2, targets=y_val, transform=None)
        val_data_loader = DataLoader(val_dataset, batch_size=batch_size, shuffle=True)

        # Train the model
        num_epochs = args.num_epochs
        logging.info(f"number of epochs for training: {num_epochs}")
        for epoch in range(num_epochs):
            model.train()
            train_loss = 0.0

            for inputs, targets in train_data_loader:
                inputs, targets = inputs.to(device), targets.to(device)
                optimizer.zero_grad()
                outputs = model(inputs)
                loss = criterion(outputs, targets)

                loss.backward()
                optimizer.step()
                train_loss += loss.item()

            model.eval()
            val_loss = 0.0
            with torch.no_grad():
                for inputs, targets in val_data_loader:
                    inputs, targets = inputs.to(device), targets.to(device)
                    outputs = model(inputs)
                    loss = criterion(outputs, targets)
                    val_loss += loss.item()

            # Update learning rate using early stopping
            early_stopping.step(val_loss)

            logging.info(
                f"Epoch [{epoch+1}/{num_epochs}], Train Loss: {train_loss:.4f}, Val Loss: {val_loss:.4f}"
            )

            print(
                f"Epoch [{epoch+1}/{num_epochs}], Train Loss: {train_loss:.4f}, Val Loss: {val_loss:.4f}"
            )

        # assess the model on test data
        x_test2 = np.asarray(x_test2)
        x_test2 = torch.tensor(x_test2, dtype=torch.float32)
        logging.info(f"converting test inputs to torch.tensor")

        predictions_test = model(x_test2)

        # round predictions
        roundedTestPreds = np.round(predictions_test.detach().numpy())

        # print out performance metrics
        print(classification_report(y_test.values, roundedTestPreds))

        logging.info(f"Training finished successfully!")

        model_file = {}
        model_file["description"] = "neural net trained for predicting multilabels"
        model_file["features"] = features_used
        model_file["labels"] = labels_used
        model_file["model"] = model
        torch.save(model_file, args.model_out)
        logging.info(f"writing model file: {args.model_out}")



    @classmethod
    def predict(cls, args: Iterable[str] = None) -> int:
        """Predict the presence or absence of select KEGG modules on bacterial
        annotation data.

        Parameters
        ----------
        args : Iterable[str], optional
            value of None, when passed to `parser.parse_args` causes the parser to
            read `sys.argv`

        Returns
        -------
        return_call : 0
            return call if the program completes successfully

        """
        parser = argparse.ArgumentParser()

        parser.add_argument(
            "--input",
            "-i",
            action = "extend",
            nargs = "+",
            dest="input",
            required=True,
            help="input file path(s) [required]",
        )
        parser.add_argument(
            "--annotation-format",
            "-a",
            dest="annotation_format",
            required=True,
            help="annotation format [kofamscan, dram, koala; default: kofamscan]",
        )
        parser.add_argument(
            "--output",
            "-o",
            dest="output",
            required=True,
            help="output file name [required; no default folder created]",
        )
        parser.add_argument(
            "--model-files",
            "-m",
            action="extend",
            nargs="+",
            dest="model_in",
            required=True,
            help="input model files location [default: MetaPathPredict folder]",
        )
        parser.add_argument(
            "--scaler-files",
            "-s",
            action="extend",
            nargs="+",
            dest="scaler_in",
            required=True,
            help="input scaler files location [default: MetaPathPredict folder]",
        )
        parser.add_argument(
            "--kegg-modules",
            "-k",
            dest="kegg_modules",
            required=False,
            default=None,
            action="extend",
            nargs="+",
            help="KEGG modules to predict [default: MetaPathPredict KEGG modules]",
        )


        args = parser.parse_args()
        
        with open(args.scaler_in[0],'rb') as f:
          commonClassifierScaler = pickle.load(f)
          
        with open(args.scaler_in[1],'rb') as g:
          rareClassifierScaler = pickle.load(g)

        models = [torch.load(args.model_in[0]), torch.load(args.model_in[1])]
        
        logging.info(f"reading model files: {args.model_in[0]}, {args.model_in[1]}")
        logging.info(f"reading scaler files: {args.scaler_in[0]}, {args.scaler_in[1]}")

        # load the input features
        files_list = InputData(files = args.input) 
        
        if args.annotation_format == "kofamscan":
          files_list.read_kofamscan_detailed_tsv()
          
        elif args.annotation_format == "dram":
          files_list.read_dram_annotation_tsv()
          
        elif args.annotation_format == "koala":
          files_list.read_koala_tsv()
        
        else:
          logging.error("""Did not recognize annotation format; use "kofamscan", "dram", or "koala""""")
          sys.exit(0)
          
          
        model_0_cols = np.ndarray.tolist(commonClassifierScaler.feature_names_in_)
        model_1_cols = np.ndarray.tolist(rareClassifierScaler.feature_names_in_)
        reqColsAll = list(set(model_0_cols).union(set(model_1_cols)))
          
        input_features = AnnotationList(
          requiredColumnsAll = reqColsAll, # add list of all required columns for model #1 and model #2
          requiredColumnsModel0 = commonClassifierScaler.feature_names_in_, # add list of all required columns for model #1
          requiredColumnsModel1 = rareClassifierScaler.feature_names_in_, # add list of all required columns for model #2
          annotations = files_list.annotations)

        input_features.create_feature_df()
        input_features.check_feature_columns()
        input_features.select_model_features()
        input_features.transform_model_features(commonClassifierScaler, rareClassifierScaler)


        predictions_list = []
        for x in range(2):
          
          # convert to pytorch.tensor
          features = torch.tensor(np.asarray(input_features.feature_df[x]), dtype=torch.float32)

          # predict
          predictions = models[x]['model'](features)

          # round predictions
          roundedPreds = np.round(predictions.detach().numpy())
          
          predsDf = pd.DataFrame(data = roundedPreds, columns = models[x]['labels']).astype(int)
  
          predictions_list.append(predsDf)

        
        out_df = pd.concat(predictions_list, axis = 1)

        if args.kegg_modules is not None:
          if all(modules in out_df.columns for module in args.kegg_modules):
            out_df = out_df[[args.kegg_modules]]
          else:
            logging.error("""Did not recognize one or more KEGG modules specified with --kegg-modules; keeping all prediction columns""")

        out_df.insert(loc = 0, column = 'file', value = args.input)
        
        logging.info(f"writing output to file: {args.output}")
        out_df.to_csv(args.output, sep='\t', index=None)

        logging.info(f"Output matrix size: {out_df.shape[0]} x {out_df.shape[1]}")
