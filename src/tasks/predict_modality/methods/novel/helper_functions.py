import torch

from torch import nn
import torch.nn.functional as F

from torch.utils.data import Dataset

from typing import Optional

import anndata
import numpy as np
import pandas as pd
import scipy.sparse
import sklearn.decomposition
import sklearn.feature_extraction.text
import sklearn.preprocessing
import sklearn.neighbors
import sklearn.utils.extmath

class tfidfTransformer():
    def __init__(self):
        self.idf = None
        self.fitted = False

    def fit(self, X):
        self.idf = X.shape[0] / X.sum(axis=0)
        self.fitted = True

    def transform(self, X):
        if not self.fitted:
            raise RuntimeError('Transformer was not fitted on any data')
        if scipy.sparse.issparse(X):
            tf = X.multiply(1 / X.sum(axis=1))
            return tf.multiply(self.idf)
        else:
            tf = X / X.sum(axis=1, keepdims=True)
            return tf * self.idf

    def fit_transform(self, X):
        self.fit(X)
        return self.transform(X)

class lsiTransformer():
    def __init__(self,
                 n_components: int = 20,
                 use_highly_variable = None
                ):
        self.n_components = n_components
        self.use_highly_variable = use_highly_variable
        self.tfidfTransformer = tfidfTransformer()
        self.normalizer =  sklearn.preprocessing.Normalizer(norm="l1")
        self.pcaTransformer = sklearn.decomposition.TruncatedSVD(n_components = self.n_components, random_state=777)
        # self.lsi_mean = None
        # self.lsi_std = None
        self.fitted = None
        
    def fit(self, adata: anndata.AnnData):
        if self.use_highly_variable is None:
            self.use_highly_variable = "hvg" in adata.var
        adata_use = adata[:, adata.var["hvg"]] if self.use_highly_variable else adata
        X = self.tfidfTransformer.fit_transform(adata_use.X)
        X_norm = self.normalizer.fit_transform(X)
        X_norm = np.log1p(X_norm * 1e4)
        X_lsi = self.pcaTransformer.fit_transform(X_norm)
        # self.lsi_mean = X_lsi.mean(axis=1, keepdims=True)
        # self.lsi_std = X_lsi.std(axis=1, ddof=1, keepdims=True)
        self.fitted = True
    
    def transform(self, adata):
        if not self.fitted:
            raise RuntimeError('Transformer was not fitted on any data')
        adata_use = adata[:, adata.var["hvg"]] if self.use_highly_variable else adata
        X = self.tfidfTransformer.transform(adata_use.X)
        X_norm = self.normalizer.transform(X)
        X_norm = np.log1p(X_norm * 1e4)
        X_lsi = self.pcaTransformer.transform(X_norm)
        X_lsi -= X_lsi.mean(axis=1, keepdims=True)
        X_lsi /= X_lsi.std(axis=1, ddof=1, keepdims=True)
        lsi_df = pd.DataFrame(X_lsi, index = adata_use.obs_names)
        return lsi_df
    
    def fit_transform(self, adata):
        self.fit(adata)
        return self.transform(adata)

class ModalityMatchingDataset(Dataset):
    def __init__(
        self, df_modality1, df_modality2, is_train=True
    ):
        super().__init__()
        self.df_modality1 = df_modality1
        self.df_modality2 = df_modality2
        self.is_train = is_train
    def __len__(self):
        return self.df_modality1.shape[0]
    
    def __getitem__(self, index: int):
        if self.is_train == True:
            x = self.df_modality1.iloc[index].values
            y = self.df_modality2.iloc[index].values
            return x, y
        else:
            x = self.df_modality1.iloc[index].values
            return x

class Swish(torch.autograd.Function):
    @staticmethod
    def forward(ctx, i):
        result = i * sigmoid(i)
        ctx.save_for_backward(i)
        return result
    @staticmethod
    def backward(ctx, grad_output):
        i = ctx.saved_variables[0]
        sigmoid_i = sigmoid(i)
        return grad_output * (sigmoid_i * (1 + i * (1 - sigmoid_i)))

class Swish_module(nn.Module):
    def forward(self, x):
        return Swish.apply(x)
    
sigmoid = torch.nn.Sigmoid()

class ModelRegressionGex2Atac(nn.Module):
    def __init__(self, dim_mod1, dim_mod2):
        super(ModelRegressionGex2Atac, self).__init__()
        #self.bn = torch.nn.BatchNorm1d(1024)
        self.input_ = nn.Linear(dim_mod1, 1024)
        self.fc = nn.Linear(1024, 256)
        self.fc1 = nn.Linear(256, 2048)
        self.dropout1 = nn.Dropout(p=0.298885630228993)
        self.dropout2 = nn.Dropout(p=0.11289717442776658)
        self.dropout3 = nn.Dropout(p=0.13523634924414762)
        self.output = nn.Linear(2048, dim_mod2)
    def forward(self, x):
        x = F.gelu(self.input_(x))
        x = self.dropout1(x)
        x = F.gelu(self.fc(x))
        x = self.dropout2(x)
        x = F.gelu(self.fc1(x))
        x = self.dropout3(x)
        x = F.gelu(self.output(x))
        return x

class ModelRegressionAtac2Gex(nn.Module): #
    def __init__(self, dim_mod1, dim_mod2):
        super(ModelRegressionAtac2Gex, self).__init__()
        self.input_ = nn.Linear(dim_mod1, 2048)
        self.fc = nn.Linear(2048, 2048)
        self.fc1 = nn.Linear(2048, 512)
        self.dropout1 = nn.Dropout(p=0.2649138776004753)
        self.dropout2 = nn.Dropout(p=0.1769628308148758)
        self.dropout3 = nn.Dropout(p=0.2516791883012817)
        self.output = nn.Linear(512, dim_mod2)
    def forward(self, x):
        x = F.gelu(self.input_(x))
        x = self.dropout1(x)
        x = F.gelu(self.fc(x))
        x = self.dropout2(x)
        x = F.gelu(self.fc1(x))
        x = self.dropout3(x)
        x = F.gelu(self.output(x))
        return x

class ModelRegressionAdt2Gex(nn.Module):
    def __init__(self, dim_mod1, dim_mod2):
        super(ModelRegressionAdt2Gex, self).__init__()
        self.input_ = nn.Linear(dim_mod1, 512)
        self.dropout1 = nn.Dropout(p=0.0)
        self.swish = Swish_module()
        self.fc = nn.Linear(512, 512)
        self.fc1 = nn.Linear(512, 512)
        self.fc2 = nn.Linear(512, 512)
        self.output = nn.Linear(512, dim_mod2)
    def forward(self, x):
        x = F.gelu(self.input_(x))
        x = F.gelu(self.fc(x))
        x = F.gelu(self.fc1(x))
        x = F.gelu(self.fc2(x))
        x = F.gelu(self.output(x))
        return x
    
class ModelRegressionGex2Adt(nn.Module):
    def __init__(self, dim_mod1, dim_mod2):
        super(ModelRegressionGex2Adt, self).__init__()
        self.input_ = nn.Linear(dim_mod1, 512)
        self.dropout1 = nn.Dropout(p=0.20335661386636347)
        self.dropout2 = nn.Dropout(p=0.15395289261127876)
        self.dropout3 = nn.Dropout(p=0.16902655078832815)
        self.fc = nn.Linear(512, 512)
        self.fc1 = nn.Linear(512, 2048)
        self.output = nn.Linear(2048, dim_mod2)
    def forward(self, x):
       # x = self.batchswap_noise(x)
        x = F.gelu(self.input_(x))
        x = self.dropout1(x)
        x = F.gelu(self.fc(x))
        x = self.dropout2(x)
        x = F.gelu(self.fc1(x))
        x = self.dropout3(x)
        x = F.gelu(self.output(x))
        return x

def rmse(y, y_pred):
    return np.sqrt(np.mean(np.square(y - y_pred)))

def train_and_valid(model, optimizer, loss_fn, dataloader_train, dataloader_test, name_model, device):
    best_score = 100000
    for i in range(100):
        train_losses = []
        test_losses = []
        model.train()

        for x, y in dataloader_train:
            optimizer.zero_grad()
            output = model(x.float().to(device))
            loss = torch.sqrt(loss_fn(output, y.float().to(device)))
            loss.backward()
            train_losses.append(loss.item())
            optimizer.step()
           
        model.eval()
        with torch.no_grad():
            for x, y in dataloader_test:
                output = model(x.float().to(device))
                output[output<0] = 0.0
                loss = torch.sqrt(loss_fn(output, y.float().to(device)))
                test_losses.append(loss.item())
        
        outputs = []
        targets = []
        model.eval()
        with torch.no_grad():
            for x, y in dataloader_test:
                output = model(x.float().to(device))
                
                outputs.append(output.detach().cpu().numpy())
                targets.append(y.float().detach().cpu().numpy())
        cat_outputs = np.concatenate(outputs)
        cat_targets = np.concatenate(targets)
        cat_outputs[cat_outputs<0.0] = 0
        
        if best_score > rmse(cat_targets,cat_outputs):
            torch.save(model.state_dict(), name_model)
            best_score = rmse(cat_targets,cat_outputs)
    print("best rmse: ", best_score)
    
