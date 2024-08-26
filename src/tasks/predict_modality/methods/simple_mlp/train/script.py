import os
import math
import logging
from pathlib import Path

import anndata as ad
import numpy as np

import torch
import pytorch_lightning as pl
from torch.utils.data import TensorDataset, DataLoader
from pytorch_lightning.callbacks import ModelCheckpoint
from pytorch_lightning.loggers import TensorBoardLogger,WandbLogger

logging.basicConfig(level=logging.INFO)

## VIASH START
par = {
    'input_train_mod1': 'resources_test/predict_modality/openproblems_neurips2021/bmmc_multiome/swap/train_mod1.h5ad',
    'input_train_mod2': 'resources_test/predict_modality/openproblems_neurips2021/bmmc_multiome/swap/train_mod2.h5ad',
    'input_test_mod1': 'resources_test/predict_modality/openproblems_neurips2021/bmmc_multiome/swap/test_mod1.h5ad',
    'output': 'output/model'
}
meta = {
    'resources_dir': 'src/tasks/predict_modality/methods/simple_mlp',
    'cpus': 10
}
## VIASH END

resources_dir = f"{meta['resources_dir']}/resources"

import sys
sys.path.append(resources_dir)
from models import MLP
import utils

def _train(X, y, Xt, yt, logger, config, num_workers):
    
    X = torch.from_numpy(X).float()
    y = torch.from_numpy(y).float()
    ymean = torch.mean(y, dim=0, keepdim=True)
    
    tr_ds = TensorDataset(X,y)
    tr_loader = DataLoader(
        tr_ds,
        batch_size=config.batch_size,
        num_workers=num_workers,
        shuffle=True,
        drop_last=True
    )
    
    Xt = torch.from_numpy(Xt).float()
    yt = torch.from_numpy(yt).float()
    te_ds = TensorDataset(Xt,yt)
    te_loader = DataLoader(
        te_ds,
        batch_size=config.batch_size,
        num_workers=num_workers,
        shuffle=False,
        drop_last=False
    )
    
    checkpoint_callback = ModelCheckpoint(
        monitor='valid_RMSE',
        dirpath=logger.save_dir,
        save_top_k=1,
    )
    
    trainer = pl.Trainer(
        devices="auto",
        enable_checkpointing=True,
        logger=logger, 
        max_epochs=config.epochs, 
        callbacks=[checkpoint_callback],
        default_root_dir=logger.save_dir,
        # progress_bar_refresh_rate=5
    )
    
    net = MLP(X.shape[1], y.shape[1], ymean, config)
    trainer.fit(net, tr_loader, te_loader)
    
    yp = trainer.predict(net, te_loader, ckpt_path='best')
    yp = torch.cat(yp, dim=0)
    
    score = ((yp-yt)**2).mean()**0.5
    print(f"VALID RMSE {score:.3f}")
    del trainer
    return score,yp.detach().numpy()



input_train_mod1 = ad.read_h5ad(par['input_train_mod1'])
input_train_mod2 = ad.read_h5ad(par['input_train_mod2'])

mod_1 = input_train_mod1.uns["modality"]
mod_2 = input_train_mod2.uns["modality"]

task = f'{mod_1}2{mod_2}'
yaml_path = f'{resources_dir}/yaml/mlp_{task}.yaml'

obs_info = utils.to_site_donor(input_train_mod1)
# TODO: if we want this method to work for other datasets, resolve dependence on site notation
sites = obs_info.site.unique()

os.makedirs(par['output'], exist_ok=True)

print('Compute ymean', flush=True)
ymean = np.asarray(input_train_mod2.layers["normalized"].mean(axis=0))
path = f"{par['output']}/{task}_ymean.npy"
np.save(path, ymean)


if task == "GEX2ATAC":
    logging.info(f"No training required for this task ({task}).")
    sys.exit(0)

if not os.path.exists(yaml_path):
    logging.error(f"No configuration file found for task '{task}'")
    sys.exit(1)

yaml_path = f'{resources_dir}/yaml/mlp_{task}.yaml'
yps = []
scores = []

msgs = {}
# TODO: if we want this method to work for other datasets, dont use hardcoded range
for fold in range(3):

    run_name = f"{task}_fold_{fold}"
    save_path = f"{par['output']}/{run_name}"
    num_workers = meta["cpus"] or 0

    Path(save_path).mkdir(parents=True, exist_ok=True)   

    X,y,Xt,yt = utils.split(input_train_mod1, input_train_mod2, fold)
    
    logger = TensorBoardLogger(save_path, name='') 
    
    config = utils.load_yaml(yaml_path)

    if config.batch_size > X.shape[0]:
        config = config._replace(batch_size=math.ceil(X.shape[0] / 2))

    score, yp = _train(X, y, Xt, yt, logger, config, num_workers)
    yps.append(yp)
    scores.append(score)
    msg = f"{task} Fold {fold} RMSE {score:.3f}"
    msgs[f'Fold {fold}'] = f'{score:.3f}'
    print(msg)

yp = np.concatenate(yps)
score = np.mean(scores)
msgs['Overall'] = f'{score:.3f}'
print('Overall', f'{score:.3f}')

