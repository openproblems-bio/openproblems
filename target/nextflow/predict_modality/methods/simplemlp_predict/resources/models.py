import torch
import pytorch_lightning as pl
import torch.nn as nn
import torch.nn.functional as F

class MLP(pl.LightningModule):
    def __init__(self,in_dim,out_dim,ymean,config):
        super(MLP, self).__init__()
        self.ymean = ymean.cuda()
        H1 = config.H1
        H2 = config.H2
        p = config.dropout
        self.config = config
        self.fc1 = nn.Linear(in_dim, H1)
        self.fc2 = nn.Linear(H1,H2)
        self.fc3 = nn.Linear(H1+H2, out_dim)
        self.dp2 = nn.Dropout(p=p)

    def forward(self, x):
        x0 = x
        x1 = F.relu(self.fc1(x))
        x1 = self.dp2(x1)
        x = F.relu(self.fc2(x1))
        x = torch.cat([x,x1],dim=1)
        x = self.fc3(x)
        x = self.apply_mask(x)
        return x
    
    def apply_mask(self,yp):
        tmp = torch.ones_like(yp).float()*self.ymean
        mask = tmp<self.config.threshold
        mask = mask.float()
        return yp*(1-mask) + tmp*mask
    
    def training_step(self, batch, batch_nb):
        x,y = batch
        yp = self(x)
        criterion = nn.MSELoss()
        loss = criterion(yp, y)
        self.log('train_loss', loss, prog_bar=True)
        return loss
    
    def validation_step(self, batch, batch_idx):
        x,y = batch
        yp = self(x)
        criterion = nn.MSELoss()
        loss = criterion(yp, y)
        self.log('valid_RMSE', loss**0.5, prog_bar=True)
        return loss
    
    def predict_step(self, batch, batch_idx):
        if len(batch) == 2:
            x,_ = batch
        else:
            x = batch
        return self(x)
    
    def configure_optimizers(self):
        lr = self.config.lr
        wd = float(self.config.wd)
        adam = torch.optim.Adam(self.parameters(), lr=lr, weight_decay=wd)
        if self.config.lr_schedule == 'adam':
            return adam
        elif self.config.lr_schedule == 'adam_cosin':
            slr = torch.optim.lr_scheduler.CosineAnnealingLR(adam, self.config.epochs)
            return [adam], [slr]
        else:
            assert 0