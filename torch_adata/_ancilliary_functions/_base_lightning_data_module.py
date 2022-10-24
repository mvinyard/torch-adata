
# -- import packages: --------------------------------------------------------------------
from pytorch_lightning import LightningDataModule
from torch.utils.data import DataLoader
import os
import anndata
from abc import ABC, abstractmethod

class BaseLightningDataModule(ABC, LightningDataModule):
    def __init__(
        self, adata: anndata.AnnData = None, batch_size=2000, num_workers=os.cpu_count(), **kwargs
    ):
        super().__init__()

        self.adata = adata
        self.batch_size = batch_size
        self.num_workers = num_workers
        self.kwargs = kwargs
        self.__configure__()
    
    @abstractmethod
    def __configure__(self):
        pass
        
    def train_dataloader(self):
        return DataLoader(
            self.train_dataset, batch_size=self.batch_size, num_workers=self.num_workers
        )

    def val_dataloader(self):
        return DataLoader(
            self.val_dataset, batch_size=self.batch_size, num_workers=self.num_workers
        )

    def test_dataloader(self):
        return DataLoader(
            self.test_dataset, batch_size=self.batch_size, num_workers=self.num_workers
        )

    def predict_dataloader(self):
        return DataLoader(
            self.dataset, batch_size=self.batch_size, num_workers=self.num_workers
        )