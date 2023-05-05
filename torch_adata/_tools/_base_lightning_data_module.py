
__module_name__ = "_base_lightning_data_module.py"
__doc__ = """Aux. module to organize AnnData/torch datasets into PyTorch-Lightning LightningDataModule."""
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu"])


# -- import native modules: ----------------------------------------------------
from abc import ABC, abstractmethod
import os


# -- import packages: ----------------------------------------------------------
from lightning import LightningDataModule
from torch.utils.data import DataLoader
import anndata


# -- main module class: --------------------------------------------------------
class BaseLightningDataModule(ABC, LightningDataModule):
    def __init__(
        self,
        adata: anndata.AnnData = None,
        batch_size: int = 2000,
        num_workers: int = os.cpu_count(),
        **kwargs
    ):
        super(BaseLightningDataModule, self).__init__()
        self.__parse__(locals())        

    def __parse__(self, kwargs, ignore=["self", "__class__"]):
        self._kwargs = {}
        for k, v in kwargs.items():
            if not k in ignore:
                setattr(self, k, v)
                self._kwargs[k] = v
                if k == "kwargs":
                    for l, w in v.items():
                        setattr(self, l, w)
                        self._kwargs[l] = w

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
