
__module_name__ = "_lightning_anndata_module.py"
__doc__ = """Module for LightningAnnDataModule."""
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu"])


# -- import python natives: --------------------------------------------------------------
import os


# -- import packages: --------------------------------------------------------------------
from lightning import LightningDataModule
from licorice_font import font_format
import pandas as pd
import numpy as np
import anndata
import torch


# -- import local dependencies: ----------------------------------------------------------
from .._AnnDataset import AnnDataset
from ._function_kwargs import function_kwargs
from ._train_val_split import TrainValSplit
from ._configure_adata import configure_adata


# -- Main module class: ------------------------------------------------------------------
class LightningAnnDataModule(LightningDataModule):
    def __init__(
        self,
        adata=None,
        h5ad_path=None,
        batch_size=2000,
        num_workers=os.cpu_count(),
        train_val_split=[0.8, 0.2],
        n_predict=2000,
        use_key="X_pca",
        obs_keys=[],
        groupby=None,
        train_key="train",
        val_key="val",
        test_key="test",
        predict_key="predict",
        shuffle=True,
        shuffle_time_labels = False,
        silent=True,
        **kwargs,
    ):
        super(LightningAnnDataModule, self).__init__()

        self.save_hyperparameters(ignore=["adata"])
        self._adata = adata

    # -- functions: ------------------------------------------------------------
    def _format_adata_obs_index(self):
        self._adata.obs.reset_index(drop=True, inplace=True)
        self._adata.obs.index = self._adata.obs.index.astype(str)

    def _configure_adata(self):
        """configures the property self.adata"""

        if not isinstance(self._adata, anndata.AnnData):
            if isinstance(self.hparams["h5ad_path"], str):
                if not os.path.exists(self.hparams["h5ad_path"]):
                    raise FileNotFoundError("Improper path to .h5ad")
                if not self.hparams["h5ad_path"].endswith(".h5ad"):
                    raise ValueError("Path does not end in .h5ad")
                self._adata = anndata.read_h5ad(self.hparams["h5ad_path"])
            elif not isinstance(self.hparams["h5ad_path"], str):
                raise ValueError("Must pass adata or h5ad_path")

        if not self.properly_formatted_index:
            self._format_adata_obs_index()

    def configure_train_val_split(self):
        
        train_val_split = TrainValSplit(
            self.adata, self.data_keys, self.hparams["train_val_split"]
        )
        self._data_keys = train_val_split.configure_validation(self.init_train_adata, force_reallocate=True)
        
    def subset_adata(self, key: str):        
        self.df = self.adata.obs.copy()
        access_key = self.data_keys[key]
        if not hasattr(self.df, access_key):
            if key == "predict":
                self.df[access_key] = True
            else:
                raise KeyError(
                    "Key: Access Key pair: {}:{} not found".format(key, access_key)
                )
        return self.adata[self.df[access_key]].copy()

    def to_dataset(self, key: str) -> torch.utils.data.Dataset:
        """key funciton to transform adata -> torch.utils.data.Dataset"""
        adata = getattr(self, "{}_adata".format(key))
        return AnnDataset(adata=adata, **self.AnnDatasetKWARGS)
    
    def shuffle_time_labels(self):
        
        df = self._adata.obs.copy()
        non_t0 = df.loc[df['t'] != 0]['t']
        
        shuffled_t = np.zeros(len(df))
        shuffled_t[non_t0.index.astype(int)] = np.random.choice(non_t0.values, len(non_t0))
        self._adata.obs["t"] = shuffled_t

    def _return_loader(self, dataset_key):

        if dataset_key in ["train", "val"]:
            # could probably add a flag to make this optional if desired
            # these happen every time the loader is called, which is useful when
            # you want to shuffle the organization of labels.
            # maybe a better way to do it.... 
            self.configure_train_val_split()
            
            if self.hparams["shuffle_time_labels"]:
                self.shuffle_time_labels()
                
        if dataset_key == "train":
            shuffle=self.hparams["shuffle"]
        else:
            shuffle = shuffle_labels = False
                    

        if dataset_key == "test":
            if not hasattr(self, "n_test_cells"):
                self.setup(stage="test")
            batch_size = self.hparams['test_batch_size']
        else:
            batch_size = self.hparams['batch_size']

        return torch.utils.data.DataLoader(getattr(self, "{}_dataset".format(dataset_key)),
                          num_workers=self.hparams["num_workers"],
                          batch_size=batch_size,
                          shuffle=shuffle,
                         )

    # -- adata properties: -----------------------------------------------------
    @property
    def properly_formatted_index(self):
        return all(
            pd.Index(range(len(self._adata))).astype(str) == self._adata.obs.index
        )

    @property
    def adata(self):
        self._configure_adata()
        return self._adata

    @property
    def cell_idx(self):
        return self.adata.obs.index

    @property
    def n_cells(self):
        return self.adata.shape[0]
    
    @property
    def n_features(self):
        return self.adata.shape[1]

    @property
    def data_keys(self):
        if not hasattr(self, "_data_keys"):
            self._data_keys = {
                attr.split("_key")[0]: val
                for attr, val in self.hparams.items()
                if attr.endswith("key")
            }
        return self._data_keys
    
    @property
    def init_train_adata(self):
        # preserve the original set of fitting data
        # we only want to set this once
        if not hasattr(self, "_init_train_adata"):
            self._init_train_adata = self.subset_adata("train")
        return self._init_train_adata

    @property
    def train_adata(self):
        return self.subset_adata("train")

    @property
    def val_adata(self):
        return self.subset_adata("val")

    @property
    def test_adata(self):
        return self.subset_adata("test")

    @property
    def predict_adata(self):
        
        return self.subset_adata("predict")

    # -- dataset properties: ---------------------------------------------------
    @property
    def AnnDatasetKWARGS(self):
        return function_kwargs(AnnDataset, self.hparams)

    @property
    def train_dataset(self):
        return self.to_dataset("train")

    @property
    def val_dataset(self):
        return self.to_dataset("val")

    @property
    def test_dataset(self):
        return self.to_dataset("test")

    @property
    def predict_dataset(self):
        return self.to_dataset("predict")
    
    def prepare_data(self):
        ...
        
    def setup(self, stage):
        if stage in ["fit"]:
            self.configure_train_val_split()

    # -- Required DataLoader methods: ------------------------------------------
    def train_dataloader(self):
        return self._return_loader("train")

    def val_dataloader(self):
        return self._return_loader("val")

    def test_dataloader(self):
        return self._return_loader("test")

    def predict_dataloader(self):        
        return self._return_loader("predict")

    def __repr__(self):
        return "⚡ {} ⚡".format(font_format("LightningAnnDataModule", ["PURPLE"]))
