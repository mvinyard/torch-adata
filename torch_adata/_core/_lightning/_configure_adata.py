
__module_name__ = "_configure_adata.py"
__doc__ = """Module to configure AnnData within LightningDataModule"""
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu"])


# -- import packages: --------------------------------------------------------------------
import pandas as pd
import anndata


# -- import python natives: --------------------------------------------------------------
import os


# -- import local dependencies: ----------------------------------------------------------
from ._auto_parse_base_class import AutoParseBase


# -- Main class: -------------------------------------------------------------------------
class AnnDataConfig(AutoParseBase):
    def __init__(self, adata=None, h5ad_path=None):
        
        self.__parse__(locals(), private="adata")
        
    @property
    def properly_formatted_index(self):
        return all(
            pd.Index(range(len(self._adata))).astype(str) == self._adata.obs.index
        )
    
    @property
    def adata_is_AnnData(self):
        return isinstance(self._adata, anndata.AnnData)
    
    @property
    def h5ad_path_is_str(self):
        return isinstance(self.h5ad_path, str)
    
    @property
    def h5ad_path_exists(self):
        return os.path.exists(self.h5ad_path)
    
    @property
    def h5ad_path_ends_properly(self):
        return self.h5ad_path.endswith(".h5ad")
    
    def format_obs_index(self):
        self._adata.obs.reset_index(drop=True, inplace=True)
        self._adata.obs.index = self._adata.obs.index.astype(str)
        
    def __call__(self):
        
        if not self.adata_is_AnnData:
            if self.h5ad_path_is_str:
                if not self.h5ad_path_exists:
                    raise FileNotFoundError("Path to AnnData does not exist.")
                if not self.h5ad_path_ends_properly:
                    raise ValueError("Path does not end in .h5ad")
                self._adata = anndata.read_h5ad(self.h5ad_path)
            else:
                raise ValueError("Must pass adata or h5ad_path")
                
        if not self.properly_formatted_index:
            self.format_obs_index()


# -- API-facing function: ----------------------------------------------------------------
def configure_adata(adata=None, h5ad_path:str = None):
    
    config = AnnDataConfig(adata=adata, h5ad_path=h5ad_path)
    return config()
