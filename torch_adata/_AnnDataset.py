
__module_name__ = "_AnnDataset.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu"])


# import packages: -------------------------------------------------------
from torch.utils.data import Dataset


# import local dependencies: ---------------------------------------------
from . import _functions as funcs


# main module class: -----------------------------------------------------
class AnnDataset(Dataset):
    """Base class"""
        
    def __init__(self, adata, use_key="X", obs_key=None):

        self._adata = adata
        funcs.do_setup(self, use_key, obs_key)
        
    def __len__(self):
        return self._len
    
    def __getitem__(self, idx):
        return self._return_item(self, idx)