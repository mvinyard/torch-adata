
__module_name__ = "_TimeResolvedAnnDataset.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu"])


# import packages: -------------------------------------------------------
import numpy as np


# import local dependencies: ---------------------------------------------
from . import _functions as funcs
from ._AnnDataset import AnnDataset


# Main module function: --------------------------------------------------
class TimeResolvedAnnDataset(AnnDataset):
    def __init__(self, adata, time_key, use_key="X_pca", obs_key=None, return_t=True):
        super().__init__(adata, use_key, obs_key)
        self.t = np.sort(self._adata.obs[time_key].unique())
        funcs.setup_time(self, time_key, return_t)