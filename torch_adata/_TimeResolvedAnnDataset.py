
__module_name__ = "_TimeResolvedAnnDataset.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu"])


# import local dependencies: ---------------------------------------------
from . import _functions as funcs
from ._AnnDataset import AnnDataset


# Main module function: --------------------------------------------------
class TimeResolvedAnnDataset(AnnDataset):
    def __init__(self, adata, time_key, use_key="X", obs_key=None):
        super().__init__(adata, use_key, obs_key)
        funcs.setup_time(self, time_key)