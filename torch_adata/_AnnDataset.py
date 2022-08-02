
__module_name__ = "_AnnDataset.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu"])


# import packages #
# --------------- #
from torch.utils.data import Dataset


# import local dependencies #
# ------------------------- #
from ._functions._use_X import _use_X


class AnnDataset(Dataset):
    """Base class"""
    def __init__(self, adata=None, use_key="X_pca"):

        self._adata = adata
        self.X = _use_X(self._adata, use_key)

    def __len__(self):
        return self.X.shape[1]

    def __getitem__(self, idx):
        return self.X[:, idx]
