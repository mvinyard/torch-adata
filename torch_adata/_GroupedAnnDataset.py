

__module_name__ = "_GroupedAnnDataset.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu"])


# import packages #
# --------------- #
import torch


# local imports #
# ------------- #
from ._BaseTorchAnnDataset import BaseTorchAnnDataset


def _index_subset_group(df, X):
    return X[df.index]


# def _grouped_getitem(X, obs, groupby):
#     return obs.groupby(groupby).apply(_index_subset_group, X).to_dict()


class GroupedAnnDataset(BaseTorchAnnDataset):
    def __init__(self, adata, groupby, use_key=None):
        super(GroupedAnnDataset, self).__init__(adata, use_key)
        self._groupby = groupby
        
        
    def _grouped_getitem(self, X, obs, groupby):
        return obs.groupby(groupby).apply(_index_subset_group, X).to_dict()

    def __getitem__(self, idx):
        """keys: groups (e.g., t); values: torch.Tensor([X])"""
        X, obs = self._X[idx], self._obs.loc[idx].reset_index(drop=True)
        return self._grouped_getitem(X, obs, self._groupby)