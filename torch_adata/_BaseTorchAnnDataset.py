
__module_name__ = "_BaseTorchAnnDataset.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu"])


# import packages #
# --------------- #
import torch
from torch.utils.data import Dataset


def _fetch_X_from_adata(adata, use_key=None):
    if not use_key:
        return torch.Tensor(adata.X)

    elif use_key in adata.obsm_keys():
        return torch.Tensor(adata.obsm[use_key])

    else:
        return "{} is not valid. Pass `None` or `use_key` should exist in `adata.obsm_keys()`.".format(
            use_key
        )


def _format_obs_idx(obs):
    """save existing idx and reset index to ensure accurate sampling"""
    return obs.reset_index().copy().rename({"index": "saved_idx"}, axis=1)


class BaseTorchAnnDataset(Dataset):
    def __init__(self, adata, use_key=None):

        self._adata = adata
        self._use_key = use_key
        self._X = _fetch_X_from_adata(self._adata, use_key=self._use_key)
        self._obs = _format_obs_idx(self._adata.obs)

    def __len__(self):
        """Returns the number of cells in the dataset."""
        return len(self._X)

    def __getitem__(self, idx):
        """fetches a cell from the dataset, given a cell idx"""
        return self._X[idx], self._obs.loc[idx]