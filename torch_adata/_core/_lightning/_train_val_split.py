
__module_name__ = "_train_val_split.py"
__doc__ = """Split train / val data."""
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu"])


# -- import python natives: --------------------------------------------------------------
from typing import Union, List


# -- import packages: --------------------------------------------------------------------
import numpy as np

# -- import local dependencies: ----------------------------------------------------------
from ._auto_parse_base_class import AutoParseBase


# -- Supporting class: -------------------------------------------------------------------
class SplitSize:
    def __init__(self, n_cells: int, n_groups: int):

        self._n_cells = n_cells
        self._n_groups = n_groups

    def _sum_norm(self, vals: Union[List, np.ndarray]) -> np.ndarray:
        return np.array(vals) / vals.sum()

    def uniformly(self):
        div, mod = divmod(self._n_cells, self._n_groups)
        return [div] * (self._n_groups - 1) + [div + mod]

    def proportioned(self, percentages=[0.8, 0.2], remainder_idx=-1):

        percentages = self._sum_norm(np.array(percentages))
        split_lengths = [int(self._n_cells * pct_i) for pct_i in percentages]
        remainder = self._n_cells - sum(split_lengths)
        split_lengths[remainder_idx] += remainder

        return split_lengths


# -- Main module class: ------------------------------------------------------------------
class TrainValSplit(AutoParseBase):
    def __init__(
        self,
        adata,
        data_keys,
        train_val_split=[0.8, 0.2],
    ):

        self.__parse__(locals())

    @property
    def obs_df(self):
        return self.adata.obs.copy()

    @property
    def n_cells(self):
        return len(self.adata)

    @property
    def obs_cols(self):
        return self.obs_df.columns.tolist()

    def _configure_training_subset(self):

        if not self.train_key in self.obs_cols:
            self.obs_df[self.train_key] = True

    @property
    def n_groups(self):
        return len(self.train_val_split)

    @property
    def uniform_split(self) -> list([int, ..., int]):
        return self.split.uniformly()

    @property
    def proportioned_split(self) -> list([int, ..., int]):
        return self.split.proportioned(percentages=self.train_val_split)

    @property
    def has_train_not_val(self):
        return (hasattr(self.obs_df, self.data_keys["train"])) and (
            not hasattr(self.obs_df, self.data_keys["val"])
        )

    @property
    def split(self):
        return SplitSize(self.n_fit, self.n_groups)

    def _set_new_idx(self, df, idx, key_added):

        tmp = np.zeros(len(df), dtype=bool)
        tmp[idx] = True
        df[key_added] = tmp.astype(bool)

    @property
    def train_cells(self):

        train_idx = np.random.choice(
            range(len(self.original_train_idx)), size=self.n_train, replace=False
        )
        train_cells = np.zeros(len(self.original_train_idx), dtype=bool)
        train_cells[train_idx] = True

        return train_cells

    def allocate_validation(self, train_adata, remainder_idx=-1):

        """
        If validation key is not found, invoke this function. Takes the train subset
        adata and breaks it into non-overlapping train and validation adata subsets.
        """
        
        self.n_fit = train_adata.shape[0]
        self.n_train = train_adata.shape[0]
        self.original_train_idx = train_adata.obs.index

        self.data_keys["train"] = "fit_train"
        self.data_keys["val"] = "fit_val"

        self.n_train, self.n_val = self.split.proportioned(
            percentages=self.train_val_split, remainder_idx=remainder_idx
        )

        fit_train_idx = self.original_train_idx[self.train_cells]
        fit_val_idx = self.original_train_idx[~self.train_cells]

        processed_df = self.obs_df.copy()

        self._set_new_idx(
            processed_df, idx=fit_train_idx.astype(int), key_added="fit_train"
        )
        self._set_new_idx(
            processed_df, idx=fit_val_idx.astype(int), key_added="fit_val"
        )

        self.adata.obs = processed_df

    def configure_validation(self, train_adata, force_reallocate=False):
        """configure train-val split if it has train but not val"""
        if self.has_train_not_val or force_reallocate:
            self.allocate_validation(train_adata)
        
        return self.data_keys
