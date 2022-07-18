
__module_name__ = "_TimeResolvedAnnDataset.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu"])


# import packages #
# --------------- #
import torch


from ._GroupedAnnDataset import GroupedAnnDataset

class _TimeResolvedAnnDataset(GroupedAnnDataset):
    def __init__(self, adata, time_key="Time point", use_key=None):
        super(TimeResolvedAnnDataset, self).__init__(adata, time_key, use_key)

        self._time_key = time_key
        self._t_init = min(adata.obs[time_key])

    def __getitem__(self, idx):
        """keys: groups (e.g., t); values: torch.Tensor([X])"""
        X, obs = self._X[idx], self._obs.loc[idx].reset_index(drop=True)

        X_dict = _grouped_getitem(X, obs, self._groupby)
        X0 = X_dict[self._t_init]
        t = torch.Tensor(list(X_dict_.keys()))

        return X0, X_dict, t, obs