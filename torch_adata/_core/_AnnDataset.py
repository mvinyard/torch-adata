
__module_name__ = "_AnnDataset.py"
__doc__ = """Main API module."""
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu"])


# import packages: -----------------------------------------------------------------------
from torch.utils.data import Dataset
import anndata


# import local dependencies: -------------------------------------------------------------
from . import _core_ancilliary as core


# supporting functions: ------------------------------------------------------------------
def count_attrs(dataset):
    setattr(dataset, "_n_attrs", len(dataset._attr_names))


# functions called directly by main class: -----------------------------------------------
def register_args_parse_adata(
    dataset, adata, groupby, use_key, obs_keys, attr_names, one_hot_encode
):
    """
    Register passed inputs to AnnDataset class.
    """

    attr_names, one_hot_encode = core.register_args(
        dataset, obs_keys, attr_names, one_hot_encode
    )

    fetch = core.Fetch(adata)

    if groupby:
        dataset._data_axis = 1
        dataset._grouped_by = groupby
        dataset.X, obs_data = fetch.grouped_adata(
            groupby=groupby,
            use_key=use_key,
            obs_keys=obs_keys,
            attr_names=attr_names,
            one_hot=one_hot_encode,
        )
        if obs_keys:
            fetch.update_obs_attrs(dataset, obs_data)

    else:
        dataset._data_axis = 0
        dataset._grouped_by = None
        dataset.X = fetch.X(use_key)
        if obs_keys:
            obs_data = fetch.multi_obs(obs_keys, attr_names, one_hot_encode)
            fetch.update_obs_attrs(dataset, obs_data)
    
    count_attrs(dataset)
    print(core.identity_msg(dataset))
    

def return_data_on_axis(dataset, idx):
    """
    Automatically return fetched / sampled batches of data along the right data axes.
    """
    if dataset._data_axis:
        return [getattr(dataset, key)[:, idx] for key in dataset._attr_names]
    return [getattr(dataset, key)[idx] for key in dataset._attr_names]


# Main module class: ---------------------------------------------------------------------
class AnnDataset(Dataset):
    """AnnDataset Module for formatting AnnData as a pytorch Dataset."""

    def __init__(
        self,
        adata: anndata.AnnData,
        use_key: str,
        groupby: str = None,
        obs_keys: list([str, ..., str]) = None,
        attr_names: list([str, ..., str]) = None,
        one_hot_encode: list([bool, ..., bool]) = False,
    ):
        super().__init__()

        register_args_parse_adata(
            self, adata, groupby, use_key, obs_keys, attr_names, one_hot_encode
        )

    def __len__(self):
        return self.X.shape[self._data_axis]

    def __getitem__(self, idx):
        return return_data_on_axis(self, idx)
    
    def __repr__(self) -> str:
        return core.identity_msg(self)
