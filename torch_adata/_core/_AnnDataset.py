
__module_name__ = "_AnnDataset.py"
__doc__ = """Main API module."""
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu"])


# -- import packages: --------------------------------------------------------------------
from torch.utils.data import Dataset
import anndata


# -- import local dependencies: ----------------------------------------------------------
from . import _core_ancilliary as core


# -- supporting functions: ---------------------------------------------------------------
def return_data_on_axis(dataset, idx):
    """
    Automatically return fetched / sampled batches of data along the right data axes.
    """
    if dataset._data_axis:
        return [getattr(dataset, key)[:, idx] for key in dataset._attr_names]
    return [getattr(dataset, key)[idx] for key in dataset._attr_names]


# -- Main module class: ------------------------------------------------------------------
class AnnDataset(Dataset):
    """AnnDataset Module for formatting AnnData as a pytorch Dataset."""

    def __init__(
        self,
        adata: anndata.AnnData,
        use_key: str,
        groupby: str = None,
        obs_keys: list([str, ..., str]) = None,
        attr_names: list([str, ..., str]) = None,
        one_hot: list([bool, ..., bool]) = False,
        aux_keys: list([str, ..., str]) = None,
        silent=False,
    )->None:
        
        super().__init__()

        core.register_init(
            self,
            adata,
            groupby,
            use_key,
            obs_keys,
            aux_keys,
            attr_names,
            one_hot,
            silent,
        )

    def __len__(self):
        return self.X.shape[self._data_axis]

    def __getitem__(self, idx):
        return return_data_on_axis(self, idx)
    
    def __repr__(self) -> str:
        return core.identity_msg(self)
