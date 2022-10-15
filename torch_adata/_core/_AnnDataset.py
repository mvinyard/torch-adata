
__module_name__ = "_AnnDataset.py"
__doc__ = """Main API module."""
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu"])


# -- import packages: --------------------------------------------------------------------
from torch.utils.data import Dataset
import anndata


# -- import local dependencies: ----------------------------------------------------------
from . import _core_ancilliary as core


# -- Main module class: ------------------------------------------------------------------
class AnnDataset(Dataset):
    """AnnDataset Module for formatting AnnData as a pytorch Dataset."""

    def __init__(
        self,
        adata: anndata.AnnData,
        use_key: str,
        groupby: str = None,
        obs_keys: list([str, ..., str]) = None,
        attr_names: dict({"obs":[str, ..., str], "aux":[str, ..., str]}) = {"obs":[], "aux":[]},
        one_hot: list([bool, ..., bool]) = False,
        aux_keys: list([str, ..., str]) = None,
        silent: bool =False,
    )->None:
        
        """
        adata
            type: anndata.AnnData

        use_key
            AnnData accession key that sets dataset.X attribute.
            type: str
        
        groupby
            Stack dataset on this axis.
            type: str
            default: None
        
        obs_keys
            type: list([str, ..., str])
            default: None
            
        attr_names
            type: dict({"obs":[str, ..., str], "aux":[str, ..., str]})
            default: {"obs":[], "aux":[]}
            
        one_hot
            type: list([bool, ..., bool])
            default: False
        
        aux_keys
            type: list([str, ..., str])
            default: None
        
        silent
            type: bool
            default: False
        """
        
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
        return core.return_data_on_axis(self, idx)
    
    def __repr__(self) -> str:
        return core.identity_msg(self)
