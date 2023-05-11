
__module_name__ = "__init__.py"
__doc__ = """Main __init__ module."""
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu"])


# -- specify package version: --------------------------------------------------
__version__ = "0.0.23rc0"


# -- import modules: -----------------------------------------------------------
from ._core._AnnDataset import AnnDataset
from ._core._lightning._lightning_anndata_module import LightningAnnDataModule
from . import _tools as tl


# -- import API-hidden core modules (for dev): ---------------------------------
from ._core import _core_ancilliary as _core
