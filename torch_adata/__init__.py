
__module_name__ = "__init__.py"
__doc__ = """Main __init__ module."""
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu"])


# -- specify package version: ------------------------------------------------------------
__version__ = "0.0.16"


# -- import modules: ---------------------------------------------------------------------
from ._core._AnnDataset import AnnDataset
from ._ancilliary_functions import *
from ._core import _core_ancilliary as dev