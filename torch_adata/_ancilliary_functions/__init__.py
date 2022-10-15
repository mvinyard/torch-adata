
__module_name__ = "__init__.py"
__doc__ = """functions __init__ module."""
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu"])


# -- specify version: --------------------------------------------------------------------
__version__ = "0.0.16"


# -- import ancilliary functions: --------------------------------------------------------
from ._split import split
from ._idx import idx
from ._dummy_batch import dummy_batch
from ._base_lightning_data_module import BaseLightningDataModule