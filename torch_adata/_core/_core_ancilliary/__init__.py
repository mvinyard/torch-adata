
__module_name__ = "__init__.py"
__doc__ = """core __init__ module."""
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu"])


# -- specify package version: ------------------------------------------------------------
__version__ = "0.0.16"


# -- import key functions: ---------------------------------------------------------------
from ._fetch_data import Fetch
from ._register_args import register_args
from ._group_and_pad_adata import group_and_pad_adata
from ._one_hot_encode import one_hot_encode
from ._identity_msg import identity_msg


# -- import data util functions: ---------------------------------------------------------
from ._data_types import is_numpy_array, is_pandas_categorical, to_np_array, tensorize