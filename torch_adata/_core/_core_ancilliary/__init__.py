
__module_name__ = "__init__.py"
__doc__ = """core __init__ module."""
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["mvinyard.ai@gmail.com"])


# -- import AnnDataset-facing functions: -------------------------------------------------
from ._register_init import register_init
from ._return_data_on_axis import return_data_on_axis

# -- import key functions: ---------------------------------------------------------------
from ._fetch_data import Fetch
from ._group_and_pad_adata import group_and_pad_adata
from ._one_hot_encode import one_hot_encode
from ._identity_msg import identity_msg


# -- import data util functions: ---------------------------------------------------------
from ._data_typing import (
    is_numpy_array,
    is_pandas_categorical,
    to_np_array,
    tensorize,
    as_list,
)