__module_name__ = "_data_typing.py"
__doc__ = """Data typing module"""
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu"])


# -- import packages: --------------------------------------------------------------------
import numpy
import pandas
import torch


# -- data type handlers: -----------------------------------------------------------------
def is_numpy_array(x_arr) -> bool:
    """
    Inspects and indicates if input is numpy.ndarray

    Parameters:
    -----------
    x_arr
        type: unknown

    Returns:
    --------
    indicator
        type: bool
    """
    return x_arr.__class__ is numpy.ndarray


def is_pandas_categorical(x_arr) -> bool:
    """
    Inspects and indicates if input is pandas.Categorical

    Parameters:
    -----------
    x_arr
        type: unknown

    Returns:
    --------
    indicator
        type: bool
    """
    return x_arr.__class__ is pandas.Categorical


def to_np_array(x_arr) -> numpy.ndarray:
    """
    Inspects input. If input is not, numpy.ndarray, it is transformed to numpy.ndarray.

    Parameters:
    -----------
    x_arr
        type: unknown

    Returns:
    --------
    x_arr_transformed
        type: numpy.ndarray
    """
    if is_numpy_array(x_arr):
        return x_arr
    if is_pandas_categorical(x_arr):
        return x_arr.to_numpy()
    return x_arr.toarray()


def tensorize(x_arr) -> torch.Tensor:
    """
    Inspects input. Ensures it is numpy.ndarray and then converts to torch.Tensor.

    Parameters:
    -----------
    x_arr
        type: unknown (np.ndarray-like)

    Returns:
    --------
    x_arr_transformed
        type: torch.Tensor
    """
    return torch.Tensor(to_np_array(x_arr))


def as_list(item):
    if item:
        if not isinstance(item, list):
            return [item]
        return item
    return []