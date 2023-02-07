
__module_name__ = "_one_hot_encode.py"
__doc__ = """sklearn.preprocessing.OneHotEncoder() wrapper module."""
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu"])


# -- import packages: --------------------------------------------------------------------
import sklearn.preprocessing
import torch
import numpy


# -- sole function: ----------------------------------------------------------------------
def one_hot_encode(y_data: numpy.ndarray) -> (sklearn.preprocessing.OneHotEncoder, torch.Tensor):

    """
    One-hot encode a y_data vector of shape y to a matrix of shape: y x n_categories.

    Parameters:
    -----------
    y_data
        numpy.ndarray

    Returns:
    --------
    encoder
        type: sklearn.preprocessing.OneHotEncoder

    y_one_hot
        type: torch.Tensor

    """

    encoder = sklearn.preprocessing.OneHotEncoder()
    y_one_hot = encoder.fit_transform(y_data).toarray()

    return encoder, torch.Tensor(y_one_hot)
