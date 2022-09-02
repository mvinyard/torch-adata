
__module_name__ = "_TimeResolvedAnnDataset.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu"])


# import local dependencies #
# ------------------------- #
from ._functions._use_X import _use_X
from ._AnnDataset import AnnDataset
from ._functions._pad_time_resolved_dataset import _pad_time_resolved_dataset

import numpy as np
import torch

class TimeResolvedAnnDataset(AnnDataset):
    def __init__(self,
                 adata,
                 time_key="Time point",
                 data_key="X_pca",
                 weight_key="fate_score",
                 fate_bias_key="X_fate_smoothed",
                ):
        super(TimeResolvedAnnDataset, self).__init__(adata, data_key)

        self._time_key = time_key
        self._data_key = data_key
        self._weight_key = weight_key
        self._fate_bias_key = fate_bias_key
        
        if not self._weight_key in adata.obs.columns:
            adata.obs[self._weight_key] = 1
        if not self._fate_bias_key in adata.uns_keys():
            adata.uns[self._fate_bias_key] = np.ones([len(adata), 11])
        
        self.F = torch.Tensor(adata.uns[self._fate_bias_key])        
        self.X, self.W = _pad_time_resolved_dataset(self._adata,
                                                    self._time_key,
                                                    self._data_key,
                                                    self._weight_key,
                                                    self._fate_bias_key,
                                                   )
        
        