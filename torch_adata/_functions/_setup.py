
__module_name__ = "_setup.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu"])


# import local dependencies: ---------------------------------------------
from ._pad_time_resolved_dataset import _pad_time_resolved_dataset
from ._fetch_labels_from_obs import _fetch_labels_from_obs
from ._use_X import _use_X


# ------------------------------------------------------------------------
def _return_X(self, idx):
    return self.X[idx]

def _return_X_and_y(self, idx):
    return self.X[idx], self.y[idx]


# ------------------------------------------------------------------------
def _setup_X(self, use_key):
    self._use_key = use_key
    self.X = _use_X(self._adata, self._use_key)
    self._X_len = self.X.shape[0]
    self._return_item = _return_X

def _setup_y(self, obs_key):

    self._obs_key = obs_key
    if self._obs_key:
        self._OneHot, self.y = _fetch_labels_from_obs(self._adata, self._obs_key)
        self._return_item = _return_X_and_y
        self._y_len = self.y.shape[0]
        assert self._X_len == self._y_len,"X and y do not have the same shape"

def _do_setup(self, use_key, obs_key):

    _setup_X(self, use_key)
    _setup_y(self, obs_key)
    self._len = self._X_len
    
    
# setup time: ------------------------------------------------------------
def _return_X_time_resolved(self, idx):
    return self.X[:, idx]

def _return_X_and_y_time_resolved(self, idx):
    return self.X[:, idx], self.y[:, idx]
    
def _setup_time(self, time_key):

    self._time_key = time_key
    
    if self._obs_key:
        self.X, self.y = _pad_time_resolved_dataset(
            self._adata, self._time_key, self._use_key, self._obs_key
        )
    else:
        self.X = _pad_time_resolved_dataset(
            self._adata, self._time_key, self._use_key, self._obs_key
        )
    self._X_len = self.X.shape[1]
    
    if self._obs_key:
        self._y_len = self.y.shape[1]
        self._return_item = _return_X_and_y_time_resolved
        assert self._X_len == self._y_len,"X and y do not have the same shape"
    else:
        self._return_item = _return_X_time_resolved
    self._len = self._X_len