
__module_name__ = "_setup.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu"])


# import local dependencies: ---------------------------------------------
from ._pad_time_resolved_dataset import _pad_time_resolved_dataset
from ._fetch_data import _fetch_y_from_obs, _use_X


# ------------------------------------------------------------------------
def _return_X(self, idx):
    return self.X[idx]

def _return_X_and_y(self, idx):
    return self.X[idx], self.y[idx]


# ------------------------------------------------------------------------
def _setup_X(self, use_key, onehot_encode_y=False):
    self._use_key = use_key
    self.X = _use_X(self._adata, self._use_key, onehot_encode_y)
    self._X_len = self.X.shape[0]
    self._return_item = _return_X

def _setup_y(self, obs_key, onehot_encode_y):

    self._obs_key = obs_key
    if self._obs_key:
        self._return_item = _return_X_and_y
        self._OneHot, self.y = _fetch_y_from_obs(self._adata, self._obs_key, onehot_encode_y)
        self._y_len = self.y.shape[0]
        assert self._X_len == self._y_len,"X and y do not have the same shape"

def _do_setup(self, use_key, obs_key, onehot_encode_y):

    _setup_X(self, use_key)
    _setup_y(self, obs_key, onehot_encode_y)
    self._len = self._X_len
    
    
# setup time: ------------------------------------------------------------
def _return_X_time_resolved(self, idx):
    return self.X[:, idx]

def _return_X_and_y_time_resolved(self, idx):
    return self.X[:, idx], self.y[:, idx]

def _return_X_time_resolved_w_time(self, idx):
    return self.X[:, idx], self.t

def _return_X_and_y_time_resolved_w_time(self, idx):
    return self.X[:, idx], self.y[:, idx], self.t
    
def _setup_time(self, time_key, return_t):

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
        assert self._X_len == self._y_len,"X and y do not have the same shape"
        
        if return_t:
            self._return_item = _return_X_and_y_time_resolved_w_time
        else:
            self._return_item = _return_X_and_y
    else:
        if return_t:
            self._return_item = _return_X_time_resolved_w_time
        else:
            self._return_item = _return_X
            
    self._len = self._X_len