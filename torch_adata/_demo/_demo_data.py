

from .. import _utils as utils
from ._pbmc_3k import PBMC3k
from ._larry import LARRY


class DemoData(utils.ABCParse):
    def __init__(self):
        """"""

    @property
    def PBMC3k(self):
        if not hasattr(self, "_pbmcs"):
            self._pbmcs = PBMC3k()
        else:
            self._pbmcs.plot()
        return self._pbmcs.adata

    @property
    def LARRY(self):
        if not hasattr(self, "_larry"):
            self._larry = LARRY()
        else:
            self._larry.plot()
        return self._larry.adata