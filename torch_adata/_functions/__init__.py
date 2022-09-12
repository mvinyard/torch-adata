
__module_name__ = "__init__.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu"])


# ------------------------------------------------------------------------
from ._fetch_labels_from_obs import _fetch_labels_from_obs as fetch_labels
from ._use_X import _use_X as fetch_data
from ._do_setup import _do_setup as do_setup