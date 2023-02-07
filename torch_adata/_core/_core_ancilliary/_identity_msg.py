
__module_name__ = "_identity_msg.py"
__doc__ = """Print the identity of the AnnDataset class."""
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu"])


# -- import packages: --------------------------------------------------------------------
import licorice_font as lf


# -- import local dependencies: ----------------------------------------------------------
from ._data_typing import as_list


# -- supporting functions: ---------------------------------------------------------------
def _header_msg(dataset):
    l = len(dataset)
    msg = "[ {}-{} ]: AnnDataset object with {} samples\n{}"
    t, a = lf.font_format("torch", ["RED", "BOLD"]), lf.font_format("adata", ["PURPLE", "BOLD"])
    u = "".join(["-"] * int(48 + len(str(l))))
    return msg.format(t, a, l, u)

def _groupby_msg(dataset):

    if not dataset._grouped_by:
        return
    msg_g1 = lf.font_format("Grouped by", ["BOLD"])
    msg_g2 = lf.font_format(dataset._grouped_by, ["RED"])
    return "{}: '{}' with attributes:".format(msg_g1, msg_g2)

def annotate_attr_size(dataset, attr_name_set):
    attr_name_set_ = []
    for attr_name in attr_name_set:
        shape = str(getattr(dataset, attr_name).shape)
        attr_fmt_name = ": ".join([attr_name, shape])
        attr_name_set_.append(attr_fmt_name)
    return as_list(attr_name_set_)
        
def _attr_msg(dataset, attr_name_set, title):
    if attr_name_set:
        attr_name_set = annotate_attr_size(dataset, attr_name_set)
        return " - {:<15}".format(lf.font_format(title + ":", ["BOLD"])) + ", ".join(attr_name_set)

def _all_attr_msg(dataset):

    msg_list = []
    aux_msg = _attr_msg(dataset, dataset._aux_attr_names, title="X_aux")
    if aux_msg:
        msg_list.append(aux_msg)
    obs_msg = _attr_msg(dataset, dataset._obs_attr_names, title="obs")
    if obs_msg:
        msg_list.append(obs_msg)
    if any([aux_msg, obs_msg]):
        return "\n".join(msg_list)

def _X_msg(dataset, use_key):
    shape = str(dataset.X.shape)
    msg = "{} (use_key = '{}' {})"
    fmt = lf.font_format(" - X", ["BOLD"])
    return msg.format(fmt, use_key, shape)


def identity_msg(dataset, do_print=False):
    header = _header_msg(dataset)
    groupby = _groupby_msg(dataset)
    if groupby:
        header = "\n".join([header, groupby])
    msg = "\n".join([header, _X_msg(dataset, dataset._use_key)])
    attrs = _all_attr_msg(dataset)
    if attrs:
        msg = "\n".join([msg, attrs])
    if not do_print:
        return msg
    print(msg)