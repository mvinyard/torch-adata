
__module_name__ = "_identity_msg.py"
__doc__ = """Print the identity of the AnnDataset class."""
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu"])


# -- import packages: --------------------------------------------------------------------
import licorice_font


from licorice_font import font_format as ff


def _header_msg(dataset):
    l = len(dataset)
    msg = "[ {}-{} ]: AnnDataset object with {} samples\n{}"
    t, a = ff("torch", ["RED", "BOLD"]), ff("adata", ["PURPLE", "BOLD"])
    u = "".join(["-"] * int(48 + len(str(l))))
    return msg.format(t, a, l, u)


def _groupby_msg(dataset):

    if not dataset._grouped_by:
        return
    msg_g1 = ff("Grouped by", ["BOLD"])
    msg_g2 = ff(dataset._grouped_by, ["RED"])
    return "{}: '{}' with attributes:".format(msg_g1, msg_g2)


def _attr_msg(attr_name_set, title):
    if attr_name_set:
        return " - {:<18}".format(ff(title + ":", ["BOLD"])) + ", ".join(attr_name_set)


def _all_attr_msg(dataset):

    aux_msg = _attr_msg(dataset._aux_attr_names, title="X_aux")
    obs_msg = _attr_msg(dataset._obs_attr_names, title="obs")
    if any([aux_msg, obs_msg]):
        return "\n".join([aux_msg, obs_msg])


def _X_msg(use_key):
    msg = "{} (use_key = '{}')"
    return msg.format(ff(" - X", ["BOLD"]), use_key)


def identity_msg(dataset, do_print=False):
    header = _header_msg(dataset)
    groupby = _groupby_msg(dataset)
    if groupby:
        header = "\n".join([header, groupby])
    msg = "\n".join([header, _X_msg(dataset._use_key)])
    attrs = _all_attr_msg(dataset)
    if attrs:
        msg = "\n".join([msg, attrs])
    if not do_print:
        return msg
    print(msg)