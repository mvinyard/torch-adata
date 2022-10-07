
__module_name__ = "_identity_msg.py"
__doc__ = """Print the identity of the AnnDataset class."""
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu"])


# -- import packages: --------------------------------------------------------------------
import licorice_font


def identity_msg(dataset):
    msg_01 = licorice_font.font_format("torch", ["RED", "BOLD"])
    msg_02 = licorice_font.font_format("adata", ["PURPLE", "BOLD"])
    msg = "[ {}-{} ]: AnnDataset object with {} samples".format(
        msg_01, msg_02, len(dataset)
    )
    underline = "".join(["-"] * int(48 + len(str(len(dataset)))))
    header = "{}\n{}".format(msg, underline)

    if dataset._data_axis:
        msg_g1 = licorice_font.font_format("Grouped by", ["BOLD"])
        msg_g2 = licorice_font.font_format(dataset._grouped_by, ["RED"])
        groupby_description = "{}: '{}' with attributes:".format(msg_g1, msg_g2)
    else:
        groupby_description = ""
    attr_print_list = []
    for attr in dataset._attr_names:
        print_attr = licorice_font.font_format(attr, ["BOLD"])
        attr_print_list.append(" - {}: {}".format(print_attr, getattr(dataset, attr).shape))
    
    return "\n".join([header, groupby_description] + attr_print_list)
