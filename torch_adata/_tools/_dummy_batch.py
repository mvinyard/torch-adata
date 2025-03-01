
__module_name__ = "_dummy_batch.py"
__doc__ = """Randomly sample a dummy batch"""
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu"])

# -- import local dependencies: ------------------------------------------------
from ._idx import idx


# -- supporting functions: -----------------------------------------------------
def _print_batch_details(batch, attr_names):
    """print the dummy batch details"""
    
    print("Batch:")
    for n, key in enumerate(attr_names):
        print(" - {}:  {}".format(key, batch[n].shape), end="\n")


# -- main module function: -----------------------------------------------------
def dummy_batch(dataset, silent=False):
    """Fetch a dummy batch"""
    
    dummy_idx = idx(dataset)
    dummy_batch = dataset[dummy_idx]
    if not silent:
        _print_batch_details(dummy_batch, dataset._attr_names)
        
    return dummy_batch