"""rpy2 bridge. Imported once so every module talks to the same R session."""
import numpy as np
import rpy2.robjects as ro
from rpy2.robjects import numpy2ri
from rpy2.robjects.packages import importr

limma = importr("limma")
NP2R = ro.default_converter + numpy2ri.converter


def to_r(arr):
    """numpy -> R. Converted explicitly (instead of inside a converter context that also wraps the
    call) because the automatic conversion of R lists changed between rpy2 versions."""
    with NP2R.context():
        return ro.conversion.get_conversion().py2rpy(np.asarray(arr, dtype=float))


def r_list_to_dict(r_out):
    """R named list -> {name: numpy array}; works across rpy2 versions."""
    return {name: np.array(r_out.rx2(name), dtype=float) for name in r_out.names}
