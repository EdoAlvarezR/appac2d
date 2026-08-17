"""
Aeropropulsive panel code, a Python wrapper of the OCTAVE source code
"""

__version__ = "0.0.1"

__all__ = [
    "__version__",
]

# ----------- IMPORT OCTAVE BINARY ---------------------------------------------
import os
_DIR = os.path.dirname(os.path.abspath(__file__))

appac_path = os.path.join(_DIR, "..")

from oct2py import octave

octave.addpath(appac_path)
octave.addpath(os.path.join(appac_path, "mesh2d"))
octave.initmsh()


# ----------- IMPORT AND FORMAT PYPLOT -----------------------------------------
import matplotlib.pyplot as plt

plt.rc("font", family="STIXGeneral")            # Text font
plt.rc("mathtext", fontset="stix")              # Math font
plt.rc("font", size=16)          # controls default text sizes
plt.rc("axes", titlesize=16)     # fontsize of the axes title
plt.rc("axes", labelsize=18)     # fontsize of the x and y labels
plt.rc("xtick", labelsize=16)    # fontsize of the tick labels
plt.rc("ytick", labelsize=16)    # fontsize of the tick labels
plt.rc("legend", fontsize=16)    # legend fontsize
plt.rc("figure", titlesize=22)   # fontsize of the figure title


# ----------- IMPORT SUBMODULES ------------------------------------------------
from appc import preprocess, postprocess, visualization, utils