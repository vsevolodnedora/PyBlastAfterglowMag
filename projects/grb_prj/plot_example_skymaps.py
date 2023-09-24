# import PyBlastAfterglowMag
import numpy as np
import h5py
from scipy.integrate import ode
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.axes3d import Axes3D
from matplotlib import gridspec
from scipy import interpolate
# import pytest
from pathlib import Path
from matplotlib.colors import Normalize, LogNorm
from matplotlib import cm
import os

# import PyBlastAfterglowMag as PBA
import package.src.PyBlastAfterglowMag as PBA
# plt.style.use('seaborn-v0_8')

# plt.style.use('dark_background')