# import PyBlastAfterglowMag
import shutil

import numpy as np
import h5py
from scipy.integrate import ode
import matplotlib.pyplot as plt
from scipy import interpolate
# import pytest
from pathlib import Path
from matplotlib.colors import Normalize, LogNorm
from matplotlib import cm
import os

from package.src.PyBlastAfterglowMag.interface import BPA_METHODS as PBA
from package.src.PyBlastAfterglowMag.interface import cgs, latex_float
from package.src.PyBlastAfterglowMag.id_maker_from_thc_ourflow import prepare_kn_ej_id_1d, prepare_kn_ej_id_2d

def main():
    pba = PBA(os.getcwd()+"/",readparfileforpaths=True)
    mag = pba.get_mag_obj()
    print(np.array(mag["n_grav"]))
    print(mag.keys())
    # plt.loglog(np.array(mag["tburst"]), np.array(mag["mdot"])/cgs.solar_m)
    plt.loglog(np.array(mag["tburst"])[1:],np.array(mag["lprop"])[1:], color='blue', label=r"$L_{\rm prop}$")
    plt.loglog(np.array(mag["tburst"])[1:],np.array(mag["ldip"])[1:], color='red', label=r"$L_{\rm dip}$")
    # plt.loglog(np.array(mag["tburst"]),-np.array(mag["n_dip"]), color='blue', label=r"$n_{\rm dip}$")
    # plt.loglog(np.array(mag["tburst"]),-np.array(mag["n_acc"]), color='red', label=r"$n_{\rm acc}$")
    # plt.loglog(np.array(mag["tburst"]),-np.array(mag["n_grav"]), color='green', label=r"$n_{\rm grav}$")
    plt.xlabel("time [s]")
    plt.ylabel("Luminosity [erg/s]")
    plt.legend()
    plt.title("Magnetar spindown")
    plt.ylim(1e47, 1e52)
    plt.show()


if __name__ == '__main__':
    main()