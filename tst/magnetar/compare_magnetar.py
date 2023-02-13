# import PyBlastAfterglowMag
import shutil

import numpy as np
import h5py
from scipy.integrate import ode
import matplotlib.pyplot as plt
from scipy import interpolate
# import pytest
import os

try:
    from PyBlastAfterglowMag.interface import modify_parfile_par_opt
    from PyBlastAfterglowMag.interface import BPA_METHODS
    from PyBlastAfterglowMag.interface import (distribute_and_run, get_str_val, set_parlists_for_pars)
    from PyBlastAfterglowMag.utils import latex_float, cgs, get_beta, get_Gamma
except ImportError:
    try:
        from package.src.PyBlastAfterglowMag.interface import modify_parfile_par_opt
        from package.src.PyBlastAfterglowMag.interface import BPA_METHODS
        from package.src.PyBlastAfterglowMag.interface import (distribute_and_run, get_str_val, set_parlists_for_pars)
        from package.src.PyBlastAfterglowMag.utils import latex_float, cgs, get_beta, get_Gamma
    except ImportError:
        raise ImportError("Cannot import PyBlastAfterglowMag")




def main():
    workdir = os.getcwd()+'/'
    parfile_name = "parfile.par"
    pba = BPA_METHODS(workingdir=workdir, readparfileforpaths=True, parfile=parfile_name)
    pba.run()
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