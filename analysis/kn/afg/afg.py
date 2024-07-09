import shutil,json,os,h5py
import pandas as pd
import numpy as np
from scipy.integrate import ode
import matplotlib.pyplot as plt
from scipy import interpolate

from pathlib import Path
from matplotlib.colors import Normalize, LogNorm
from matplotlib import cm

try:
    import package.src.PyBlastAfterglowMag as PBA
except ImportError:
    try:
        import PyBlastAfterglowMag as PBA
    except:
        raise ImportError("Cannot import PyBlastAfterglowMag")

DATA_PATH = str(__file__).split("analysis/kn/")[0] + "analysis/kn/kenta_data/"
get_ej_data = lambda name : DATA_PATH+name+'/'+"ej_collated.h5"

# load the metadata
with open(DATA_PATH+"metadata.json") as json_file:
    json_data = json.load(json_file)
    print(json_data)
SIMS = pd.DataFrame.from_dict(json_data).T
SIMS.set_index("name")
# select only new simulations
df = SIMS[SIMS["given_time"] == "new"]

def main(sim_dict:pd.Series):

    pba = PBA.wrappers.run_kn(
        working_dir=os.getcwd()+'/'+'working_dir/',
        struct=dict(struct="numeric",
                    n_layers_pw=None,
                    corr_fpath_kenta=get_ej_data(sim_dict['name']),
                    text=25,
                    t0=1e3
                    ),
        P=dict(main=dict(n_ism=1e-1,d_l=100e6 * PBA.utils.cgs.pc,z=0.001,theta_obs=0.,integrator="DOP853"),
               kn=dict(p_fs=2.5,eps_e_fs=1.e-1,eps_b_fs=1.e-1,
                       method_synchrotron_fs="Joh06",method_ele_fs="analytic",
                       use_1d_id="no",do_skymap="no", do_lc="yes",
                       ebl_tbl_fpath="none",method_spread="None",method_nonrel_dist_fs="use_Sironi",
                       method_ne_fs="usenprime",method_comp_mode="observFlux"
                       )),
        run=False)

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(4.6, 3.2))
    ax.plot(pba.KN.get_lc_times() / PBA.utils.cgs.day, pba.KN.get_lc(freq=3.e9) * 1e3,
            color=sim_dict["color"],ls=sim_dict["ls"],lw=1,label=sim_dict["label"])
    for ishell in [10,20,30,40,50,60,70,80,90,98]:
        ax.plot(pba.KN.get_lc_times() / PBA.utils.cgs.day, pba.KN.get_lc(freq=3.e9,ishell=ishell) * 1e3,
                color='gray',ls='-',lw=0.7)
    pba.clear()

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.minorticks_on()
    ax.set_xlabel(r"$t_{\rm obs}$ [day]")
    ax.set_ylabel(r"$F_{\nu}$ [$\mu$Jy]")
    ax.set_xlim(1e-1, 1e5)
    ax.set_ylim(1e-1, 1e3)
    ax.legend()

    plt.show()

if __name__ == '__main__':
    main(df.loc["SFHo_13_14_res150"])