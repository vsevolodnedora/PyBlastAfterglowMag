import numpy as np
import matplotlib.pyplot as plt
import h5py
import os

from mpl_toolkits.axes_grid1.inset_locator import (inset_axes, InsetPosition, mark_inset)
from matplotlib.colors import Normalize, LogNorm
from matplotlib import ticker, cm, rc, rcParams
import matplotlib.colors as colors


import os
import math
import mpl_toolkits.axisartist.floating_axes as floating_axes
from matplotlib.projections import PolarAxes
from mpl_toolkits.axisartist.grid_finder import FixedLocator, MaxNLocator, DictFormatter
import random

from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle

from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator

rc('text', usetex=True) # $ sudo apt-get install cm-super
rc('font', family='serif')
rcParams['font.size'] = 10

curdir = os.getcwd() + '/' #"/home/vsevolod/Work/GIT/GitHub/PyBlastAfterglow_dev/PyBlastAfterglow/src/PyBlastAfterglow/tests/dyn/"

# paperdir = "/home/vsevolod/Work/GIT/overleaf/61bd9127bbc4c5790d6db8c2/figs/"

class cgs:

    pi = 3.141592653589793

    tmp = 1221461.4847847277

    c = 2.9979e10
    mp = 1.6726e-24
    me = 9.1094e-28
    # e = 1.602176634e-19 # Si
    # h = 6.62607015e-34 # ???? Si
    h = 6.6260755e-27 # erg s
    mpe = mp + me
    mue = me / mp
    hcgs = 6.6260755e-27  # Planck constant in cgs
    # si_h = 6.62607015e-34
    kB = 1.380658e-16
    sigmaT = 6.6524e-25
    qe = 4.803204e-10
    # si_qe = 1.602176634e-19
    sigma_B = 5.6704e-5  ### Stephan Boltzann constant in cgs
    lambda_c = (h / (me * c)) # Compton wavelength
    mecc2MeV = 0.511
    mec2 = 8.187105649650028e-07  # erg # electron mass in erg, mass_energy equivalence
    gamma_c_w_fac = 6 * np.pi * me * c / sigmaT
    rad_const = 4 * sigma_B / c   #### Radiation constant
    mppme = mp + me

    pc = 3.0857e18 # cm
    year= 3.154e+7 # sec
    day = 86400

    solar_m = 1.989e+33

    ns_rho = 1.6191004634e-5
    time_constant = 0.004925794970773136  # to to to ms
    energy_constant = 1787.5521500932314
    volume_constant = 2048

    sTy = 365. * 24. * 60. * 60.    # seconds to years conversion factor
    sTd = 24. * 60. * 60.           # seconds to days conversion factor
    rad2mas = 206264806.247

# def rho(r,mej,rej,delta=1):
#     ''' Eq. 2 in doi:10.1093/mnras/stt1922 '''
#     return (3.-delta) / 4.*np.pi * mej / (rej**3) * (r/rej)**(-delta)
#
# x_arr = np.array([ 2.99792e+09, 5.96587e+09, 8.93382e+09, 1.19018e+10, 1.48697e+10, 1.78377e+10, 2.08056e+10, 2.37735e+10, 2.67415e+10, 2.97094e+10, 3.26774e+10, 3.56453e+10, 3.86133e+10, 4.15812e+10, 4.45492e+10, 4.75171e+10, 5.04851e+10, 5.3453e+10, 5.64209e+10, 5.93889e+10, 6.23568e+10, 6.53248e+10, 6.82927e+10, 7.12607e+10, 7.42286e+10, 7.71966e+10, 8.01645e+10, 8.31324e+10, 8.61004e+10, 8.90683e+10, 9.20363e+10, 9.50042e+10, 9.79722e+10, 1.0094e+11, 1.03908e+11, 1.06876e+11, 1.09844e+11, 1.12812e+11, 1.1578e+11, 1.18748e+11, 1.21716e+11, 1.24684e+11, 1.27652e+11, 1.3062e+11, 1.33588e+11, 1.36555e+11, 1.39523e+11, 1.42491e+11, 1.45459e+11, 1.48427e+11, 1.51395e+11, 1.54363e+11, 1.57331e+11, 1.60299e+11, 1.63267e+11, 1.66235e+11, 1.69203e+11, 1.72171e+11, 1.75139e+11, 1.78107e+11, 1.81075e+11, 1.84043e+11, 1.87011e+11, 1.89978e+11, 1.92946e+11, 1.95914e+11, 1.98882e+11, 2.0185e+11, 2.04818e+11, 2.07786e+11, 2.10754e+11, 2.13722e+11, 2.1669e+11, 2.19658e+11, 2.22626e+11, 2.25594e+11, 2.28562e+11, 2.3153e+11, 2.34498e+11, 2.37466e+11, 2.40434e+11, 2.43401e+11, 2.46369e+11, 2.49337e+11, 2.52305e+11, 2.55273e+11, 2.58241e+11, 2.61209e+11, 2.64177e+11, 2.67145e+11, 2.70113e+11, 2.73081e+11, 2.76049e+11, 2.79017e+11, 2.81985e+11, 2.84953e+11, 2.87921e+11, 2.90889e+11])
# y_arr = np.array([ 0.0338132, 0.0239693, 0.0181677, 0.014473, 0.0121081, 0.0103214, 0.00912469, 0.00786751, 0.0069139, 0.00603652, 0.00497203, 0.00423341, 0.00387051, 0.00359267, 0.00325797, 0.0029195, 0.00263854, 0.00224722, 0.00201092, 0.00182701, 0.00164119, 0.0013476, 0.00114731, 0.000959042, 0.000804705, 0.000716644, 0.000669936, 0.000618169, 0.000531163, 0.000462963, 0.000401329, 0.000342046, 0.000296272, 0.000239405, 0.000198068, 0.000170822, 0.000149442, 0.000129222, 0.00011397, 0.000103007, 9.35766e-05, 8.54766e-05, 7.75451e-05, 7.10069e-05, 6.34217e-05, 5.66388e-05, 5.18301e-05, 4.73526e-05, 4.28353e-05, 3.78702e-05, 3.39341e-05, 2.99516e-05, 2.70452e-05, 2.43674e-05, 2.15786e-05, 1.7089e-05, 1.40623e-05, 1.03363e-05, 7.36981e-06, 5.40693e-06, 3.55046e-06, 2.16996e-06, 1.36508e-06, 1.07368e-06, 9.22994e-07, 8.45615e-07, 7.65885e-07, 6.94481e-07, 6.09351e-07, 5.36406e-07, 4.74121e-07, 4.12352e-07, 3.71102e-07, 3.2306e-07, 2.81433e-07, 2.33967e-07, 2.0114e-07, 1.70503e-07, 1.42518e-07, 1.2425e-07, 1.0948e-07, 9.00745e-08, 7.45861e-08, 5.83542e-08, 4.72846e-08, 3.91133e-08, 2.94786e-08, 2.79353e-08, 2.81956e-08, 1.78446e-08, 8.40481e-09, 3.66657e-09, 1.39995e-09, 6.68156e-10, 4.16988e-10, 2.70648e-10, 8.30977e-11, 0])
# plt.loglog(x_arr, y_arr, ls='--', label='SFHo q=1.00 t=0')
# plt.loglog(x_arr, rho(x_arr,7.8616e+29,1e11), ls='--', label='Kasen \& Bildsten 2010')
# plt.axhline(y=6.86255e-05)
# plt.xlabel(r"Radius [cm]")
# plt.ylabel(r"Density $[{\rm g/cm }^3]$")
# plt.title("Comparison between analytic and numeric density profile")
# plt.legend()
# plt.show()


def Beta(Gamma):
    return np.sqrt(1. - np.power(Gamma, -2))
def load_data(fname):
    with open(fname) as f:
        first_line = f.readline()
    return (first_line.split()[1:], np.loadtxt(fname))

def plot_jet_layers():
    layers = ["layer=0", "layer=20", "layer=30", "layer=40", "layer=50", "layer=60"]
    # v_ns = ["Gamma"]

    dfile = h5py.File(curdir+"jet_dynamics_layers.h5", "r")
    print(dfile.keys())
    print(dfile["layer=0"]["M2"])

    fid, axes = plt.subplots(ncols=4, nrows=1, figsize=(8,2))
    norm = Normalize(vmin=0, vmax=dfile.attrs["nlayers"])
    cmap = cm.viridis

    v_n_x = "R"
    v_n_ys = ["Gamma", "M2", "tt", "theta"]
    for iv_n, v_n in enumerate(v_n_ys):


        for il, layer in enumerate(layers):
            x_arr = np.array(dfile[layer][v_n_x])
            y_arr = np.array(dfile[layer][v_n])
            axes[iv_n].plot(x_arr, y_arr, ls='-', color=cmap(norm(int(layer.split("=")[-1]))))

        axes[iv_n].set_xlabel(v_n_x)
        axes[iv_n].set_ylabel(v_n)
        axes[iv_n].set_xscale("log")
        axes[iv_n].set_yscale("log")
        # axes[iv_n].legend()
        axes[iv_n].grid()

    plt.show()
# plot_jet_layers()

# layer by layer for selected shells
def plot_ejecta_layers(ishells=(1,), ilayers=(0,10,22), v_n_x = "R", v_n_ys = ("rho", "mom"), colors_by="layers",legend=False):
    layers = []
    for i in ishells:
        for j in ilayers:
            layers.append("shell={} layer={}".format(i,j))

    # v_ns = ["Gamma"]

    dfile = h5py.File(curdir+"magnetar_driven_ej.h5", "r")
    print(dfile.keys())
    # print(dfile["layer=0"]["M2"])

    fid, axes = plt.subplots(ncols=1, nrows=len(v_n_ys), figsize=(6,6),sharex="all")
    # norm = Normalize(vmin=0, vmax=dfile.attrs["nlayers"])
    cmap = cm.viridis
    mynorm = Normalize(vmin=0,vmax=len(ishells)*len(ilayers))#norm(len(ishells)*len(ilayers))

    for iv_n, v_n in enumerate(v_n_ys):
        i = 0
        ax = axes[iv_n] if len(v_n_ys) > 1 else axes
        for il, layer in enumerate(layers):
            x_arr = np.array(dfile[layer][v_n_x])
            y_arr = np.array(dfile[layer][v_n])
            y_arr = y_arr[x_arr > 0]
            x_arr = x_arr[x_arr > 0]
            # if (v_n == "R"):
                # y_arr = y_arr/y_arr.max()
            if (colors_by=="layers"): color=cmap(mynorm(int(i)))#color=cmap(norm(int(layer.split("layer=")[-1])))
            else: color=cmap(mynorm(int(i)))#color=cmap(norm(int(layer.split("shell=")[-1].split("layer=")[0])))
            if (v_n_x == "tburst"): x_arr /=cgs.day;
            ax.plot(x_arr, y_arr, ls='-', color=color, label=layer)
            i=i+1
        ax.set_xlabel(v_n_x)
        if (v_n_x == "tburst"): ax.set_xlabel(v_n_x + " [day]")
        ax.set_ylabel(v_n)
        ax.set_xscale("log")
        ax.set_yscale("log")
        # axes[iv_n].legend()
        ax.grid()
        ax.set_xscale("log")
        ax.set_yscale("log")
        # ax.set_xlim(5e-1,4e3)
    if legend: plt.legend()
    plt.savefig("./many_bw_dynamics.png", dpi=256)
    plt.show()
# plot_ejecta_layers(ishells=(90,91,92,93,94,95,95,96,97), ilayers=(0,), v_n_x = "tburst", v_n_ys = ("Eint2", "mom", "R", "M2"), colors_by="shell")
# plot_ejecta_layers(ishells=([i for i in range(98)]), ilayers=(0,), v_n_x = "tburst", v_n_ys = ("Eint2", "mom"), colors_by="shell")
# plot_ejecta_layers(ishells=([i for i in range(98)]), ilayers=(0,), v_n_x = "tburst", v_n_ys = (["R"]), colors_by="shell")

def plot_ej_and_magnetar_layers(ishells=(1,), ilayers=(0,10,22),
                                v_n_x = "R", v_n_ys = ("rho"),
                                pwn_v_n_ys = ("Epwn"),
                                colors_by="layers",legend=False,
                                scale = False, same_ax = False,
                                figname="./pwn.png"):

    layers = []
    for j in ilayers:
        for i in ishells:
            layers.append("shell={} layer={}".format(i,j))

    dfile = h5py.File(curdir+"magnetar_driven_ej.h5", "r")
    # print(dfile[list(dfile.keys())[0]].keys())
    dfile_pwn = h5py.File(curdir+"pwn.h5","r")
    # print(dfile_pwn[list(dfile_pwn.keys())[0]].keys())
    print(dfile_pwn.keys())
    # print(dfile["layer=0"]["M2"])

    if (same_ax): nrows = 1
    else: nrows = len(v_n_ys)+len(pwn_v_n_ys)
    fid, axes = plt.subplots(ncols=1, nrows=nrows, figsize=(6,6),sharex="all")
    if not hasattr(axes,"__len__"): axes = [axes]
    # norm = Normalize(vmin=0, vmax=dfile.attrs["nlayers"])
    cmap = cm.Reds_r
    mynorm = Normalize(vmin=0,vmax=len(ishells)*len(ilayers))#norm(len(ishells)*len(ilayers))

    for iv_n, v_n in enumerate(v_n_ys):
        i = 0
        ax = axes[iv_n] if (len(v_n_ys)+len(pwn_v_n_ys)) >= 1 else axes
        for il, layer in enumerate(layers):
            x_arr = np.array(dfile[layer+f" key={v_n_x}"])
            y_arr = np.array(dfile[layer+f" key={v_n}"])
            # if (v_n == "R"):
            # y_arr = y_arr/y_arr.max()
            if (colors_by=="layers"): color=cmap(mynorm(int(i)))#color=cmap(norm(int(layer.split("layer=")[-1])))
            else: color=cmap(mynorm(int(i)))#color=cmap(norm(int(layer.split("shell=")[-1].split("layer=")[0])))
            if (v_n_x == "tburst"): x_arr /=cgs.day;
            if (v_n_x == "tburst" and v_n == "R" and scale):
                y_arr_ = np.array(dfile[layers[-1]+f" key={v_n}"])
                y_arr = y_arr/y_arr_
            # if (v_n_x == "tburst" and v_n == "mom"):
            #     y_arr_ = np.array(dfile[layer][v_n])[0]
            #     y_arr = y_arr/y_arr_
            y_arr = y_arr[x_arr > 0]
            x_arr = x_arr[x_arr > 0]
            if (len(x_arr)>0 and len(y_arr)>0):
                ax.plot(x_arr, y_arr, ls='-', color=color, label="Ej. "+layer)
            else:
                print("Empty shell?")
            # ------------
            i=i+1

        ax.set_xlabel(v_n_x)
        if (v_n_x == "tburst"): ax.set_xlabel(v_n_x + " [day]")
        ax.set_ylabel(v_n)
        ax.set_xscale("log")
        ax.set_yscale("log")
        # axes[iv_n].legend()
        ax.grid()
        ax.set_xscale("log")
        ax.set_yscale("log")
        # ax.set_xlim(1e-4,1e-1)
        # ax.set_ylim(1,1.012)

    cmap = cm.Blues_r
    layers = []
    for j in ilayers:
        for i in [0]:
            layers.append("shell={} layer={}".format(i,j))

    if same_ax: offset = 0
    else: offset = len(v_n_ys)
    for iv_n, v_n in enumerate(pwn_v_n_ys):
        i = 0
        ax = axes[offset + iv_n] if len(v_n_ys) + len(pwn_v_n_ys) > 1 else axes
        for il, layer in enumerate(layers):
            x_arr = np.array(dfile_pwn[layer+f" key={v_n_x}"])
            y_arr = np.array(dfile_pwn[layer+f" key={v_n}"])
            y_arr = y_arr[x_arr > 0]
            x_arr = x_arr[x_arr > 0]
            if (colors_by=="layers"): color=cmap(mynorm(int(i)))#color=cmap(norm(int(layer.split("layer=")[-1])))
            else: color=cmap(mynorm(int(i)))#color=cmap(norm(int(layer.split("shell=")[-1].split("layer=")[0])))
            if (v_n_x == "tburst"): x_arr /=cgs.day;
            if (len(x_arr)>0 and len(y_arr)>0):
                ax.plot(x_arr, y_arr, ls='-', color=color, label="PWN "+layer)
            # ------------
            i=i+1
        ax.set_xlabel(v_n_x)
        if (v_n_x == "tburst"): ax.set_xlabel(v_n_x + " [day]")
        ax.set_ylabel(v_n)
        ax.set_xscale("log")
        ax.set_yscale("log")
        # axes[iv_n].legend()
        ax.grid()
        ax.set_xscale("log")
        ax.set_yscale("log")

    if legend: plt.legend()
    plt.savefig(figname, dpi=256)
    plt.show()
# plot_ej_and_magnetar_layers(ishells=([0]), ilayers=(0,), v_n_x = "tburst",
#                             v_n_ys = (["mom"]), pwn_v_n_ys=([]), same_ax=True, legend=True,
#                             colors_by="shell",figname="./pwn.png")
# plot_ej_and_magnetar_layers(ishells=([0]), ilayers=(9,), v_n_x = "tburst",
#                             v_n_ys = (["mom"]), pwn_v_n_ys=(["mom"]), same_ax=True, legend=True,
#                             colors_by="shell",figname="./pwn.png")
plot_ej_and_magnetar_layers(ishells=([i for i in range(60)]), ilayers=(0,), v_n_x = "tburst",
                            v_n_ys = (["mom","Eint2"]), pwn_v_n_ys=([]),
                            colors_by="shell",figname="./pwn_driv_ejecta.png")


def plot_dynamics_layers(jet_layers=(0,20,40,60,69), v_n_x = "R", v_n_ys = ("rho", "mom", "tburst"),
                         ej_shell=(0), ej_layers=(0,2,4,6,8,10,12,14,16,18,20,22)):

    # -------- JET ---------
    layers = []
    for il in jet_layers:
        layers.append("layer={}".format(il))
    # v_ns = ["Gamma"]

    dfile = h5py.File(curdir+"jet_dynamics_layers.h5", "r")
    print(dfile.keys())
    print(dfile["layer=0"].keys())
    print(dfile["layer=0"]["U_e"][0])
    print("TT=({},{})".format(dfile["layer=0"]["tt"][0],dfile["layer=0"]["tt"][-1]))

    fid, axes = plt.subplots(ncols=1, nrows=3, figsize=(8,5),sharex="all")
    norm = Normalize(vmin=0, vmax=len(dfile.keys()))
    cmap = cm.Blues

    for iv_n, v_n in enumerate(v_n_ys):


        for il, layer in enumerate(layers):
            x_arr = np.array(dfile[layer][v_n_x])
            y_arr = np.array(dfile[layer][v_n])
            axes[iv_n].plot(x_arr, y_arr, '.', color=cmap(norm(int(layer.split("=")[-1]))))

        axes[iv_n].set_xlabel(v_n_x)
        axes[iv_n].set_ylabel(v_n)
        axes[iv_n].set_xscale("log")
        axes[iv_n].set_yscale("log")
        # axes[iv_n].legend()
        axes[iv_n].grid()


    # ------- EJECTA ----------
    layers = []
    for i in ej_shell:
        for j in ej_layers:
            layers.append("shell={} layer={}".format(i,j))

    dfile = h5py.File(curdir+"ejecta_dynamics_shells_layers.h5", "r")
    print(dfile.keys())

    ishell = 20
    ilayer = 23
    print("Beta0={}".format( dfile["shell={} layer={}".format(ishell,ilayer)]["beta"][2]))
    print("M2={}".format( dfile["shell={} layer={}".format(ishell,ilayer)]["M2"][2]))
    print("rho2={}".format( dfile["shell={} layer={}".format(ishell,ilayer)]["rho2"][2]))
    print("U_e={}".format( dfile["shell={} layer={}".format(ishell,ilayer)]["U_e"][2]))
    # print("B={}".format( dfile["shell={} layer={}".format(ishell,ilayer)]["B"][0]))




    # exit(1)
    # print(dfile["layer=0"]["M2"])

    # fid, axes = plt.subplots(ncols=4, nrows=1, figsize=(8,2))
    norm = Normalize(vmin=0, vmax=dfile.attrs["nlayers"])
    cmap = cm.Reds

    #
    # v_n_x = "R"
    # v_n_ys = ["rho", "Gamma", "tburst"]
    for iv_n, v_n in enumerate(v_n_ys):

        for il, layer in enumerate(layers):
            x_arr = np.array(dfile[layer][v_n_x])
            y_arr = np.array(dfile[layer][v_n])
            axes[iv_n].plot(x_arr, y_arr, '.', color=cmap(norm(int(layer.split("layer=")[-1]))))

        axes[iv_n].set_xlabel(v_n_x)
        axes[iv_n].set_ylabel(v_n.replace("_","\_"))
        axes[iv_n].set_xscale("log")
        axes[iv_n].set_yscale("log")
        # axes[iv_n].legend()
        axes[iv_n].grid()

    plt.tight_layout()
    plt.show()
# plot_dynamics_layers(jet_layers=(0,20,40,60,69), v_n_x = "R", v_n_ys = ("rho", "mom", "tburst"),
#                      ej_shell=(1,20),ej_layers=(0,10,20))


def plot_dynamics_pretty_plot():
    # -------- JET ---------

    v_n_x = "tburst"
    v_n_y = "R"
    def plot_jet(ax, xlim, ylim):
        dfile = h5py.File(curdir+"jet_dynamics_layers.h5", "r")
        print(dfile.keys())
        layers = ["layer=0",
                  "layer=19",
                  "layer=49",
                  "layer=69"]
        dfile = h5py.File(curdir+"jet_dynamics_layers.h5", "r")
        norm = Normalize(vmin=0, vmax=len( dfile.keys() ))
        cmap = cm.Blues
        for il, layer in enumerate(layers):
            x_arr = np.array(dfile[layer][v_n_x])
            if (v_n_y == "mom"):
                y_arr = np.array(dfile[layer]["Gamma"]) * np.array(dfile[layer]["beta"])
            else:
                y_arr = np.array(dfile[layer][v_n_y])
            if (len(xlim)>0):
                y_arr = y_arr[(x_arr > xlim[0]) & (x_arr < xlim[-1])]
                x_arr = x_arr[(x_arr > xlim[0]) & (x_arr < xlim[-1])]
            if (len(ylim)>0):
                x_arr = x_arr[(y_arr > ylim[0]) & (y_arr < ylim[-1])]
                y_arr = y_arr[(y_arr > ylim[0]) & (y_arr < ylim[-1])]
            ax.plot(x_arr, y_arr, ls='-', color=cmap(norm(int(layer.split("=")[-1]))))

    def plot_ejecta(ax, xlim, ylim):
        ishell = 10
        layers = ["shell={} layer=0".format(ishell),
                  # "shell={} layer=2".format(ishell),
                  "shell={} layer=4".format(ishell),
                  # "shell={} layer=6".format(ishell),
                  "shell={} layer=8".format(ishell),
                  # "shell={} layer=10".format(ishell),
                  "shell={} layer=12".format(ishell),
                  # "shell={} layer=14".format(ishell),
                  "shell={} layer=16".format(ishell),
                  # "shell={} layer=18".format(ishell),
                  "shell={} layer=20".format(ishell),
                  # "shell={} layer=22".format(ishell),
                  "shell={} layer=23".format(ishell)
                  ]
        ishell = 15
        layers2 = ["shell={} layer=0".format(ishell),
                   # "shell={} layer=2".format(ishell),
                   "shell={} layer=4".format(ishell),
                   # "shell={} layer=6".format(ishell),
                   "shell={} layer=8".format(ishell),
                   # "shell={} layer=10".format(ishell),
                   "shell={} layer=12".format(ishell),
                   # "shell={} layer=14".format(ishell),
                   "shell={} layer=16".format(ishell),
                   # "shell={} layer=18".format(ishell),
                   "shell={} layer=20".format(ishell),
                   # "shell={} layer=22".format(ishell),
                   "shell={} layer=23".format(ishell)
                   ]
        layers+=layers2
        dfile = h5py.File(curdir+"ejecta_dynamics_shells_layers.h5", "r")
        norm = Normalize(vmin=0, vmax=dfile.attrs["nlayers"])
        cmap = cm.Reds_r
        for il, layer in enumerate(layers):
            x_arr = np.array(dfile[layer][v_n_x])
            if (v_n_y == "mom"):
                y_arr = np.array(dfile[layer]["Gamma"]) * np.array(dfile[layer]["beta"])
            else:
                y_arr = np.array(dfile[layer][v_n_y])
            if (len(xlim)>0):
                y_arr = y_arr[(x_arr > xlim[0]) & (x_arr < xlim[-1])]
                x_arr = x_arr[(x_arr > xlim[0]) & (x_arr < xlim[-1])]
            if (len(ylim)>0):
                x_arr = x_arr[(y_arr > ylim[0]) & (y_arr < ylim[-1])]
                y_arr = y_arr[(y_arr > ylim[0]) & (y_arr < ylim[-1])]
            ax.plot(x_arr, y_arr, ls='-', color=cmap(norm(int(layer.split("=")[-1]))))


    # T, Cp = np.loadtxt('Ta-Cp.txt', unpack=True)
    # T_E, CV_E = np.loadtxt('Ta-CV_Einstein.txt', unpack=True)
    # T_D, CV_D = np.loadtxt('Ta-CV_Debye.txt', unpack=True)

    fig, ax1 = plt.subplots(figsize=(4,3))
    # T_E = np.arange(1,max(T)+1,1)
    # The data.
    plot_jet(ax1, (), ())
    plot_ejecta(ax1, (), ())
    # ax1.plot(T, Cp, 'x', c='b', mew=2, alpha=0.8, label='Experiment')
    # The Einstein fit.
    # ax1.plot(T_E, CV_E, c='m', lw=2, alpha=0.5, label='Einstein model')
    ax1.set_xlabel(r'Radius [cm]')#(r'$T\,/\mathrm{K}$')
    ax1.set_ylabel(r'Momentum $\Gamma\beta$')#(r'$C_p\,/\mathrm{J\,K^{-1}\,mol^{-1}}$')
    ax1.legend(loc=0)

    # Create a set of inset Axes: these should fill the bounding box allocated to
    # them.
    ax2 = plt.axes([0,0,2,2])
    # Manually set the position and relative size of the inset axes within ax1
    ip = InsetPosition(ax1, [0.20,0.1,0.4,0.4]) # [0.4,0.2,0.5,0.5]
    ax2.set_axes_locator(ip)
    # Mark the region corresponding to the inset axes on ax1 and draw lines
    # in grey linking the two axes.
    mark_inset(ax1, ax2, loc1=2, loc2=4, fc="none", ec='0.5')

    # The data: only display for low temperature in the inset figure.
    # Tmax = max(T_D)

    plot_jet(ax2, (1.5e18,3.8e18), (1e-1,4e-1))
    plot_ejecta(ax2, (1.5e18,3.8e18), (1e-1,4e-1))

    # ax2.plot(T[T<=Tmax], Cp[T<=Tmax], 'x', c='b', mew=2, alpha=0.8, label='Experiment')
    # The Einstein fit (not very good at low T).
    # ax2.plot(T_E[T_E<=Tmax], CV_E[T_E<=Tmax], c='m', lw=2, alpha=0.5, label='Einstein model')
    # The Debye fit.
    # ax2.plot(T_D, CV_D, c='r', lw=2, alpha=0.5, label='Debye model')


    ax1.plot([1e15,1e16],[1e-4,1e-4],color='blue',lw=2,label='GRB ejecta layers')
    ax1.plot([1e15,1e16],[1e-4,1e-4],color='red',lw=2,label='kN ejecta layers')
    ax1.legend(loc="upper right")
    # Some ad hoc tweaks.
    # ax1.set_ylim(-2,2)
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.set_xlim(1e16,1e20)
    ax2.set_xscale("log")
    ax2.set_yscale("log")
    ax1.set_title("GRB and kN ejecta dynamics")
    # ax2.set_xticklabels([])
    ax1.tick_params(axis='both', which='both', labelleft=True,
                    labelright=False, tick1On=True, tick2On=True,
                    # labelsize=plotdic["fontsize"],
                    direction='in',
                    bottom=True, top=True, left=True, right=True)
    ax2.tick_params(axis='both', which='both', labelleft=True,
                    labelright=False, tick1On=True, tick2On=True,
                    # labelsize=plotdic["fontsize"],
                    direction='in',
                    bottom=True, top=True, left=True, right=True)

    # ax2.set_xticks(np.array([3e18,3e18]))
    ax2.set_xticklabels(ax2.get_xticks(), backgroundcolor='w')
    ax2.set_yticklabels(ax2.get_yticks(), backgroundcolor='w')
    # ax2.tick_params(axis='x', which='major', pad=8)
    # ax2.tick_params(axis='y', which='major', pad=8)
    plt.tight_layout()
    plt.savefig(curdir+"plot_pretty_dynamics.png",dpi=256)
    plt.show()
# plot_dynamics_pretty_plot()
def plot_dynamics_pretty_plot2():

    dfile_jet = h5py.File(curdir+"jet_dynamics_layers.h5", "r")
    dfile_ej = h5py.File(curdir+"ejecta_dynamics_shells_layers.h5", "r")

    v_n_x = "tburst"
    v_n_y = "mom"
    v_n_y2= "R"
    def plot_jet(ax, xlim, ylim):
        print(dfile_jet.keys())
        layers = ["layer=0",
                  "layer=19",
                  "layer=49",
                  "layer=69"]
        norm = Normalize(vmin=0, vmax=len( dfile_jet.keys() ))
        cmap = cm.Blues
        for il, layer in enumerate(layers):
            x_arr = np.array(dfile_jet[layer][v_n_x])
            if (v_n_y == "mom"):
                y_arr = np.array(dfile_jet[layer]["Gamma"]) * np.array(dfile_jet[layer]["beta"])
                # y_arr = y_arr/y_arr[0]
            else:
                y_arr = np.array(dfile_jet[layer][v_n_y])
            if (len(xlim)>0):
                y_arr = y_arr[(x_arr > xlim[0]) & (x_arr < xlim[-1])]
                x_arr = x_arr[(x_arr > xlim[0]) & (x_arr < xlim[-1])]
            if (len(ylim)>0):
                x_arr = x_arr[(y_arr > ylim[0]) & (y_arr < ylim[-1])]
                y_arr = y_arr[(y_arr > ylim[0]) & (y_arr < ylim[-1])]
            ax.plot(x_arr, y_arr, ls='-', color=cmap(norm(int(layer.split("=")[-1]))))

    def plot_ejecta(ax, xlim, ylim):
        ishell = 2
        layers = ["shell={} layer=0".format(ishell),
                  # "shell={} layer=2".format(ishell),
                  # "shell={} layer=4".format(ishell),
                  # "shell={} layer=6".format(ishell),
                  "shell={} layer=8".format(ishell),
                  # "shell={} layer=10".format(ishell),
                  # "shell={} layer=12".format(ishell),
                  # "shell={} layer=14".format(ishell),
                  "shell={} layer=16".format(ishell),
                  # "shell={} layer=18".format(ishell),
                  # "shell={} layer=20".format(ishell),
                  # "shell={} layer=22".format(ishell),
                  "shell={} layer=23".format(ishell)
                  ]
        ishell = 25
        layers2 = ["shell={} layer=0".format(ishell),
                   # "shell={} layer=2".format(ishell),
                   # "shell={} layer=4".format(ishell),
                   # "shell={} layer=6".format(ishell),
                   "shell={} layer=8".format(ishell),
                   # "shell={} layer=10".format(ishell),
                   # "shell={} layer=12".format(ishell),
                   # "shell={} layer=14".format(ishell),
                   "shell={} layer=16".format(ishell),
                   # "shell={} layer=18".format(ishell),
                   # "shell={} layer=20".format(ishell),
                   # "shell={} layer=22".format(ishell),
                   "shell={} layer=23".format(ishell)
                   ]
        layers+=layers2
        norm = Normalize(vmin=0, vmax=dfile_ej.attrs["nlayers"])
        cmap = cm.Reds_r
        for il, layer in enumerate(layers):

            ijl = "layer=69" # the one to vompare with
            beta_jet = np.array(dfile_jet[ijl]["beta"])
            R_jet = np.array(dfile_jet[ijl]["R"])
            beta_ej = np.array(dfile_ej[layer]["beta"])
            R_ej = np.array(dfile_ej[layer]["R"])
            beta_jet = np.array(dfile_ej[layer][v_n_x])
            beta12 = (beta_ej - beta_jet) / (1. - beta_ej * beta_jet)
            coll_idx = np.argmax( R_jet < R_ej )
            cs_jet2 = 5. * beta_jet ** 2 / 9.

            x_arr = np.array(dfile_ej[layer][v_n_x])
            if (v_n_y == "mom"):
                y_arr = np.array(dfile_ej[layer]["Gamma"]) * np.array(dfile_ej[layer]["beta"])
                # y_arr = y_arr/y_arr[0]
            else:
                y_arr = np.array(dfile_ej[layer][v_n_y])
            if (len(xlim)>0):
                y_arr = y_arr[(x_arr > xlim[0]) & (x_arr < xlim[-1])]
                x_arr = x_arr[(x_arr > xlim[0]) & (x_arr < xlim[-1])]
            if (len(ylim)>0):
                x_arr = x_arr[(y_arr > ylim[0]) & (y_arr < ylim[-1])]
                y_arr = y_arr[(y_arr > ylim[0]) & (y_arr < ylim[-1])]
            ax.plot(x_arr, y_arr, ls='-', color=cmap(norm(int(layer.split("=")[-1]))))

    fig, ax1 = plt.subplots(figsize=(4,3))
    # T_E = np.arange(1,max(T)+1,1)
    # The data.
    plot_jet(ax1, (), ())
    plot_ejecta(ax1, (), ())
    # ax1.plot(T, Cp, 'x', c='b', mew=2, alpha=0.8, label='Experiment')
    # The Einstein fit.
    # ax1.plot(T_E, CV_E, c='m', lw=2, alpha=0.5, label='Einstein model')
    # ax1.set_xlabel(r'Radius [cm]')#(r'$T\,/\mathrm{K}$')
    # ax1.set_ylabel(r'Momentum $\Gamma\beta$')#(r'$C_p\,/\mathrm{J\,K^{-1}\,mol^{-1}}$')

    ax1.set_ylabel(r"momentum $\Gamma \beta$")
    ax1.set_xlabel(r"time=$\int dr \beta c$ [s]")
    ax1.legend(loc=0)

    # Create a set of inset Axes: these should fill the bounding box allocated to
    # them.
    ax2 = plt.axes([0,0,2,2])
    # Manually set the position and relative size of the inset axes within ax1
    ip = InsetPosition(ax1, [0.18,0.1,0.8,0.4]) # [0.4,0.2,0.5,0.5]
    ax2.set_axes_locator(ip)
    # Mark the region corresponding to the inset axes on ax1 and draw lines
    # in grey linking the two axes.
    mark_inset(ax1, ax2, loc1=2, loc2=4, fc="none", ec='0.5')

    # The data: only display for low temperature in the inset figure.
    # Tmax = max(T_D)

    plot_jet(ax2, (5e7,5e8), (6e-2,8e-1))
    plot_ejecta(ax2, (5e7,5e8), (6e-2,8e-1))

    # ax2.plot(T[T<=Tmax], Cp[T<=Tmax], 'x', c='b', mew=2, alpha=0.8, label='Experiment')
    # The Einstein fit (not very good at low T).
    # ax2.plot(T_E[T_E<=Tmax], CV_E[T_E<=Tmax], c='m', lw=2, alpha=0.5, label='Einstein model')
    # The Debye fit.
    # ax2.plot(T_D, CV_D, c='r', lw=2, alpha=0.5, label='Debye model')


    ax1.plot([1e15,1e16],[1e-4,1e-4],color='blue',lw=2,label='GRB ejecta layers')
    ax1.plot([1e15,1e16],[1e-4,1e-4],color='red',lw=2,label='kN ejecta layers')
    ax1.legend(loc="upper right")
    # Some ad hoc tweaks.
    # ax1.set_ylim(-2,2)
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.set_xlim(4e7,1e11)
    ax1.set_ylim(1e-5,1e0)
    ax2.set_xscale("log")
    ax2.set_yscale("log")
    ax1.set_title("GRB and kN ejecta dynamics")
    # ax2.set_xticklabels([])
    ax1.tick_params(axis='both', which='both', labelleft=True,
                    labelright=False, tick1On=True, tick2On=True,
                    # labelsize=plotdic["fontsize"],
                    direction='in',
                    bottom=True, top=True, left=True, right=True)
    ax2.tick_params(axis='both', which='both', labelleft=True,
                    labelright=False, tick1On=True, tick2On=True,
                    # labelsize=plotdic["fontsize"],
                    direction='in',
                    bottom=True, top=True, left=True, right=True)



    # ax2.set_xticks(np.array([3e18,3e18]))
    # ax2.set_xticklabels(ax2.get_xticks(), backgroundcolor='w')
    # ax2.set_yticklabels(ax2.get_yticks(), backgroundcolor='w')
    # ax2.tick_params(axis='x', which='major', pad=8)
    # ax2.tick_params(axis='y', which='major', pad=8)
    plt.tight_layout()
    plt.savefig(curdir+"plot_pretty_dynamics.png",dpi=256)
    plt.show()

    # ax = axes[0]
    # ax.plot(dyn2t2.get("tburst"),dyn2t2.get("beta"),ls='-', color='blue',label='Jet BW')
    # ax.plot(dyn2t2.get("tburst"), dyn2t2.get("beta_ej"),ls='-', color='red',label='kN BW')
    # ax.axvline(x=dyn2t2.get("tburst")[np.argmax( dyn2t2.get("beta") < dyn2t2.get("beta_ej") )], color="gray", linestyle="-",label=r"$R$ at $\beta_{\rm ej}=\beta_{\rm jet}$")
    # ax.axvline(x=dyn2t2.get("tburst")[np.argmax( dyn2t2.get("R") < dyn2t2.get("R_ej") )], color="gray", linestyle="--",
    #            label=r"$R$ at $R_{\rm ej}=R_{\rm jet}$, "
    #                  r"$\beta_{12}=$"+"{:.2f}".format(beta12[coll_idx])
    #                  +" $\Gamma_{12}=$"+"{:.2f}".format(get_Gamma(beta12[coll_idx]))
    #                  +" $c_s=$"+"{:.2f}".format(np.sqrt(cs_jet2[coll_idx]))
    #                  +" $\mathcal{M}_s=$"+"{:.2f}".format(beta12[coll_idx]/np.sqrt(cs_jet2[coll_idx])))
    # ax.set_xlabel(r"time=$\int dr \beta c$ [s]")
    # ax.set_ylabel(r"velocity $\beta$ [c] (solid lines)")
    # ax.set_xscale("log")
    # ax.set_yscale("log")
    # # ax.set_xlim(3e17,1e19)
    # ax.set_xlim(2e6,2e9)
    # ax.set_ylim(8e-2,1.1e0)
    # ax.legend(loc="best")
    # ax.tick_params(axis='both', which='both', labelleft=True,
    #                labelright=False, tick1On=True, tick2On=True,
    #                labelsize=12,
    #                direction='in',
    #                bottom=True, top=True, left=True, right=True)
    # ax2 = ax.twinx()
    # ax2.plot(dyn2t2.get("tburst"), dyn2t2.get("R_ej"), ls='--', color='lime')
    # ax2.plot(dyn2t2.get("tburst"), np.abs(dyn2t2.get("R") - dyn2t2.get("R_ej")), ls='--', color='lime', label=r"$| \Delta R |$")
    # ax2.plot(dyn2t2.get("tburst"), dyn2t2.get("R_ej"), ls='--', color='red')
    # ax2.legend(loc="best")
    # ax2.set_ylabel("Radius [cm] (dashed lines)")
    # ax2.set_xscale("log")
    # ax2.set_yscale("log")
    #
    # ctheta_j = dyn2t2.kwargs["ctheta_j"] + 0.5 * (2.0 * dyn2t2.get("theta") - 2.0 * dyn2t2.kwargs["theta_w_j"])
    # ctheta_ej = dyn2t2.kwargs["ctheta_ej"] + 0.5 * (2.0 * dyn2t2.get("theta_ej") - 2.0 * dyn2t2.kwargs["theta_w_ej"])
    #
    # plt.tight_layout()
    # plt.savefig("time_intersection2.png",dpi=256)
    # plt.show()
# plot_dynamics_pretty_plot2()


def plot_dynamics_pretty_plot2Accel():

    dfile_jet = h5py.File(curdir+"jet_dynamics_layers.h5", "r")
    dfile_ej = h5py.File(curdir+"ejecta_dynamics_shells_layers.h5", "r")

    v_n_x = "tburst"
    v_n_y = "mom"
    v_n_y2= "R"
    def plot_jet(ax, xlim, ylim):
        print(dfile_jet.keys())
        layers = ["layer=0",
                  "layer=19",
                  "layer=49",
                  "layer=69"]
        norm = Normalize(vmin=0, vmax=len( dfile_jet.keys() ))
        cmap = cm.Blues
        for il, layer in enumerate(layers):
            x_arr = np.array(dfile_jet[layer][v_n_x])
            if (v_n_y == "mom"):
                y_arr = np.array(dfile_jet[layer]["Gamma"]) * np.array(dfile_jet[layer]["beta"])
                y_arr = y_arr/y_arr[0]
            else:
                y_arr = np.array(dfile_jet[layer][v_n_y])
            if (len(xlim)>0):
                y_arr = y_arr[(x_arr > xlim[0]) & (x_arr < xlim[-1])]
                x_arr = x_arr[(x_arr > xlim[0]) & (x_arr < xlim[-1])]
            if (len(ylim)>0):
                x_arr = x_arr[(y_arr > ylim[0]) & (y_arr < ylim[-1])]
                y_arr = y_arr[(y_arr > ylim[0]) & (y_arr < ylim[-1])]
            ax.plot(x_arr, y_arr, ls='-', color=cmap(norm(int(layer.split("=")[-1]))))

    def plot_ejecta(ax1, ax2, xlim, ylim):
        ishell = 2
        layers = ["shell={} layer=0".format(ishell),
                  # "shell={} layer=2".format(ishell),
                  # "shell={} layer=4".format(ishell),
                  # "shell={} layer=6".format(ishell),
                  "shell={} layer=8".format(ishell),
                  # "shell={} layer=10".format(ishell),
                  # "shell={} layer=12".format(ishell),
                  # "shell={} layer=14".format(ishell),
                  "shell={} layer=16".format(ishell),
                  # "shell={} layer=18".format(ishell),
                  # "shell={} layer=20".format(ishell),
                  # "shell={} layer=22".format(ishell),
                  "shell={} layer=23".format(ishell)
                  ]
        norm = Normalize(vmin=0, vmax=dfile_ej.attrs["nlayers"])

        for il, layer in enumerate(layers):

            ijl = "layer=69" # the one to vompare with
            beta_jet = np.array(dfile_jet[ijl]["beta"])
            R_jet = np.array(dfile_jet[ijl]["R"])
            beta_ej = np.array(dfile_ej[layer]["beta"])
            R_ej = np.array(dfile_ej[layer]["R"])
            beta_jet = np.array(dfile_ej[layer][v_n_x])
            beta12 = (beta_ej - beta_jet) / (1. - beta_ej * beta_jet)
            coll_idx = np.argmax( R_jet < R_ej )
            cs_jet2 = 5. * beta_jet ** 2 / 9.

            cmap = cm.Greys_r
            x_arr = np.array(dfile_ej[layer][v_n_x])
            y_arr1 = Beta( np.array(dfile_ej[layer]["GammaRho"]) ) * np.array(dfile_ej[layer]["GammaRho"])
            # y_arr = y_arr / y_arr[0]
            if (len(xlim)>0):
                y_arr1 = y_arr1[(x_arr > xlim[0]) & (x_arr < xlim[-1])]
                x_arr = x_arr[(x_arr > xlim[0]) & (x_arr < xlim[-1])]
            if (len(ylim)>0):
                x_arr = x_arr[(y_arr1 > ylim[0]) & (y_arr1 < ylim[-1])]
                y_arr1 = y_arr1[(y_arr1 > ylim[0]) & (y_arr1 < ylim[-1])]
            ax1.plot(x_arr, y_arr1, ls='--', color=cmap(norm(int(layer.split("=")[-1]))))

            tmp = np.array(dfile_ej[layer]["GammaRho"])
            y_arr1 = Beta( np.array(dfile_ej[layer]["GammaRho"]) ) * np.array(dfile_ej[layer]["GammaRho"])
            y_arr2 = Beta( np.array(dfile_ej[layer]["Gamma"]) ) * np.array(dfile_ej[layer]["Gamma"])
            print("Mom CBM[{}, {}]".format(y_arr1.min(),y_arr1.max()))
            print("MOm [{}, {}]".format(y_arr2.min(),y_arr2.max()))

            mask1 = y_arr1 > y_arr2
            if (len(y_arr2[mask1]) > 1):
                cmap = cm.Reds_r
                y_arr = y_arr2[mask1]
                x_arr = np.array(dfile_ej[layer][v_n_x])[mask1]
                y_arr = y_arr / y_arr[0]
                if (len(xlim)>0):
                    y_arr = y_arr[(x_arr > xlim[0]) & (x_arr < xlim[-1])]
                    x_arr = x_arr[(x_arr > xlim[0]) & (x_arr < xlim[-1])]
                if (len(ylim)>0):
                    x_arr = x_arr[(y_arr > ylim[0]) & (y_arr < ylim[-1])]
                    y_arr = y_arr[(y_arr > ylim[0]) & (y_arr < ylim[-1])]
                ax2.plot(x_arr, y_arr, ls='--', lw=2, color=cmap(norm(int(layer.split("=")[-1]))))

            mask2 = y_arr1 < y_arr2
            if (len(y_arr2[mask2])>1):
                cmap = cm.Greens_r
                y_arr = y_arr2[mask2]
                x_arr = np.array(dfile_ej[layer][v_n_x])[mask2]
                y_arr = y_arr / y_arr[0]
                if (len(xlim)>0):
                    y_arr = y_arr[(x_arr > xlim[0]) & (x_arr < xlim[-1])]
                    x_arr = x_arr[(x_arr > xlim[0]) & (x_arr < xlim[-1])]
                if (len(ylim)>0):
                    x_arr = x_arr[(y_arr > ylim[0]) & (y_arr < ylim[-1])]
                    y_arr = y_arr[(y_arr > ylim[0]) & (y_arr < ylim[-1])]
                ax2.plot(x_arr, y_arr, ls='--', color=cmap(norm(int(layer.split("=")[-1]))))

            # ax.plot(x_arr, y_arr, ls='-', color=cmap(norm(int(layer.split("=")[-1]))))
        # # -----------------------------------------------------
        # ishell = 25
        # layers2 = ["shell={} layer=0".format(ishell),
        #            # "shell={} layer=2".format(ishell),
        #            # "shell={} layer=4".format(ishell),
        #            # "shell={} layer=6".format(ishell),
        #            "shell={} layer=8".format(ishell),
        #            # "shell={} layer=10".format(ishell),
        #            # "shell={} layer=12".format(ishell),
        #            # "shell={} layer=14".format(ishell),
        #            "shell={} layer=16".format(ishell),
        #            # "shell={} layer=18".format(ishell),
        #            # "shell={} layer=20".format(ishell),
        #            # "shell={} layer=22".format(ishell),
        #            "shell={} layer=23".format(ishell)
        #            ]
        # layers=layers2
        # norm = Normalize(vmin=0, vmax=dfile_ej.attrs["nlayers"])
        # cmap = cm.Greens_r
        # for il, layer in enumerate(layers):
        #
        #     ijl = "layer=69" # the one to vompare with
        #     beta_jet = np.array(dfile_jet[ijl]["beta"])
        #     R_jet = np.array(dfile_jet[ijl]["R"])
        #     beta_ej = np.array(dfile_ej[layer]["beta"])
        #     R_ej = np.array(dfile_ej[layer]["R"])
        #     beta_jet = np.array(dfile_ej[layer][v_n_x])
        #     beta12 = (beta_ej - beta_jet) / (1. - beta_ej * beta_jet)
        #     coll_idx = np.argmax( R_jet < R_ej )
        #     cs_jet2 = 5. * beta_jet ** 2 / 9.
        #
        #     x_arr = np.array(dfile_ej[layer][v_n_x])
        #     y_arr = Beta( np.array(dfile_ej[layer]["GammaCBM"]) ) * np.array(dfile_ej[layer]["GammaCBM"])
        #     # y_arr = y_arr / y_arr[0]
        #     if (len(xlim)>0):
        #         y_arr = y_arr[(x_arr > xlim[0]) & (x_arr < xlim[-1])]
        #         x_arr = x_arr[(x_arr > xlim[0]) & (x_arr < xlim[-1])]
        #     if (len(ylim)>0):
        #         x_arr = x_arr[(y_arr > ylim[0]) & (y_arr < ylim[-1])]
        #         y_arr = y_arr[(y_arr > ylim[0]) & (y_arr < ylim[-1])]
        #     ax1.plot(x_arr, y_arr, ls='-', color=cmap(norm(int(layer.split("=")[-1]))))
        #
        #     x_arr = np.array(dfile_ej[layer][v_n_x])
        #     y_arr = Beta( np.array(dfile_ej[layer]["Gamma"]) ) * np.array(dfile_ej[layer]["Gamma"])
        #     y_arr = y_arr / y_arr[0]
        #     if (len(xlim)>0):
        #         y_arr = y_arr[(x_arr > xlim[0]) & (x_arr < xlim[-1])]
        #         x_arr = x_arr[(x_arr > xlim[0]) & (x_arr < xlim[-1])]
        #     if (len(ylim)>0):
        #         x_arr = x_arr[(y_arr > ylim[0]) & (y_arr < ylim[-1])]
        #         y_arr = y_arr[(y_arr > ylim[0]) & (y_arr < ylim[-1])]
        #     ax2.plot(x_arr, y_arr, ls=':', color=cmap(norm(int(layer.split("=")[-1]))))

    def plot_ejecta_dlnrhodr(ax1, ax2, xlim, ylim):
        ishell = 2
        layers = ["shell={} layer=0".format(ishell),
                  # "shell={} layer=2".format(ishell),
                  # "shell={} layer=4".format(ishell),
                  # "shell={} layer=6".format(ishell),
                  "shell={} layer=8".format(ishell),
                  # "shell={} layer=10".format(ishell),
                  # "shell={} layer=12".format(ishell),
                  # "shell={} layer=14".format(ishell),
                  "shell={} layer=16".format(ishell),
                  # "shell={} layer=18".format(ishell),
                  # "shell={} layer=20".format(ishell),
                  # "shell={} layer=22".format(ishell),
                  "shell={} layer=23".format(ishell)
                  ]
        norm = Normalize(vmin=0, vmax=dfile_ej.attrs["nlayers"])

        for il, layer in enumerate(layers):

            ijl = "layer=69" # the one to vompare with
            beta_jet = np.array(dfile_jet[ijl]["beta"])
            R_jet = np.array(dfile_jet[ijl]["R"])
            beta_ej = np.array(dfile_ej[layer]["beta"])
            R_ej = np.array(dfile_ej[layer]["R"])
            beta_jet = np.array(dfile_ej[layer][v_n_x])
            beta12 = (beta_ej - beta_jet) / (1. - beta_ej * beta_jet)
            coll_idx = np.argmax( R_jet < R_ej )
            cs_jet2 = 5. * beta_jet ** 2 / 9.

            cmap = cm.Greys_r
            x_arr = np.array(dfile_ej[layer][v_n_x])
            y_arr1 = Beta( np.array(dfile_ej[layer]["GammaCBM"]) ) * np.array(dfile_ej[layer]["GammaCBM"])
            # y_arr = y_arr / y_arr[0]
            if (len(xlim)>0):
                y_arr1 = y_arr1[(x_arr > xlim[0]) & (x_arr < xlim[-1])]
                x_arr = x_arr[(x_arr > xlim[0]) & (x_arr < xlim[-1])]
            if (len(ylim)>0):
                x_arr = x_arr[(y_arr1 > ylim[0]) & (y_arr1 < ylim[-1])]
                y_arr1 = y_arr1[(y_arr1 > ylim[0]) & (y_arr1 < ylim[-1])]
            ax1.plot(x_arr, y_arr1, ls='--', color=cmap(norm(int(layer.split("=")[-1]))))


            # y_arr1 = Beta( np.array(dfile_ej[layer]["GammaCBM"]) ) * np.array(dfile_ej[layer]["GammaCBM"])
            y_arr =np.array(dfile_ej[layer]["dlnrhodr"])
            print("Mom CBM[{}, {}]".format(y_arr1.min(),y_arr1.max()))
            # print("MOm [{}, {}]".format(y_arr2.min(),y_arr2.max()))

            cmap = cm.Reds_r
            x_arr = np.array(dfile_ej[layer][v_n_x])
            y_arr = y_arr / y_arr[0]
            if (len(xlim)>0):
                y_arr = y_arr[(x_arr > xlim[0]) & (x_arr < xlim[-1])]
                x_arr = x_arr[(x_arr > xlim[0]) & (x_arr < xlim[-1])]
            if (len(ylim)>0):
                x_arr = x_arr[(y_arr > ylim[0]) & (y_arr < ylim[-1])]
                y_arr = y_arr[(y_arr > ylim[0]) & (y_arr < ylim[-1])]
            ax2.plot(x_arr, y_arr, ls='--', lw=2, color=cmap(norm(int(layer.split("=")[-1]))))

    fig, ax1 = plt.subplots(figsize=(4,3))
    ax1_ = ax1.twinx()
    # T_E = np.arange(1,max(T)+1,1)
    # The data.
    # plot_jet(ax1, (), ())
    plot_ejecta(ax1_, ax1, (), ())
    # ax1.plot(T, Cp, 'x', c='b', mew=2, alpha=0.8, label='Experiment')
    # The Einstein fit.
    # ax1.plot(T_E, CV_E, c='m', lw=2, alpha=0.5, label='Einstein model')
    # ax1.set_xlabel(r'Radius [cm]')#(r'$T\,/\mathrm{K}$')
    # ax1.set_ylabel(r'Momentum $\Gamma\beta$')#(r'$C_p\,/\mathrm{J\,K^{-1}\,mol^{-1}}$')

    ax1.set_ylabel(r"momentum $\Gamma \beta$")
    ax1.set_xlabel(r"time=$\int dr \beta c$ [s]")
    ax1.legend(loc=0)

    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.set_xlim(1e7, 1e10)#(7e6,5e10)
    ax1.set_ylim(0.994, 1.002)#(1e-1,1.1e0)
    # ax1.set_title("GRB and kN ejecta dynamics")
    ax1.tick_params(axis='both', which='both', labelleft=True,
                    labelright=False, tick1On=True, tick2On=False,
                    # labelsize=plotdic["fontsize"],
                    direction='in',
                    bottom=True, top=True, left=True, right=False)

    ax1_.set_ylabel(r"momentum $\Gamma \beta$")
    ax1_.set_xlabel(r"time=$\int dr \beta c$ [s]")
    ax1_.legend(loc=0)

    ax1_.set_yscale("log")
    ax1_.set_ylim(8e-3, 0.994)#(1e-1,1.1e0)
    # ax1.set_title("GRB and kN ejecta dynamics")
    ax1_.tick_params(axis='both', which='both', labelleft=False,
                     labelright=True, tick1On=False, tick2On=True,
                     # labelsize=plotdic["fontsize"],
                     direction='in',
                     bottom=True, top=True, left=False, right=True)
    # ax1.plot([1e15,1e16],[1e-4,1e-4],color='blue',lw=2,label='GRB ejecta layers')
    # ax1.plot([1e15,1e16],[1e-4,1e-4],color='red',lw=2,label='kN ejecta layers')

    # Create a set of inset Axes: these should fill the bounding box allocated to
    # them.
    '''
    ax2 = plt.axes([0,0,2,2])
    # Manually set the position and relative size of the inset axes within ax1
    ip = InsetPosition(ax1, [0.20,0.1,0.75,0.4]) # [0.4,0.2,0.5,0.5]
    ax2.set_axes_locator(ip)
    # Mark the region corresponding to the inset axes on ax1 and draw lines
    # in grey linking the two axes.
    mark_inset(ax1, ax2, loc1=2, loc2=4, fc="none", ec='0.5')

    # The data: only display for low temperature in the inset figure.
    # Tmax = max(T_D)

    # plot_jet(ax2, (5e7,5e8), (6e-2,8e-1))
    # plot_ejecta(ax2, (1e7,1e10), (0.994,1.002))

    # ax2.plot(T[T<=Tmax], Cp[T<=Tmax], 'x', c='b', mew=2, alpha=0.8, label='Experiment')
    # The Einstein fit (not very good at low T).
    # ax2.plot(T_E[T_E<=Tmax], CV_E[T_E<=Tmax], c='m', lw=2, alpha=0.5, label='Einstein model')
    # The Debye fit.
    # ax2.plot(T_D, CV_D, c='r', lw=2, alpha=0.5, label='Debye model')

    ax2.set_xscale("log")
    ax2.set_yscale("log")
    # ax2.set_xticklabels([])

    ax2.tick_params(axis='both', which='both', labelleft=True,
                    labelright=False, tick1On=True, tick2On=True,
                    # labelsize=plotdic["fontsize"],
                    direction='in',
                    bottom=True, top=True, left=True, right=True)
    '''


    # ax2.set_xticks(np.array([3e18,3e18]))
    # ax2.set_xticklabels(ax2.get_xticks(), backgroundcolor='w')
    # ax2.set_yticklabels(ax2.get_yticks(), backgroundcolor='w')
    # ax2.tick_params(axis='x', which='major', pad=8)
    # ax2.tick_params(axis='y', which='major', pad=8)
    plt.tight_layout()
    plt.savefig(curdir+"plot_pretty_dynamics.png",dpi=256)
    plt.show()

    # ax = axes[0]
    # ax.plot(dyn2t2.get("tburst"),dyn2t2.get("beta"),ls='-', color='blue',label='Jet BW')
    # ax.plot(dyn2t2.get("tburst"), dyn2t2.get("beta_ej"),ls='-', color='red',label='kN BW')
    # ax.axvline(x=dyn2t2.get("tburst")[np.argmax( dyn2t2.get("beta") < dyn2t2.get("beta_ej") )], color="gray", linestyle="-",label=r"$R$ at $\beta_{\rm ej}=\beta_{\rm jet}$")
    # ax.axvline(x=dyn2t2.get("tburst")[np.argmax( dyn2t2.get("R") < dyn2t2.get("R_ej") )], color="gray", linestyle="--",
    #            label=r"$R$ at $R_{\rm ej}=R_{\rm jet}$, "
    #                  r"$\beta_{12}=$"+"{:.2f}".format(beta12[coll_idx])
    #                  +" $\Gamma_{12}=$"+"{:.2f}".format(get_Gamma(beta12[coll_idx]))
    #                  +" $c_s=$"+"{:.2f}".format(np.sqrt(cs_jet2[coll_idx]))
    #                  +" $\mathcal{M}_s=$"+"{:.2f}".format(beta12[coll_idx]/np.sqrt(cs_jet2[coll_idx])))
    # ax.set_xlabel(r"time=$\int dr \beta c$ [s]")
    # ax.set_ylabel(r"velocity $\beta$ [c] (solid lines)")
    # ax.set_xscale("log")
    # ax.set_yscale("log")
    # # ax.set_xlim(3e17,1e19)
    # ax.set_xlim(2e6,2e9)
    # ax.set_ylim(8e-2,1.1e0)
    # ax.legend(loc="best")
    # ax.tick_params(axis='both', which='both', labelleft=True,
    #                labelright=False, tick1On=True, tick2On=True,
    #                labelsize=12,
    #                direction='in',
    #                bottom=True, top=True, left=True, right=True)
    # ax2 = ax.twinx()
    # ax2.plot(dyn2t2.get("tburst"), dyn2t2.get("R_ej"), ls='--', color='lime')
    # ax2.plot(dyn2t2.get("tburst"), np.abs(dyn2t2.get("R") - dyn2t2.get("R_ej")), ls='--', color='lime', label=r"$| \Delta R |$")
    # ax2.plot(dyn2t2.get("tburst"), dyn2t2.get("R_ej"), ls='--', color='red')
    # ax2.legend(loc="best")
    # ax2.set_ylabel("Radius [cm] (dashed lines)")
    # ax2.set_xscale("log")
    # ax2.set_yscale("log")
    #
    # ctheta_j = dyn2t2.kwargs["ctheta_j"] + 0.5 * (2.0 * dyn2t2.get("theta") - 2.0 * dyn2t2.kwargs["theta_w_j"])
    # ctheta_ej = dyn2t2.kwargs["ctheta_ej"] + 0.5 * (2.0 * dyn2t2.get("theta_ej") - 2.0 * dyn2t2.kwargs["theta_w_ej"])
    #
    # plt.tight_layout()
    # plt.savefig("time_intersection2.png",dpi=256)
    # plt.show()
# plot_dynamics_pretty_plot2Accel()

def plot_dynamics_pretty_plot2_rho(jet_layers=(0,19,59,69),
                                   ej_shells=(1,), ej_layers=(0,8,16,23), v_n_x = "tburst", v_n_y = "nism"):

    dfile_jet = h5py.File(curdir+"jet_dynamics_layers.h5", "r")
    dfile_ej = h5py.File(curdir+"ejecta_dynamics_shells_layers.h5", "r")


    v_n_y2= "mom"
    def plot_jet(ax, xlim, ylim):
        print(dfile_jet.keys())
        layers = []
        for i in jet_layers:
            layers.append("layer={}".format(i))

        norm = Normalize(vmin=0, vmax=len( dfile_jet.keys() ))
        cmap = cm.Blues
        for il, layer in enumerate(layers):
            x_arr = np.array(dfile_jet[layer][v_n_x])
            if (v_n_y == "nism"):
                y_arr = np.array(dfile_jet[layer]["rho"]) / cgs.mp
            else:
                y_arr = np.array(dfile_jet[layer][v_n_y])
            if (len(xlim)>0):
                y_arr = y_arr[(x_arr > xlim[0]) & (x_arr < xlim[-1])]
                x_arr = x_arr[(x_arr > xlim[0]) & (x_arr < xlim[-1])]
            if (len(ylim)>0):
                x_arr = x_arr[(y_arr > ylim[0]) & (y_arr < ylim[-1])]
                y_arr = y_arr[(y_arr > ylim[0]) & (y_arr < ylim[-1])]
            ax.plot(x_arr, y_arr, ls='-', color=cmap(norm(int(layer.split("=")[-1]))))

    def plot_ejecta(ax, xlim, ylim):
        layers = []
        for i in ej_shells:
            for j in ej_layers:
                layers.append("shell={} layer={}".format(i,j))
        # layers+=layers2
        norm = Normalize(vmin=0, vmax=dfile_ej.attrs["nlayers"])
        cmap = cm.Reds_r
        for il, layer in enumerate(layers):

            ijl = "layer=69" # the one to vompare with
            beta_jet = np.array(dfile_jet[ijl]["beta"])
            R_jet = np.array(dfile_jet[ijl]["R"])
            beta_ej = np.array(dfile_ej[layer]["beta"])
            R_ej = np.array(dfile_ej[layer]["R"])
            beta_jet = np.array(dfile_ej[layer][v_n_x])
            beta12 = (beta_ej - beta_jet) / (1. - beta_ej * beta_jet)
            coll_idx = np.argmax( R_jet < R_ej )
            cs_jet2 = 5. * beta_jet ** 2 / 9.

            x_arr = np.array(dfile_ej[layer][v_n_x])
            if (v_n_y == "nism"):
                y_arr = np.array(dfile_ej[layer]["rho"]) / cgs.mp
            else:
                y_arr = np.array(dfile_ej[layer][v_n_y])
            if (len(xlim)>0):
                y_arr = y_arr[(x_arr > xlim[0]) & (x_arr < xlim[-1])]
                x_arr = x_arr[(x_arr > xlim[0]) & (x_arr < xlim[-1])]
            if (len(ylim)>0):
                x_arr = x_arr[(y_arr > ylim[0]) & (y_arr < ylim[-1])]
                y_arr = y_arr[(y_arr > ylim[0]) & (y_arr < ylim[-1])]
            ax.plot(x_arr, y_arr, ls='-', color=cmap(norm(int(layer.split("=")[-1]))))

    fig, ax1 = plt.subplots(figsize=(4,3))
    # T_E = np.arange(1,max(T)+1,1)
    # The data.
    plot_jet(ax1, (), ())
    plot_ejecta(ax1, (), ())
    # ax1.plot(T, Cp, 'x', c='b', mew=2, alpha=0.8, label='Experiment')
    # The Einstein fit.
    # ax1.plot(T_E, CV_E, c='m', lw=2, alpha=0.5, label='Einstein model')
    # ax1.set_xlabel(r'Radius [cm]')#(r'$T\,/\mathrm{K}$')
    # ax1.set_ylabel(r'Momentum $\Gamma\beta$')#(r'$C_p\,/\mathrm{J\,K^{-1}\,mol^{-1}}$')

    ax1.set_ylabel(r"number density $n$ [cm$^{-3}$]")
    ax1.set_xlabel(r"time=$\int dr \beta c$ [s]")
    ax1.legend(loc=0)

    # Create a set of inset Axes: these should fill the bounding box allocated to
    # them.
    ax2 = plt.axes([0,0,2,2])
    # Manually set the position and relative size of the inset axes within ax1
    ip = InsetPosition(ax1, [0.18,0.1,0.8,0.4]) # [0.4,0.2,0.5,0.5]
    ax2.set_axes_locator(ip)
    # Mark the region corresponding to the inset axes on ax1 and draw lines
    # in grey linking the two axes.
    mark_inset(ax1, ax2, loc1=2, loc2=4, fc="none", ec='0.5')

    # The data: only display for low temperature in the inset figure.
    # Tmax = max(T_D)

    plot_jet(ax2, (5e9,6e10), (4e-5,5e-4))
    plot_ejecta(ax2, (5e9,6e10), (4e-5,5e-4))

    # ax2.plot(T[T<=Tmax], Cp[T<=Tmax], 'x', c='b', mew=2, alpha=0.8, label='Experiment')
    # The Einstein fit (not very good at low T).
    # ax2.plot(T_E[T_E<=Tmax], CV_E[T_E<=Tmax], c='m', lw=2, alpha=0.5, label='Einstein model')
    # The Debye fit.
    # ax2.plot(T_D, CV_D, c='r', lw=2, alpha=0.5, label='Debye model')

    ax1.plot([1e15,1e16],[1e-4,1e-4],color='blue',lw=2,label='GRB ejecta layers')
    ax1.plot([1e15,1e16],[1e-4,1e-4],color='red',lw=2,label='kN ejecta layers')
    ax1.legend(loc="best")
    # Some ad hoc tweaks.
    # ax1.set_ylim(-2,2)
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.set_xlim(4e7,1e11)
    ax1.set_ylim(2e-11,4e-4)
    ax2.set_xscale("log")
    ax2.set_yscale("log")
    ax1.set_title("GRB and kN ejecta dynamics")
    # ax2.set_xticklabels([])
    ax1.tick_params(axis='both', which='both', labelleft=True,
                    labelright=False, tick1On=True, tick2On=True,
                    # labelsize=plotdic["fontsize"],
                    direction='in',
                    bottom=True, top=True, left=True, right=True)
    ax2.tick_params(axis='both', which='both', labelleft=True,
                    labelright=False, tick1On=True, tick2On=True,
                    # labelsize=plotdic["fontsize"],
                    direction='in',
                    bottom=True, top=True, left=True, right=True)



    # ax2.set_xticks(np.array([3e18,3e18]))
    # ax2.set_xticklabels(ax2.get_xticks(), backgroundcolor='w')
    # ax2.set_yticklabels(ax2.get_yticks(), backgroundcolor='w')
    # ax2.tick_params(axis='x', which='major', pad=8)
    # ax2.tick_params(axis='y', which='major', pad=8)
    plt.tight_layout()
    plt.savefig(curdir+"plot_pretty_dynamics_rho.png",dpi=256)
    plt.show()

    # ax = axes[0]
    # ax.plot(dyn2t2.get("tburst"),dyn2t2.get("beta"),ls='-', color='blue',label='Jet BW')
    # ax.plot(dyn2t2.get("tburst"), dyn2t2.get("beta_ej"),ls='-', color='red',label='kN BW')
    # ax.axvline(x=dyn2t2.get("tburst")[np.argmax( dyn2t2.get("beta") < dyn2t2.get("beta_ej") )], color="gray", linestyle="-",label=r"$R$ at $\beta_{\rm ej}=\beta_{\rm jet}$")
    # ax.axvline(x=dyn2t2.get("tburst")[np.argmax( dyn2t2.get("R") < dyn2t2.get("R_ej") )], color="gray", linestyle="--",
    #            label=r"$R$ at $R_{\rm ej}=R_{\rm jet}$, "
    #                  r"$\beta_{12}=$"+"{:.2f}".format(beta12[coll_idx])
    #                  +" $\Gamma_{12}=$"+"{:.2f}".format(get_Gamma(beta12[coll_idx]))
    #                  +" $c_s=$"+"{:.2f}".format(np.sqrt(cs_jet2[coll_idx]))
    #                  +" $\mathcal{M}_s=$"+"{:.2f}".format(beta12[coll_idx]/np.sqrt(cs_jet2[coll_idx])))
    # ax.set_xlabel(r"time=$\int dr \beta c$ [s]")
    # ax.set_ylabel(r"velocity $\beta$ [c] (solid lines)")
    # ax.set_xscale("log")
    # ax.set_yscale("log")
    # # ax.set_xlim(3e17,1e19)
    # ax.set_xlim(2e6,2e9)
    # ax.set_ylim(8e-2,1.1e0)
    # ax.legend(loc="best")
    # ax.tick_params(axis='both', which='both', labelleft=True,
    #                labelright=False, tick1On=True, tick2On=True,
    #                labelsize=12,
    #                direction='in',
    #                bottom=True, top=True, left=True, right=True)
    # ax2 = ax.twinx()
    # ax2.plot(dyn2t2.get("tburst"), dyn2t2.get("R_ej"), ls='--', color='lime')
    # ax2.plot(dyn2t2.get("tburst"), np.abs(dyn2t2.get("R") - dyn2t2.get("R_ej")), ls='--', color='lime', label=r"$| \Delta R |$")
    # ax2.plot(dyn2t2.get("tburst"), dyn2t2.get("R_ej"), ls='--', color='red')
    # ax2.legend(loc="best")
    # ax2.set_ylabel("Radius [cm] (dashed lines)")
    # ax2.set_xscale("log")
    # ax2.set_yscale("log")
    #
    # ctheta_j = dyn2t2.kwargs["ctheta_j"] + 0.5 * (2.0 * dyn2t2.get("theta") - 2.0 * dyn2t2.kwargs["theta_w_j"])
    # ctheta_ej = dyn2t2.kwargs["ctheta_ej"] + 0.5 * (2.0 * dyn2t2.get("theta_ej") - 2.0 * dyn2t2.kwargs["theta_w_ej"])
    #
    # plt.tight_layout()
    # plt.savefig("time_intersection2.png",dpi=256)
    # plt.show()
# plot_dynamics_pretty_plot2_rho()
from matplotlib.transforms import  Transform

class LogPolarTransform(PolarAxes.PolarTransform):
    input_dims = 2
    output_dims = 2
    is_separable = False

    def __init__(self, axis=None, use_rmin=True):
        Transform.__init__(self)
        self._axis = axis
        self._use_rmin = use_rmin
    def transform_non_affine(self, tr):
        xy = np.empty(tr.shape, np.float_)
        if self._axis is not None:
            if self._use_rmin:
                rmin = self._axis.viewLim.ymin
            else:
                rmin = 0
            theta_offset = self._axis.get_theta_offset()
            theta_direction = self._axis.get_theta_direction()
        else:
            rmin = 0
            theta_offset = 0
            theta_direction = 1

        t = tr[:, 0:1]
        r = tr[:, 1:2]
        x = xy[:, 0:1]
        y = xy[:, 1:2]

        t *= theta_direction
        t += theta_offset

        r = r - rmin
        mask = r < 0
        x[:] = np.where(mask, np.nan, np.log(r) * np.cos(t))
        y[:] = np.where(mask, np.nan, np.log(r) * np.sin(t))

        return xy


    def inverted(self):
        return InvertedLogPolarTransform(self._axis, self._use_rmin)
    inverted.__doc__ = Transform.inverted.__doc__

class InvertedLogPolarTransform(Transform):
    """
    The inverse of the polar transform, mapping Cartesian
    coordinate space *x* and *y* back to *theta* and *r*.
    """
    input_dims = 2
    output_dims = 2
    is_separable = False

    def __init__(self, axis=None, use_rmin=True):
        Transform.__init__(self)
        self._axis = axis
        self._use_rmin = use_rmin

    def transform_non_affine(self, xy):
        if self._axis is not None:
            if self._use_rmin:
                rmin = self._axis.viewLim.ymin
            else:
                rmin = 0
            theta_offset = self._axis.get_theta_offset()
            theta_direction = self._axis.get_theta_direction()
        else:
            rmin = 0
            theta_offset = 0
            theta_direction = 1

        x = xy[:, 0:1]
        y = xy[:, 1:]
        r = np.exp(np.sqrt(x*x + y*y))
        with np.errstate(invalid='ignore'):
            # At x=y=r=0 this will raise an
            # invalid value warning when doing 0/0
            # Divide by zero warnings are only raised when
            # the numerator is different from 0. That
            # should not happen here.
            theta = np.arccos(x / r)
        theta = np.where(y < 0, 2 * np.pi - theta, theta)

        theta -= theta_offset
        theta *= theta_direction
        theta %= 2 * np.pi

        r += rmin

        return np.concatenate((theta, r), 1)

    def inverted(self):
        return PolarAxes.LogPolarTransform(self._axis, self._use_rmin)


def setup_arc_radial_axes(fig, rect, angle_ticks, radius_ticks, min_rad, max_rad):

    tr = PolarAxes.PolarTransform()
    # tr = LogPolarTransform()
    # tr = InvertedLogPolarTransform()

    pi = np.pi

    grid_locator1 = FixedLocator([v for v, s in angle_ticks])
    tmp = dict(angle_ticks)
    # tmp["rotation"] = 40
    tick_formatter1 = DictFormatter(tmp)


    grid_locator2 = FixedLocator([a for a, b in radius_ticks])
    tick_formatter2 = DictFormatter(dict(radius_ticks))

    grid_helper = floating_axes.GridHelperCurveLinear(tr,
                                                      # extremes=((370.0*(pi/180.0)), (170.0*(pi/180.0)), max_rad, min_rad),
                                                      extremes=((90.0*(pi/180.0)), (00.0*(pi/180.0)), max_rad, min_rad),
                                                      # extremes=((170.0*(pi/180.0)), (370.0*(pi/180.0)), max_rad, min_rad),
                                                      grid_locator1=grid_locator1,
                                                      grid_locator2=grid_locator2,
                                                      tick_formatter1=tick_formatter1,
                                                      tick_formatter2=tick_formatter2,
                                                      )

    ax1 = floating_axes.FloatingSubplot(fig, rect, grid_helper=grid_helper)
    fig.add_subplot(ax1)


    ax1.grid(True)

    # create a parasite axes whose transData in RA, cz
    aux_ax = ax1.get_aux_axes(tr)

    aux_ax.patch = ax1.patch
    ax1.patch.zorder=0.9



    # ax1.axis["left"].set_ticklabel_direction("+")
    ax1.axis["bottom"].set_axis_direction("top")
    ax1.axis["bottom"].set_ticklabel_direction("+")
    ax1.axis["bottom"].label.set_rotation(0)
    ax1.axis["bottom"].label.set_pad(-20)
    # ax1.axis["bottom"].toggle(ticklabels=True, label=True)
    # ax1.axis["bottom"].major_ticklabels.set_axis_direction("left")
    # ax1.axis["bottom"].invert_yaxis()
    # plt.gca().invert_xaxis()
    # plt.gca().invert_yaxis()




    # ax = aux_ax
    # labels = []
    # for i in range(len(angles)) :
    #     angles[i] = 0
    # for label, angle in zip(ax.get_xticklabels(), angles):
    #     x,y = label.get_position()
    #     lab = ax.text(x,y, label.get_text(), transform=label.get_transform(),
    #                   ha=label.get_ha(), va=label.get_va())
    #     lab.set_rotation(angle)
    #     labels.append(lab)
    # ax.set_xticklabels([])
    # ax1.axis["bottom"].set_axis_direction("bottom")
    # ax1.axis["bottom"].set_tick_label_direction("bottom")
    # ax1.axis["bottom"].set_axis_direction("left")
    # ax1.axis["bottom"].set_rotation(0)
    # ax1.axis["bottom"].set_axis_direction("top")

    return ax1, aux_ax
def example():
    # write angle values to the plotting array
    angles = []
    for mic_num in range(38):
        angle = float(mic_num)*(180.0/36.0)*(math.pi/180.0)+math.pi
        angles.append(angle)
    # angles = angles[::-1]
    # ------------------------------------ #
    ### these are merely the ticks that appear on the plot axis
    ### these don't actually get plotted

    angle_ticks = range(0,100,10)
    angle_ticks_rads = [a * math.pi / 180.0 for a in angle_ticks]
    angle_ticks_rads_plus_offset = [a for a in angle_ticks_rads]
    angle_ticks_for_plot = []
    for i in range(len(angle_ticks)):
        angle_ticks_for_plot.append((angle_ticks_rads_plus_offset[i],r"$"+str(angle_ticks[i])+"$"))

    # ------------------------------------ #

    plot_real_min = 30.0
    plot_real_max = 100.0

    plot_fake_min = 0.0
    plot_fake_max = 5000.0

    rad_tick_increment = 500.0

    radius_ticks = []
    for i in range(int(plot_fake_min),int(plot_fake_max)+int(rad_tick_increment),int(rad_tick_increment)):
        plot_fake_val = ((i-plot_fake_min)/(plot_fake_max-plot_fake_min))*(plot_real_max-plot_real_min)+plot_real_min
        radius_ticks.append((plot_fake_val, r"$"+str(i)+"$"))

    # ------------------------------------ #

    azimuths = np.radians(np.linspace(0, 90, 91))
    azimuths_adjusted = [ (x) for x in azimuths ]
    zeniths = np.arange(0, 5050, 50)
    zeniths_adjusted = [((x-plot_fake_min)/(plot_fake_max-plot_fake_min))*(plot_real_max-plot_real_min)+plot_real_min for x in zeniths]

    # ------------------------------------ #

    r, theta = np.meshgrid(zeniths_adjusted, azimuths_adjusted)
    values = 90.0+5.0*np.random.random((len(azimuths), len(zeniths)))

    # ------------------------------------ #

    scale = 1.0
    aspect = 1.50
    height = 4.0
    fig = plt.figure(1, figsize=(height*aspect*scale, height*scale))
    fig.subplots_adjust(wspace=0.3, left=0.05, right=0.95, top=0.84)
    fig.subplots_adjust()


    ax2, aux_ax2 = setup_arc_radial_axes(fig, 111, angle_ticks_for_plot, radius_ticks, plot_real_min, plot_real_max)





    # r, theta = np.meshgrid(np.log10(Rs),cthetas)
    # values=val
    aux_ax2.contourf(theta, r, values)

    cbar = plt.colorbar(aux_ax2.contourf(theta, r, values), orientation='vertical')
    cbar.ax.set_ylabel('Contour Value [Unit]', fontsize = 16)

    # ax2.axis["bottom"].axes.set_xlabel('Angle [deg]', fontsize=20, weight="bold")

    plt.suptitle('Plot Title ', fontsize = 24, weight="bold")
    plt.legend(loc=3,prop={'size':20})
    # plt.xlabel('Angle [deg]', fontsize=20, weight="bold", rotation=30)
    plt.ylabel('Frequency [Hz]', fontsize=20, weight="bold")

    # plt.show()
    plt.savefig('plot.png', dpi=100)
    plt.show()
    plt.close()

def plot_quarter_plot_ejects_rho(ishell=90, v_n_x="tburst", logx="log10", offset_x_ticks=(0,2),
                                 v_n="Eint2", norm_val=1e30, log="log10", vcenter=1):
    # dfile_jet = h5py.File(curdir+"jet_dynamics_layers.h5", "r")
    dfile_ej = h5py.File(curdir+"magnetar_driven_ej_1.h5", "r")
    print(dfile_ej.keys())
    print(dfile_ej[list(dfile_ej.keys())[0]].keys())
    print(dfile_ej.attrs.keys())

    all_shells = dfile_ej.attrs["nshells"]
    all_layers = dfile_ej.attrs["nlayers"]
    print("all_shells={} all_layers={}".format(all_shells, all_shells))

    cthetas = []
    Rs = []
    val = []
    for ilayer in range(int(all_layers)):
        key = "shell={} layer={}".format(ishell, ilayer)
        cthetas.append( np.array(dfile_ej[key]["ctheta"])[0] )
        Rs = np.array( dfile_ej[key][v_n_x] )
        if (log=="log10"):
            val.append( np.log10( np.array( dfile_ej[key][v_n] ) / norm_val ))
        elif (log=="log2"):
            val.append( np.log2( np.array( dfile_ej[key][v_n] ) / norm_val ))
        else:
            val.append( np.array( dfile_ej[key][v_n]  / norm_val ))

        # val.append( np.array( dfile_ej[key]["beta"] ) * np.array( dfile_ej[key]["Gamma"] ) )
        # val[ilayer] /= val[ilayer][0]
        # val.append( np.array( dfile_ej[key]["GammaRho"] ) * np.array( dfile_ej[key]["GammaRho"] ) )
        # val.append( np.log10( np.array( dfile_ej[key]["rho"] ) / np.array( dfile_ej[key]["rho"] )[-1] ))#* ( np.array( dfile_ej[key]["rho"][-1] ) ) ) )
        # val.append( np.log10( np.array( dfile_ej[key][v_n] ) / norm ))#* ( np.array( dfile_ej[key]["rho"][-1] ) ) ) )
        # val.append( np.log10(np.array( dfile_ej[key]["rho2"] ) / ( np.array( dfile_ej[key]["rho2"][0] ) ) ) )
        # plt.semilogx(Rs, val[-1])
        # plt.show()
    cthetas = np.array(cthetas)#[::-1] # INVERT DATA ----------------------------------------------------------

    print("ctheta: {}".format(cthetas.shape))
    print("{} {}".format(v_n_x, Rs.shape))

    val = np.reshape(np.array(val), (len(cthetas), len(Rs)))
    print("{} {}".format(v_n, val.shape))

    print(cthetas * 180 / np.pi)
    print(Rs)

    # angle_ticks = range(0,100,10)
    angle_ticks = cthetas * 180 / np.pi
    angle_ticks_rads = [a * math.pi / 180.0 for a in angle_ticks]#[::-1] # INVERT TICKS -------------------
    angle_ticks_rads_plus_offset = [a for a in angle_ticks_rads]
    # angle_ticks_rads_plus_offset = angle_ticks_rads_plus_offset[::-1] # Polar Angle
    angle_ticks_for_plot = []
    for i in range(len(angle_ticks)):
        angle_ticks_for_plot.append((angle_ticks_rads_plus_offset[i],r"$"+"{:.0f}".format(angle_ticks[i])+"$"))

    # ------------------------------------ #

    if (logx=="log10"):
        lRs = np.log10(Rs)#/Rs[0])
    elif (logx=="log2"):
        lRs = np.log2(Rs)#/Rs[0])
    else:
        lRs = Rs

    radius_ticks = range(int(lRs[np.isfinite(lRs)].min()+offset_x_ticks[0]),
                         int(lRs[np.isfinite(lRs)].max()+offset_x_ticks[1]),
                         1)
    radius_ticks_for_plot=[]
    for i in range(len(radius_ticks)):
        radius_ticks_for_plot.append((radius_ticks[i], r"$"+str(radius_ticks[i])+"$"))
    #
    # plot_real_min = 30
    # plot_real_max = 100
    #
    # plot_fake_min = 0.0
    # plot_fake_max = 5000.0
    #
    # rad_tick_increment = 500.0
    #
    # radius_ticks = []
    # for i in range(int(plot_fake_min),int(plot_fake_max)+int(rad_tick_increment),int(rad_tick_increment)):
    #     plot_fake_val = ( (i -plot_fake_min) / (plot_fake_max - plot_fake_min)) * (plot_real_max - plot_real_min) + plot_real_min
    #     radius_ticks.append((plot_fake_val, r"$"+str(i)+"$"))
    #
    # # ------------------------------------ #
    #
    # azimuths = np.radians(np.linspace(0, 90, 91))
    # azimuths_adjusted = [ (x) for x in azimuths ]
    # zeniths = np.arange(0, 5050, 50)
    # zeniths_adjusted = [((x-plot_fake_min)/(plot_fake_max-plot_fake_min))*(plot_real_max-plot_real_min)+plot_real_min for x in zeniths]

    # ------------------------------------ #

    # r, theta = np.meshgrid(zeniths_adjusted, azimuths_adjusted)
    # values = 90.0+5.0*np.random.random((len(azimuths), len(zeniths)))

    # ------------------------------------ #

    scale = 1.2
    aspect = 1.20
    height = 3.0
    fig = plt.figure(1, figsize=(height*aspect*scale, height*scale))
    fig.subplots_adjust(wspace=0.01, left=0.1, right=0.95, top=0.93)
    fig.subplots_adjust()


    ax2, aux_ax2 = setup_arc_radial_axes(fig, 111, angle_ticks_for_plot, radius_ticks_for_plot, radius_ticks[0], radius_ticks[-1])

    # ax2.axis["left"].set_lim(50)


    # ax2.set_ylim(10,40)
    r, theta = np.meshgrid(lRs,cthetas)
    values=val
    values[~np.isfinite(values)] = -1
    r[~np.isfinite(r)] = r[np.isfinite(r)].max()
    # r = np.ma.masked_array(r, np.isfinite(r))
    # values = np.ma.masked_array(values, np.isfinite(values))

    levels = MaxNLocator(nbins=40).tick_values(values.min(), values.max())
    cmap = plt.get_cmap('seismic')
    norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
    norm=colors.TwoSlopeNorm(vmin=values.min(), vcenter=vcenter, vmax=values.max())
    # im = aux_ax2.contourf(theta, r, values, levels=np.logspace(np.log10(1e-4),np.log10(1e0),20), cmap=cmap, norm=norm)
    # im = aux_ax2.contourf(theta, r, values, levels=np.linspace(val.min(),val.max(),20), cmap=cmap, norm=norm)
    im = aux_ax2.pcolormesh(theta, r, values, cmap=cmap, norm=norm)
    # norm=colors.LogNorm(vmin=1e-3, vmax=1e0),
    # norm=colors.Normalize(vmin=1e-3, vmax=1e0),
    # cmap='inferno_r', shading='auto')
    # ax2.plot([-1,-1],[-1,-1], marker='d', color='blue', label=r'$\rho<\rho_{\rm j;shock}$')
    # ax2.plot([-1,-1],[-1,-1], marker='d', color='red', label=r'$\rho>\rho_{\rm j;shock}$')

    # Create patch collection with specified colour/alpha
    rect1 = Rectangle((-200,-100), 0, 0, color='blue', label=r'$\rho<\rho_{\rm j;sh}$')
    rect2 = Rectangle((-200,-100), 0, 0, color='red', label=r'$\rho>\rho_{\rm j;sh}$')
    ax2.add_patch(rect1)
    ax2.add_patch(rect2)

    cbar = plt.colorbar(im, orientation='vertical')
    cbar.ax.set_ylabel(r'$\log_{10}( \rho/\rho_{\rm ISM} )$', fontsize = 10)

    # ax2.axis["bottom"].axes.set_xlabel('Angle [deg]', fontsize=20, weight="bold")

    # ax2.axis["left"].label.set_text(r"cz [km$^{-1}$]")
    ax2.set_xlabel('Angle [deg]', fontsize=20, weight="bold", rotation=30)
    ax2.set_ylabel(r'$\log_{10}(t_{\rm b})$', fontsize=20, weight="bold", rotation=30)

    plt.suptitle(' Density profile seen by ejecta ', fontsize = 14, weight="bold")
    plt.legend(**{"fancybox": False, "loc": 'lower left',
                  #"bbox_to_anchor": (0.5, 1.2),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
                  "shadow": "False", "ncol": 1, "fontsize": 10,
                  "framealpha": 0., "borderaxespad": 0., "frameon": False})
    # plt.xlabel('Angle [deg]', fontsize=20, weight="bold", rotation=30)
    # plt.ylabel(r'$\log_{10}(t_{\rm b})$', fontsize=50, weight="bold")

    # plt.show()
    plt.savefig('plot.png', dpi=256)
    plt.show()
    plt.close()

    # # exit(1)
    #
    # fig, ax = plt.subplots(figsize=(4,3),subplot_kw={'projection' :'polar'},ncols=1, nrows=1)
    # ax.set_rscale('log')
    # plt.show()
    #
    #
    # fig = plt.figure()
    fig, ax = plt.subplots(figsize=(4,3),subplot_kw={'projection' :'polar'},ncols=1, nrows=1)
    # ax = plt.subplot(111, polar=True)
    levels = MaxNLocator(nbins=15).tick_values(val.min(), val.max())
    cmap = plt.get_cmap('inferno')
    norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
    im = ax.pcolor(cthetas, np.log10(Rs), val.T, cmap=cmap, norm=norm)
    fig.colorbar(im, ax=ax)
    # ax.set_title('pcolormesh with levels')
    print(ax.get_rmin(), ax.get_rmax())
    # ax.set_rmax(10)
    # ax.set_rmin(20)
    print(ax.get_rmin(), ax.get_rmax())

    max_theta = 90
    ax.set_thetamax(max_theta)
    # ax.set_rlim(10,20)

    # ax.set_rlim(Rs.min(), Rs.max())
    # ax.set_rscale("log")
    # ax.set_rscale('log')
    #
    # ticklabels = ax.get_yticklabels()
    # labels = range(80, 0, -10)
    # for i in range(0, len(labels)):
    #     ticklabels[i] = str(labels[i])
    # ax.set_yticklabels(ticklabels)

    plt.show()
plot_quarter_plot_ejects_rho(ishell=1, v_n_x="tburst", logx="log10", offset_x_ticks=(0,2),
                             v_n="rho", norm_val=pow(10, -4.1) * cgs.mp, log="log10", vcenter=0)


def plot_dynamics_pretty_plot2_gammaCBM(v_n_x = "tburst", v_n_y = "GammaCBM"):

    dfile_jet = h5py.File(curdir+"jet_dynamics_layers.h5", "r")
    dfile_ej = h5py.File(curdir+"ejecta_dynamics_shells_layers.h5", "r")

    def plot_ejecta_rho(ax, xlim, ylim):
        ishell = 1
        layers = [###"shell={} layer=0".format(ishell),
            # "shell={} layer=2".format(ishell),
            # "shell={} layer=4".format(ishell),
            # "shell={} layer=6".format(ishell),
            "shell={} layer=8".format(ishell),
            # "shell={} layer=10".format(ishell),
            # "shell={} layer=12".format(ishell),
            # "shell={} layer=14".format(ishell),
            ###"shell={} layer=16".format(ishell),
            # "shell={} layer=18".format(ishell),
            # "shell={} layer=20".format(ishell),
            # "shell={} layer=22".format(ishell),
            ###"shell={} layer=23".format(ishell)
        ]
        norm = Normalize(vmin=0, vmax=dfile_ej.attrs["nlayers"])
        cmap = cm.Reds_r
        for il, layer in enumerate(layers):

            ijl = "layer=69" # the one to vompare with
            beta_jet = np.array(dfile_jet[ijl]["beta"])
            R_jet = np.array(dfile_jet[ijl]["R"])
            beta_ej = np.array(dfile_ej[layer]["beta"])
            R_ej = np.array(dfile_ej[layer]["R"])
            beta_jet = np.array(dfile_ej[layer][v_n_x])
            beta12 = (beta_ej - beta_jet) / (1. - beta_ej * beta_jet)
            coll_idx = np.argmax( R_jet < R_ej )
            cs_jet2 = 5. * beta_jet ** 2 / 9.

            x_arr = np.array(dfile_ej[layer][v_n_x])
            y_arr = np.array(dfile_ej[layer]["rho"])
            # y_arr = y_arr / y_arr[0]
            if (len(xlim)>0):
                y_arr = y_arr[(x_arr > xlim[0]) & (x_arr < xlim[-1])]
                x_arr = x_arr[(x_arr > xlim[0]) & (x_arr < xlim[-1])]
            if (len(ylim)>0):
                x_arr = x_arr[(y_arr > ylim[0]) & (y_arr < ylim[-1])]
                y_arr = y_arr[(y_arr > ylim[0]) & (y_arr < ylim[-1])]
            ax.plot(x_arr, y_arr, marker='.', color=cmap(norm(int(layer.split("=")[-1]))))

            # x_arr = np.array(dfile_ej[layer][v_n_x])
            # y_arr = Beta( np.array(dfile_ej[layer]["Gamma"]) ) * np.array(dfile_ej[layer]["Gamma"])
            # y_arr = y_arr / y_arr[0]
            # if (len(xlim)>0):
            #     y_arr = y_arr[(x_arr > xlim[0]) & (x_arr < xlim[-1])]
            #     x_arr = x_arr[(x_arr > xlim[0]) & (x_arr < xlim[-1])]
            # if (len(ylim)>0):
            #     x_arr = x_arr[(y_arr > ylim[0]) & (y_arr < ylim[-1])]
            #     y_arr = y_arr[(y_arr > ylim[0]) & (y_arr < ylim[-1])]
            # ax.plot(x_arr, y_arr, ls=':', color=cmap(norm(int(layer.split("=")[-1]))))
            # ax.plot(x_arr, y_arr, ls='-', color=cmap(norm(int(layer.split("=")[-1]))))
    def plot_ejecta(ax, xlim, ylim):
        ishell = 2
        layers = [###"shell={} layer=0".format(ishell),
            # "shell={} layer=2".format(ishell),
            # "shell={} layer=4".format(ishell),
            # "shell={} layer=6".format(ishell),
            "shell={} layer=8".format(ishell),
            # "shell={} layer=10".format(ishell),
            # "shell={} layer=12".format(ishell),
            # "shell={} layer=14".format(ishell),
            ###"shell={} layer=16".format(ishell),
            # "shell={} layer=18".format(ishell),
            # "shell={} layer=20".format(ishell),
            # "shell={} layer=22".format(ishell),
            ###"shell={} layer=23".format(ishell)
        ]
        norm = Normalize(vmin=0, vmax=dfile_ej.attrs["nlayers"])
        cmap = cm.Reds_r
        for il, layer in enumerate(layers):

            ijl = "layer=69" # the one to vompare with
            beta_jet = np.array(dfile_jet[ijl]["beta"])
            R_jet = np.array(dfile_jet[ijl]["R"])
            beta_ej = np.array(dfile_ej[layer]["beta"])
            R_ej = np.array(dfile_ej[layer]["R"])
            beta_jet = np.array(dfile_ej[layer][v_n_x])
            beta12 = (beta_ej - beta_jet) / (1. - beta_ej * beta_jet)
            coll_idx = np.argmax( R_jet < R_ej )
            cs_jet2 = 5. * beta_jet ** 2 / 9.

            x_arr = np.array(dfile_ej[layer][v_n_x])
            y_arr = Beta( np.array(dfile_ej[layer]["GammaRho"]) ) * np.array(dfile_ej[layer]["GammaRho"])
            y_arr = y_arr / y_arr[0]
            if (len(xlim)>0):
                y_arr = y_arr[(x_arr > xlim[0]) & (x_arr < xlim[-1])]
                x_arr = x_arr[(x_arr > xlim[0]) & (x_arr < xlim[-1])]
            if (len(ylim)>0):
                x_arr = x_arr[(y_arr > ylim[0]) & (y_arr < ylim[-1])]
                y_arr = y_arr[(y_arr > ylim[0]) & (y_arr < ylim[-1])]
            ax.plot(x_arr, y_arr, ls='-', color=cmap(norm(int(layer.split("=")[-1]))))

            x_arr = np.array(dfile_ej[layer][v_n_x])
            y_arr = Beta( np.array(dfile_ej[layer]["Gamma"]) ) * np.array(dfile_ej[layer]["Gamma"])
            y_arr = y_arr / y_arr[0]
            if (len(xlim)>0):
                y_arr = y_arr[(x_arr > xlim[0]) & (x_arr < xlim[-1])]
                x_arr = x_arr[(x_arr > xlim[0]) & (x_arr < xlim[-1])]
            if (len(ylim)>0):
                x_arr = x_arr[(y_arr > ylim[0]) & (y_arr < ylim[-1])]
                y_arr = y_arr[(y_arr > ylim[0]) & (y_arr < ylim[-1])]
            ax.plot(x_arr, y_arr, ls=':', color=cmap(norm(int(layer.split("=")[-1]))))
            # ax.plot(x_arr, y_arr, ls='-', color=cmap(norm(int(layer.split("=")[-1]))))
        # -----------------------------------------------------
        ishell = 25
        layers2 = [###"shell={} layer=0".format(ishell),
            # "shell={} layer=2".format(ishell),
            # "shell={} layer=4".format(ishell),
            # "shell={} layer=6".format(ishell),
            "shell={} layer=8".format(ishell),
            # "shell={} layer=10".format(ishell),
            # "shell={} layer=12".format(ishell),
            # "shell={} layer=14".format(ishell),
            ###"shell={} layer=16".format(ishell),
            # "shell={} layer=18".format(ishell),
            # "shell={} layer=20".format(ishell),
            # "shell={} layer=22".format(ishell),
            ###"shell={} layer=23".format(ishell)
        ]
        layers=layers2
        norm = Normalize(vmin=0, vmax=dfile_ej.attrs["nlayers"])
        cmap = cm.Greens_r
        for il, layer in enumerate(layers):

            ijl = "layer=69" # the one to vompare with
            beta_jet = np.array(dfile_jet[ijl]["beta"])
            R_jet = np.array(dfile_jet[ijl]["R"])
            beta_ej = np.array(dfile_ej[layer]["beta"])
            R_ej = np.array(dfile_ej[layer]["R"])
            beta_jet = np.array(dfile_ej[layer][v_n_x])
            beta12 = (beta_ej - beta_jet) / (1. - beta_ej * beta_jet)
            coll_idx = np.argmax( R_jet < R_ej )
            cs_jet2 = 5. * beta_jet ** 2 / 9.

            x_arr = np.array(dfile_ej[layer][v_n_x])
            y_arr = Beta( np.array(dfile_ej[layer]["GammaRho"]) ) * np.array(dfile_ej[layer]["GammaRho"])
            if (len(xlim)>0):
                y_arr = y_arr[(x_arr > xlim[0]) & (x_arr < xlim[-1])]
                x_arr = x_arr[(x_arr > xlim[0]) & (x_arr < xlim[-1])]
            if (len(ylim)>0):
                x_arr = x_arr[(y_arr > ylim[0]) & (y_arr < ylim[-1])]
                y_arr = y_arr[(y_arr > ylim[0]) & (y_arr < ylim[-1])]
            ax.plot(x_arr, y_arr, ls='-', color=cmap(norm(int(layer.split("=")[-1]))))

            x_arr = np.array(dfile_ej[layer][v_n_x])
            y_arr = Beta( np.array(dfile_ej[layer]["Gamma"]) ) * np.array(dfile_ej[layer]["Gamma"])
            if (len(xlim)>0):
                y_arr = y_arr[(x_arr > xlim[0]) & (x_arr < xlim[-1])]
                x_arr = x_arr[(x_arr > xlim[0]) & (x_arr < xlim[-1])]
            if (len(ylim)>0):
                x_arr = x_arr[(y_arr > ylim[0]) & (y_arr < ylim[-1])]
                y_arr = y_arr[(y_arr > ylim[0]) & (y_arr < ylim[-1])]
            ax.plot(x_arr, y_arr, ls=':', color=cmap(norm(int(layer.split("=")[-1]))))
    def plot_ejecta_dlnrhodr1(ax, xlim, ylim):
        # ishell = 2
        # layers = ["shell={} layer=0".format(ishell),
        #           # "shell={} layer=2".format(ishell),
        #           # "shell={} layer=4".format(ishell),
        #           # "shell={} layer=6".format(ishell),
        #           "shell={} layer=8".format(ishell),
        #           # "shell={} layer=10".format(ishell),
        #           # "shell={} layer=12".format(ishell),
        #           # "shell={} layer=14".format(ishell),
        #           "shell={} layer=16".format(ishell),
        #           # "shell={} layer=18".format(ishell),
        #           # "shell={} layer=20".format(ishell),
        #           # "shell={} layer=22".format(ishell),
        #           "shell={} layer=23".format(ishell)
        #           ]
        # norm = Normalize(vmin=0, vmax=dfile_ej.attrs["nlayers"])
        # cmap = cm.Blues_r
        # for il, layer in enumerate(layers):
        #
        #     ijl = "layer=69" # the one to vompare with
        #     beta_jet = np.array(dfile_jet[ijl]["beta"])
        #     R_jet = np.array(dfile_jet[ijl]["R"])
        #     beta_ej = np.array(dfile_ej[layer]["beta"])
        #     R_ej = np.array(dfile_ej[layer]["R"])
        #     beta_jet = np.array(dfile_ej[layer][v_n_x])
        #     beta12 = (beta_ej - beta_jet) / (1. - beta_ej * beta_jet)
        #     coll_idx = np.argmax( R_jet < R_ej )
        #     cs_jet2 = 5. * beta_jet ** 2 / 9.
        #
        #     x_arr = np.array(dfile_ej[layer][v_n_x])
        #     y_arr =  np.array(dfile_ej[layer]["dlnrhodr"])
        #     if (len(xlim)>0):
        #         y_arr = y_arr[(x_arr > xlim[0]) & (x_arr < xlim[-1])]
        #         x_arr = x_arr[(x_arr > xlim[0]) & (x_arr < xlim[-1])]
        #     if (len(ylim)>0):
        #         x_arr = x_arr[(y_arr > ylim[0]) & (y_arr < ylim[-1])]
        #         y_arr = y_arr[(y_arr > ylim[0]) & (y_arr < ylim[-1])]
        #     ax.plot(x_arr, y_arr, ls='-', color=cmap(norm(int(layer.split("=")[-1]))))

        # x_arr = np.array(dfile_ej[layer][v_n_x])
        # y_arr =  np.array(dfile_ej[layer]["dlnrhodr"])
        # if (len(xlim)>0):
        #     y_arr = y_arr[(x_arr > xlim[0]) & (x_arr < xlim[-1])]
        #     x_arr = x_arr[(x_arr > xlim[0]) & (x_arr < xlim[-1])]
        # if (len(ylim)>0):
        #     x_arr = x_arr[(y_arr > ylim[0]) & (y_arr < ylim[-1])]
        #     y_arr = y_arr[(y_arr > ylim[0]) & (y_arr < ylim[-1])]
        # ax.plot(x_arr, y_arr, ls=':', color=cmap(norm(int(layer.split("=")[-1]))))

        # ax.plot(x_arr, y_arr, ls='-', color=cmap(norm(int(layer.split("=")[-1]))))
        # -----------------------------------------------------
        ishell = 2
        layers2 = [###"shell={} layer=0".format(ishell),
            # "shell={} layer=2".format(ishell),
            # "shell={} layer=4".format(ishell),
            # "shell={} layer=6".format(ishell),
            "shell={} layer=8".format(ishell),
            # "shell={} layer=10".format(ishell),
            # "shell={} layer=12".format(ishell),
            # "shell={} layer=14".format(ishell),
            ###"shell={} layer=16".format(ishell),
            # "shell={} layer=18".format(ishell),
            # "shell={} layer=20".format(ishell),
            # "shell={} layer=22".format(ishell),
            ###"shell={} layer=23".format(ishell)
        ]
        layers=layers2
        norm = Normalize(vmin=0, vmax=dfile_ej.attrs["nlayers"])
        cmap = cm.Greens_r
        for il, layer in enumerate(layers):

            ijl = "layer=69" # the one to vompare with
            beta_jet = np.array(dfile_jet[ijl]["beta"])
            R_jet = np.array(dfile_jet[ijl]["R"])
            beta_ej = np.array(dfile_ej[layer]["beta"])
            R_ej = np.array(dfile_ej[layer]["R"])
            beta_jet = np.array(dfile_ej[layer][v_n_x])
            beta12 = (beta_ej - beta_jet) / (1. - beta_ej * beta_jet)
            coll_idx = np.argmax( R_jet < R_ej )
            cs_jet2 = 5. * beta_jet ** 2 / 9.

            y_arr1 = np.array(dfile_ej[layer]["drhodr"])
            y_arr2 = np.copy(y_arr1)
            print("Mom CBM[{}, {}]".format(y_arr1.min(),y_arr1.max()))
            print("MOm [{}, {}]".format(y_arr2.min(),y_arr2.max()))

            mask0 = y_arr1 == 0.
            if (len(y_arr2[mask0]) > 1):
                cmap = cm.Greys_r
                y_arr = y_arr2[mask0]
                x_arr = np.array(dfile_ej[layer][v_n_x])[mask0]
                y_arr = y_arr + 1e-48# / y_arr[0]
                if (len(xlim)>0):
                    y_arr = y_arr[(x_arr > xlim[0]) & (x_arr < xlim[-1])]
                    x_arr = x_arr[(x_arr > xlim[0]) & (x_arr < xlim[-1])]
                if (len(ylim)>0):
                    x_arr = x_arr[(y_arr > ylim[0]) & (y_arr < ylim[-1])]
                    y_arr = y_arr[(y_arr > ylim[0]) & (y_arr < ylim[-1])]
                ax2.plot(x_arr, y_arr, marker='.', ls='none', color=cmap(norm(int(layer.split("=")[-1]))))


            mask1 = y_arr1 > 0.
            if (len(y_arr2[mask1]) > 1):
                cmap = cm.Greens_r
                y_arr = y_arr2[mask1]
                x_arr = np.array(dfile_ej[layer][v_n_x])[mask1]
                y_arr = y_arr# / y_arr[0]
                if (len(xlim)>0):
                    y_arr = y_arr[(x_arr > xlim[0]) & (x_arr < xlim[-1])]
                    x_arr = x_arr[(x_arr > xlim[0]) & (x_arr < xlim[-1])]
                if (len(ylim)>0):
                    x_arr = x_arr[(y_arr > ylim[0]) & (y_arr < ylim[-1])]
                    y_arr = y_arr[(y_arr > ylim[0]) & (y_arr < ylim[-1])]
                ax2.plot(x_arr, y_arr, marker='.', ls='none', color=cmap(norm(int(layer.split("=")[-1]))))

            mask2 = y_arr1 < 0.
            if (len(y_arr2[mask2])>1):
                cmap = cm.Reds_r
                y_arr = y_arr2[mask2]
                x_arr = np.array(dfile_ej[layer][v_n_x])[mask2]
                y_arr = -y_arr# / y_arr[0]
                if (len(xlim)>0):
                    y_arr = y_arr[(x_arr > xlim[0]) & (x_arr < xlim[-1])]
                    x_arr = x_arr[(x_arr > xlim[0]) & (x_arr < xlim[-1])]
                if (len(ylim)>0):
                    x_arr = x_arr[(y_arr > ylim[0]) & (y_arr < ylim[-1])]
                    y_arr = y_arr[(y_arr > ylim[0]) & (y_arr < ylim[-1])]
                ax2.plot(x_arr, y_arr, marker='.', ls='none', color=cmap(norm(int(layer.split("=")[-1]))))

    fig, axes = plt.subplots(figsize=(4,5), nrows=2, ncols=1, sharex="all")


    ax1 = axes[0]
    ax2 = ax1.twinx()
    # T_E = np.arange(1,max(T)+1,1)
    # The data.
    # plot_jet(ax1, (), ())
    plot_ejecta(ax1, (), ())
    # ax1.plot(T, Cp, 'x', c='b', mew=2, alpha=0.8, label='Experiment')
    # The Einstein fit.
    # ax1.plot(T_E, CV_E, c='m', lw=2, alpha=0.5, label='Einstein model')
    # ax1.set_xlabel(r'Radius [cm]')#(r'$T\,/\mathrm{K}$')
    # ax1.set_ylabel(r'Momentum $\Gamma\beta$')#(r'$C_p\,/\mathrm{J\,K^{-1}\,mol^{-1}}$')

    ax1.set_ylabel(r"momentum $\Gamma \beta$")
    ax1.set_xlabel(r"time=$\int dr \beta c$ [s]")
    ax1.legend(loc=0)

    ax1.plot([1e15,1e16],[1e-4,1e-4],color='green',lw=2,label=r'dlnrhodr $ > $0')
    ax1.plot([1e15,1e16],[1e-4,1e-4],color='red',lw=2,label=r'dlnrhodr $ < $0')
    ax1.plot([1e15,1e16],[1e-4,1e-4],color='gray',lw=2,label=r'dlnrhodr $ = $0')
    ax1.legend(loc="lower left")
    # Some ad hoc tweaks.
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.set_xlim(1e7, 1e10)#(7e6,5e10)
    ax1.set_ylim(0.994, 1.002)#(1e-1,1.1e0)
    ax1.set_title("GRB and kN ejecta dynamics")
    ax1.tick_params(axis='both', which='both', labelleft=True,
                    labelright=False, tick1On=True, tick2On=True,
                    # labelsize=plotdic["fontsize"],
                    direction='in',
                    bottom=True, top=True, left=True, right=True)

    plot_ejecta_dlnrhodr1(ax2, (), ())
    ax2.set_xscale("log")
    ax2.set_yscale("log")

    ax1 = axes[1]
    plot_ejecta_rho(ax1, (), ())
    ax1.set_xscale("log")
    ax1.set_yscale("log")

    # # Create a set of inset Axes: these should fill the bounding box allocated to them.
    # ax2 = plt.axes([0,0,2,2])
    # # Manually set the position and relative size of the inset axes within ax1
    # ip = InsetPosition(ax1, [0.18,0.1,0.8,0.4]) # [0.4,0.2,0.5,0.5]
    # ax2.set_axes_locator(ip)
    # # Mark the region corresponding to the inset axes on ax1 and draw lines in grey linking the two axes.
    # mark_inset(ax1, ax2, loc1=2, loc2=4, fc="none", ec='0.5')
    #
    # # The data: only display for low temperature in the inset figure.
    # # plot_jet(ax2, (5e7,5e8), (6e-2,8e-1))
    # plot_ejecta(ax2, (1e8, 1e11), (0.,0.055))
    #
    #
    # ax2.set_xscale("log")
    # ax2.set_yscale("log")
    #
    # ax2.tick_params(axis='both', which='both', labelleft=True,
    #                 labelright=False, tick1On=True, tick2On=True,
    #                 # labelsize=plotdic["fontsize"],
    #                 direction='in',
    #                 bottom=True, top=True, left=True, right=True)



    # ax2.set_xticks(np.array([3e18,3e18]))
    # ax2.set_xticklabels(ax2.get_xticks(), backgroundcolor='w')
    # ax2.set_yticklabels(ax2.get_yticks(), backgroundcolor='w')
    # ax2.tick_params(axis='x', which='major', pad=8)
    # ax2.tick_params(axis='y', which='major', pad=8)
    plt.tight_layout()
    plt.savefig(curdir+"plot_pretty_dynamicsCBM.png",dpi=256)
    plt.show()

    # ax = axes[0]
    # ax.plot(dyn2t2.get("tburst"),dyn2t2.get("beta"),ls='-', color='blue',label='Jet BW')
    # ax.plot(dyn2t2.get("tburst"), dyn2t2.get("beta_ej"),ls='-', color='red',label='kN BW')
    # ax.axvline(x=dyn2t2.get("tburst")[np.argmax( dyn2t2.get("beta") < dyn2t2.get("beta_ej") )], color="gray", linestyle="-",label=r"$R$ at $\beta_{\rm ej}=\beta_{\rm jet}$")
    # ax.axvline(x=dyn2t2.get("tburst")[np.argmax( dyn2t2.get("R") < dyn2t2.get("R_ej") )], color="gray", linestyle="--",
    #            label=r"$R$ at $R_{\rm ej}=R_{\rm jet}$, "
    #                  r"$\beta_{12}=$"+"{:.2f}".format(beta12[coll_idx])
    #                  +" $\Gamma_{12}=$"+"{:.2f}".format(get_Gamma(beta12[coll_idx]))
    #                  +" $c_s=$"+"{:.2f}".format(np.sqrt(cs_jet2[coll_idx]))
    #                  +" $\mathcal{M}_s=$"+"{:.2f}".format(beta12[coll_idx]/np.sqrt(cs_jet2[coll_idx])))
    # ax.set_xlabel(r"time=$\int dr \beta c$ [s]")
    # ax.set_ylabel(r"velocity $\beta$ [c] (solid lines)")
    # ax.set_xscale("log")
    # ax.set_yscale("log")
    # # ax.set_xlim(3e17,1e19)
    # ax.set_xlim(2e6,2e9)
    # ax.set_ylim(8e-2,1.1e0)
    # ax.legend(loc="best")
    # ax.tick_params(axis='both', which='both', labelleft=True,
    #                labelright=False, tick1On=True, tick2On=True,
    #                labelsize=12,
    #                direction='in',
    #                bottom=True, top=True, left=True, right=True)
    # ax2 = ax.twinx()
    # ax2.plot(dyn2t2.get("tburst"), dyn2t2.get("R_ej"), ls='--', color='lime')
    # ax2.plot(dyn2t2.get("tburst"), np.abs(dyn2t2.get("R") - dyn2t2.get("R_ej")), ls='--', color='lime', label=r"$| \Delta R |$")
    # ax2.plot(dyn2t2.get("tburst"), dyn2t2.get("R_ej"), ls='--', color='red')
    # ax2.legend(loc="best")
    # ax2.set_ylabel("Radius [cm] (dashed lines)")
    # ax2.set_xscale("log")
    # ax2.set_yscale("log")
    #
    # ctheta_j = dyn2t2.kwargs["ctheta_j"] + 0.5 * (2.0 * dyn2t2.get("theta") - 2.0 * dyn2t2.kwargs["theta_w_j"])
    # ctheta_ej = dyn2t2.kwargs["ctheta_ej"] + 0.5 * (2.0 * dyn2t2.get("theta_ej") - 2.0 * dyn2t2.kwargs["theta_w_ej"])
    #
    # plt.tight_layout()
    # plt.savefig("time_intersection2.png",dpi=256)
    # plt.show()
plot_dynamics_pretty_plot2_gammaCBM()


def plot_compare_gauss_dyn(fname = "dyn_gauss"):
    layers = [9, 19, 29, 39, 49]

    fid, axes = plt.subplots(ncols=4, nrows=1)

    for il in layers:
        names, data = load_data( curdir + fname + "pw_layer_" + str(il)+ ".txt" )
        rnames, rdata = load_data( curdir + fname + "_layer_" + str(il)+ "_ref.txt" )

        i = 0
        axes[i].plot(data[:,names.index("R")], data[:,names.index("Gamma")]*Beta(data[:,names.index("Gamma")]), ls='-',color='black', label='PBA')
        axes[i].plot(rdata[:,rnames.index("R")], rdata[:,rnames.index("Gamma")]*Beta(rdata[:,rnames.index("Gamma")]), ls='--',color='gray', label='ref')
        axes[i].set_xlabel("R")
        axes[i].set_ylabel("GammaBeta")
        axes[i].set_xscale("log")
        axes[i].set_yscale("log")
        axes[i].legend()
        axes[i].grid()

        i = 1
        axes[i].plot(data[:,names.index("R")], data[:,names.index("M2")], ls='-',color='black', label='PBA')
        axes[i].plot(rdata[:,rnames.index("R")], rdata[:,rnames.index("M2")], ls='--',color='gray', label='ref')
        axes[i].set_xlabel("R")
        axes[i].set_ylabel("M2")
        axes[i].set_xscale("log")
        axes[i].set_yscale("log")
        axes[i].legend()
        axes[i].grid()

        i = 2
        axes[i].plot(data[:,names.index("R")], data[:,names.index("TT")], ls='-',color='black', label='PBA')
        axes[i].plot(rdata[:,rnames.index("R")], rdata[:,rnames.index("TT")], ls='--',color='gray', label='Old')
        axes[i].set_xlabel("R")
        axes[i].set_ylabel("TT")
        axes[i].set_xscale("log")
        axes[i].set_yscale("log")
        axes[i].legend()
        axes[i].grid()

        i = 3
        axes[i].plot(data[:,names.index("R")], data[:,names.index("theta")], ls='-',color='black', label='PBA')
        axes[i].plot(rdata[:,rnames.index("R")], rdata[:,rnames.index("theta")], ls='--',color='gray', label='Old')
        axes[i].set_xlabel("R")
        axes[i].set_ylabel("theta")
        axes[i].set_xscale("log")
        axes[i].set_yscale("log")
        axes[i].legend()
        axes[i].grid()

    plt.show()

    fid, axes = plt.subplots(ncols=2, nrows=1)

    names, data = load_data( curdir + fname + "_layers.txt" )
    rnames, rdata = load_data( curdir + fname + "_layers_ref.txt" )

    i = 0
    axes[i].plot(rdata[:,rnames.index("cthetas")], np.log10(rdata[:,rnames.index("E0")]), marker='o',color='blue', label='E0',ms=6.)
    axes[i].plot(rdata[:,rnames.index("cthetas")], np.log10(rdata[:,rnames.index("G0")]), marker='o',color='red', label='G0',ms=6.)
    axes[i].plot(rdata[:,rnames.index("cthetas")], np.log10(rdata[:,rnames.index("M0")]), marker='o',color='green', label='M0',ms=6.)


    axes[i].plot(data[:,names.index("cthetas")], np.log10(data[:,names.index("E0")]), marker='x',color='blue',ms=12.)
    axes[i].plot(data[:,names.index("cthetas")], np.log10(data[:,names.index("G0")]), marker='x',color='red',ms=12.)
    axes[i].plot(data[:,names.index("cthetas")], np.log10(data[:,names.index("M0")]), marker='x',color='green',ms=12.)

    axes[i].set_xlabel("Data")
    axes[i].set_ylabel("cthetas")
    axes[i].set_xscale("linear")
    axes[i].set_yscale("linear")
    axes[i].legend()
    axes[i].grid()

    plt.show()
###plot_compare_gauss_dyn()

def plot_compare_gauss_dyn_a(fname = "dyn_gauss"):
    layers = [0, 2, 4, 6, 9]

    fid, axes = plt.subplots(ncols=2, nrows=1)

    for il in layers:
        names, data = load_data( curdir + fname + "a_layer_" + str(il)+ ".txt" )
        rnames, rdata = load_data( curdir + fname + "a_layer_" + str(il)+ "_ref.txt" )

        i = 0
        axes[i].plot(data[:,names.index("R")], data[:,names.index("Gamma")]*Beta(data[:,names.index("Gamma")]), ls='-',color='black', label='PBA')
        axes[i].plot(rdata[:,rnames.index("R")], rdata[:,rnames.index("Gamma")]*Beta(rdata[:,rnames.index("Gamma")]), ls='--',color='gray', label='ref')
        axes[i].set_xlabel("R")
        axes[i].set_ylabel("GammaBeta")
        axes[i].set_xscale("log")
        axes[i].set_yscale("log")
        axes[i].legend()
        axes[i].grid()

        i = 1
        # axes[i].plot(data[:,names.index("R")], data[:,names.index("M2")], ls='-',color='black', label='PBA')
        # axes[i].plot(rdata[:,rnames.index("R")], rdata[:,rnames.index("M2")], ls='--',color='gray', label='ref')
        # axes[i].set_xlabel("R")
        # axes[i].set_ylabel("M2")
        # axes[i].set_xscale("log")
        # axes[i].set_yscale("log")
        # axes[i].legend()
        # axes[i].grid()

        i = 2
        # axes[i].plot(data[:,names.index("R")], data[:,names.index("TT")], ls='-',color='black', label='PBA')
        # axes[i].plot(rdata[:,rnames.index("R")], rdata[:,rnames.index("TT")], ls='--',color='gray', label='Old')
        # axes[i].set_xlabel("R")
        # axes[i].set_ylabel("TT")
        # axes[i].set_xscale("log")
        # axes[i].set_yscale("log")
        # axes[i].legend()
        # axes[i].grid()

        i = 1
        axes[i].plot(data[:,names.index("R")], data[:,names.index("theta")], ls='-',color='black', label='PBA')
        axes[i].plot(rdata[:,rnames.index("R")], rdata[:,rnames.index("theta")], ls='--',color='gray', label='Old')
        axes[i].set_xlabel("R")
        axes[i].set_ylabel("theta")
        axes[i].set_xscale("log")
        axes[i].set_yscale("log")
        axes[i].legend()
        axes[i].grid()

    plt.show()

    fid, axes = plt.subplots(ncols=2, nrows=1)

    names, data = load_data( curdir + fname + "_layers.txt" )
    rnames, rdata = load_data( curdir + fname + "_layers_ref.txt" )

    i = 0
    axes[i].plot(rdata[:,rnames.index("cthetas")], np.log10(rdata[:,rnames.index("E0")]), marker='o',color='blue', label='E0',ms=6.)
    axes[i].plot(rdata[:,rnames.index("cthetas")], np.log10(rdata[:,rnames.index("G0")]), marker='o',color='red', label='G0',ms=6.)
    axes[i].plot(rdata[:,rnames.index("cthetas")], np.log10(rdata[:,rnames.index("M0")]), marker='o',color='green', label='M0',ms=6.)


    axes[i].plot(data[:,names.index("cthetas")], np.log10(data[:,names.index("E0")]), marker='x',color='blue',ms=12.)
    axes[i].plot(data[:,names.index("cthetas")], np.log10(data[:,names.index("G0")]), marker='x',color='red',ms=12.)
    axes[i].plot(data[:,names.index("cthetas")], np.log10(data[:,names.index("M0")]), marker='x',color='green',ms=12.)

    axes[i].set_xlabel("Data")
    axes[i].set_ylabel("cthetas")
    axes[i].set_xscale("linear")
    axes[i].set_yscale("linear")
    axes[i].legend()
    axes[i].grid()

    plt.show()
### plot_compare_gauss_dyn_a()

class quarter_plots():
    def __init__(self):
        pass

    # axis
    def setup_arc_radial_axes(self, fig, rect, angle_ticks, radius_ticks, min_rad, max_rad):

        tr = PolarAxes.PolarTransform()
        # tr = LogPolarTransform()
        # tr = InvertedLogPolarTransform()

        pi = np.pi

        grid_locator1 = FixedLocator([v for v, s in angle_ticks])
        tmp = dict(angle_ticks)
        # tmp["rotation"] = 40
        tick_formatter1 = DictFormatter(tmp)


        grid_locator2 = FixedLocator([a for a, b in radius_ticks])
        tick_formatter2 = DictFormatter(dict(radius_ticks))

        grid_helper = floating_axes.GridHelperCurveLinear(tr,
                                                          # extremes=((370.0*(pi/180.0)), (170.0*(pi/180.0)), max_rad, min_rad),
                                                          extremes=((90.0*(pi/180.0)), (00.0*(pi/180.0)), max_rad, min_rad),
                                                          # extremes=((170.0*(pi/180.0)), (370.0*(pi/180.0)), max_rad, min_rad),
                                                          grid_locator1=grid_locator1,
                                                          grid_locator2=grid_locator2,
                                                          tick_formatter1=tick_formatter1,
                                                          tick_formatter2=tick_formatter2,
                                                          )

        ax1 = floating_axes.FloatingSubplot(fig, rect, grid_helper=grid_helper)
        fig.add_subplot(ax1)

        # ax1.grid(True)

        # create a parasite axes whose transData in RA, cz
        aux_ax = ax1.get_aux_axes(tr)

        aux_ax.patch = ax1.patch
        ax1.patch.zorder=0.9

        ax1.set_xlabel('Angle [deg]', fontsize=20, weight="bold", rotation=30)

        # ax1.axis["left"].set_ticklabel_direction("+")
        ax1.axis["bottom"].set_axis_direction("top")
        ax1.axis["bottom"].set_ticklabel_direction("+")
        ax1.axis["bottom"].label.set_rotation(0)
        ax1.axis["bottom"].label.set_pad(-20)
        # ax1.axis["bottom"].toggle(ticklabels=True, label=True)
        # ax1.axis["bottom"].major_ticklabels.set_axis_direction("left")





        # ax = aux_ax
        # labels = []
        # for i in range(len(angles)) :
        #     angles[i] = 0
        # for label, angle in zip(ax.get_xticklabels(), angles):
        #     x,y = label.get_position()
        #     lab = ax.text(x,y, label.get_text(), transform=label.get_transform(),
        #                   ha=label.get_ha(), va=label.get_va())
        #     lab.set_rotation(angle)
        #     labels.append(lab)
        # ax.set_xticklabels([])
        # ax1.axis["bottom"].set_axis_direction("bottom")
        # ax1.axis["bottom"].set_tick_label_direction("bottom")
        # ax1.axis["bottom"].set_axis_direction("left")
        # ax1.axis["bottom"].set_rotation(0)
        # ax1.axis["bottom"].set_axis_direction("top")

        return ax1, aux_ax

    # plot full evolution jet
    def plot_full_evol_jet(self):
        dfile_jet = h5py.File(curdir+"jet_dynamics_layers.h5", "r")
        nlayers = dfile_jet.keys()
        print("layers={}".format(nlayers))

        # cthetas = []
        # cthetas0 = []
        # val = []
        # Rs0 = []
        # Rs = []

        r_grid = np.array( dfile_jet["layer={}".format(0)]["R"] )

        cthetas = np.empty_like(r_grid)
        rs = np.empty_like(r_grid)
        val = np.empty_like(r_grid)

        for ilayer in range(len(nlayers)):
            key = "layer={}".format(ilayer)
            cthetas = np.vstack(( cthetas, np.array(dfile_jet[key]["ctheta"]) ))
            rs = np.vstack(( rs, np.array(dfile_jet[key]["R"]) ))
            i_val = np.array( dfile_jet[key]["beta"] ) * np.array( dfile_jet[key]["Gamma"] )
            val = np.vstack(( val, i_val ))

            # cthetas.append( np.array(dfile_jet[key]["ctheta0"]) )
            # cthetas0.append( np.array(dfile_jet[key]["ctheta0"][0] ) )
            # Rs0 = np.array( dfile_jet[key]["R"] )
            # Rs.append( np.array( dfile_jet[key]["R"] ) )
            # val.append( np.array( dfile_jet[key]["beta"] ) * np.array( dfile_jet[key]["Gamma"] ) )
            # val[ilayer] /= val[ilayer][0]

        print(cthetas[1,0]*180/np.pi, cthetas[1,-1]*180/np.pi)

        rs = rs[1:, :]
        cthetas = cthetas[1:, :]
        val = val[1:, :]

        # rs = rs[::-1, :]
        # cthetas = cthetas[::-1, :]
        # val = val[::-1, :]

        print("Jet Data loaded")

        # cthetas = np.array(cthetas)#[::-1] # INVERT DATA ----------------------------------------------------------
        # cthetas0 = np.array(cthetas0)#[::-1] # INVERT DATA ----------------------------------------------------------
        # cthetas = np.array(cthetas)[::-1,:] # INVERT DATA ----------------------------------------------------------
        #
        # print("ctheta: {}".format(cthetas0.shape))
        # print("Rs: {}".format(Rs0.shape))
        #
        # val = np.reshape(np.array(val), (len(cthetas0), len(Rs0)))
        # ctheta = np.reshape(np.array(cthetas), (len(cthetas0), len(Rs0)))
        # Rs = np.reshape(np.array(Rs), (len(cthetas0), len(Rs0)))
        # print("Val: {}".format(val.shape))
        #
        # print("ctheta={}, {}".format(cthetas[0] * 180 / np.pi, cthetas[-1] * 180 / np.pi))
        # print("Rs={}, {}".format(Rs0[0],Rs0[-1]))

        # ---------------------------------

        angle_ticks = range(0,100,10)
        angle_ticks_rads = [a * math.pi / 180.0 for a in angle_ticks]#[::-1] # INVERT TICKS -------------------
        angle_ticks_rads_plus_offset = [a for a in angle_ticks_rads]
        # angle_ticks_rads_plus_offset = angle_ticks_rads_plus_offset[::-1] # Polar Angle
        angle_ticks_for_plot = []
        for i in range(len(angle_ticks)):
            angle_ticks_for_plot.append((angle_ticks_rads_plus_offset[i],r"$"+str(angle_ticks[i])+"$"))

        print("Angle ticks prepared")

        # lRs = np.log10(Rs/Rs0[0])
        # lRs0 = np.log10(Rs0/Rs0[0])

        rs = rs/rs.min()
        lrs = np.log10(rs)
        radius_ticks = range(int(lrs[0, 0]), int(lrs[0, -1]), 1)
        radius_ticks_for_plot=[]
        for i in range(len(radius_ticks)):
            radius_ticks_for_plot.append((radius_ticks[i], r"$"+str(radius_ticks[i])+"$"))

        print("Radial ticks prepared")

        # ---------------------------------------

        scale = 1.5
        aspect = 1.20
        height = 3.0
        fig = plt.figure(1, figsize=(height*aspect*scale, height*scale))
        fig.subplots_adjust(wspace=0.1, left=0.05, right=0.95, top=0.94)
        fig.subplots_adjust()

        ax2, aux_ax2 = setup_arc_radial_axes(fig, 111, angle_ticks_for_plot, radius_ticks_for_plot, 1, radius_ticks[-1])

        # r, theta = np.meshgrid(lRs,cthetas)
        values=val

        levels = ticker.LogLocator(base=10,numticks=100).tick_values(val.min(), val.max())
        # levels = MaxNLocator(nbins=20).tick_values(val.min(), val.max())
        cmap = plt.get_cmap('inferno') # seismic
        norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
        # norm=colors.TwoSlopeNorm(vmin=val.min(), vcenter=0.9, vmax=val.max())
        # [:,100:950]
        im = aux_ax2.pcolormesh(cthetas[1:4,:], lrs[1:4,:], values[1:4,:], cmap=cmap, norm=norm, alpha=0.7)

        cbar = plt.colorbar(im, orientation='vertical')
        cbar.ax.set_ylabel(r'$\log_{10}( \rho/\rho_{\rm ISM} )$', fontsize = 12)

        # ax2.axis["bottom"].axes.set_xlabel('Angle [deg]', fontsize=20, weight="bold")

        plt.suptitle(' Jet layer dynamics ', fontsize = 14, weight="bold")
        plt.legend(loc=3,prop={'size':14})
        # plt.xlabel('Angle [deg]', fontsize=20, weight="bold", rotation=30)
        plt.ylabel(r'Radius $[R/R_0]$', fontsize=14, weight="bold")

        # plt.show()
        plt.savefig('plot.png', dpi=256)
        plt.show()
        plt.close()

        fig, ax = plt.subplots(figsize=(4,3),subplot_kw={'projection' :'polar'},ncols=1, nrows=1)
        # ax = plt.subplot(111, polar=True)
        levels = MaxNLocator(nbins=15).tick_values(val.min(), val.max())
        levels = ticker.LogLocator(base=10,numticks=100).tick_values(val.min(), val.max())

        cmap = plt.get_cmap('inferno')
        norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
        im = ax.pcolormesh(cthetas, lrs, val, cmap=cmap, norm=norm)
        fig.colorbar(im, ax=ax)
        # ax.set_title('pcolormesh with levels')
        print(ax.get_rmin(), ax.get_rmax())
        # ax.set_rmax(10)
        # ax.set_rmin(20)
        print(ax.get_rmin(), ax.get_rmax())

        max_theta = 90
        # ax.set_theta_zero_location("U")
        ax.set_theta_offset(np.pi/2)



        ax.set_thetamax(90)
        # ax.set_thetamin(0)
        # ax.set_rlim(10,20)

        # ax.set_rlim(Rs.min(), Rs.max())
        # ax.set_rscale("log")
        # ax.set_rscale('log')
        #
        # ticklabels = ax.get_yticklabels()
        # labels = range(80, 0, -10)
        # for i in range(0, len(labels)):
        #     ticklabels[i] = str(labels[i])
        # ax.set_yticklabels(ticklabels)

        plt.show()

    def plot_timestep_ejecta(self):
        dfile_ej = h5py.File(curdir+"ejecta_dynamics_shells_layers.h5", "r")
        print(dfile_ej.keys())
        print(dfile_ej[list(dfile_ej.keys())[0]].keys())
        print(dfile_ej.attrs.keys())


        ishell = 20
        all_shells = int(dfile_ej.attrs["nshells"])
        all_layers = int(dfile_ej.attrs["nlayers"])
        print("all_shells={} all_layers={}".format(all_shells, all_shells))

        r_grid = np.array( dfile_ej["shell={} layer={}".format(0,0)]["R"] )

        n = all_layers * all_layers

        cthetas = np.zeros((all_layers,all_shells))
        rs = np.zeros((all_layers,all_shells))
        val = np.zeros((all_layers,all_shells))

        ir = 100


        for ilayer in range(int(all_layers)):
            for ishell in range(int(all_shells)):
                key = "shell={} layer={}".format(ishell, ilayer)
                cthetas[ilayer,ishell] = np.array(dfile_ej[key]["ctheta"])[ir]
                rs[ilayer,ishell] = np.array(dfile_ej[key]["R"])[ir]
                i_val = np.array( dfile_ej[key]["beta"] ) * np.array( dfile_ej[key]["Gamma"] )
                val[ilayer,ishell] = i_val[ir]

                # cthetas = np.vstack(( cthetas, np.array(dfile_ej[key]["ctheta0"]) ))
                # rs = np.vstack(( rs, np.array(dfile_ej[key]["R"]) ))
                # i_val = np.array( dfile_ej[key]["beta"] ) * np.array( dfile_ej[key]["Gamma"] )
                # val = np.vstack(( val, i_val ))

        # lrs = np.log10(rs)


        # val.append( np.array( dfile_ej[key]["GammaRho"] ) * np.array( dfile_ej[key]["GammaRho"] ) )
        # val.append( np.log10( np.array( dfile_ej[key]["rho"] )/np.array( dfile_ej[key]["rho"] )[-1] ))#* ( np.array( dfile_ej[key]["rho"][-1] ) ) ) )
        # val.append( np.log10(np.array( dfile_ej[key]["rho2"] ) / ( np.array( dfile_ej[key]["rho2"][0] ) ) ) )
        # plt.semilogx(Rs, val[-1])
        # plt.show()
        cthetas = np.array(cthetas)#[::-1] # INVERT DATA ----------------------------------------------------------

        print("ctheta: {}".format(cthetas.shape))
        print("Rs: {}".format(rs.shape))


        print("Jet Data loaded")

        # cthetas = np.array(cthetas)#[::-1] # INVERT DATA ----------------------------------------------------------
        # cthetas0 = np.array(cthetas0)#[::-1] # INVERT DATA ----------------------------------------------------------
        # cthetas = np.array(cthetas)[::-1,:] # INVERT DATA ----------------------------------------------------------
        #
        # print("ctheta: {}".format(cthetas0.shape))
        # print("Rs: {}".format(Rs0.shape))
        #
        # val = np.reshape(np.array(val), (len(cthetas0), len(Rs0)))
        # ctheta = np.reshape(np.array(cthetas), (len(cthetas0), len(Rs0)))
        # Rs = np.reshape(np.array(Rs), (len(cthetas0), len(Rs0)))
        # print("Val: {}".format(val.shape))
        #
        # print("ctheta={}, {}".format(cthetas[0] * 180 / np.pi, cthetas[-1] * 180 / np.pi))
        # print("Rs={}, {}".format(Rs0[0],Rs0[-1]))

        # ---------------------------------

        angle_ticks = range(0,100,10)
        angle_ticks_rads = [a * math.pi / 180.0 for a in angle_ticks]#[::-1] # INVERT TICKS -------------------
        angle_ticks_rads_plus_offset = [a for a in angle_ticks_rads]
        # angle_ticks_rads_plus_offset = angle_ticks_rads_plus_offset[::-1] # Polar Angle
        angle_ticks_for_plot = []
        for i in range(len(angle_ticks)):
            angle_ticks_for_plot.append((angle_ticks_rads_plus_offset[i],r"$"+str(angle_ticks[i])+"$"))

        print("Angle ticks prepared")

        # lRs = np.log10(Rs/Rs0[0])
        # lRs0 = np.log10(Rs0/Rs0[0])

        # rs = rs/rs.min()
        lrs = np.log2(rs)- 0.9*np.log2(rs)[0,0]
        # lrs /= lrs.min()
        radius_ticks = range(int(lrs[0, 0]), int(lrs[0, -1]), 1)
        radius_ticks_for_plot=[]
        for i in range(len(radius_ticks)):
            radius_ticks_for_plot.append((radius_ticks[i], r"$"+str(radius_ticks[i])+"$"))

        print("Radial ticks prepared")

        # ---------------------------------------

        scale = 1.5
        aspect = 1.20
        height = 3.0
        fig = plt.figure(1, figsize=(height*aspect*scale, height*scale))
        fig.subplots_adjust(wspace=0.1, left=0.05, right=0.95, top=0.94)
        fig.subplots_adjust()

        ax2, aux_ax2 = self.setup_arc_radial_axes(fig, 111, angle_ticks_for_plot, radius_ticks_for_plot, 1, radius_ticks[-1])

        # r, theta = np.meshgrid(lRs,cthetas)
        values=val

        # levels = ticker.LogLocator(base=10,numticks=1000).tick_values(1e-4, 0.8)
        levels = MaxNLocator(nbins=20).tick_values(val.min(), val.max())
        cmap = plt.get_cmap('viridis') # seismic
        norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
        norm=colors.TwoSlopeNorm(vmin=val.min(), vcenter=0.1, vmax=val.max())

        # im = aux_ax2.pcolormesh(cthetas, lrs, values, cmap=cmap, norm=norm, alpha=0.7)
        im = aux_ax2.pcolormesh(cthetas, lrs, values, cmap=cmap,  alpha=0.7, norm=LogNorm(vmin=1e-3,vmax=5e-1), rasterized=True)

        cbar = plt.colorbar(im, orientation='vertical')
        cbar.ax.set_ylabel(r'$\log_{10}( \rho/\rho_{\rm ISM} )$', fontsize = 12)

        # ax2.axis["bottom"].axes.set_xlabel('Angle [deg]', fontsize=20, weight="bold")

        plt.suptitle(' Jet layer dynamics ', fontsize = 14, weight="bold")
        plt.legend(loc=3,prop={'size':14})
        # plt.xlabel('Angle [deg]', fontsize=20, weight="bold", rotation=30)
        plt.ylabel(r'Radius $[R/R_0]$', fontsize=14, weight="bold")

        # plt.show()
        plt.savefig('plot.png', dpi=256)
        plt.show()
        plt.close()


        fig, ax = plt.subplots(figsize=(4,3),subplot_kw={'projection' :'polar'},ncols=1, nrows=1)
        # ax = plt.subplot(111, polar=True)
        # levels = MaxNLocator(nbins=25).tick_values(val.min(), val.max())
        levels = ticker.LogLocator(base=2,numticks=100).tick_values(val.min(), val.max())

        cmap = plt.get_cmap('inferno')
        norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
        im = ax.pcolormesh(cthetas, lrs, val, cmap=cmap, norm=norm)
        fig.colorbar(im, ax=ax)
        # ax.set_title('pcolormesh with levels')
        print(ax.get_rmin(), ax.get_rmax())
        # ax.set_rmax(10)
        # ax.set_rmin(20)
        print(ax.get_rmin(), ax.get_rmax())

        max_theta = 90
        # ax.set_theta_zero_location("U")
        # ax.set_theta_offset(np.pi/2)



        ax.set_thetamax(90)
        # ax.set_thetamin(0)
        # ax.set_rlim(10,20)

        # ax.set_rlim(Rs.min(), Rs.max())
        # ax.set_rscale("log")
        # ax.set_rscale('log')
        #
        # ticklabels = ax.get_yticklabels()
        # labels = range(80, 0, -10)
        # for i in range(0, len(labels)):
        #     ticklabels[i] = str(labels[i])
        # ax.set_yticklabels(ticklabels)

        plt.show()


q = quarter_plots()
# q.plot_full_evol_jet()
q.plot_timestep_ejecta()