# import PyBlastAfterglowMag
import copy
import shutil
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import os
import time as timer

from matplotlib.colors import Normalize, LogNorm
from matplotlib import cm, rc, rcParams
import matplotlib.colors as colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1 import ImageGrid
from mpl_toolkits.axisartist.grid_finder import MaxNLocator
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
import matplotlib.gridspec as gridspec

import package.src.PyBlastAfterglowMag as PBA
from package.src.PyBlastAfterglowMag.utils import cgs

# temprary dir to save PyBlastAfterglow output into
working_dir = os.getcwd() + '/working_dirs/'
fig_dir = os.getcwd() + '/figs/'

class GRB170817A(object):

    def __init__(self):

        # load data
        times, freqs, data, errs = self._load_data()

        # unique freqs in the data
        ufreqs = np.unique(freqs)

        # concatenate the data for unique frequencies
        # data_ord = np.array([])
        # err_ord = np.array([])
        #
        # for frequency in freqs:
        #     # print(frequency)
        #     data_ord = np.concatenate([data_ord, data[freq == frequency]])
        #     err_ord = np.concatenate([err_ord, err[freq == frequency]])
        # assert np.shape(err_ord) == np.shape(data_ord)

        # print("--- Observations for GRB170817 ---")
        # print("Total amount of time: {}".format(len(times)))
        # print("Total amount of data: {}".format(len(data)))
        # print("Total amount of errs: {}".format(len(errs)))
        # print("Unique frequencies   ({})".format(len(ufreqs)))
        # print(ufreqs)
        # print("Data per frequency:")
        # for ifreq, freq in enumerate(ufreqs):
        #     print("freq={:.2e} N={:d}".format(freq, len(data[freq==freqs])))

        self.times = times
        self.data = data
        self.errs = errs
        self.freqs = freqs
        self.ufreqs = ufreqs

    def __call__(self, freq):
        return self.get(freq)

    def _load_data(self):
        """
            Chandra21 :: https://ui.adsabs.harvard.edu/abs/2020GCN.29055....1H/abstract
            Hajela:2020
        :return:
        """

        # Updated optical and X-ray data from Fong ea 2019, and Hajela ea 2019
        time = np.array([9.2, 14.9, 16.4, 17.4, 18.3, 18.7,
                         19.4, 21.4, 22.4, 23.4, 24.2, 31.3,
                         35.3, 39.2, 46.3, 53.3, 54.3, 57.2,
                         65.9, 66.6, 67.2, 72.2, 75.5, 75.5,
                         77.6, 79.2, 80.1, 92.4, 93.1, 93.1,
                         93.2, 97.1, 107., 107., 109., 109.,
                         111., 112., 115., 115., 115., 125.,
                         125., 126., 133., 137., 149., 150.,
                         152., 158., 161., 163., 163., 163.,
                         163., 165., 167., 170., 172., 183.,
                         197., 197., 207., 209., 216., 217.,
                         217., 217., 217., 218., 218., 222.,
                         229., 252., 257., 259., 261., 267.,
                         267., 273., 273., 289., 289., 294.,
                         297., 298., 320., 324., 328., 357.,
                         359., 362., 380., 489., 545., 580.,
                         581., 741., 767., 938., 1211 # Chandra21
                         ])  # time in days

        # flux in mJy
        data = np.array(
            [5.66e-04, 6.58e-04, 1.87e+01, 1.51e+01, 1.45e+01, 1.54e+01,
             1.59e+01, 1.36e+01, 2.25e+01, 2.00e+01, 2.56e+01, 3.40e+01,
             4.40e+01, 2.28e+01, 4.40e+01, 3.20e+01, 4.80e+01, 6.10e+01,
             1.48e+02, 9.80e+01, 4.26e+01, 5.80e+01, 3.59e+01, 3.96e+01,
             7.70e+01, 4.50e+01, 4.17e+01, 3.17e+01, 9.80e+01, 7.00e+01,
             2.60e+01, 1.99e+02, 1.27e+02, 5.32e+01, 2.96e-03, 1.09e-01,
             1.11e-01, 6.29e+01, 9.62e+01, 5.12e+01, 4.12e+01, 5.82e+01,
             1.28e+02, 2.21e+02, 3.37e-03, 8.40e-02, 6.06e+01, 9.00e+01,
             1.84e+02, 3.03e-03, 2.66e-03, 9.73e+01, 6.73e+01, 4.74e+01,
             3.96e+01, 9.10e-02, 5.79e+01, 1.13e-01, 8.50e-02, 2.11e+02,
             7.59e+01, 8.93e+01, 4.20e+01, 8.20e-02, 3.63e+01, 6.05e+01,
             4.17e+01, 3.26e+01, 2.47e+01, 6.47e+01, 6.30e-02, 3.97e+01,
             4.80e+01, 7.13e+01, 4.32e+01, 1.55e-03, 6.26e+01, 2.50e+01,
             4.03e+01, 3.48e+01, 2.72e+01, 3.63e+01, 2.70e+01, 3.12e+01,
             4.40e-02, 2.34e+01, 2.31e+01, 4.72e+01, 3.40e-02, 9.70e-04,
             1.55e+01, 2.70e-02, 3.79e+01, 1.48e+01, 5.90e+00, 1.80e+01,
             3.54e-04, 2.68e-04, 4.90e+00, 1.95e-04, 3.46e-04 # Chandra21 3.46(+1.06 -1.31) e-15 erg/cm2/s
             ])

        # frequency of obs.
        freq = np.array([2.41e+17, 2.41e+17, 3.00e+09, 3.00e+09, 3.00e+09, 7.25e+09,
                         6.20e+09, 6.20e+09, 3.00e+09, 6.00e+09, 3.00e+09, 3.00e+09,
                         1.50e+09, 6.00e+09, 3.00e+09, 6.00e+09, 3.00e+09, 3.00e+09,
                         6.70e+08, 1.30e+09, 6.00e+09, 4.50e+09, 7.35e+09, 7.35e+09,
                         1.40e+09, 4.50e+09, 6.00e+09, 7.25e+09, 1.50e+09, 3.00e+09,
                         1.50e+10, 6.70e+08, 1.30e+09, 1.30e+09, 2.41e+17, 3.80e+14,
                         5.06e+14, 6.00e+09, 3.00e+09, 1.00e+10, 1.50e+10, 7.25e+09,
                         1.30e+09, 6.70e+08, 2.41e+17, 5.06e+14, 7.25e+09, 5.10e+09,
                         1.30e+09, 2.41e+17, 2.41e+17, 3.00e+09, 6.00e+09, 1.00e+10,
                         1.50e+10, 5.06e+14, 7.25e+09, 3.80e+14, 5.06e+14, 6.50e+08,
                         3.00e+09, 1.30e+09, 5.00e+09, 5.06e+14, 1.00e+10, 3.00e+09,
                         6.00e+09, 1.00e+10, 1.50e+10, 3.00e+09, 5.06e+14, 7.25e+09,
                         4.50e+09, 1.30e+09, 3.00e+09, 2.41e+17, 1.30e+09, 7.25e+09,
                         3.00e+09, 3.00e+09, 6.00e+09, 3.00e+09, 6.00e+09, 3.00e+09,
                         5.06e+14, 7.25e+09, 7.25e+09, 1.30e+09, 5.06e+14, 2.41e+17,
                         7.25e+09, 5.06e+14, 1.30e+09, 3.00e+09, 6.00e+09, 7.25e+09,
                         2.41e+17, 2.41e+17, 3.00e+09, 2.41e+17, 2.41e+17 # Chandra21
                         ])

        # error on flux
        err = np.array([1.70e-04, 1.30e-04, 6.30e+00, 3.90e+00, 3.70e+00, 4.80e+00,
                        5.50e+00, 2.90e+00, 3.40e+00, 3.10e+00, 2.90e+00, 3.60e+00,
                        1.00e+01, 2.60e+00, 4.00e+00, 4.00e+00, 6.00e+00, 9.00e+00,
                        2.20e+01, 2.00e+01, 4.10e+00, 5.00e+00, 4.30e+00, 7.00e+00,
                        1.90e+01, 7.00e+00, 4.70e+00, 4.30e+00, 1.40e+01, 5.70e+00,
                        4.40e+00, 1.60e+01, 1.80e+01, 4.50e+00, 2.60e-04, 1.70e-02,
                        1.90e-02, 3.20e+00, 8.00e+00, 3.40e+00, 1.90e+00, 5.00e+00,
                        2.10e+01, 1.90e+01, 4.00e-04, 1.80e-02, 4.30e+00, 3.00e+01,
                        1.90e+01, 2.60e-04, 2.70e-04, 1.13e+01, 4.10e+00, 3.60e+00,
                        2.00e+00, 1.60e-02, 6.90e+00, 1.90e-02, 1.70e-02, 3.40e+01,
                        5.20e+00, 1.39e+01, 1.20e+01, 2.00e-02, 3.60e+00, 7.50e+00,
                        7.50e+00, 4.00e+00, 3.10e+00, 2.70e+00, 1.80e-02, 7.20e+00,
                        6.00e+00, 6.70e+00, 5.80e+00, 1.90e-04, 7.00e+00, 4.10e+00,
                        2.70e+00, 4.90e+00, 2.10e+00, 3.90e+00, 2.80e+00, 3.60e+00,
                        1.40e-02, 4.20e+00, 4.00e+00, 1.28e+01, 1.10e-02, 1.90e-04,
                        5.00e+00, 7.00e-03, 1.18e+01, 2.90e+00, 1.90e+00, 4.20e+00,
                        9.00e-05, 9.00e-05, 1.80e+00, 7.00e-05, 1.3e-04 # Chandra21 3.46(+1.06 -1.31) e-15 erg/cm2/s
                        ])

        assert np.shape(time) == np.shape(freq)
        assert np.shape(data) == np.shape(err)

        data *= 1.e-29 # muJy -> ergs
        err *= 1.e-29 # muJy -> ergs

        return(time, freq, data, err) # [days, Hz, ergs, ergs]

    def get(self, freq):
        if not freq in self.ufreqs:
            raise NameError("freq: {} is not in data. Available are: {}".format(freq, self.ufreqs))
        mask = freq == self.freqs
        if len(mask) < 1:
            raise ValueError("Failed to get data for freq:{}".format(freq))
        return (self.times[mask], self.data[mask], self.errs[mask], np.zeros_like(self.times[mask]))

    def get_chandra(self):
        # Time since Merger (days) Fluxd (mJy) Err (mJy)
        data = np.array(
            [[2.19000, 1.40000e-07,  0.00000,     0.00000],
             [9.19679,  2.10000e-07,  8.72000e-08, -9.02000e-08],
             [15.3900,  6.44000e-07,  1.72000e-07, -1.20000e-07],
             [108.386,  2.21000e-06,  2.14000e-07, -2.02000e-07],
             [157.755,  2.41000e-06,  2.62000e-07, -2.11000e-07],
             [259.665,  1.07000e-06,  1.64000e-07, -1.55000e-07],
             [358.609,  9.12000e-07,  2.10000e-07, -1.73000e-07],
             [581.818,  2.14000e-07,  1.09000e-07, -7.72000e-08],
             [741.478,  1.26000e-07,  6.64000e-08, -5.39000e-08],
             [939.310,  1.54977e-07,  8.46022e-08, -6.26098e-08],
             [1234.11,  2.16000e-07,  5.45000e-08, -7.95000e-08]
             ])
        # return [day] [erg] [erg] [erg]
        errs = np.abs(data[:, 2:][:,::-1]).T
        # errs[:, 0] += data[:, 1]
        # errs[:, 1] += data[:, 1]
        uplims = np.zeros_like(data[:, 0], dtype=bool)
        return (data[:, 0], data[:, 1] / 1e3 / 1e23, errs / 1e3 / 1e23, uplims)##,  data[:, 3] * 1e3 * 1e23)

    def get_vla_3ggz(self):

        data = np.array(
            [
                # Time since Merger (days) Fluxd (mJy) Err (mJy)
                [3.34000,     0.032000002,       0.0000000],
                [16.4200,     0.018700000,    0.0063000000],
                [17.3900,     0.015100000,    0.0038999999],
                [18.3300,     0.014500000,    0.0037000000],
                [22.3600,     0.022500001,    0.0034000000],
                [24.2600,     0.025599999,    0.0029000000],
                [31.3200,     0.034000002,    0.0035999999],
                [46.2600,     0.044000000,    0.0040000002],
                [54.2700,     0.048000000,    0.0060000001],
                [57.2200,     0.061000001,    0.0089999996],
                [93.1300,     0.070000000,    0.0057000001],
                ##[196.790     0.081793795,    0.0081793800],
                [196.790,     0.073213301,    0.0066756521],
                [115.200,     0.075690837,     0.019037671],
                [115.200,      0.10307770,     0.011835645],
                ##[163.070,     0.096567894,     0.020509181],
                [163.070,     0.098128255,     0.018721838],
                [216.910,     0.068999998,     0.015000000],
                [220.000,     0.064700000,    0.0027000001],
                [256.760,     0.055000000,     0.012300000],
                [267.000,     0.040300000,    0.0027000001],
                [272.670,     0.043900002,     0.010500000],
                [288.610,     0.046399999,     0.011400000],
                [294.000,     0.031199999,    0.0035999999],
                [489.000,     0.014800000,    0.0029000000],
                [767.000,    0.0049000001,    0.0018000000],
                # [1221.50,    0.0054000001,      0.00000000],
                [1243.000,   2.86*1e-3,          0.99*1e-3], # https://arxiv.org/pdf/2103.04821.pdf
                [4.6*cgs.year/cgs.day,   4.5*1e-3,          1.1*1e-3] # https://arxiv.org/pdf/2205.14788.pdf
            ]
        )
        uplims = np.zeros_like(data[:, 0], dtype=bool)
        # data[:, 2] = 0.00200000
        uplims[-1] = True
        return (data[:, 0], data[:, 1] / 1e3 / 1e23, data[:, 2] / 1e3 / 1e23, uplims)
def plot_grb170817a_lc(ax, freq, color):
    o_obs = GRB170817A()
    if (freq == 3e9):
        x_obs, y_obs, yerr_obs, uplims = o_obs.get_vla_3ggz()
        lbl = r"GRB170817A VLA 3\,GHz"
    elif (freq == 2.41e+17):
        x_obs, y_obs, yerr_obs, uplims = o_obs.get_chandra()
        lbl = r"GRB170817A Chandra 1\,keV"
    elif freq in o_obs.freqs:
        x_obs, y_obs, yerr_obs, uplims = o_obs.get(freq)
        lbl = "GRB170817A"
    else:
        return
    # plotdic = {}
    plotdic = {
        "x_vals": x_obs, "y_vals": y_obs * 1e23 * 1e3, "yerr_vals": yerr_obs * 1e23 * 1e3,
        "plot": {'ecolor': 'gray', 'mec': color,  "uplims":uplims,
                 'marker': 'o', 'markersize': 8,
                 "capsize": 2, "linestyle": "None", 'mfc': "white",
                 "label": lbl, "zorder": 100},
        "legend": {"fancybox": False, "loc": 'lower left',
                   # "bbox_to_anchor":(1.0, 0.0),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
                   "shadow": "False", "ncol": 1, "fontsize": 12,
                   "framealpha": 0., "borderaxespad": 0., "frameon": False}
    }
    tmp = copy.deepcopy(plotdic)
    obs_times, obs_fluxes, obs_errs = tmp["x_vals"], tmp["y_vals"], tmp["yerr_vals"]
    tmptmp = copy.deepcopy(plotdic["plot"])
    if "label" in tmptmp.keys():
        lbl = tmptmp["label"]
        del tmptmp["label"]
    _l = ax.errorbar(obs_times, obs_fluxes, yerr=obs_errs, **tmptmp)
    # tmptmp["marker"] = "v"
    # for ulim in uplims:
    #     ax.plot(ulim)
    if "label" in plotdic["plot"].keys():
        obs_legend = copy.deepcopy(tmp["legend"])
        obs_legend["loc"] = "lower right"
        leg2 = ax.legend([_l], [lbl], **obs_legend)
        ax.add_artist(leg2)

def get_vla_and_hst_data():
    """
    Make a dataset of observed flux centroid positions for 8, 75, 206 and 230 days after merger.
    Here both, VLA and HST, HSA data is given.
    Radio data is expected to be given at 4.5e9 Hz; while optical data at 8 days is presumably at 4.6e14 Hz
    Resulting dataset has the following structure

    # time[days] freq[Hz] xc[mas] yc[mas] xc_err_left xc_err_right yc_err_low yc_err_up
    # ...        ...      ...     ...     ...         ...          ...        ...
    :return:
    """
    # extracted from https://www.nature.com/articles/s41586-022-05145-7 Fig.1
    data = np.array([
        [6.92063492063492,	1.68247422680412],
        [6.91609977324263,	1.28041237113402],
        [6.91609977324263,	2.07835051546392],
        [7.08390022675737,	1.68247422680412],
        [6.75283446712018,	1.68247422680412],
        [5.83673469387755,	1.52783505154639],
        [5.83673469387755,	1.29278350515464],
        [5.84126984126984,	1.77525773195876],
        [6.04988662131519,	1.5340206185567],
        [5.63265306122449,	1.54020618556701],
        [3.99546485260771,	1.50309278350515],
        [4.0,	1.10103092783505],
        [4.0,	1.90515463917526],
        [4.11337868480726,	1.50309278350515],
        [3.87755102040816,	1.50309278350515],
        [1.37868480725624,	1.58350515463918],
        [1.38321995464853,	1.40412371134021],
        [1.38321995464853,	1.76907216494845],
        [1.70068027210884,	1.58969072164948],
        [1.06575963718821,	1.58969072164948],
    ])
    times = np.array([8., 75., 206., 230.])
    freqs = np.array([4.6e14, 4.5e9, 4.5e9, 4.5e9])

    xs = data[0::5,0]
    left = data[3::5,0]
    right = data[4::5,0]

    ys = data[0::5,1]
    low = data[1::5,1]
    up = data[2::5,1]

    print('x', xs)
    print('left', left)
    print('right', right)

    print('ys',ys)
    print('lower',low)
    print('upper',up)
    print('---------------------')

    # noramize
    low-=ys[2]; up-=ys[2]; ys-=ys[2]# set so that data at 75 days is at [0,0]
    left-=xs[2]; right-=xs[2];xs-=xs[2] # set so that data at 75 days is at [0,0]

    print('x', xs)
    print('left', left)
    print('right', right)

    print('ys',ys)
    print('lower',low)
    print('upper',up)

    left = left - xs
    right = xs - right

    low = ys - low
    up = up - ys

    res = np.column_stack(
        (times, freqs, xs, ys, left, right, low, up)
    )
    return (res)

    # fig,ax = plt.subplots(ncols=1,nrows=1,figsize=(7,3),layout='constrained')
    #
    # ax.errorbar(xs,ys,
    #             xerr=np.column_stack((left,right)).T,
    #             yerr=np.column_stack((low,up)).T,
    #             **{'ecolor': 'gray', 'mec': 'gray',
    #                'marker': 'o', 'markersize': 8,
    #                "capsize": 2, "linestyle": "None", 'mfc': "white",
    #                "label": "VLBI GRB170817a", "zorder": 100})
    #
    # ax.set_xlim(4,-4)
    # ax.set_ylim(-1.5,1.5)
    # ax.grid(ls=':')
    # ax.set_aspect(aspect=1/2)
    # plt.show()
def compare_lcs_and_skymap_props(pp:dict, plot:dict,working_dir:str,run:bool=True,process_skymaps:bool=True):
    start_time = timer.perf_counter()
    pba = PBA.wrappers.run_grb(
        working_dir=working_dir,
        P=pp,
        run=run,
        process_skymaps=process_skymaps,
        loglevel="err",#"err"
    )
    finish_time = timer.perf_counter()
    print(f"PyBlastAfterglow finished in {finish_time-start_time:.2f} seconds. ({(finish_time-start_time)/60.:.2f} min.) ")

    # ---------

    freqs=(3e9, 3.80e+14, 2.41e+17)
    colors=('blue','green','orange')
    offsets=(1,1,1)
    times = pba.GRB.get_lc_times()
    fig, axes = plt.subplots(nrows=4,ncols=1,figsize=(5,7),layout='constrained')

    i=0
    for freq,color,offset in zip(freqs,colors,offsets):
        lc = pba.GRB.get_lc(freq=freq)
        axes[i].plot(times/cgs.day, lc,color=color,ls='-', label="${}$".format(PBA.utils.latex_float(freq)))
        plot_grb170817a_lc(ax=axes[i], freq=freq, color=color)
    # axes[i].set_xscale('log')
    axes[i].set_yscale('log')
    # axes[i].set_xlabel("time [day]")
    # axes[i].grid(ls=':')
    axes[i].set_ylabel("Flux density [mJy]")
    # ax.legend(fancybox=False, loc='lower right', columnspacing=0.8,
    #           # bbox_to_anchor=(0.5, 0.5),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
    #           shadow=False, ncol=1, fontsize=12, labelcolor='black',
    #           framealpha=0.0, borderaxespad=0.)
    # axes[i].set_xlim(1,2e3)
    axes[i].set_ylim(1e-10,4e-1)
    # plt.show()

    # ----------------

    freq = 4.5e9
    times = pba.GRB.get_skymap_times()
    skymaps = [pba.GRB.get_skymap(time=t, freq=freq) for t in times]
    xcs = [skymap.xc for skymap in skymaps]
    delta_x = [abs(skymap.x2-skymap.x1) for skymap in skymaps]
    fluxes = [skymap.flux for skymap in skymaps]
    beta_app = np.diff(xcs)*cgs.rad2mas / times[1:] # / cgs.c

    # fig,axes = plt.subplots(ncols=1,nrows=2,sharex='all',layout='constrained',figsize=(6,4))
    i+=1
    axes[i].plot(times/cgs.day, xcs, color='black', ls='-', marker='.')
    axes[i].set_ylabel("x_c [mas]")
    axes[i].set_ylim(0,7)
    # axes[i].grid(ls=':')
    # axes[i].set_xlim(1,1e4)
    # axes[i].set_xscale('log')

    axes[i].errorbar([74.5719, 205.4533, 229.2961],
                     [2.4054, 4.089, 5.0562],
                     xerr=np.array([2, 3, 2, 3, 2, 3]).reshape(-1,3),
                     yerr=np.array([2.4054-2.0239, 2.7802-2.4054, 4.089-3.6593, 4.4974-4.089, 5.0562-4.6610, 5.4719-5.0562]).reshape(-1,3),
                     **{'ecolor': 'gray', 'mec': 'gray',
                        'marker': 'o', 'markersize': 8,
                        "capsize": 2, "linestyle": "None", 'mfc': "white",
                        "label": "VLBI GRB170817a", "zorder": 100})
    sublum_dxcs = np.zeros_like(times)
    for j in range(len(times)):
        beta_app_sub = 1.
        sublum_dxcs[j] = beta_app_sub / cgs.rad2mas * times[j]
    sublum_xcs = np.cumsum(sublum_dxcs)#[:-1] + sublum_dxcs[1:]
    axes[i].plot(times[1:]/cgs.day, sublum_xcs[1:], color='gray', ls=':')
    axes[i].fill_between(times[1:]/cgs.day, sublum_xcs[1:], where=sublum_xcs[1:]>=0, interpolate=True, color='gray')

    i+=1
    axes[i].plot(times[1:]/cgs.day, beta_app, color='black', ls='-', marker='.')
    axes[i].set_ylabel("beta_app [c]")
    axes[i].set_ylim(-1,8)
    # axes[i].grid(ls=':')
    # axes[i].set_xlim(1,1e4)
    # axes[i].set_xscale('log')
    #4.1±0.5 # Mooley et al 2018
    axes[i].errorbar([(75-8)/2, (206-8)/2, (230-8)/2],
                     [7.6, 4.7, 5.2],
                     xerr=np.array([(75-8)/2-8, 75-(75-8)/2, (206-8)/2-8, 206-(206-8)/2, (230-8)/2-8, 230-(230-8)/2]).reshape(-1,3),
                     yerr=np.array([1.3,1.3, 0.6,0.6, 0.5,0.5]).reshape(-1,3),
                     **{'ecolor': 'gray', 'mec': 'gray',
                        'marker': 'o', 'markersize': 8,
                        "capsize": 2, "linestyle": "None", 'mfc': "white",
                        "label": "VLBI GRB170817a", "zorder": 100})

    axes[i].legend(**{"fancybox": False, "loc": 'lower left',
                      # "bbox_to_anchor":(1.0, 0.0),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
                      "shadow": "False", "ncol": 1, "fontsize": 12,
                      "framealpha": 0., "borderaxespad": 0., "frameon": False})
    for ax in axes:
        ax.set_xlim(1,2e3)
        ax.grid(ls=':')
        ax.set_xscale('log')
        ax.tick_params(axis='both', which='both', labelleft=True,
                       labelright=False, tick1On=True, tick2On=True,
                       labelsize=12, direction='in',
                       bottom=True, top=True, left=True, right=True)
    axes[0].set_xticklabels([])
    axes[1].set_xticklabels([])
    axes[2].set_xlabel("time [day]",fontsize=12)
    # plt.show()


    # -----------
    ax =axes[3]

    obs = get_vla_and_hst_data()
    obs_times = obs[:,0]
    obs_freqs = obs[:,1]


    # fig,ax = plt.subplots(ncols=1,nrows=1,figsize=(5,2),layout='constrained')

    ref_skymap = pba.GRB.get_skymap(time=obs_times[1]*cgs.day, freq=obs_freqs[1])
    offset = ref_skymap.xc
    for t,nu in zip(obs_times,obs_freqs):
        skymap=pba.GRB.get_skymap(time=t*cgs.day, freq=nu)
        int_zz = skymap.im_hist; int_x=skymap.grid_x-offset; int_y=skymap.grid_y
        xc = skymap.xc-offset
        yc = skymap.yc
        for i,x in enumerate(int_x):
            for j,y in enumerate(int_y):
                if (abs(x - xc) > 1. or abs(y - yc) > 1.):
                    int_zz[i,j] = np.nan

        norm = Normalize(int_zz[np.isfinite(int_zz)].max()*1e-2, int_zz[np.isfinite(int_zz)].max())
        cmap = cm.get_cmap('Greys')
        cmap.set_under('white')
        cmap.set_over('black')
        im = ax.pcolormesh(int_x, int_y, int_zz.T, norm=norm, cmap=cmap, alpha=0.7)
        # ax.set_facecolor(cmap(float(0.)))
        ax.plot(xc, yc, **dict(marker='s',color='black',markersize=12,fillstyle='none'))
        # ax.errorbar([skymap.xc, skymap.xc], [skymap.y1, skymap.y2],
        #             xerr=[0.1, 0.1], **dict(color='black',ls="-"))
        # ax.errorbar([skymap.x1, skymap.x2], [skymap.yc, skymap.yc],
        #             yerr=[0.1, 0.1], **dict(color='black',ls="-"))
    for t,x in zip(obs_times,[0.85,0.5,0.3,0.15]):
        bbox = dict(boxstyle='round', fc='white', ec='black', alpha=1.)
        ax.text(x, 0.80, f"{int(t)}", fontsize=12, bbox=bbox,
                     transform=ax.transAxes, horizontalalignment='right')

    # observational data of grb170817a
    obs = get_vla_and_hst_data()
    ax.errorbar(obs[:,2],obs[:,3],
                xerr=np.column_stack((obs[:,4],obs[:,5])).T,
                yerr=np.column_stack((obs[:,6],obs[:,7])).T,
                **{'ecolor': 'gray', 'mec': 'gray',
                   'marker': 'o', 'markersize': 8,
                   "capsize": 2, "linestyle": "None", 'mfc': "white",
                   "label": "HST, HSA and gVLBI", "zorder": 100})


    # im.set_rasteraized(True)


    ax.set_yscale("linear")
    ax.set_xscale("linear")
    if "xlim" in plot.keys() and len(plot["xlim"]) == 2: ax.set_xlim(plot["xlim"])
    else: ax.set_xlim(skymap.grid_x.min() * 1.1, skymap.grid_x.max() * 1.1)
    if "ylim" in plot.keys() and len(plot["ylim"]) == 2: ax.set_ylim(plot["ylim"])
    else: ax.set_ylim(skymap.grid_y.min() * 1.1, skymap.grid_y.max() * 1.1)
    ax.tick_params(axis='both', which='both', labelleft=True,
                        labelright=False, tick1On=True, tick2On=True,
                        labelsize=12, direction='in',
                        bottom=True, top=True, left=True, right=True)
    ax.minorticks_on()
    ax.axhline(y=0, linestyle=':', linewidth=0.4)
    ax.axvline(x=0, linestyle=':', linewidth=0.4)
    ax.grid(ls=':')
    ax.set_xlabel("$X$ [mas]",fontsize=12)
    ax.set_ylabel("$Z$ [mas]",fontsize=12)
    ax.set_aspect(aspect=1/2)
    ax.legend(fancybox=False, loc='lower right', columnspacing=0.8,
           # bbox_to_anchor=(0.5, 0.5),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
           shadow=False, ncol=1, fontsize=12, labelcolor='black',
           framealpha=0.0, borderaxespad=0.)
    if "figname" in plot.keys():
        plt.savefig(fig_dir+plot["figname"]+'.png',dpi=256)
    plt.savefig(fig_dir+plot["figname"]+'.pdf')
    if plot["show"]: plt.show()



if __name__ == '__main__':
    # https://arxiv.org/pdf/1808.00469
    # https://www.nature.com/articles/s41586-022-05145-7
    struct = dict(struct="gaussian",Eiso_c=np.power(10, 54.1), Gamma0c= 300., M0c= -1.,
                  theta_c= np.deg2rad(3.50), theta_w= np.deg2rad(25.))
    skymap_conf=dict(nx=128, ny=64, extend_grid=2, fwhm_fac=0.5, lat_dist_method="integ",
                     intp_filter=dict(type="gaussian", size=2, sigma=1.5, mode='reflect'),  # "gaussian"
                     hist_filter=dict(type="gaussian", size=2, sigma=1.5, mode='reflect'))

    compare_lcs_and_skymap_props(
        # struct = dict(struct="gaussian",Eiso_c=1.e52, Gamma0c= 300., M0c= -1., theta_c= 0.085, theta_w= 0.2618),
        pp = dict(main=dict(lc_freqs='array logspace 1e8 1e29 96',
                            lc_times='array logspace 3e3 1e10 128',
                            tb0=3e2,tb1=1e14,ntb=1000,
                            skymap_freqs="array 4.5e9 4.6e14", # Radio of VLA and Optical of HST
                            skymap_times="array logspace 1e3 1e8 100",
                            n_ism = np.power(10, -1.60),theta_obs=np.deg2rad(20.8),d_l = 1.27e+26, z = 0.0099),
                  grb=dict(structure=struct,eats_type='a',do_rs='no',bw_type='fs',do_lc='yes',do_skymap='yes',
                           eps_e_fs = np.power(10,-3.42), eps_b_fs =np.power(10,-4.02), p_fs = 2.10,
                           eps_e_rs = np.power(10,-3.42), eps_b_rs =np.power(10,-4.02), p_rs = 2.10,
                           method_ssc_fs='none',method_pp_fs='none',use_ssa_fs='no',
                           method_ssc_rs='none',method_pp_rs='none',use_ssa_rs='no',
                           ebl_tbl_fpath='none',
                           skymap_conf=skymap_conf)),
        plot = dict(xlim=(4, -4), ylim=(-2, 2), show=True, figname="lcs_170817a_gauss"),
        working_dir=working_dir+"tmp/",
        run=False,process_skymaps=False
    )

    # struct = dict(struct="gaussian",Eiso_c=np.power(10, 54.2), Gamma0c= 300., M0c= -1.,
    #               theta_c= np.deg2rad(4.0), theta_w= np.deg2rad(40.))
    # skymap_conf=dict(nx=128, ny=64, extend_grid=2, fwhm_fac=0.5, lat_dist_method="integ",
    #                  intp_filter=dict(type="gaussian", size=2, sigma=1.5, mode='reflect'),  # "gaussian"
    #                  hist_filter=dict(type="gaussian", size=2, sigma=1.5, mode='reflect'))
    # compare_lcs_and_skymap_props(
    #     # struct = dict(struct="gaussian",Eiso_c=1.e52, Gamma0c= 300., M0c= -1., theta_c= 0.085, theta_w= 0.2618),
    #     pp = dict(main=dict(lc_freqs='array logspace 1e8 1e29 96',
    #                         lc_times='array logspace 3e3 1e10 128',
    #                         tb0=3e1,tb1=1e15,ntb=2000,iout=2,
    #                         skymap_freqs="array 4.5e9 4.6e14", # Radio of VLA and Optical of HST
    #                         skymap_times="array logspace 1e3 1e8 100",
    #                         n_ism = np.power(10, -1.95),theta_obs=np.deg2rad(20.),d_l = 1.27e+26, z = 0.0099),
    #               grb=dict(structure=struct,eats_type='a',do_rs='no',bw_type='fs',do_lc='yes',do_skymap='yes',
    #                        eps_e_fs = np.power(10,-3.42), eps_b_fs =np.power(10,-4.02), p_fs = 2.10,
    #                        eps_e_rs = np.power(10,-3.42), eps_b_rs =np.power(10,-4.02), p_rs = 2.10,
    #                        method_ssc_fs='none',method_pp_fs='none',use_ssa_fs='no',
    #                        method_ssc_rs='none',method_pp_rs='none',use_ssa_rs='no',
    #                        ebl_tbl_fpath='none',do_mphys_in_situ='no',do_mphys_in_ppr='yes',
    #                        skymap_conf=skymap_conf)),
    #     plot = dict(xlim=(4, -4), ylim=(-2, 2), show=True, figname="lcs_170817a_gauss"),
    #     working_dir=working_dir+"tmp/",
    #     run=True,process_skymaps=True
    # )