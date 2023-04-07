import copy
import os.path
import numpy as np
import h5py
import hashlib
from scipy import interpolate

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, LogNorm
from matplotlib import cm, rc, rcParams
import matplotlib.colors as colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1 import ImageGrid
from mpl_toolkits.axisartist.grid_finder import MaxNLocator
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator

rc('text', usetex=True) # $ sudo apt-get install cm-super
rc('font', family='serif')
rcParams['font.size'] = 8

# from data.observations import GRB170817A, _plot_obs
# from prj_kenta.kenta_simulations import *
# from interfaces.afterglow import *
# from prj_skymaps.run_postprocess import make_hash, tmp_for_file2


from .utils import cgs, make_hash, latex_float, find_nearest_index
from .interface import PyBlastAfterglow, tmp_for_file2, combine_images, get_skymap_lat_dist, get_skymap_fwhm, \
    smooth_interpolated_skymap_with_gaussian_kernel

def precompute_skymaps(tasks_to_plot, settings):
    settings["precompute"] = True
    ix = -1
    iy = -1
    print("len={}".format(len(tasks_to_plot)))
    print("len={}".format(len(tasks_to_plot[0])))

    tmp_tile = h5py.File(settings["precompute_fpath"], "w")
    hashes = []
    i = 0
    for tmp in tasks_to_plot:
        ix += 1
        for task in tmp:
            iy += 1

            plot_dict = tasks_to_plot[ix][iy]
            time = plot_dict["time"]
            freq = plot_dict["freq"]

            hash, s  = make_hash(plot_dict, time, freq)
            print("ix={}/{} iy={}/{} i={}/{} time={} freq={:.1e} hash={}"
                  .format(ix,len(tasks_to_plot), iy, len(tmp), i, len(tasks_to_plot)*len(tmp),time,freq, hash))
            for _key in s.split("--"):
                print("\t{}".format(_key))
            if not hash in hashes:
                hashes.append(hash)
            else:
                raise ValueError("Hash:{} already in hashes={}".format(hash, hashes))
            grp = tmp_tile.create_group(hash)
            for key in ["res_dir_kn", "ejecta_prefix1", "ejecta_prefix2", "res_dir_grb", "jet_prefix"]:
                grp.attrs.create(key, task[key])
            grp.attrs.create("time", time)
            grp.attrs.create("freq", freq)
            grp.attrs.create("hash", s)

            xc, yc, x1, x2, y1, y2, int_x, int_y, int_zz = \
                tmp_for_file2(time=time, freq=freq, grp=grp, plot_dict=plot_dict, settings=settings)
            i = i + 1
        iy = -1
    for key in tmp_tile.keys():
        print(" group {} = {}".format(key, tmp_tile[key]))
        for _key in tmp_tile[key].keys():
            print("   data {} = {}".format(_key, np.array(tmp_tile[key][_key]).shape))
        for atr in tmp_tile[key].attrs.keys():
            print("   atr {} = {}".format(atr, tmp_tile[key].attrs[atr]))
    # for atr in tmp_tile.attrs.keys():
    #     print(" atr {} = {} ".format(atr, tmp_tile[atr]))
    # tmp_tile.close()
    print("Done precomputing. File saved:\n{}".format(settings["precompute_fpath"]))


def plot_text2(ax, dic):
    """'text':{'x':0.5, 'y':0.5, 'text':'my_text', 'fs':12, 'color':'black',
               'horizontalalignment':'center', 'transform':True}"""
    # print("---------")
    xcorr = dic['x']
    ycorr = dic['y']
    text = dic['text']
    fs = dic['fs']
    color = dic['color']
    horal = dic['horizontalalignment']
    if 'transform' in dic.keys() and dic['transform']:
        ax.text(xcorr, ycorr, text, color=color, fontsize=fs, horizontalalignment=horal, transform=ax.transAxes)
    else:
        ax.text(xcorr, ycorr, text, color=color, fontsize=fs, horizontalalignment=horal)
    return 0

def plot_pcolomesh(ax, task_dict, int_x, int_y, int_zz, outfname=None):
    norm_min = task_dict["norm"][1]
    norm_max = task_dict["norm"][2]

    if not hasattr(norm_min, "__contains__"):
        pass
    else:
        if norm_min.__contains__("min"):
            norm_min = float(norm_min.replace("min", "")) * int_zz[np.isfinite(int_zz)].min()
        elif norm_min.__contains__("max"):
            norm_min = float(norm_min.replace("max", "")) * int_zz[np.isfinite(int_zz)].max()
        else:
            raise KeyError("norm_min does not contain 'min' or 'max' ")

    if not hasattr(norm_max, "__contains__"):
        pass
    else:
        if norm_max.__contains__("max"):
            norm_max = float(norm_max.replace("max", "")) * int_zz[np.isfinite(int_zz)].max()
        elif norm_max.__contains__("min"):
            norm_max = float(norm_max.replace("min", "")) * int_zz[np.isfinite(int_zz)].min()
        else:
            raise KeyError("norm_max does not contain 'min' or 'max' ")

    if task_dict["norm"][0] == "linear":
        norm = Normalize(norm_min, norm_max)
    elif task_dict["norm"][0] == "log":
        norm = LogNorm(norm_min, norm_max)
    elif task_dict["norm"][0] == "twoslope":
        norm = colors.TwoSlopeNorm(vmin=task_dict["norm"][1], vcenter=task_dict["norm"][2], vmax=task_dict["norm"][3])
    else:
        raise KeyError("norm is wrong")

    # levels = MaxNLocator(nbins=40).tick_values(int_ration.min(), int_ration.max())
    # # levels = MaxNLocator(nbins=40).tick_values(-5, 1)
    # cmap = plt.get_cmap(tmp["pcolormesh"]["cmap"])
    # norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

    cmap = cm.get_cmap(task_dict["cmap"])
    if not task_dict["set_under"] is None: cmap.set_under(task_dict["set_under"])
    if not task_dict["set_over"] is None: cmap.set_over(task_dict["set_over"])
    im = ax.pcolormesh(int_x, int_y, int_zz, norm=norm, cmap=cmap, alpha=float(task_dict["alpha"]))
    # im = ax_main.contourf(int_x, int_y, int_zz, norm=norm, cmap=cmap, alpha=1.0)
    # im.set_rasteraized(True)
    if not (task_dict["facecolor"] is None):
        ax.set_facecolor(cmap(float(task_dict["facecolor"])))

    ''' --------- SAVE FILE --------- '''
    # if not outfname is None:
    #     arr = copy.deepcopy(int_zz)
    #     arr[~np.isfinite(arr)] = 0.0
    #     print("max={} min={}".format(arr.max(), arr.min()))
    #     arr = arr / arr.max()
    #     # arr = arr.T
    #     out = np.zeros((len(int_x) * len(int_y), 3))
    #     i = 0
    #     for ix in range(len(int_x)):
    #         for iy in range(len(int_y)):
    #             out[i][0] = int_x[ix]
    #             out[i][1] = int_y[iy]
    #             out[i][2] = arr[iy, ix]
    #             i = i+1
    #     print(out.shape)
    #     np.savetxt(outfname + ".txt", out, fmt=("%.4f", "%.4f", "%.4e"), header="X[mas] Z[mas] I[normalized]",
    #                footer="# X has {} points Z has {} points".format(len(int_x), len(int_y)))
    #     # arr = int_zz/int_zz.max()
    #
    #     dfile = h5py.File(outfname + ".h5", "w")
    #     dfile.create_dataset(name="X[mas]", data=int_x)
    #     dfile.create_dataset(name="Z[mas]", data=int_y)
    #     dfile.create_dataset(name="intensity[normalized]", data=arr)
    #     dfile.close()
    #     exit(1)

    return im

def read_precompute_file_get_gals(settings,plot_dict,hashes=None):
    try:
        tmp_tile = h5py.File(plot_dict["precompute_fpath"], "r")
        print("Loading: {}".format(plot_dict["precompute_fpath"]))
    except OSError:
        raise OSError("failed reading file: \n{} \n Try recomputing pre-compute file".format(plot_dict["precompute_fpath"]))
    times = []
    freqs = []
    print(tmp_tile.keys())
    for key in tmp_tile.keys():
        grp = tmp_tile[key]
        if not "time" in grp.attrs.keys():
            print("Available keys: {}".format(grp.attrs.keys()))
            raise KeyError("not found 'time' in grp.attrs.keys() : \n{}".format(grp.attrs.keys()))
        times.append(float(grp.attrs["time"]))
        freqs.append(float(grp.attrs["freq"]))
    times = set(times)
    print("Times located: {}".format(times))
    req_freq = plot_dict["freq"]
    if not req_freq in freqs:
        raise ValueError("Req. freq. {} not found in pre-computed file. \n{} \n Available:{}"
                         .format(req_freq, plot_dict["precompute_fpath"], freqs))
    # vals = {
    #     "xc":[],
    #     "ysize":[],
    #     "xsize":[],
    # }
    _vals = {
        "kn_skymap":{
            "xc": [],
            "ysize": [],
            "xsize": []
        },
        "grb_skymap":{
            "xc": [],
            "ysize": [],
            "xsize": []
        },
        "kn_grb_skymap":{
            "xc": [],
            "ysize": [],
            "xsize": []
        },
        "kn_w_skymap":{
            "xc": [],
            "ysize": [],
            "xsize": []
        }
    }
    vals = copy.deepcopy(_vals)
    for _time in times:

        hash, s = make_hash(plot_dict, _time, req_freq)
        if (not hashes is None) and (not hash in hashes):
            hashes.append({ hash : s })
        elif (not hashes is None):
            raise ValueError("Hash:{} already in hashes={}".format(hash, hashes))
        else:
            pass
        if (not hash in tmp_tile.keys()):
            print("ERROR! Hash not found: {}".format(hash))
            print("--------------- Searching: ---------------------- ")
            for _key in s.split("--"):
                print("\t{}".format(_key))
            print("--------------- Available: ---------------------- ")
            for dataset_hash in tmp_tile.keys():
                dataset_str = tmp_tile[dataset_hash].attrs["hash"]
                print(" Hash={}".format(dataset_hash))
                for _key in dataset_str.split("--"):
                    print("\t{}".format(_key))
            # for key in hash.split(" "):
            #     print(key)
            print("FIle: {}".format(plot_dict["precompute_fpath"]))
            raise KeyError("Exiting...")

        grp = tmp_tile[hash]
        _plot_dict = copy.deepcopy(plot_dict)
        if "precompute_fpath" in _plot_dict.keys():  _plot_dict.pop("precompute_fpath")
        if "line" in _plot_dict.keys():  _plot_dict.pop("line")

        ''' ---- plot size/cm of the ejecta ---- '''
        pb = BPA_METHODS()
        pb.res_dir_kn = plot_dict["res_dir_kn"]
        pb.ejecta_prefix = plot_dict["ejecta_prefix1"]
        pb.res_dir_grb = plot_dict["res_dir_grb"]
        pb.jet_prefix = plot_dict["jet_prefix"]
        if len(settings["kn_skymap"].keys()) > 0:
            tmp = copy.deepcopy(settings["kn_skymap"])
            xc, yc, x1, x2, y1, y2, int_x, int_y, int_zz = tmp_for_file2(time=_time, freq=req_freq, grp=grp,
                                                                         plot_dict=plot_dict,
                                                                         settings=settings, kn_skymap=True)
            if (("xc" in tmp.keys()) and (len(tmp["xc"].keys()) > 0)) or ("dxcdt" in tmp.keys()): vals["kn_skymap"]["xc"].append(xc)
            if ("ysize" in tmp.keys()) and (len(tmp["ysize"].keys()) > 0): vals["kn_skymap"]["ysize"].append(y2 - y1)
            if ("xsize" in tmp.keys()) and (len(tmp["xsize"].keys()) > 0): vals["kn_skymap"]["xsize"].append(x2 - x1)

        ''' ---- plot size/cm of the grb ---- '''
        if (len(settings["grb_skymap"].keys()) > 0):
            tmp = copy.deepcopy(settings["grb_skymap"])
            xc, yc, x1, x2, y1, y2, int_x, int_y, int_zz = tmp_for_file2(time=_time, freq=req_freq, grp=grp,
                                                                         plot_dict=plot_dict,
                                                                         settings=settings, grb_skymap=True)
            if (("xc" in tmp.keys()) and (len(tmp["xc"].keys()) > 0)) or ("dxcdt" in tmp.keys()): vals["grb_skymap"]["xc"].append(xc)
            if ("ysize" in tmp.keys()) and (len(tmp["ysize"].keys()) > 0): vals["grb_skymap"]["ysize"].append(y2 - y1)
            if ("xsize" in tmp.keys()) and (len(tmp["xsize"].keys()) > 0): vals["grb_skymap"]["xsize"].append(x2 - x1)

        ''' --- plot size/cm of the kn with grb (using previous two) ---'''
        if (len(settings["kn_grb_skymap"].keys()) > 0):
            tmp = copy.deepcopy(settings["kn_grb_skymap"])
            xc, yc, x1, x2, y1, y2, int_x, int_y, int_zz = tmp_for_file2(time=_time, freq=req_freq, grp=grp,
                                                                         plot_dict=plot_dict,
                                                                         settings=settings, kn_grb_skymap=True)
            if (("xc" in tmp.keys()) and (len(tmp["xc"].keys()) > 0)) or ("dxcdt" in tmp.keys()): vals["kn_grb_skymap"]["xc"].append(xc)
            if ("ysize" in tmp.keys()) and (len(tmp["ysize"].keys()) > 0): vals["kn_grb_skymap"]["ysize"].append(y2 - y1)
            if ("xsize" in tmp.keys()) and (len(tmp["xsize"].keys()) > 0): vals["kn_grb_skymap"]["xsize"].append(x2 - x1)

        ''' ---- plot size/cm of the ejecta with interaction ---- '''
        pb_w = BPA_METHODS()
        pb_w.res_dir_kn = plot_dict["res_dir_kn"]
        pb_w.ejecta_prefix = plot_dict["ejecta_prefix2"]
        pb_w.res_dir_grb = plot_dict["res_dir_grb"]
        pb_w.jet_prefix = plot_dict["jet_prefix"]
        if (len(settings["kn_w_skymap"].keys()) > 0):
            tmp = copy.deepcopy(settings["kn_w_skymap"])
            xc, yc, x1, x2, y1, y2, int_x, int_y, int_zz = tmp_for_file2(time=_time, freq=req_freq, grp=grp,
                                                                         plot_dict=plot_dict,
                                                                         settings=settings, kn_w_skymap=True)
            if (("xc" in tmp.keys()) and (len(tmp["xc"].keys()) > 0)) or ("dxcdt" in tmp.keys()): vals["kn_w_skymap"]["xc"].append(xc)
            if ("ysize" in tmp.keys()) and (len(tmp["ysize"].keys()) > 0): vals["kn_w_skymap"]["ysize"].append(y2 - y1)
            if ("xsize" in tmp.keys()) and (len(tmp["xsize"].keys()) > 0): vals["kn_w_skymap"]["xsize"].append(x2 - x1)



    # _vals = copy.deepcopy(_vals)
    for mode in vals.keys():
        for key in vals[mode].keys():
            if (len(vals[mode][key]) > 1):
                _times, _vals_arr = zip(*sorted(zip(times, vals[mode][key])))
                _vals[mode][key] = np.array(_vals_arr)

    # times, vals = zip(*sorted(zip(times, vals)))
    times = np.array(_times)
    # vals = np.array(vals)

    keys = ["kn_skymap", "grb_skymap", "kn_grb_skymap"]
    for mode in vals.keys():
        if len(settings[mode].keys()) > 0:
            if ("dxcdt" in settings[mode].keys()):
                if not len(_vals[mode]["xc"]) > 0:
                    raise KeyError("In order to add 'dxcdt' the 'xc' must be computed first for {}".format(mode))
                # _vals[mode]["dxcdt"] = np.diff(_vals[mode]["xc"]*cgs.rad2mas) / np.diff(times*cgs.day)
                _vals[mode]["dxcdt"] = np.diff(_vals[mode]["xc"]) / np.diff(times) * 236.0 # ms->rad? days->sec?... Conversion factor from Fernandez table

                t1 = 75.
                t2 = 230.
                xc1 = _vals[mode]["xc"][ find_nearest_index(times,t1) ]
                xc2 = _vals[mode]["xc"][ find_nearest_index(times,t2) ]
                print( (xc2-xc1) )
                print( (xc2-xc1)/(t2-t1) )
                print( (xc2-xc1)/(t2-t1)*236.0 )

                t=800.
                print("{} beta_c={} at {} days".format(mode, _vals[mode]["dxcdt"][find_nearest_index(times, t)], t))

                dx = _vals[mode]["dxcdt"]
                _vals[mode]["dxcdt"] = np.insert(_vals[mode]["dxcdt"], -1, _vals[mode]["dxcdt"][-1])



    for mode in vals.keys():
        for key in vals[mode].keys():
            if len(_vals[mode][key]) > 0:
                assert len(times) == len(_vals[mode][key])

    return (times, _vals)

    # ax.plot(times, vals, **plot_dict["line"])

def plot_skymaps(tasks_to_plot, settings):
    # # tasks_to_plot = [
    # #     [{"res_dir1": "/media/vsevolod/Data/KentaData/SFHo_135_135_150m_11//afterglow/testdata/", "time1": 1000,
    # #       "ejecta_prefix1": "ejecta_text26_theta0_nism005_p205_epse01_epsb001_epst1_",
    # #       "res_dir2": "/media/vsevolod/Data/KentaData/SFHo_135_135_150m_11//afterglow/testdata/", "time2": 1000,
    # #       "ejecta_prefix2": "ejecta_text26_theta0_nism005_p205_epse01_epsb001_epst1_"}
    # #         ,
    # #          {"res_dir1": "/media/vsevolod/Data/KentaData/SFHo_135_135_150m_11//afterglow/testdata/", "time1": 1000,
    # #           "ejecta_prefix1": "ejecta_text26_theta45_nism005_p205_epse01_epsb001_epst1_",
    # #           "res_dir2": "/media/vsevolod/Data/KentaData/SFHo_135_135_150m_11//afterglow/testdata/", "time2": 1000,
    # #           "ejecta_prefix2": "ejecta_text26_theta45_nism005_p205_epse01_epsb001_epst1_"
    # #           }]
    # #     ,
    # #     [{"res_dir1": "/media/vsevolod/Data/KentaData/SFHo_135_135_150m_11//afterglow/testdata/", "time1": 1000,
    # #       "ejecta_prefix1": "ejecta_text26_theta75_nism005_p205_epse01_epsb001_epst1_",
    # #       "res_dir2": "/media/vsevolod/Data/KentaData/SFHo_135_135_150m_11//afterglow/testdata/", "time2": 1000,
    # #       "ejecta_prefix2": "ejecta_text26_theta75_nism005_p205_epse01_epsb001_epst1_"
    # #       }
    # #         ,
    # #          {"res_dir1": "/media/vsevolod/Data/KentaData/SFHo_135_135_150m_11//afterglow/testdata/", "time1": 1000,
    # #           "ejecta_prefix1": "ejecta_text26_theta90_nism005_p205_epse01_epsb001_epst1_",
    # #           "res_dir2": "/media/vsevolod/Data/KentaData/SFHo_135_135_150m_11//afterglow/testdata/", "time2": 1000,
    # #           "ejecta_prefix2": "ejecta_text26_theta90_nism005_p205_epse01_epsb001_epst1_"
    # #           }]
    # # ]
    #
    # # tasks_to_plot = [
    # #     [{"res_dir_kn": "/media/vsevolod/Data/KentaData/SFHo_135_135_150m_11//afterglow/testdata/", "time_kn": 220,
    # #       "ejecta_prefix": "ejecta_text26_theta0_nism005_p205_epse01_epsb001_epst1_",
    # #       "res_dir_grb": "/media/vsevolod/Data/KentaData/SFHo_135_135_150m_11//afterglow/testdata/", "time_jet": 220,
    # #       "jet_prefix": "jet_G_nism005_thetaobs0_epse007079_epsb000525_p216_epst00_nlayers180_"}
    # #         ,
    # #          {"res_dir_kn": "/media/vsevolod/Data/KentaData/SFHo_135_135_150m_11//afterglow/testdata/", "time_kn": 220,
    # #           "ejecta_prefix": "ejecta_text26_theta25_nism005_p205_epse01_epsb001_epst1_",
    # #           "res_dir_grb": "/media/vsevolod/Data/KentaData/SFHo_135_135_150m_11//afterglow/testdata/", "time_jet": 220,
    # #           "jet_prefix": "jet_G_nism005_thetaobs25_epse007079_epsb000525_p216_epst00_nlayers180_"
    # #           }]
    # #     ,
    # #     [{"res_dir_kn": "/media/vsevolod/Data/KentaData/SFHo_135_135_150m_11//afterglow/testdata/", "time_kn": 220,
    # #       "ejecta_prefix": "ejecta_text26_theta75_nism005_p205_epse01_epsb001_epst1_",
    # #       "res_dir_grb": "/media/vsevolod/Data/KentaData/SFHo_135_135_150m_11//afterglow/testdata/", "time_jet": 220,
    # #       "jet_prefix": "jet_G_nism005_thetaobs75_epse007079_epsb000525_p216_epst00_nlayers180_"
    # #       }
    # #         ,
    # #          {"res_dir_kn": "/media/vsevolod/Data/KentaData/SFHo_135_135_150m_11//afterglow/testdata/", "time_kn": 220,
    # #           "ejecta_prefix": "ejecta_text26_theta90_nism005_p205_epse01_epsb001_epst1_",
    # #           "res_dir_grb": "/media/vsevolod/Data/KentaData/SFHo_135_135_150m_11//afterglow/testdata/", "time_jet": 220,
    # #           "jet_prefix": "jet_G_nism005_thetaobs90_epse007079_epsb000525_p216_epst00_nlayers180_"
    # #           }
    # #     ]
    # # ]
    # tasks_to_plot = [
    #     [
    #      # {"res_dir_kn": "/media/vsevolod/Data/KentaData/SFHo_135_135_150m_11//afterglow/testdata/", "time_kn": 220,
    #      #  "ejecta_prefix": "ejecta_text26_theta0_nism000031_p205_epse01_epsb001_epst1_",
    #      #  "res_dir_grb": "/media/vsevolod/Data/KentaData/SFHo_135_135_150m_11//afterglow/testdata/", "time_jet": 220,
    #      #  "jet_prefix": "jet_G_nism000031_thetaobs0_epse007079_epsb000525_p216_epst00_nlayers180_"}
    #      #    ,
    #      #     {"res_dir_kn": "/media/vsevolod/Data/KentaData/SFHo_135_135_150m_11//afterglow/testdata/", "time_kn": 140,
    #      #      "ejecta_prefix": "ejecta_text26_theta25_nism000031_p205_epse01_epsb001_epst1_",
    #      #      "res_dir_grb": "/media/vsevolod/Data/KentaData/SFHo_135_135_150m_11//afterglow/testdata/", "time_jet": 140,
    #      #      "jet_prefix": "jet_G_nism000031_thetaobs25_epse007079_epsb000525_p216_epst00_nlayers180_"
    #      #      }
    #     # ]
    #     # ,
    #     # [
    #       # {"res_dir_kn": "/media/vsevolod/Data/KentaData/SFHo_135_135_150m_11//afterglow/testdata/", "time_kn": 220,
    #       # "ejecta_prefix": "ejecta_text26_theta75_nism000031_p205_epse01_epsb001_epst1_",
    #       # "res_dir_grb": "/media/vsevolod/Data/KentaData/SFHo_135_135_150m_11//afterglow/testdata/", "time_jet": 220,
    #       # "jet_prefix": "jet_G_nism000031_thetaobs75_epse007079_epsb000525_p216_epst00_nlayers180_"
    #       # }
    #       #   ,
    #       #    {"res_dir_kn": "/media/vsevolod/Data/KentaData/SFHo_135_135_150m_11//afterglow/testdata/", "time_kn": 220,
    #       #     "ejecta_prefix": "ejecta_text26_theta90_nism000031_p205_epse01_epsb001_epst1_",
    #       #     "res_dir_grb": "/media/vsevolod/Data/KentaData/SFHo_135_135_150m_11//afterglow/testdata/", "time_jet": 220,
    #       #     "jet_prefix": "jet_G_nism000031_thetaobs90_epse007079_epsb000525_p216_epst00_nlayers180_"
    #       #     }
    #     ]
    # ]
    # tasks_to_plot = [
    #     [
    #         {"res_dir_kn": "/media/vsevolod/Data/KentaData/SFHo_135_135_150m_11//afterglow/testdata/", "time_kn": 120,
    #          "ejecta_prefix1": "ejecta_text26_theta45_nism000031_p205_epse01_epsb001_epst1_",
    #          "ejecta_prefix2": "ejecta_WG_text26_theta45_nism000031_p205_epse01_epsb001_epst1_",
    #          "res_dir_grb": "/media/vsevolod/Data/KentaData/SFHo_135_135_150m_11//afterglow/testdata/", "time_jet": 120,
    #          "jet_prefix": "jet_G_nism000031_thetaobs45_epse007079_epsb000525_p216_epst00_nlayers180_"
    #          }
    #     ]
    # ]
    # theta="75"
    # nism="000031"
    # tasks_to_plot_ = [
    #     [
    #         {"res_dir_kn": "/media/vsevolod/Data/KentaData/SFHo_135_135_150m_11//afterglow/testdata/", "time_kn": 20,
    #          "ejecta_prefix1": "ejecta_text26_theta{}_nism{}_p205_epse01_epsb001_epst1_".format(theta,nism),
    #          "ejecta_prefix2": "ejecta_WG_text26_theta{}_nism{}_p205_epse01_epsb001_epst1_".format(theta,nism),
    #          "res_dir_grb": "/media/vsevolod/Data/KentaData/SFHo_135_135_150m_11//afterglow/testdata/", "time_jet": 20,
    #          "jet_prefix": "jet_G_nism{}_thetaobs{}_epse007079_epsb000525_p216_epst00_nlayers180_".format(nism,theta)
    #          },
    #         {"res_dir_kn": "/media/vsevolod/Data/KentaData/SFHo_135_135_150m_11//afterglow/testdata/", "time_kn": 40,
    #          "ejecta_prefix1": "ejecta_text26_theta{}_nism{}_p205_epse01_epsb001_epst1_".format(theta,nism),
    #          "ejecta_prefix2": "ejecta_WG_text26_theta{}_nism{}_p205_epse01_epsb001_epst1_".format(theta,nism),
    #          "res_dir_grb": "/media/vsevolod/Data/KentaData/SFHo_135_135_150m_11//afterglow/testdata/", "time_jet": 40,
    #          "jet_prefix": "jet_G_nism{}_thetaobs{}_epse007079_epsb000525_p216_epst00_nlayers180_".format(nism,theta)
    #          },
    #         {"res_dir_kn": "/media/vsevolod/Data/KentaData/SFHo_135_135_150m_11//afterglow/testdata/", "time_kn": 80,
    #          "ejecta_prefix1": "ejecta_text26_theta{}_nism{}_p205_epse01_epsb001_epst1_".format(theta,nism),
    #          "ejecta_prefix2": "ejecta_WG_text26_theta{}_nism{}_p205_epse01_epsb001_epst1_".format(theta,nism),
    #          "res_dir_grb": "/media/vsevolod/Data/KentaData/SFHo_135_135_150m_11//afterglow/testdata/", "time_jet": 80,
    #          "jet_prefix": "jet_G_nism{}_thetaobs{}_epse007079_epsb000525_p216_epst00_nlayers180_".format(nism,theta)
    #          },
    #     ],
    #     [
    #         {"res_dir_kn": "/media/vsevolod/Data/KentaData/SFHo_135_135_150m_11//afterglow/testdata/", "time_kn": 120,
    #          "ejecta_prefix1": "ejecta_text26_theta{}_nism{}_p205_epse01_epsb001_epst1_".format(theta,nism),
    #          "ejecta_prefix2": "ejecta_WG_text26_theta{}_nism{}_p205_epse01_epsb001_epst1_".format(theta,nism),
    #          "res_dir_grb": "/media/vsevolod/Data/KentaData/SFHo_135_135_150m_11//afterglow/testdata/", "time_jet": 120,
    #          "jet_prefix": "jet_G_nism{}_thetaobs{}_epse007079_epsb000525_p216_epst00_nlayers180_".format(nism,theta)
    #          },
    #         {"res_dir_kn": "/media/vsevolod/Data/KentaData/SFHo_135_135_150m_11//afterglow/testdata/", "time_kn": 160,
    #          "ejecta_prefix1": "ejecta_text26_theta{}_nism{}_p205_epse01_epsb001_epst1_".format(theta,nism),
    #          "ejecta_prefix2": "ejecta_WG_text26_theta{}_nism{}_p205_epse01_epsb001_epst1_".format(theta,nism),
    #          "res_dir_grb": "/media/vsevolod/Data/KentaData/SFHo_135_135_150m_11//afterglow/testdata/", "time_jet": 160,
    #          "jet_prefix": "jet_G_nism{}_thetaobs{}_epse007079_epsb000525_p216_epst00_nlayers180_".format(nism,theta)
    #          },
    #         {"res_dir_kn": "/media/vsevolod/Data/KentaData/SFHo_135_135_150m_11//afterglow/testdata/", "time_kn": 220,
    #          "ejecta_prefix1": "ejecta_text26_theta{}_nism{}_p205_epse01_epsb001_epst1_".format(theta,nism),
    #          "ejecta_prefix2": "ejecta_WG_text26_theta{}_nism{}_p205_epse01_epsb001_epst1_".format(theta,nism),
    #          "res_dir_grb": "/media/vsevolod/Data/KentaData/SFHo_135_135_150m_11//afterglow/testdata/", "time_jet": 220,
    #          "jet_prefix": "jet_G_nism{}_thetaobs{}_epse007079_epsb000525_p216_epst00_nlayers180_".format(nism,theta)
    #          }
    #     ]
    #     # [
    #     #     {"res_dir_kn": "/media/vsevolod/Data/KentaData/SFHo_135_135_150m_11//afterglow/testdata/", "time_kn": 260,
    #     #      "ejecta_prefix1": "ejecta_text26_theta45_nism000031_p205_epse01_epsb001_epst1_",
    #     #      "ejecta_prefix2": "ejecta_WG_text26_theta45_nism000031_p205_epse01_epsb001_epst1_",
    #     #      "res_dir_grb": "/media/vsevolod/Data/KentaData/SFHo_135_135_150m_11//afterglow/testdata/", "time_jet": 260,
    #     #      "jet_prefix": "jet_G_nism000031_thetaobs45_epse007079_epsb000525_p216_epst00_nlayers180_"
    #     #      },
    #     #     {"res_dir_kn": "/media/vsevolod/Data/KentaData/SFHo_135_135_150m_11//afterglow/testdata/", "time_kn": 320,
    #     #      "ejecta_prefix1": "ejecta_text26_theta45_nism000031_p205_epse01_epsb001_epst1_",
    #     #      "ejecta_prefix2": "ejecta_WG_text26_theta45_nism000031_p205_epse01_epsb001_epst1_",
    #     #      "res_dir_grb": "/media/vsevolod/Data/KentaData/SFHo_135_135_150m_11//afterglow/testdata/", "time_jet": 320,
    #     #      "jet_prefix": "jet_G_nism000031_thetaobs45_epse007079_epsb000525_p216_epst00_nlayers180_"
    #     #      },
    #     #     {"res_dir_kn": "/media/vsevolod/Data/KentaData/SFHo_135_135_150m_11//afterglow/testdata/", "time_kn": 400,
    #     #      "ejecta_prefix1": "ejecta_text26_theta45_nism000031_p205_epse01_epsb001_epst1_",
    #     #      "ejecta_prefix2": "ejecta_WG_text26_theta45_nism000031_p205_epse01_epsb001_epst1_",
    #     #      "res_dir_grb": "/media/vsevolod/Data/KentaData/SFHo_135_135_150m_11//afterglow/testdata/", "time_jet": 400,
    #     #      "jet_prefix": "jet_G_nism000031_thetaobs45_epse007079_epsb000525_p216_epst00_nlayers180_"
    #     #      }
    #     # ],
    #     # [
    #     #     {"res_dir_kn": "/media/vsevolod/Data/KentaData/SFHo_135_135_150m_11//afterglow/testdata/", "time_kn": 500,
    #     #      "ejecta_prefix1": "ejecta_text26_theta45_nism000031_p205_epse01_epsb001_epst1_",
    #     #      "ejecta_prefix2": "ejecta_WG_text26_theta45_nism000031_p205_epse01_epsb001_epst1_",
    #     #      "res_dir_grb": "/media/vsevolod/Data/KentaData/SFHo_135_135_150m_11//afterglow/testdata/", "time_jet": 500,
    #     #      "jet_prefix": "jet_G_nism000031_thetaobs45_epse007079_epsb000525_p216_epst00_nlayers180_"
    #     #      },
    #     #     {"res_dir_kn": "/media/vsevolod/Data/KentaData/SFHo_135_135_150m_11//afterglow/testdata/", "time_kn": 600,
    #     #      "ejecta_prefix1": "ejecta_text26_theta45_nism000031_p205_epse01_epsb001_epst1_",
    #     #      "ejecta_prefix2": "ejecta_WG_text26_theta45_nism000031_p205_epse01_epsb001_epst1_",
    #     #      "res_dir_grb": "/media/vsevolod/Data/KentaData/SFHo_135_135_150m_11//afterglow/testdata/", "time_jet": 600,
    #     #      "jet_prefix": "jet_G_nism000031_thetaobs45_epse007079_epsb000525_p216_epst00_nlayers180_"
    #     #      },
    #     #     {"res_dir_kn": "/media/vsevolod/Data/KentaData/SFHo_135_135_150m_11//afterglow/testdata/", "time_kn": 800,
    #     #      "ejecta_prefix1": "ejecta_text26_theta45_nism000031_p205_epse01_epsb001_epst1_",
    #     #      "ejecta_prefix2": "ejecta_WG_text26_theta45_nism000031_p205_epse01_epsb001_epst1_",
    #     #      "res_dir_grb": "/media/vsevolod/Data/KentaData/SFHo_135_135_150m_11//afterglow/testdata/", "time_jet": 800,
    #     #      "jet_prefix": "jet_G_nism000031_thetaobs45_epse007079_epsb000525_p216_epst00_nlayers180_"
    #     #      }
    #     # ],
    #     # [
    #     #     {"res_dir_kn": "/media/vsevolod/Data/KentaData/SFHo_135_135_150m_11//afterglow/testdata/", "time_kn": 1000,
    #     #      "ejecta_prefix1": "ejecta_text26_theta45_nism000031_p205_epse01_epsb001_epst1_",
    #     #      "ejecta_prefix2": "ejecta_WG_text26_theta45_nism000031_p205_epse01_epsb001_epst1_",
    #     #      "res_dir_grb": "/media/vsevolod/Data/KentaData/SFHo_135_135_150m_11//afterglow/testdata/", "time_jet": 1000,
    #     #      "jet_prefix": "jet_G_nism000031_thetaobs45_epse007079_epsb000525_p216_epst00_nlayers180_"
    #     #      },
    #     #     {"res_dir_kn": "/media/vsevolod/Data/KentaData/SFHo_135_135_150m_11//afterglow/testdata/", "time_kn": 1200,
    #     #      "ejecta_prefix1": "ejecta_text26_theta45_nism000031_p205_epse01_epsb001_epst1_",
    #     #      "ejecta_prefix2": "ejecta_WG_text26_theta45_nism000031_p205_epse01_epsb001_epst1_",
    #     #      "res_dir_grb": "/media/vsevolod/Data/KentaData/SFHo_135_135_150m_11//afterglow/testdata/", "time_jet": 1200,
    #     #      "jet_prefix": "jet_G_nism000031_thetaobs45_epse007079_epsb000525_p216_epst00_nlayers180_"
    #     #      },
    #     #     {"res_dir_kn": "/media/vsevolod/Data/KentaData/SFHo_135_135_150m_11//afterglow/testdata/", "time_kn": 1400,
    #     #      "ejecta_prefix1": "ejecta_text26_theta45_nism000031_p205_epse01_epsb001_epst1_",
    #     #      "ejecta_prefix2": "ejecta_WG_text26_theta45_nism000031_p205_epse01_epsb001_epst1_",
    #     #      "res_dir_grb": "/media/vsevolod/Data/KentaData/SFHo_135_135_150m_11//afterglow/testdata/", "time_jet": 1400,
    #     #      "jet_prefix": "jet_G_nism000031_thetaobs45_epse007079_epsb000525_p216_epst00_nlayers180_"
    #     #      }
    #     # ]
    # ]
    # # tasks_to_plot = [
    # #     # [
    # #     #     {"res_dir_kn": "/media/vsevolod/Data/KentaData/SFHo_135_135_150m_11//afterglow/testdata/", "time_kn": 40,
    # #     #      "ejecta_prefix1": "ejecta_text26_nlayers150_theta{}_nism{}_p205_epse01_epsb001_epst1_".format(theta,nism),
    # #     #      "ejecta_prefix2": "ejecta_WG_text26_nlayers150_theta{}_nism{}_p205_epse01_epsb001_epst1_".format(theta,nism),
    # #     #      "res_dir_grb": "/media/vsevolod/Data/KentaData/SFHo_135_135_150m_11//afterglow/testdata/", "time_jet": 40,
    # #     #      "jet_prefix": "jet_G_nism{}_thetaobs{}_epse007079_epsb000525_p216_epst00_nlayers180_".format(nism,theta)
    # #     #      },
    # #     #     {"res_dir_kn": "/media/vsevolod/Data/KentaData/SFHo_135_135_150m_11//afterglow/testdata/", "time_kn": 80,
    # #     #      "ejecta_prefix1": "ejecta_text26_nlayers150_theta{}_nism{}_p205_epse01_epsb001_epst1_".format(theta,nism),
    # #     #      "ejecta_prefix2": "ejecta_WG_text26_nlayers150_theta{}_nism{}_p205_epse01_epsb001_epst1_".format(theta,nism),
    # #     #      "res_dir_grb": "/media/vsevolod/Data/KentaData/SFHo_135_135_150m_11//afterglow/testdata/", "time_jet": 80,
    # #     #      "jet_prefix": "jet_G_nism{}_thetaobs{}_epse007079_epsb000525_p216_epst00_nlayers180_".format(nism,theta)
    # #     #      },
    # #     #     {"res_dir_kn": "/media/vsevolod/Data/KentaData/SFHo_135_135_150m_11//afterglow/testdata/", "time_kn": 120,
    # #     #      "ejecta_prefix1": "ejecta_text26_nlayers150_theta{}_nism{}_p205_epse01_epsb001_epst1_".format(theta,nism),
    # #     #      "ejecta_prefix2": "ejecta_WG_text26_nlayers150_theta{}_nism{}_p205_epse01_epsb001_epst1_".format(theta,nism),
    # #     #      "res_dir_grb": "/media/vsevolod/Data/KentaData/SFHo_135_135_150m_11//afterglow/testdata/", "time_jet": 120,
    # #     #      "jet_prefix": "jet_G_nism{}_thetaobs{}_epse007079_epsb000525_p216_epst00_nlayers180_".format(nism,theta)
    # #     #      },
    # #     # ],
    # #     [
    # #         {"res_dir_kn": "/media/vsevolod/Data/KentaData/SFHo_135_135_150m_11//afterglow/testdata/", "time_kn": 160,
    # #          "ejecta_prefix1": "ejecta_text26_nlayers150_theta{}_nism{}_p205_epse01_epsb001_epst1_".format(theta,nism),
    # #          "ejecta_prefix2": "ejecta_WG_text26_nlayers150_theta{}_nism{}_p205_epse01_epsb001_epst1_".format(theta,nism),
    # #          "res_dir_grb": "/media/vsevolod/Data/KentaData/SFHo_135_135_150m_11//afterglow/testdata/", "time_jet": 160,
    # #          "jet_prefix": "jet_G_nism{}_thetaobs{}_epse007079_epsb000525_p216_epst00_nlayers180_".format(nism,theta)
    # #          },
    # #         {"res_dir_kn": "/media/vsevolod/Data/KentaData/SFHo_135_135_150m_11//afterglow/testdata/", "time_kn": 180,
    # #          "ejecta_prefix1": "ejecta_text26_nlayers150_theta{}_nism{}_p205_epse01_epsb001_epst1_".format(theta,nism),
    # #          "ejecta_prefix2": "ejecta_WG_text26_nlayers150_theta{}_nism{}_p205_epse01_epsb001_epst1_".format(theta,nism),
    # #          "res_dir_grb": "/media/vsevolod/Data/KentaData/SFHo_135_135_150m_11//afterglow/testdata/", "time_jet": 180,
    # #          "jet_prefix": "jet_G_nism{}_thetaobs{}_epse007079_epsb000525_p216_epst00_nlayers180_".format(nism,theta)
    # #          },
    # #         {"res_dir_kn": "/media/vsevolod/Data/KentaData/SFHo_135_135_150m_11//afterglow/testdata/", "time_kn": 220,
    # #          "ejecta_prefix1": "ejecta_text26_nlayers150_theta{}_nism{}_p205_epse01_epsb001_epst1_".format(theta,nism),
    # #          "ejecta_prefix2": "ejecta_WG_text26_nlayers150_theta{}_nism{}_p205_epse01_epsb001_epst1_".format(theta,nism),
    # #          "res_dir_grb": "/media/vsevolod/Data/KentaData/SFHo_135_135_150m_11//afterglow/testdata/", "time_jet": 220,
    # #          "jet_prefix": "jet_G_nism{}_thetaobs{}_epse007079_epsb000525_p216_epst00_nlayers180_".format(nism,theta)
    # #          }
    # #     ],
    # #     [
    # #         {"res_dir_kn": "/media/vsevolod/Data/KentaData/SFHo_135_135_150m_11//afterglow/testdata/", "time_kn": 240,
    # #          "ejecta_prefix1": "ejecta_text26_nlayers150_theta{}_nism{}_p205_epse01_epsb001_epst1_".format(theta, nism),
    # #          "ejecta_prefix2": "ejecta_WG_text26_nlayers150_theta{}_nism{}_p205_epse01_epsb001_epst1_".format(theta, nism),
    # #          "res_dir_grb": "/media/vsevolod/Data/KentaData/SFHo_135_135_150m_11//afterglow/testdata/", "time_jet": 240,
    # #          "jet_prefix": "jet_G_nism{}_thetaobs{}_epse007079_epsb000525_p216_epst00_nlayers180_".format(nism, theta)
    # #          },
    # #         {"res_dir_kn": "/media/vsevolod/Data/KentaData/SFHo_135_135_150m_11//afterglow/testdata/", "time_kn": 260,
    # #          "ejecta_prefix1": "ejecta_text26_nlayers150_theta{}_nism{}_p205_epse01_epsb001_epst1_".format(theta, nism),
    # #          "ejecta_prefix2": "ejecta_WG_text26_nlayers150_theta{}_nism{}_p205_epse01_epsb001_epst1_".format(theta, nism),
    # #          "res_dir_grb": "/media/vsevolod/Data/KentaData/SFHo_135_135_150m_11//afterglow/testdata/", "time_jet": 260,
    # #          "jet_prefix": "jet_G_nism{}_thetaobs{}_epse007079_epsb000525_p216_epst00_nlayers180_".format(nism, theta)
    # #          },
    # #         {"res_dir_kn": "/media/vsevolod/Data/KentaData/SFHo_135_135_150m_11//afterglow/testdata/", "time_kn": 290,
    # #          "ejecta_prefix1": "ejecta_text26_nlayers150_theta{}_nism{}_p205_epse01_epsb001_epst1_".format(theta, nism),
    # #          "ejecta_prefix2": "ejecta_WG_text26_nlayers150_theta{}_nism{}_p205_epse01_epsb001_epst1_".format(theta, nism),
    # #          "res_dir_grb": "/media/vsevolod/Data/KentaData/SFHo_135_135_150m_11//afterglow/testdata/", "time_jet": 290,
    # #          "jet_prefix": "jet_G_nism{}_thetaobs{}_epse007079_epsb000525_p216_epst00_nlayers180_".format(nism, theta)
    # #          }
    # #     ],
    # #     [
    # #         {"res_dir_kn": "/media/vsevolod/Data/KentaData/SFHo_135_135_150m_11//afterglow/testdata/", "time_kn": 320,
    # #          "ejecta_prefix1": "ejecta_text26_nlayers150_theta{}_nism{}_p205_epse01_epsb001_epst1_".format(theta, nism),
    # #          "ejecta_prefix2": "ejecta_WG_text26_nlayers150_theta{}_nism{}_p205_epse01_epsb001_epst1_".format(theta, nism),
    # #          "res_dir_grb": "/media/vsevolod/Data/KentaData/SFHo_135_135_150m_11//afterglow/testdata/", "time_jet": 320,
    # #          "jet_prefix": "jet_G_nism{}_thetaobs{}_epse007079_epsb000525_p216_epst00_nlayers180_".format(nism, theta)
    # #          },
    # #         {"res_dir_kn": "/media/vsevolod/Data/KentaData/SFHo_135_135_150m_11//afterglow/testdata/", "time_kn": 360,
    # #          "ejecta_prefix1": "ejecta_text26_nlayers150_theta{}_nism{}_p205_epse01_epsb001_epst1_".format(theta, nism),
    # #          "ejecta_prefix2": "ejecta_WG_text26_nlayers150_theta{}_nism{}_p205_epse01_epsb001_epst1_".format(theta, nism),
    # #          "res_dir_grb": "/media/vsevolod/Data/KentaData/SFHo_135_135_150m_11//afterglow/testdata/", "time_jet": 360,
    # #          "jet_prefix": "jet_G_nism{}_thetaobs{}_epse007079_epsb000525_p216_epst00_nlayers180_".format(nism, theta)
    # #          },
    # #         {"res_dir_kn": "/media/vsevolod/Data/KentaData/SFHo_135_135_150m_11//afterglow/testdata/", "time_kn": 400,
    # #          "ejecta_prefix1": "ejecta_text26_nlayers150_theta{}_nism{}_p205_epse01_epsb001_epst1_".format(theta, nism),
    # #          "ejecta_prefix2": "ejecta_WG_text26_nlayers150_theta{}_nism{}_p205_epse01_epsb001_epst1_".format(theta, nism),
    # #          "res_dir_grb": "/media/vsevolod/Data/KentaData/SFHo_135_135_150m_11//afterglow/testdata/", "time_jet": 400,
    # #          "jet_prefix": "jet_G_nism{}_thetaobs{}_epse007079_epsb000525_p216_epst00_nlayers180_".format(nism, theta)
    # #          }
    # #     ]
    # # ]


    settings["precompute"] = False

    if settings["grid"]:
        fig = plt.figure(figsize=settings["figsize"])#(figsize=(2 * 4.6, 4.8))
        axes = ImageGrid(fig, 111,  # as in plt.subplot(111)
                         nrows_ncols=(len(tasks_to_plot[0]), len(tasks_to_plot)),
                         **settings["imagegrid"]
                         )
        print("ax[]={}".format(len(axes)))
    else:
        fig, axes = plt.subplots(figsize=settings["figsize"],#(figsize=(2 * 4.6, 4.8),
                                 ncols=len(tasks_to_plot[0]),
                                 nrows=len(tasks_to_plot),
                                 sharex="all",
                                 sharey="all")
        if not hasattr(axes, "__len__"):
            axes = [[axes]]
        if not hasattr(axes[0], "__len__"):
            axes = [[ax] for ax in axes]
        # print("ax[]={} ax[][]={}".format(len(axes), len(axes[0])))


    ix = -1
    iy = -1
    print("len={}".format(len(tasks_to_plot)))
    print("len={}".format(len(tasks_to_plot[0])))

    if settings["precompute"]: tmp_tile = h5py.File(settings["precompute_fpath"], "w")
    else: tmp_tile = h5py.File(settings["precompute_fpath"], "r")

    hashes = []

    for tmp in tasks_to_plot:
        ix += 1
        for task in tmp:
            iy += 1
            print("ix={} iy={}".format(ix, iy))

            if settings["grid"]: ax = axes[iy * len(tasks_to_plot) + ix]
            else: ax = axes[ix][iy]
            plot_dict = tasks_to_plot[ix][iy]

            freq = plot_dict["freq"]
            time = plot_dict["time"]

            hash, s = make_hash(plot_dict, time, freq)            # hash = "ix={} iy={}".format(ix, iy)
            if not hash in hashes: hashes.append(hash)
            else: raise ValueError("Hash:{} already in hashes={}".format(hash,hashes))

            print("ix={} iy={} time={} freq={:.2e} hash={}".format(ix, iy, time, freq, hash))

            if not hash in tmp_tile.keys():
                raise KeyError("hash={} is not in the file \n Missing: {}".format(hash, plot_dict))

            # --- saving ------
            if settings["precompute"]:
                grp = tmp_tile.create_group(hash)
                for key in task.keys():
                    grp.attrs.create(key, task[key])
            else:
                grp = tmp_tile[hash]
                for key in task.keys():
                    print("Warning! task[key] != grp.attrs[key]")

            ''' ------------------------------------ '''
            if settings["precompute"]:
                xc, yc, x1, x2, y1, y2, int_x, int_y, int_zz = \
                    tmp_for_file2(time=time, freq=freq, grp=grp, plot_dict=plot_dict, settings=settings)

            ''' ---- plot size/cm of the ejecta ---- '''
            pb = BPA_METHODS()
            pb.res_dir_kn = plot_dict["res_dir_kn"]
            pb.ejecta_prefix = plot_dict["ejecta_prefix1"]
            pb.res_dir_grb = plot_dict["res_dir_grb"]
            pb.jet_prefix = plot_dict["jet_prefix"]
            if (len(settings["kn_skymap"].keys()) > 0):
                if (not settings["precompute"]):
                    xc, yc, x1, x2, y1, y2, int_x, int_y, int_zz = \
                        tmp_for_file2(time=time, freq=freq, grp=grp, plot_dict=plot_dict, settings=settings, kn_skymap=True)
                tmp = copy.deepcopy(settings["kn_skymap"])
                # xc, yc, x1, x2, y1, y2, int_x, int_y, int_zz = tmp_for_file2(time=time, freq=freq, grp=grp,
                #                                                              plot_dict=plot_dict,
                #                                                              settings=settings)
                # assert x2 > x1
                if "cm" in tmp.keys(): ax.plot(xc, yc, **tmp["cm"])# color='red', marker="o")
                if "ysize" in tmp.keys(): ax.errorbar([xc, xc], [y1, y2], xerr=[int_x.max() / 10, int_x.max() / 10], **tmp["ysize"])
                if "xsize" in tmp.keys(): ax.errorbar([x1, x2], [yc, yc], yerr=[int_y.max() / 10, int_y.max() / 10], **tmp["xsize"])
                # --------------------
                if "pcolormesh" in tmp.keys():
                    # levels = np.geomspace(int_ration[~np.isnan(int_ration)].min(), int_ration[~np.isnan(int_ration)].max(), 50)
                    int_zz = int_zz.T
                    int_zz[~np.isfinite(int_zz)] = float(tmp["pcolormesh"]["isnan"])
                    if (len(tmp["smooth"].keys()) > 0):
                        int_zz = pb.smooth_interpolated_skymap_with_gaussian_kernel(i_zz=int_zz,  type=tmp["smooth"]["type"])
                    im = plot_pcolomesh(ax=ax, task_dict=tmp["pcolormesh"], int_x=int_x, int_y=int_y, int_zz=int_zz,
                                        outfname=settings["figpath"] + (task["res_dir_kn"].replace("/", "_") + task["ejecta_prefix1"]).lower() )
                    if tmp["pcolormesh"]["set_rasterized"]: ax.set_rasterized(im)

            ''' ---- plot size/cm of the grb ---- '''
            if (len(settings["grb_skymap"].keys()) > 0):
                if (not settings["precompute"]):
                    xc, yc, x1, x2, y1, y2, int_x, int_y, int_zz = \
                        tmp_for_file2(time=time, freq=freq, grp=grp, plot_dict=plot_dict, settings=settings, grb_skymap=True)
                tmp = copy.deepcopy(settings["grb_skymap"])
                # xc, yc, x1, x2, y1, y2, int_x, int_y, int_zz = tmp_for_file2(time=time, freq=freq, grp=grp,
                #                                                              plot_dict=plot_dict,
                #                                                              settings=settings)
                if "cm" in tmp.keys():  ax.plot(xc, yc, **tmp["cm"])
                if "ysize" in tmp.keys(): ax.errorbar([xc, xc], [y1, y2],
                                                      xerr=[int_x.max() / 10, int_x.max() / 10], **tmp["ysize"])
                if "xsize" in tmp.keys(): ax.errorbar([x1, x2], [yc, yc],
                                                      yerr=[int_y.max() / 10, int_y.max() / 10], **tmp["xsize"])
                # --------------------
                if "pcolormesh" in tmp.keys():
                    # int_zz = int_zz.T
                    # int_zz[~np.isfinite(int_zz)] = 0.
                    # im = plot_pcolomesh(ax=ax, task_dict=tmp["pcolormesh"], int_x=int_x, int_y=int_y, int_zz=int_zz)
                    # if tmp["pcolormesh"]["set_rasterized"]: ax.set_rasterized(im)
                    int_zz = int_zz.T
                    int_zz[~np.isfinite(int_zz)] = float(tmp["pcolormesh"]["isnan"])
                    if (len(tmp["smooth"].keys()) > 0):
                        int_zz = pb.smooth_interpolated_skymap_with_gaussian_kernel(i_zz=int_zz,
                                                                                    type=tmp["smooth"]["type"])
                    im = plot_pcolomesh(ax=ax, task_dict=tmp["pcolormesh"], int_x=int_x, int_y=int_y, int_zz=int_zz,
                                        outfname=settings["figpath"] + (task["res_dir_kn"].replace("/", "_") + task["ejecta_prefix1"]).lower())
                    if tmp["pcolormesh"]["set_rasterized"]: ax.set_rasterized(im)

            ''' --- plot size/cm of the kn with grb (using previous two) ---'''
            if len(settings["kn_grb_skymap"].keys()) > 0:
                if (not settings["precompute"]):
                    xc, yc, x1, x2, y1, y2, int_x, int_y, int_zz = \
                        tmp_for_file2(time=time, freq=freq, grp=grp, plot_dict=plot_dict, settings=settings, kn_grb_skymap=True)
                # assert (len(settings["grb_skymap"].keys()) > 0)
                # assert (len(settings["kn_skymap"].keys()) > 0)
                tmp = copy.deepcopy(settings["kn_grb_skymap"])
                # xc, yc, x1, x2, y1, y2, int_x, int_y, int_zz = \
                #     tmp_for_file2(time=time, freq=freq, grp=grp, plot_dict=plot_dict, settings=settings)
                if "cm" in tmp.keys(): ax.plot(xc, yc, **tmp["cm"])
                if "ysize" in tmp.keys(): ax.errorbar([xc, xc], [y1, y2],
                                                      xerr=[int_x.max() / 10, int_x.max() / 10], **tmp["ysize"])
                if "xsize" in tmp.keys(): ax.errorbar([x1, x2], [yc, yc],
                                                      yerr=[int_y.max() / 10, int_y.max() / 10], **tmp["xsize"])
                # --------------------
                if "pcolormesh" in tmp.keys():
                    # levels = np.geomspace(int_ration[~np.isnan(int_ration)].min(), int_ration[~np.isnan(int_ration)].max(), 50)
                    # int_zz_pj = int_zz.T
                    # int_zz_pj[~np.isfinite(int_zz_pj)] = 0.
                    # im = plot_pcolomesh(ax=ax, task_dict=tmp["pcolormesh"], int_x=int_x, int_y=int_y, int_zz=int_zz_pj)
                    # if tmp["pcolormesh"]["set_rasterized"]: ax.set_rasterized(im)
                    int_zz = int_zz.T
                    int_zz[~np.isfinite(int_zz)] = float(tmp["pcolormesh"]["isnan"])
                    if (len(tmp["smooth"].keys()) > 0):
                        int_zz = pb.smooth_interpolated_skymap_with_gaussian_kernel(i_zz=int_zz,
                                                                                    type=tmp["smooth"]["type"])
                    im = plot_pcolomesh(ax=ax, task_dict=tmp["pcolormesh"], int_x=int_x, int_y=int_y, int_zz=int_zz,
                                        outfname=settings["figpath"] + (task["res_dir_kn"].replace("/", "_") + task["ejecta_prefix1"]).lower())
                    if tmp["pcolormesh"]["set_rasterized"]: ax.set_rasterized(im)
            ''' ---- plot size/cm of the ejecta with interaction ---- '''
            pb_w = BPA_METHODS()
            pb_w.res_dir_kn = plot_dict["res_dir_kn"]
            pb_w.ejecta_prefix = plot_dict["ejecta_prefix2"]
            pb_w.res_dir_grb = plot_dict["res_dir_grb"]
            pb_w.jet_prefix = plot_dict["jet_prefix"]
            if len(settings["kn_w_skymap"].keys()) > 0:
                if (not settings["precompute"]):
                    xc, yc, x1, x2, y1, y2, int_x, int_y, int_zz = \
                        tmp_for_file2(time=time, freq=freq, grp=grp, plot_dict=plot_dict, settings=settings, kn_w_skymap=True)
                tmp = copy.deepcopy(settings["kn_w_skymap"])
                # xc, yc, x1, x2, y1, y2, int_x, int_y, int_zz = \
                #     tmp_for_file2(time=time, freq=freq, grp=grp, plot_dict=plot_dict, settings=settings)
                if "cm" in tmp.keys(): ax.plot(xc, yc, **tmp["cm"])
                if "ysize" in tmp.keys(): ax.errorbar([xc, xc], [y1, y2], xerr=[int_x.max() / 10, int_x.max() / 10], **tmp["ysize"])
                if "xsize" in tmp.keys(): ax.errorbar([x1, x2], [yc, yc], yerr=[int_y.max() / 10, int_y.max() / 10], **tmp["xsize"])
                # --------------------
                if "pcolormesh" in tmp.keys():
                    # int_zz = int_zz.T
                    # int_zz[~np.isfinite(int_zz)] = 0.
                    # im = plot_pcolomesh(ax=ax, task_dict=tmp["pcolormesh"], int_x=int_x, int_y=int_y, int_zz=int_zz)
                    # if tmp["pcolormesh"]["set_rasterized"]: ax.set_rasterized(im)
                    int_zz = int_zz.T
                    int_zz[~np.isfinite(int_zz)] = float(tmp["pcolormesh"]["isnan"])
                    if (len(tmp["smooth"].keys()) > 0):
                        int_zz = pb.smooth_interpolated_skymap_with_gaussian_kernel(i_zz=int_zz,
                                                                                    type=tmp["smooth"]["type"])
                    im = plot_pcolomesh(ax=ax, task_dict=tmp["pcolormesh"], int_x=int_x, int_y=int_y, int_zz=int_zz,
                                        outfname=settings["figpath"] + (task["res_dir_kn"].replace("/", "_") + task["ejecta_prefix1"]).lower())
                    if tmp["pcolormesh"]["set_rasterized"]: ax.set_rasterized(im)
            ''' --- plot relation --- '''
            if len(settings["kn_skymap_ratio"].keys()) > 0:
                if (not settings["precompute"]):
                    xc, yc, x1, x2, y1, y2, int_x, int_y, int_zz = \
                        tmp_for_file2(time=time, freq=freq, grp=grp, plot_dict=plot_dict, settings=settings, kn_skymap_ratio=True)
                assert (len(settings["kn_skymap"].keys()) > 0)
                assert (len(settings["kn_w_skymap"].keys()) > 0)
                tmp = copy.deepcopy(settings["kn_skymap_ratio"])
                # _, _, _, _, _, _, int_x, int_y, int_zz = tmp_for_file2(time=time, freq=freq, grp=grp,
                #                                                              plot_dict=plot_dict,
                #                                                              settings=settings)

                # norm = Normalize(int_ration[np.isfinite(int_ration)].min(), int_ration[np.isfinite(int_ration)].max())
                # # levels = np.geomspace(int_ration[~np.isnan(int_ration)].min(), int_ration[~np.isnan(int_ration)].max(), 50)
                # cmap = cm.get_cmap('viridis')
                # int_ration = int_ration.T
                # int_ration[~np.isfinite(int_ration)] = 0.
                # im = ax.pcolormesh(int_x_w, int_y_w, int_ration, norm=norm, cmap=cmap, alpha=1.0)
                # im.set_rasterized(True)
                # ax.set_facecolor(cmap(0.0))

                # vmin = tmp["pcolormesh"]["vmin"]
                # vmax = tmp["pcolormesh"]["vmax"]
                # vcenter = tmp["pcolormesh"]["vcenter"]
                # levels = MaxNLocator(nbins=40).tick_values(int_ration.min(), int_ration.max())
                # # levels = MaxNLocator(nbins=40).tick_values(-5, 1)
                # cmap = plt.get_cmap(tmp["pcolormesh"]["cmap"])
                # norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
                # norm = colors.TwoSlopeNorm(vmin=vmin, vcenter=vcenter, vmax=vmax)
                # # im = aux_ax2.contourf(theta, r, values, levels=np.logspace(np.log10(1e-4),np.log10(1e0),20), cmap=cmap, norm=norm)
                # # im = aux_ax2.contourf(theta, r, values, levels=np.linspace(val.min(),val.max(),20), cmap=cmap, norm=norm)
                # im = ax.pcolormesh(int_x_w, int_y_w, int_ration, cmap=cmap, norm=norm, alpha=1.0)
                # im.set_rasterized(True)

                # im = plot_pcolomesh(ax=ax, task_dict=tmp["pcolormesh"], int_x=int_x, int_y=int_y, int_zz=int_zz)
                # if tmp["pcolormesh"]["set_rasterized"]: ax.set_rasterized(im)
                # int_zz = int_zz.T
                int_zz[~np.isfinite(int_zz)] = float(tmp["pcolormesh"]["isnan"])
                if (len(tmp["smooth"].keys()) > 0):
                    int_zz = pb.smooth_interpolated_skymap_with_gaussian_kernel(i_zz=int_zz, type=tmp["smooth"]["type"])
                im = plot_pcolomesh(ax=ax, task_dict=tmp["pcolormesh"], int_x=int_x, int_y=int_y, int_zz=int_zz,
                                    outfname=settings["figpath"] + (task["res_dir_kn"].replace("/", "_") + task["ejecta_prefix1"]).lower())
                if tmp["pcolormesh"]["set_rasterized"]: ax.set_rasterized(im)



            ax.set_yscale("linear")
            ax.set_xscale("linear")
            ax.minorticks_on()
            ax.tick_params(  # axis='both',
                **settings["tick_params"]
            )
            if "xlim" in settings.keys() and len(settings["xlim"]) == 2:
                ax.set_xlim(settings["xlim"])
            else:
                ax.set_xlim(int_x.min() * settings["xlim_fac"], int_x.max() * settings["xlim_fac"])
            if "ylim" in settings.keys() and len(settings["ylim"]) == 2:
                ax.set_ylim(settings["ylim"])
            else:
                ax.set_ylim(int_y.min() * settings["ylim_fac"], int_y.max() * settings["ylim_fac"])
            # if settings["grid"]:
            #     ax.set_xlim(min(settings["xlim"]) * 1.1, max(settings["xlim"]) * 1.1)
            #     ax.set_ylim(min(settings["ylim"]) * 1.1, max(settings["ylim"]) * 1.1)

            # ax.tick_params(axis='both', which='both', labelleft=True,
            #                labelright=False, tick1On=True, tick2On=True,
            #                labelsize=11, color="white",
            #                direction='in',
            #                bottom=True, top=True, left=True, right=True)
            # ax.set_title(r"$t_{\rm obs}=$" + r"{:d} day".format(time))
            # ax.set_title(line["label"])

            ax.axhline(y=0, linestyle=':', linewidth=0.4, color="black")
            ax.axvline(x=0, linestyle=':', linewidth=0.4, color="black")






            def plot_title(settings, ax, time):
                if settings["title"]["title"] == "time_fluxratio":
                    times = pb.get_ej_skymap_times()
                    fluxes = pb.get_ej_skymap_totfluxes(freq=freq)
                    fluxes_w = pb_w.get_ej_skymap_totfluxes(freq=freq)
                    fluxes_ration = fluxes_w / fluxes
                    idx = find_nearest_index(times, time * cgs.day)
                    # ax.set_title("t={:.0f} d. F/F={}".format(time, fluxes_ration[idx]))
                    ax.set_title("t={:.0f} d. ".format(time) +
                                 r"$F_{\nu}^{\rm w}/F_{\nu}^{\rm w/o}=$" +
                                 "{:.2f}".format(fluxes_ration[idx]))
                elif settings["title"]["title"] == "time":
                    ax.set_title("t={:.0f} day ".format(time))
                elif settings["title"]["title"] == "tasktitle":
                    ax.set_title(task["title"])
                else:
                    ax.set_title(task["title"])
            # if (len(settings["title"].keys()) > 0):
            #     if ("where" in settings["title"].keys()):
            #         if settings["title"]["where"] == "first_row":
            #             if ix==0:
            #                 plot_title(settings, ax)
            #     else:
            #         plot_title(settings, ax)

            # def plot_txt(settings, ax):
            def plot_text(settings, ax, time):
                __text_dict = copy.deepcopy(_text_dict)
                if __text_dict["text"] == "<time>":
                    __text_dict["text"] = "{:.0f} d.".format(time)
                elif __text_dict["text"] == "<tasktext>":
                    __text_dict["text"] = task["text"]["text"]
                plot_text2(ax, __text_dict)


            if ("title" in settings.keys()) and (len(settings["title"].keys()) > 0):
                if ("where" in settings["title"].keys()):
                    if settings["title"]["where"] == "first_row":
                        if iy == 0: plot_title(settings, ax=ax, time=time)
                    if settings["title"]["where"] == "first_column":
                        if ix == 0: plot_title(settings, ax=ax, time=time)
                    else:
                        plot_title(settings, ax, time=time)
                else:
                    plot_title(settings, ax, time=time)
            if "text" in settings.keys():
                plot_text2(ax, settings["text"])
            if "texts" in settings.keys():
                for _text_dict in settings["texts"]:
                    if ("where" in _text_dict.keys()):
                        if _text_dict["where"] == "first_row":
                            if iy == 0:
                                plot_text(_text_dict, ax, time=time)
                        if _text_dict["where"] == "first_column":
                            if ix == 0:
                                plot_text(_text_dict, ax, time=time)
                    else:
                        plot_text(_text_dict, ax, time=time)

                        # __text_dict = copy.deepcopy(_text_dict)
                        # if __text_dict["text"] == "<time>":
                        #     __text_dict["text"] = "{:.0f} d.".format(time)
                        # elif __text_dict["text"] == "<tasktext>":
                        #     __text_dict["text"] = task["text"]["text"]
                        # if __text_dict["where"] == "first_row":
                        #     if ix==0: plot_text2(ax, __text_dict)
                        # if __text_dict["where"] == "first_column":
                        #     if iy == 0: plot_text2(ax, __text_dict)
                    #
                    # if settings["title"]["where"] == "first_row":
                    #     if ix==0:
                    #         for _text_dict in settings["texts"]:
                    #             plot_text2(ax, _text_dict)
                    # elif settings["title"]["where"] == "first_column":
                    #     if iy==0:
                    #         for _text_dict in settings["texts"]:
                    #             plot_text2(ax, _text_dict)
                # else:
                #     for _text_dict in settings["texts"]:
                #         plot_text2(ax, _text_dict)

                # assert len(settings["texts"]) > 0
                # for _text_dict in settings["texts"]:
                #     plot_text2(ax, _text_dict)


            # x = 1
        iy = -1
        # fwhm = np.array(fwhm)
        # xcs = np.array(xcs)
        #
        # axes[0].plot(times, fwhm,  # drawstyle='steps',
        #              color=line["color"], lw=line["lw"], label=line["label"], ls=line["ls"])
        # axes[1].plot(times, xcs,  # drawstyle='steps',
        #              color=line["color"], lw=line["lw"], label=line["label"], ls=line["ls"])

    tmp_tile.close()

    if settings["grid"] and len(settings["cbar"].keys()) > 0:
        cbar = ax.cax.colorbar(im)
        ax.cax.toggle_label(True)
        ax.cax.set_title(settings["cbar"]["title"])#(r'$I_{\nu}^{\rm w}/I_{\nu}^{\rm w/o}$')
        # cbar.ax.set_ylabel(r'$I_{\nu}^{\rm w}/I_{\nu}^{\rm w/o}$', fontsize = 10)
        ax.cax.tick_params(  # axis='both',
            which='both',  # labelleft=True,
            # labelright=False,
            tick1On=True, tick2On=True,
            labelsize=12,
            direction='in',
            color="white"
        )

    plt.subplots_adjust(**settings["subplots_adjust"])
    fig.add_subplot(111, frame_on=False)
    plt.tick_params(labelcolor="none", bottom=False, left=False)
    plt.xlabel(settings["xlabel"], fontsize=settings["fontsize"])
    plt.ylabel(settings["ylabel"], fontsize=settings["fontsize"])
    figname = settings["figname"]  # make_prefix_for_pars(pars)
    if settings["save_figs"]:
        print(r"Saving:\n{}".format(settings["figpath"] + figname + ".png"))
        plt.savefig(settings["figpath"] + figname + ".png", dpi=256)
    if settings["save_pdf"]:
        print(r"Saving:\n{}".format(settings["paperpath"] + figname + ".pdf"))
        plt.savefig(settings["paperpath"] + figname + ".pdf")
    if settings["show_figs"]: plt.show()

def plot_skymap_properties_evolution(tasks_to_plot, settings):

    settings["precompute"] = False

    if settings["grid"]:
        fig = plt.figure(figsize=settings["figsize"])  # (figsize=(2 * 4.6, 4.8))
        axes = ImageGrid(fig, 111,  # as in plt.subplot(111)
                         nrows_ncols=(len(tasks_to_plot[0]), len(tasks_to_plot)),
                         **settings["imagegrid"]
                         )
        print("ax[]={}".format(len(axes)))
    else:
        fig, axes = plt.subplots(figsize=settings["figsize"],  # (figsize=(2 * 4.6, 4.8),
                                 ncols=len(tasks_to_plot[0]),
                                 nrows=len(tasks_to_plot),
                                 sharex="all",
                                 sharey="all")
        if not hasattr(axes, "__len__"):
            axes = [[axes]]
        if not hasattr(axes[0], "__len__"):
            axes = [[ax] for ax in axes]
        # print("ax[]={} ax[][]={}".format(len(axes), len(axes[0])))
    ix = -1
    iy = -1
    print("len={}".format(len(tasks_to_plot)))
    print("len={}".format(len(tasks_to_plot[0])))

    hashes = [{}]
    hashes2 = {}
    for tmp in tasks_to_plot:
        ix += 1
        for sim_name in tmp:
            iy += 1
            print("ix={} iy={}".format(ix, iy))

            if settings["grid"]: ax = axes[iy * len(tasks_to_plot) + ix]
            else: ax = axes[ix][iy]
            plot_dicts = tasks_to_plot[ix][iy]
            if not hasattr(plot_dicts, "__len__"):
                raise IOError()

            inset_axes = None
            if (len(settings["show_subplot"].keys()) > 0):
                from mpl_toolkits.axes_grid.inset_locator import inset_axes
                inset_axes = inset_axes(ax, **settings["show_subplot"]["inset_axes"])
                inset_axes.minorticks_on()
                inset_axes.tick_params( **settings["tick_params"] )
                inset_axes.set_xlim(*settings["show_subplot"]["xlim"])
                inset_axes.set_ylim(*settings["show_subplot"]["ylim"])
                inset_axes.set_xscale(settings["show_subplot"]["xscale"])
                inset_axes.set_yscale(settings["show_subplot"]["yscale"])

            # ============================================================================
            for plot_dict in plot_dicts:

                pb = BPA_METHODS()
                pb.res_dir_kn = plot_dict["res_dir_kn"]
                pb.ejecta_prefix = plot_dict["ejecta_prefix1"]
                pb.res_dir_grb = plot_dict["res_dir_grb"]
                pb.jet_prefix = plot_dict["jet_prefix"]
                freq = float(plot_dict["freq"])
                vals = None
                if settings["kn_specidx"]:
                    times = pb.get_ej_lc_times() / cgs.day
                    vals = pb.get_ej_lc_spec_idx(freq=freq)
                else:
                    # --------------------------------------------------
                    mode = "kn_skymap"
                    if (len(settings[mode].keys()) > 3):
                        times, _vals = read_precompute_file_get_gals(settings, plot_dict, hashes=hashes)
                        for key in _vals[mode].keys():
                            if (len(_vals[mode][key]) > 0):
                                vals = _vals[mode][key]
                        # if settings["plot_actual_values"]:
                        #     ax.plot(times, vals, **plot_dict["line"])
                    # --------------------------------------------------
                    mode = "grb_skymap"
                    if( len(settings[mode].keys()) > 3 ):
                        times, _vals = read_precompute_file_get_gals(settings, plot_dict, hashes=hashes)
                        for key in _vals[mode].keys():
                            if (len(_vals[mode][key]) > 0):
                                vals = _vals[mode][key]
                        # if settings["plot_actual_values_grb"]:
                        #     ax.plot(times, vals, color="black", ls="-")
                        # if settings["plot_actual_values_grb"]:
                        #     inset_axes.plot(times, vals, color="black", ls="-")
                    # --------------------------------------------------
                    mode = "kn_grb_skymap"
                    if( len(settings[mode].keys()) > 3 ):
                        times, _vals = read_precompute_file_get_gals(settings, plot_dict, hashes=hashes)
                        for key in _vals[mode].keys():
                            if (len(_vals[mode][key]) > 0):
                                vals = _vals[mode][key]
                                if settings[mode]["substruct_from_grb_val"]:
                                    vals = (_vals["grb_skymap"][key] - _vals[mode][key])
                                if settings[mode]["substruct_from_grb_val_norm"]:
                                    vals = (_vals["grb_skymap"][key] - _vals[mode][key]) / _vals["grb_skymap"][key]

                        # if settings["plot_actual_values"]:
                        #     ax.plot(times, vals, **plot_dict["line"])
                    # --------------------------------------------------

                if settings["plot_actual_values_grb"]:
                    ax.plot(times, vals, color="black", ls="-")
                if ((len(settings["show_subplot"].keys()) > 0) and settings["plot_actual_values_grb"]):
                    inset_axes.plot(times, vals, color="black", ls="-")
                if settings["plot_actual_values"]:
                    ax.plot(times, vals, **plot_dict["line"])
                if( settings["plot_actual_values"] and (len(settings["show_subplot"].keys()) > 0)):
                    inset_axes.plot(times, vals, **plot_dict["line"])

                if (len(settings["show_lc_peak1"].keys()) > 0):
                    tmp = settings["show_lc_peak1"]

                    if "marker" in plot_dict.keys():
                        marker = plot_dict["marker"]
                    else:
                        marker = tmp["marker"]["marker"]
                    # sim = str(plot_dict["sim"])

                    fmax = pb.get_ej_lc_totalflux(freq).max()
                    tmax = pb.get_ej_lc_times()[find_nearest_index(pb.get_ej_lc_totalflux(freq), fmax)]
                    tmax /= cgs.day

                    # tmax_sky = pb.get_ej_skymap_times()[find_nearest_index(pb.get_ej_skymap_times(), tmax)]

                    if tmax > times.max():
                        raise ValueError("{} tmax={} > times.max()={}".format(plot_dict["res_dir_kn"], tmax, times.max()))
                    if tmax < times.min():
                        raise ValueError("{} times={} < times.min()={}".format(plot_dict["res_dir_kn"], tmax, times.min()))

                    val_int = interpolate.interp1d(times, vals, kind="linear")(tmax)

                    # if "substruct_from_grb_val" in tmp.keys() and tmp["substruct_from_grb_val"]:
                    #     val_int_grb = interpolate.interp1d(pb.get_jet_lc_times(), pb.get_jet_lc_totalflux(freq), kind="linear")(tmax*cgs.day)
                    #     val_int = val_int_grb - val_int

                    ax.plot(tmax, val_int, ls="none",marker=marker,color=plot_dict["line"]["color"], fillstyle=plot_dict["fillstyle"])

                    # ax2 = ax.twinx()
                    # ax2.plot(pb.get_ej_lc_times()/cgs.day,pb.get_ej_lc_spec_idx(freq=freq),**plot_dict["line"])

                    if (len(settings["show_subplot"].keys()) > 0):
                        # if settings["plot_actual_values"]:
                        #     inset_axes.plot(times, vals, **plot_dict["line"])
                        # inset_axes.plot(times, vals)
                        # _t =  pb.get_ej_lc_times()/cgs.day
                        # _val = pb.get_ej_lc_spec_idx(freq=freq)
                        # inset_axes.plot(_t, _val)
                        ax.plot(tmax, val_int, ls="none", marker=marker, color=plot_dict["line"]["color"],fillstyle=plot_dict["fillstyle"])
                        inset_axes.plot(tmax, val_int, ls="none", marker=marker, color=plot_dict["line"]["color"],fillstyle=plot_dict["fillstyle"])

                if (len(settings["show_grb_lc_peak_time"].keys()) > 0):
                    tmp = settings["show_grb_lc_peak_time"]
                    fmax = pb.get_jet_lc_totalflux(freq).max()
                    tmax = pb.get_jet_lc_times()[find_nearest_index(pb.get_jet_lc_totalflux(freq), fmax)]
                    tmax /= cgs.day
                    ax.axvline(x=tmax, color=tmp["color"], linestyle=tmp["linestyle"])
                    if (len(settings["show_subplot"].keys()) > 0):
                        inset_axes.axvline(x=tmax, color=tmp["color"], linestyle=tmp["linestyle"])

                if (len(settings["show_specidx_min1"].keys())):
                    tmp = settings["show_specidx_min1"]
                    if "marker" in plot_dict.keys():
                        marker = plot_dict["marker"]
                    else:
                        marker = tmp["marker"]["marker"]
                    _times = pb.get_ej_lc_times()
                    _values = pb.get_ej_lc_spec_idx(freq=freq)
                    vmin = _values.min()
                    tmin = pb.get_ej_lc_times()[find_nearest_index(_values, vmin)]
                    tmin /= cgs.day
                    if tmin > times.max():
                        # raise ValueError("{} tmin={} > times.max()={}".format(plot_dict["res_dir_kn"], tmin, times.max()))
                        print("ERROR! {} tmin={} > times.max()={}".format(plot_dict["res_dir_kn"], tmin, times.max()))
                        tmin = times.max()
                    if tmin < times.min():
                        # raise ValueError("{} tmin={} < times.min()={}".format(plot_dict["res_dir_kn"], tmin, times.min()))
                        print("ERROR! {} tmin={} < times.min()={}".format(plot_dict["res_dir_kn"], tmin, times.min()))
                        tmin = times.min()
                    val_int = interpolate.interp1d(times, vals, kind="linear")(tmin)
                    ax.plot(tmin, val_int, ls="none", marker=marker, color=plot_dict["line"]["color"], fillstyle=plot_dict["fillstyle"])
                    if (len(settings["show_subplot"].keys()) > 0):
                        if settings["plot_actual_values"]: inset_axes.plot(times, vals, **plot_dict["line"])
                        # inset_axes.plot(times, vals)
                        # _t =  pb.get_ej_lc_times()/cgs.day
                        # _val = pb.get_ej_lc_spec_idx(freq=freq)
                        # inset_axes.plot(_t, _val)
                        ax.plot(tmin, val_int, ls="none", marker=marker, color=plot_dict["line"]["color"],fillstyle=plot_dict["fillstyle"])
                        inset_axes.plot(tmin, val_int, ls="none", marker=marker, color=plot_dict["line"]["color"],fillstyle=plot_dict["fillstyle"])
            # ===============================================================================

            ax.set_yscale("linear")
            ax.set_xscale("linear")
            ax.minorticks_on()
            ax.tick_params(  # axis='both',
                **settings["tick_params"]
            )
            if "xlim" in settings.keys() and len(settings["xlim"]) == 2:
                ax.set_xlim(settings["xlim"])
            else:
                ax.set_xlim(times.min() * settings["xlim_fac"], times.max() * settings["xlim_fac"])
            if "ylim" in settings.keys() and len(settings["ylim"]) == 2:
                ax.set_ylim(settings["ylim"])
            else:
                ax.set_ylim(vals.min() * settings["ylim_fac"], vals.max() * settings["ylim_fac"])

            if "text" in settings.keys():
                plot_text2(ax, settings["text"])
            if "texts" in settings.keys():
                # assert len(settings["texts"]) > 0
                for _text_dict in settings["texts"]:
                    plot_text2(ax, _text_dict)

            if (len(settings["show_subplot"].keys()) > 0):
                inset_axes.minorticks_on()

        iy = -1

    ax = axes[0]
    if hasattr(ax, "__len__"):
        ax = ax[0]
    if not (settings["legend"] is None) and \
            len(settings["legend"].keys()) > 0:

        multilegend = settings["legend"]

        if "lines" in multilegend.keys():
            _tasks = multilegend["lines"]
            for line in _tasks:
                ax.plot([0, 0], [0, 0], **line)
        elif "fill_between" in multilegend.keys():
            _tasks = multilegend["fill_between"]
            for _dict in _tasks:
                ax.fill_between(**_dict)
        else:
            raise NameError("not recognized")
        #     if "label" in line.keys():
        #         _lbs.append(task["line"]["label"])
        #         del task["line"]["label"]
        #
        #     _lls.append(_l)

        han, lab = ax.get_legend_handles_labels()

        # if "obs" in plotdic.keys():
        #     obs_legend = copy.deepcopy(plotdic["legend"])
        #     obs_legend["loc"] = "upper left"
        #     ax.add_artist(ax.legend([han[0]], [lab[0]], **obs_legend))
        #     han, lab = han[1:], lab[1:]

        # ax.add_artist(ax.legend(han[:-1 * len(_tasks)], lab[:-1 * len(_tasks)], **plotdic["legend"]))
        #
        ax.add_artist(ax.legend(han[len(han) - len(_tasks):], lab[len(lab) - len(_tasks):], **multilegend["legend"]))

    multilegend = {
        "lines": [
            {"color": 'green', "ls": '-', "label": "BLh"},  # ,, "zorder":-1},
            {"color": 'blue', "ls": '-', "label": "DD2"},  # ,, "zorder":-1},
            {"color": 'red', "ls": '-', "label": "LS220"},  # ,, "zorder":-1},
            {"color": 'orange', "ls": '-', "label": "SFHo"},  # ,, "zorder":-1},
            {"color": 'magenta', "ls": '-', "label": "SLy4"},  # ,, "zorder":-1},
            # {"color": 'black', "ls": '-', "label": "GRB"}  # ,, "zorder":-1},

        ],
        "legend": {"fancybox": False, "loc": 'upper right',
                   # "bbox_to_anchor": (0.5, 0.7),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
                   "shadow": "False", "ncol": 1, "fontsize": 11,
                   "framealpha": 0., "borderaxespad": 0., "frameon": False}
    }
    # ax = axes[0]
    if not (settings["multilegend1"] is None) and \
            len(settings["multilegend1"].keys()) > 0:
        multilegend = settings["multilegend1"]
        if "lines" in multilegend.keys():
            _tasks = multilegend["lines"]
            for line in _tasks:
                ax.plot([0, 0], [0, 0], **line)
        elif "fill_between" in multilegend.keys():
            _tasks = multilegend["fill_between"]
            for _dict in _tasks:
                ax.fill_between(**_dict)
        else:
            raise NameError("not recognized")
        #     if "label" in line.keys():
        #         _lbs.append(task["line"]["label"])
        #         del task["line"]["label"]
        #
        #     _lls.append(_l)

        han, lab = ax.get_legend_handles_labels()

        # if "obs" in plotdic.keys():
        #     obs_legend = copy.deepcopy(plotdic["legend"])
        #     obs_legend["loc"] = "upper left"
        #     ax.add_artist(ax.legend([han[0]], [lab[0]], **obs_legend))
        #     han, lab = han[1:], lab[1:]

        # ax.add_artist(ax.legend(han[:-1 * len(_tasks)], lab[:-1 * len(_tasks)], **plotdic["legend"]))
        #
        ax.add_artist(ax.legend(han[len(han) - len(_tasks):], lab[len(lab) - len(_tasks):], **multilegend["legend"]))

    multilegend = {
        "lines": [
            {"color": 'black', "ls": '-', "label": "GRB"}  # ,, "zorder":-1},

        ],
        "legend": {"fancybox": False, "loc": 'center left',
                   # "bbox_to_anchor": (0.5, 0.7),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
                   "shadow": "False", "ncol": 1, "fontsize": 11,
                   "framealpha": 0., "borderaxespad": 0., "frameon": False}
    }
    # ax = axes[0]
    if not (settings["multilegend2"] is None) and \
            len(settings["multilegend2"].keys()) > 0:
        multilegend = settings["multilegend2"]
        if "lines" in multilegend.keys():
            _tasks = multilegend["lines"]
            for line in _tasks:
                ax.plot([0, 0], [0, 0], **line)
        elif "fill_between" in multilegend.keys():
            _tasks = multilegend["fill_between"]
            for _dict in _tasks:
                ax.fill_between(**_dict)
        else:
            raise NameError("not recognized")
        #     if "label" in line.keys():
        #         _lbs.append(task["line"]["label"])
        #         del task["line"]["label"]
        #
        #     _lls.append(_l)

        han, lab = ax.get_legend_handles_labels()

        # if "obs" in plotdic.keys():
        #     obs_legend = copy.deepcopy(plotdic["legend"])
        #     obs_legend["loc"] = "upper left"
        #     ax.add_artist(ax.legend([han[0]], [lab[0]], **obs_legend))
        #     han, lab = han[1:], lab[1:]

        # ax.add_artist(ax.legend(han[:-1 * len(_tasks)], lab[:-1 * len(_tasks)], **plotdic["legend"]))
        #
        ax.add_artist(ax.legend(han[len(han) - len(_tasks):], lab[len(lab) - len(_tasks):], **multilegend["legend"]))

    multilegend = {
        "lines": [
            {"color": 'black', "ls": '-', "label": "GRB"}  # ,, "zorder":-1},

        ],
        "legend": {"fancybox": False, "loc": 'center left',
                   # "bbox_to_anchor": (0.5, 0.7),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
                   "shadow": "False", "ncol": 1, "fontsize": 11,
                   "framealpha": 0., "borderaxespad": 0., "frameon": False}
    }
    # ax = axes[0]
    if not (settings["multilegend3"] is None) and \
            len(settings["multilegend3"].keys()) > 0:
        multilegend = settings["multilegend3"]
        if "lines" in multilegend.keys():
            _tasks = multilegend["lines"]
            for line in _tasks:
                ax.plot([0, 0], [0, 0], **line)
        elif "fill_between" in multilegend.keys():
            _tasks = multilegend["fill_between"]
            for _dict in _tasks:
                ax.fill_between(**_dict)
        else:
            raise NameError("not recognized")
        #     if "label" in line.keys():
        #         _lbs.append(task["line"]["label"])
        #         del task["line"]["label"]
        #
        #     _lls.append(_l)

        han, lab = ax.get_legend_handles_labels()

        # if "obs" in plotdic.keys():
        #     obs_legend = copy.deepcopy(plotdic["legend"])
        #     obs_legend["loc"] = "upper left"
        #     ax.add_artist(ax.legend([han[0]], [lab[0]], **obs_legend))
        #     han, lab = han[1:], lab[1:]

        # ax.add_artist(ax.legend(han[:-1 * len(_tasks)], lab[:-1 * len(_tasks)], **plotdic["legend"]))
        #
        ax.add_artist(ax.legend(han[len(han) - len(_tasks):], lab[len(lab) - len(_tasks):], **multilegend["legend"]))

    if ("title" in settings.keys()) and len(settings["title"].keys()) > 0:
        ax.set_title(settings["title"]["title"], fontsize=12)

    plt.subplots_adjust(**settings["subplots_adjust"])
    fig.add_subplot(111, frame_on=False)
    plt.tick_params(labelcolor="none", bottom=False, left=False)
    plt.xlabel(settings["xlabel"], fontsize=settings["fontsize"])
    plt.ylabel(settings["ylabel"], fontsize=settings["fontsize"])
    figname = settings["figname"]  # make_prefix_for_pars(pars)
    if settings["save_figs"]:
        print(r"Saving:\n{}".format(settings["figpath"] + figname + ".png"))
        plt.savefig(settings["figpath"] + figname + ".png", dpi=256)
    if settings["save_pdf"]:
        print(r"Saving:\n{}".format(settings["paperpath"] + figname + ".pdf"))
        plt.savefig(settings["paperpath"] + figname + ".pdf")
    if settings["show_figs"]: plt.show()


def _plot_skymap_with_hists(ax_main, ax_histx, ax_histy, pb, tmp,
                            grid_x, grid_y,
                            int_x, int_y,
                            i_zz_x, i_zz_y, int_zz,
                            xc, yc):
    # assert x2 > x1
    # ax.axvline(x=x1, color='gray', linestyle='dotted')
    # ax.axvline(x=x2, color='gray', linestyle='dotted')
    # ax.axvline(x=xcs_m, color='gray', linestyle='solid')
    if ("cm" in tmp.keys()) and len(tmp["cm"].keys()) > 0:
        ax_main.plot(xc, yc, **tmp["cm"])
        ax_histx.axvline(x=xc, color=tmp["cm"]["color"], linestyle='dashed')
        ax_histy.axhline(y=yc, color=tmp["cm"]["color"], linestyle='dashed')
    if ("ysize" in tmp.keys()) and len(tmp["ysize"].keys()) > 0:
        if (len(tmp["smooth"].keys()) > 0):
            i_zz_y = smooth_interpolated_skymap_with_gaussian_kernel(i_zz=i_zz_y, type=tmp["smooth"]["type"],
                                                                        sigma=tmp["smooth"]["sigma"])
        y1, y2 = get_skymap_fwhm(grid_y, i_zz_y, yc)
        ax_main.errorbar([xc, xc], [y1, y2], xerr=[int_x.max() / 10, int_x.max() / 10], **tmp["ysize"])
        ax_histy.plot(i_zz_y * 1e3, grid_y, lw=1., ls='-', drawstyle='steps', color=tmp["ysize"]["color"])
        ax_histy.axhline(y=y2, color=tmp["ysize"]["color"], linestyle='dotted')
        ax_histy.axhline(y=y1, color=tmp["ysize"]["color"], linestyle='dotted')
    if ("xsize" in tmp.keys()) and len(tmp["xsize"].keys()) > 0:
        if (len(tmp["smooth"].keys()) > 0):
            i_zz_x = smooth_interpolated_skymap_with_gaussian_kernel(i_zz=i_zz_x, type=tmp["smooth"]["type"],
                                                                        sigma=tmp["smooth"]["sigma"])
        x1, x2 = get_skymap_fwhm(grid_x, i_zz_x, xc)
        ax_main.errorbar([x1, x2], [yc, yc], yerr=[int_y.max() / 10, int_y.max() / 10], **tmp["xsize"])
        ax_histx.plot(grid_x, i_zz_x * 1e3, lw=1., ls='-', drawstyle='steps',
                      color=tmp["xsize"]["color"])  # , color='blue')
        ax_histx.axvline(x=x2, color=tmp["xsize"]["color"], linestyle='dotted')
        ax_histx.axvline(x=x1, color=tmp["xsize"]["color"], linestyle='dotted')
    # --------------------
    # if len(tmp["pcolormesh"].keys()) > 0:
    #     # levels = np.geomspace(int_ration[~np.isnan(int_ration)].min(), int_ration[~np.isnan(int_ration)].max(), 50)
    #     int_zz = int_zz.T
    #     int_zz[~np.isfinite(int_zz)] = 0.
    #     if (len(tmp["smooth"].keys()) > 0):
    #         int_zz = pb.smooth_interpolated_skymap_with_gaussian_kernel(i_zz=int_zz, type=tmp["smooth"]["type"], sigma=tmp["smooth"]["sigma"])
    #     plot_pcolomesh(ax=ax_main, task_dict=tmp["pcolormesh"], int_x=int_x, int_y=int_y, int_zz=int_zz)
    if ("pcolormesh" in tmp.keys()) and "pcolormesh" in tmp.keys():
        # levels = np.geomspace(int_ration[~np.isnan(int_ration)].min(), int_ration[~np.isnan(int_ration)].max(), 50)
        int_zz = int_zz.T
        int_zz[~np.isfinite(int_zz)] = float(tmp["pcolormesh"]["isnan"])
        if (len(tmp["smooth"].keys()) > 0):
            int_zz = smooth_interpolated_skymap_with_gaussian_kernel(i_zz=int_zz, type=tmp["smooth"]["type"])
        im = plot_pcolomesh(ax=ax_main, task_dict=tmp["pcolormesh"], int_x=int_x, int_y=int_y, int_zz=int_zz, outfname=None)
        if tmp["pcolormesh"]["set_rasterized"]: ax_main.set_rasterized(im)
        return im
def plot_one_skymap_with_dists(task_to_plot, settings):

    plot_dict = task_to_plot

    fig = plt.figure(figsize=settings["figsize"])#(figsize=(5, 5))

    # Add a gridspec with two rows and two columns and a ratio of 2 to 7 between
    # the size of the marginal axes and the main axes in both directions.
    # Also adjust the subplot parameters for a square plot.
    gs = fig.add_gridspec(2, 2,**settings["gridspec"])

    ax_main = fig.add_subplot(gs[1, 0])
    ax_histx = fig.add_subplot(gs[0, 0], sharex=ax_main)
    ax_histy = fig.add_subplot(gs[1, 1], sharey=ax_main)
    ax_cbar = fig.add_subplot(gs[0, 1], sharey=ax_main)

    if settings["precompute"]:
        print("Initializing precompute file:\n\t{}".format(settings["precompute_fpath"]))
        if os.path.isfile(settings["precompute_fpath"]):
            print("Deleting... {}".format(settings["precompute_fpath"]))
            os.remove(settings["precompute_fpath"])
        tmp_tile = h5py.File(settings["precompute_fpath"], "w")
    else:
        if not os.path.isfile(settings["precompute_fpath"]):
            raise IOError("Precompute file not found\n\t{}".format(settings["precompute_fpath"]))
        tmp_tile = h5py.File(settings["precompute_fpath"], "r")

    time = plot_dict["time"]
    freq = plot_dict["freq"]

    # --- make hash name from skymap names ---
    # pb_grb = BPA_METHODS(plot_dict["workingdir"], readparfileforpaths=True, parfile=plot_dict["grb_parfile"])
    pb1 = PyBlastAfterglow(plot_dict["workingdir"], readparfileforpaths=True, parfile=plot_dict["parfile1"])
    pb2 = PyBlastAfterglow(plot_dict["workingdir"], readparfileforpaths=True, parfile=plot_dict["parfile2"])
    hash, s = make_hash(plot_dict["workingdir"],
                        pb1.GRB.fpath_sky_map,
                        pb1.KN.fpath_sky_map,
                        pb2.KN.fpath_sky_map,
                        time, freq)

    if settings["rerun"]: pb1.run()
    # if settings["rerun"]: pb2.run()

    print("time={} freq={:.2e} hash={}".format(time, freq, hash))

    # --- saving ------
    if settings["precompute"]:
        grp = tmp_tile.create_group(hash)
        for key in plot_dict.keys():
            grp.attrs.create(key, plot_dict[key])
    else:
        grp = tmp_tile[hash]
        for key in plot_dict.keys():
            assert plot_dict[key] == grp.attrs[key]

    ''' ---- plot size/cm of the ejecta ---- '''
    # pb = BPA_METHODS(settings["workingdir"])
    # pb.res_dir_kn = plot_dict["res_dir_kn"]
    # pb.ejecta_prefix = plot_dict["ejecta_prefix1"]
    # pb.res_dir_grb = plot_dict["res_dir_grb"]
    # pb.jet_prefix = plot_dict["jet_prefix"]
    # tot_fluxes = np.zeros_like(time)
    # tot_fluxes_kn = np.zeros_like(time)
    # tot_fluxes_grb = np.zeros_like(time)



    # if (len(settings["kn_skymap"].keys()) > 0) or (len(settings["kn_grb_skymap"].keys()) > 0):
    #     tot_fluxes_kn = pb1.get_ej_skymap_totfluxes(freq=freq)
    #
    # if (len(settings["grb_skymap"].keys()) > 0) or (len(settings["kn_grb_skymap"].keys()) > 0):
    #     tot_fluxes_grb = pb1.get_jet_skymap_totfluxes(freq=freq)

    # plt.delaxes(ax_main)
    # plt.delaxes(ax_histx)
    # plt.delaxes(ax_histy)
    # plt.loglog(pb.get_jet_skymap_times()/cgs.day, tot_fluxes_grb, color="blue", label="grb")
    # plt.loglog(pb.get_ej_skymap_times()/cgs.day, tot_fluxes_kn, color="red", label="kn")
    # plt.loglog(pb.get_ej_skymap_times()/cgs.day, tot_fluxes_kn+tot_fluxes_grb, color="black", label="total")
    # plt.axvline(x=time)
    # plt.ylim(1e-3,1e-1)
    # plt.legend()
    # plt.show()

    # tot_flux = tot_fluxes[find_nearest_index(pb.get_ej_skymap_times(), time)]

    if (len(settings["kn_skymap"].keys()) > 0):
        tmp = copy.deepcopy(settings["kn_skymap"])
        if settings["precompute"]:
            all_x, all_y, all_fluxes \
                = pb1.KN.get_skymap(time=time * cgs.day, freq=freq, verbose=True, remove_mu=True, renormalize=True)

            # plt.delaxes(ax_main)
            # plt.delaxes(ax_histy)
            # plt.delaxes(ax_histx)
            # plt.loglog(pb.get_jet_skymap_times()/cgs.day, pb.get_jet_skymap_totfluxes(freq=freq), color="blue")
            # plt.loglog(pb.get_jet_lc_times()/cgs.day, pb.get_jet_lc_totalflux(freq=freq), ls="--", color="blue")
            # plt.loglog(pb.get_ej_skymap_times()/cgs.day, pb.get_ej_skymap_totfluxes(freq=freq), color="red")
            # plt.loglog(pb.get_ej_lc_times()/cgs.day, pb.get_ej_lc_totalflux(freq=freq), ls="--", color="red")
            # plt.loglog(pb.get_ej_skymap_times()/cgs.day, fnus_tot, ls=":", color="red")
            # plt.show()

            int_x, int_y, int_zz = combine_images(all_x, all_y, all_fluxes, hist_or_int="hist", shells=True, nx=tmp["hist_nx"], ny=tmp["hist_ny"], extend=2)
            grid_y, _, i_zz_y, _ = get_skymap_lat_dist(all_x, all_y, all_fluxes, collapse_axis="x", nx=tmp["hist_nx"], ny=tmp["hist_ny"])
            grid_x, _, i_zz_x, _ = get_skymap_lat_dist(all_x, all_y, all_fluxes, collapse_axis="y", nx=tmp["hist_nx"], ny=tmp["hist_ny"])
            xc, yc = pb1.KN.get_skymap_cm(all_x, all_y, all_fluxes)


            print("all_fluxes=[{:.2e}, {:.2e}]".format(np.array(all_fluxes).min(), np.array(all_fluxes).max()))
            # print("i_zz_y_j=[{:.2e}, {:.2e}]".format(i_zz_y_j.min(),i_zz_y_j.max()))
            # print("i_zz_x_j=[{:.2e}, {:.2e}]".format(i_zz_x_j.min(),i_zz_x_j.max()))
            print("int_zz_j=[{:.2e}, {:.2e}]".format(int_zz.min(), int_zz.max()))
            print("i_zz_x=[{:.2e}, {:.2e}]".format(i_zz_x.min(), i_zz_x.max()))
            print("i_zz_y=[{:.2e}, {:.2e}]".format(i_zz_y.min(), i_zz_y.max()))
            print("i_zz_y=[{:.2e}, {:.2e}]".format(i_zz_y.min(), i_zz_y.max()))

            grp_kn = grp.create_group("kn nx={} ny={}".format(tmp["hist_nx"], tmp["hist_ny"]))
            grp_kn.create_dataset("int_x", data=int_x)
            grp_kn.create_dataset("int_y", data=int_y)
            grp_kn.create_dataset("int_zz", data=int_zz)
            grp_kn.create_dataset("grid_y", data=grid_y)
            grp_kn.create_dataset("i_zz_y", data=i_zz_y)
            grp_kn.create_dataset("grid_x", data=grid_x)
            grp_kn.create_dataset("i_zz_x", data=i_zz_x)
            grp_kn.attrs.create("xc", data=xc)
            grp_kn.attrs.create("yc", data=yc)
        else:
            grp_kn = grp["kn nx={} ny={}".format(tmp["hist_nx"], tmp["hist_ny"])]
            int_x = np.array(np.array(grp_kn["int_x"]))
            int_y = np.array(np.array(grp_kn["int_y"]))
            int_zz = np.array(np.array(grp_kn["int_zz"]))  # plotting
            grid_y = np.array(np.array(grp_kn["grid_y"]))
            i_zz_y = np.array(np.array(grp_kn["i_zz_y"]))
            grid_x = np.array(np.array(grp_kn["grid_x"]))
            i_zz_x = np.array(np.array(grp_kn["i_zz_x"]))
            xc = float(grp_kn.attrs["xc"])
            yc = float(grp_kn.attrs["yc"])
        print("Total flux kn={:.2e} grb={:.2e} x_c={}".format(
                pb1.KN.get_skymap_totfluxes(freq=freq, time=time*cgs.day),
                0.,
                xc
       ))
        # i_zz_x /= i_zz_x.max()
        # i_zz_y /= i_zz_y.max()
        im = _plot_skymap_with_hists(ax_main, ax_histx, ax_histy, pb1, tmp,
                                     grid_x, grid_y, int_x, int_y,
                                     i_zz_x, i_zz_y, int_zz, xc, yc)

        # y1, y2 = pb.get_skymap_fwhm(grid_y, i_zz_y, yc)
        # x1, x2 = pb.get_skymap_fwhm(grid_x, i_zz_x, xc)
        # assert x2 > x1
        # # ax.axvline(x=x1, color='gray', linestyle='dotted')
        # # ax.axvline(x=x2, color='gray', linestyle='dotted')
        # # ax.axvline(x=xcs_m, color='gray', linestyle='solid')
        # if ("cm" in tmp.keys()) and (len(tmp["cm"].keys()) > 0):
        #     ax_main.plot(xc, yc, **tmp["cm"])
        #     ax_histx.axvline(x=xc, color=tmp["cm"]["color"], linestyle='dashed')
        #     ax_histy.axhline(y=yc, color=tmp["cm"]["color"], linestyle='dashed')
        # if ("ysize" in tmp.keys()) and (len(tmp["ysize"].keys()) > 0):
        #     ax_main.errorbar([xc, xc], [y1, y2], xerr=[int_x.max() / 10, int_x.max() / 10], **tmp["ysize"])
        #     if (len(tmp["smooth"].keys()) > 0):
        #         i_zz_y = pb.smooth_interpolated_skymap_with_gaussian_kernel(i_zz=i_zz_y, type=tmp["smooth"]["type"], sigma=tmp["smooth"]["sigma"])
        #     ax_histy.plot(i_zz_y * 1e3, grid_y,  lw=1., ls='-', drawstyle='steps', color=tmp["ysize"]["color"])
        #     ax_histy.axhline(y=y2, color=tmp["ysize"]["color"], linestyle='dotted')
        #     ax_histy.axhline(y=y1, color=tmp["ysize"]["color"], linestyle='dotted')
        # if ("xsize" in tmp.keys()) and (len(tmp["xsize"].keys()) > 0):
        #     ax_main.errorbar([x1, x2], [yc, yc], yerr=[int_y.max() / 10, int_y.max() / 10], **tmp["xsize"])
        #     if (len(tmp["smooth"].keys()) > 0):
        #         i_zz_x = pb.smooth_interpolated_skymap_with_gaussian_kernel(i_zz=i_zz_x, type=tmp["smooth"]["type"], sigma=tmp["smooth"]["sigma"])
        #     ax_histx.plot(grid_x, i_zz_x * 1e3,  lw=1., ls='-', drawstyle='steps', color=tmp["xsize"]["color"])  # , color='blue')
        #     ax_histx.axvline(x=x2, color=tmp["xsize"]["color"], linestyle='dotted')
        #     ax_histx.axvline(x=x1, color=tmp["xsize"]["color"], linestyle='dotted')
        # # --------------------
        # # if len(tmp["pcolormesh"].keys()) > 0:
        # #     # levels = np.geomspace(int_ration[~np.isnan(int_ration)].min(), int_ration[~np.isnan(int_ration)].max(), 50)
        # #     int_zz = int_zz.T
        # #     int_zz[~np.isfinite(int_zz)] = 0.
        # #     if (len(tmp["smooth"].keys()) > 0):
        # #         int_zz = pb.smooth_interpolated_skymap_with_gaussian_kernel(i_zz=int_zz, type=tmp["smooth"]["type"], sigma=tmp["smooth"]["sigma"])
        # #     plot_pcolomesh(ax=ax_main, task_dict=tmp["pcolormesh"], int_x=int_x, int_y=int_y, int_zz=int_zz)
        # if "pcolormesh" in tmp.keys():
        #     # levels = np.geomspace(int_ration[~np.isnan(int_ration)].min(), int_ration[~np.isnan(int_ration)].max(), 50)
        #     int_zz = int_zz.T
        #     int_zz[~np.isfinite(int_zz)] = float(tmp["pcolormesh"]["isnan"])
        #     if (len(tmp["smooth"].keys()) > 0):
        #         int_zz = pb.smooth_interpolated_skymap_with_gaussian_kernel(i_zz=int_zz, type=tmp["smooth"]["type"])
        #     im = plot_pcolomesh(ax=ax_main, task_dict=tmp["pcolormesh"], int_x=int_x, int_y=int_y, int_zz=int_zz)
        #     if tmp["pcolormesh"]["set_rasterized"]: ax_main.set_rasterized(im)

    ''' ---- plot size/cm of the grb ---- '''
    if (len(settings["grb_skymap"].keys()) > 0):
        tmp = copy.deepcopy(settings["grb_skymap"])
        if settings["precompute"]:

            all_x_jet, all_y_jet, all_fluxes_jet \
                = pb1.GRB.get_skymap(time=time * cgs.day, freq=freq, verbose=False, remove_mu=True)

            int_x_j, int_y_j, int_zz_j = combine_images(all_x_jet, all_y_jet, all_fluxes_jet,
                                                           hist_or_int="hist", shells=False, nx=tmp["hist_nx"], ny=tmp["hist_ny"], extend=2)
            grid_y_j, _i_zz_y_j, i_zz_y_j, _ = get_skymap_lat_dist(all_x_jet, all_y_jet, all_fluxes_jet,
                                                                      collapse_axis="x", nx=tmp["hist_nx"], ny=tmp["hist_ny"])
            grid_x_j, _i_zz_x_j, i_zz_x_j, _ = get_skymap_lat_dist(all_x_jet, all_y_jet, all_fluxes_jet,
                                                                      collapse_axis="y", nx=tmp["hist_nx"], ny=tmp["hist_ny"])

            # print("all_fluxes_jet=[{:.2e}, {:.2e}]".format(all_fluxes_jet.min(),all_fluxes_jet.max()))
            # # print("i_zz_y_j=[{:.2e}, {:.2e}]".format(i_zz_y_j.min(),i_zz_y_j.max()))
            # # print("i_zz_x_j=[{:.2e}, {:.2e}]".format(i_zz_x_j.min(),i_zz_x_j.max()))
            # print("int_zz_j=[{:.2e}, {:.2e}]".format(int_zz_j.min(),int_zz_j.max()))
            # print("i_zz_y_j=[{:.2e}, {:.2e}]".format(i_zz_y_j.min(),i_zz_y_j.max()))
            # print("i_zz_x_j=[{:.2e}, {:.2e}]".format(i_zz_x_j.min(),i_zz_x_j.max()))
            # print("MAX i_zz_x_j = {} | i_zz_x_j = {}".format(i_zz_x_j.max(), i_zz_y_j.max()))
            # exit(1)
            xc_m_j, yc_m_j = pb1.GRB.get_skymap_cm(all_x_jet, all_y_jet, all_fluxes_jet)
            #
            grp_j = grp.create_group("grb nx={} ny={}".format(tmp["hist_nx"], tmp["hist_ny"]))
            grp_j.create_dataset("int_x", data=int_x_j)
            grp_j.create_dataset("int_y", data=int_y_j)
            grp_j.create_dataset("int_zz", data=int_zz_j)
            grp_j.create_dataset("grid_y", data=grid_y_j)
            grp_j.create_dataset("i_zz_y", data=i_zz_y_j)
            grp_j.create_dataset("grid_x", data=grid_x_j)
            grp_j.create_dataset("i_zz_x", data=i_zz_x_j)
            grp_j.attrs.create("xc", data=xc_m_j)
            grp_j.attrs.create("yc", data=yc_m_j)
        else:
            grp_j = grp["grb nx={} ny={}".format(tmp["hist_nx"], tmp["hist_ny"])]
            int_x_j = np.array(np.array(grp_j["int_x"]))
            int_y_j = np.array(np.array(grp_j["int_y"]))
            int_zz_j = np.array(np.array(grp_j["int_zz"]))  # plotting
            grid_y_j = np.array(np.array(grp_j["grid_y"]))
            i_zz_y_j = np.array(np.array(grp_j["i_zz_y"]))
            grid_x_j = np.array(np.array(grp_j["grid_x"]))
            i_zz_x_j = np.array(np.array(grp_j["i_zz_x"]))
            xc_m_j = float(grp_j.attrs["xc"])
            yc_m_j = float(grp_j.attrs["yc"])
        print("Total flux kn={:.2e} grb={:.2e} x_c={}".format(
            0.,# pb1.get_ej_skymap_totfluxes(freq=freq, time=time*cgs.day),
            pb1.GRB.get_skymap_totfluxes(freq=freq, time=time*cgs.day),
            xc_m_j
        ))
        im = _plot_skymap_with_hists(ax_main, ax_histx, ax_histy, pb1, tmp,
                                     grid_x_j, grid_y_j, int_x_j, int_y_j, i_zz_x_j, i_zz_y_j, int_zz_j, xc_m_j, yc_m_j)

        # x1_j, x2_j = pb.get_skymap_fwhm(grid_x_j, i_zz_x_j, xc_m_j)
        # y1_j, y2_j = pb.get_skymap_fwhm(grid_y_j, i_zz_y_j, yc_m_j)
        # assert x2_j > x1_j
        # if len(tmp["cm"].keys()) > 0:
        #     ax.plot(xc_m_j, yc_m_j, **tmp["cm"])  # color='red', marker="o")
        # if len(tmp["ysize"].keys()) > 0:
        #     ax.errorbar([xc_m_j, xc_m_j], [y1_j, y2_j], xerr=[int_x_j.max() / 10, int_x_j.max() / 10], **tmp["ysize"])
        # if len(tmp["xsize"].keys()) > 0:
        #     ax.errorbar([x1_j, x2_j], [yc_m_j, yc_m_j], yerr=[int_y_j.max() / 10, int_y_j.max() / 10], **tmp["xsize"])
        # # --------------------
        # if len(tmp["pcolormesh"].keys()) > 0:
        #     norm = Normalize(int_zz_j[np.isfinite(int_zz_j)].min(), int_zz_j[np.isfinite(int_zz_j)].max())
        #     # levels = np.geomspace(int_ration[~np.isnan(int_ration)].min(), int_ration[~np.isnan(int_ration)].max(), 50)
        #     cmap = cm.get_cmap(tmp["pcolormesh"]["cmap"])
        #     int_zz_j = int_zz_j.T
        #     int_zz_j[~np.isfinite(int_zz_j)] = 0.
        #     im = ax.pcolormesh(int_x_j, int_y_j, int_zz_j, norm=norm, cmap=cmap, alpha=1.0)
        #     im.set_rasterized(True)
        #     ax.set_facecolor(cmap(0.0))

    ''' --- plot size/cm of the kn with grb (using previous two) ---'''
    if (len(settings["kn_grb_skymap"].keys()) > 0):
        assert (len(settings["grb_skymap"].keys()) > 0)
        assert (len(settings["kn_skymap"].keys()) > 0)
        tmp = copy.deepcopy(settings["kn_grb_skymap"])
        if settings["precompute"]:
            all_x_pj = all_x_jet + all_x
            all_y_pj = all_y_jet + all_y
            all_fluxes_pj = all_fluxes_jet + all_fluxes
            int_x_pj, int_y_pj, int_zz_pj = \
                combine_images(all_x_pj, all_y_pj, all_fluxes_pj, hist_or_int="hist", shells=False, nx=tmp["hist_nx"], ny=tmp["hist_ny"], extend=2)
            # _, _, int_zz_wjet_nojet = pb.combine_images(all_x_pj, all_y_pj, all_fluxes_pj, hist_or_int="hist", shells=False, nx=225, ny=175)
            xc_m_pj, yc_m_pj = pb1.GRB.get_skymap_cm(all_x_pj, all_y_pj, all_fluxes_pj)
            grid_y_pj, _i_zz_y_pj, i_zz_y_pj, _ = get_skymap_lat_dist(all_x_pj, all_y_pj, all_fluxes_pj,
                                                                         collapse_axis="x", nx=tmp["hist_nx"], ny=tmp["hist_ny"])
            grid_x_pj, _i_zz_x_pj, i_zz_x_pj, _ = get_skymap_lat_dist(all_x_pj, all_y_pj, all_fluxes_pj,
                                                                         collapse_axis="y", nx=tmp["hist_nx"], ny=tmp["hist_ny"])

            fnu_jet = pb1.GRB.get_skymap_totfluxes(freq=freq)[find_nearest_index(pb1.GRB.get_skymap_times(),time*cgs.day)]
            fnu_ej = pb1.KN.get_skymap_totfluxes(freq=freq)[find_nearest_index(pb1.KN.get_skymap_times(),time*cgs.day)]


            print("JET    ", pb1.GRB.get_skymap_totfluxes(freq=freq))
            print("EJECTA ", pb1.KN.get_skymap_totfluxes(freq=freq))
            print("Time={} Fnu jet={:.2e} ejecta={:.2e}".format(time, fnu_jet, fnu_ej))

            grp_kn_plus_grb = grp.create_group("kn_and_grb nx={} ny={}".format(tmp["hist_nx"], tmp["hist_ny"]))
            grp_kn_plus_grb.create_dataset("int_x", data=int_x_pj)
            grp_kn_plus_grb.create_dataset("int_y", data=int_y_pj)
            grp_kn_plus_grb.create_dataset("int_zz", data=int_zz_pj)
            grp_kn_plus_grb.create_dataset("grid_y", data=grid_y_pj)
            grp_kn_plus_grb.create_dataset("i_zz_y", data=i_zz_y_pj)
            grp_kn_plus_grb.create_dataset("grid_x", data=grid_x_pj)
            grp_kn_plus_grb.create_dataset("i_zz_x", data=i_zz_x_pj)
            grp_kn_plus_grb.attrs.create("xc", data=xc_m_pj)
            grp_kn_plus_grb.attrs.create("yc", data=yc_m_pj)
        else:
            grp_kn_plus_grb = grp["kn_and_grb nx={} ny={}".format(tmp["hist_nx"], tmp["hist_ny"])]
            int_x_pj = np.array(np.array(grp_kn_plus_grb["int_x"]))
            int_y_pj = np.array(np.array(grp_kn_plus_grb["int_y"]))
            int_zz_pj = np.array(np.array(grp_kn_plus_grb["int_zz"]))  # plotting
            grid_y_pj = np.array(np.array(grp_kn_plus_grb["grid_y"]))
            i_zz_y_pj = np.array(np.array(grp_kn_plus_grb["i_zz_y"]))
            grid_x_pj = np.array(np.array(grp_kn_plus_grb["grid_x"]))
            i_zz_x_pj = np.array(np.array(grp_kn_plus_grb["i_zz_x"]))
            xc_m_pj = float(grp_kn_plus_grb.attrs["xc"])
            yc_m_pj = float(grp_kn_plus_grb.attrs["yc"])

        im = _plot_skymap_with_hists(ax_main, ax_histx, ax_histy, pb1, tmp,
                                     grid_x_pj, grid_y_pj, int_x_pj, int_y_pj, i_zz_x_pj, i_zz_y_pj, int_zz_pj, xc_m_pj, yc_m_pj)

        # x1_pj, x2_pj = pb.get_skymap_fwhm(grid_x_pj, i_zz_x_pj, xc_m_pj)
        # y1_pj, y2_pj = pb.get_skymap_fwhm(grid_y_pj, i_zz_y_pj, yc_m_pj)
        # assert x2_pj > x1_pj
        # if len(tmp["cm"].keys()) > 0: ax.plot(xc_m_pj, yc_m_pj, **tmp["cm"])
        # if len(tmp["ysize"].keys()) > 0: ax.errorbar([xc_m_pj, xc_m_pj], [y1_pj, y2_pj],
        #                                              xerr=[int_x_pj.max() / 10, int_x_pj.max() / 10], **tmp["ysize"])
        # if len(tmp["xsize"].keys()) > 0: ax.errorbar([x1_pj, x2_pj], [yc_m_pj, yc_m_pj],
        #                                              yerr=[int_y_pj.max() / 10, int_y_pj.max() / 10], **tmp["xsize"])
        # # --------------------
        # if len(tmp["pcolormesh"].keys()) > 0:
        #     norm = Normalize(int_zz_pj[np.isfinite(int_zz_pj)].min(), int_zz_pj[np.isfinite(int_zz_pj)].max())
        #     # levels = np.geomspace(int_ration[~np.isnan(int_ration)].min(), int_ration[~np.isnan(int_ration)].max(), 50)
        #     cmap = cm.get_cmap(tmp["pcolormesh"]["cmap"])
        #     int_zz_pj = int_zz_pj.T
        #     int_zz_pj[~np.isfinite(int_zz_pj)] = 0.
        #     im = ax.pcolormesh(int_x_pj, int_y_pj, int_zz_pj, norm=norm, cmap=cmap, alpha=1.0)
        #     im.set_rasterized(True)
        #     ax.set_facecolor(cmap(0.0))

    ''' ---- plot size/cm of the ejecta with interaction ---- '''
    # pb_w = BPA_METHODS()
    # pb_w.res_dir_kn = plot_dict["res_dir_kn"]
    # pb_w.ejecta_prefix = plot_dict["ejecta_prefix2"]
    # pb_w.res_dir_grb = plot_dict["res_dir_grb"]
    # pb_w.jet_prefix = plot_dict["jet_prefix"]
    if len(settings["kn_w_skymap"].keys()) > 0:
        tmp = copy.deepcopy(settings["kn_w_skymap"])
        if settings["precompute"]:
            all_x_w, all_y_w, all_fluxes_w \
                = pb2.KN.get_skymap(time=time * cgs.day, freq=freq, verbose=False, remove_mu=True)
            int_x_w, int_y_w, int_zz_w = combine_images(all_x_w, all_y_w, all_fluxes_w, hist_or_int="hist",
                                                             shells=True, nx=tmp["hist_nx"], ny=tmp["hist_ny"], extend=2)
            grid_y_w, i_zz_y_w, _ = get_skymap_lat_dist(all_x_w, all_y_w, all_fluxes_w, collapse_axis="x",
                                                             nx=tmp["hist_nx"], ny=tmp["hist_ny"])
            grid_x_w, i_zz_x_w, _ = get_skymap_lat_dist(all_x_w, all_y_w, all_fluxes_w, collapse_axis="y",
                                                             nx=tmp["hist_nx"], ny=tmp["hist_ny"])
            xc_w, yc_w = pb2.KN.get_skymap_cm(all_x_w, all_y_w, all_fluxes_w)
            #
            grp_kn_w = grp.create_group("kn_w nx={} ny={}".format(tmp["hist_nx"], tmp["hist_ny"]))
            grp_kn_w.create_dataset("int_x", data=int_x_w)
            grp_kn_w.create_dataset("int_y", data=int_y_w)
            grp_kn_w.create_dataset("int_zz", data=int_zz_w)
            grp_kn_w.create_dataset("grid_y", data=grid_y_w)
            grp_kn_w.create_dataset("i_zz_y", data=i_zz_y_w)
            grp_kn_w.create_dataset("grid_x", data=grid_x_w)
            grp_kn_w.create_dataset("i_zz_x", data=i_zz_x_w)
            grp_kn_w.attrs.create("xc", data=xc_w)
            grp_kn_w.attrs.create("yc", data=yc_w)
        else:
            grp_kn_w = grp["kn_w nx={} ny={}".format(tmp["hist_nx"], tmp["hist_ny"])]
            int_x_w = np.array(np.array(grp_kn_w["int_x"]))
            int_y_w = np.array(np.array(grp_kn_w["int_y"]))
            int_zz_w = np.array(np.array(grp_kn_w["int_zz"]))  # plotting
            grid_y_w = np.array(np.array(grp_kn_w["grid_y"]))
            i_zz_y_w = np.array(np.array(grp_kn_w["i_zz_y"]))
            grid_x_w = np.array(np.array(grp_kn_w["grid_x"]))
            i_zz_x_w = np.array(np.array(grp_kn_w["i_zz_x"]))
            xc_w = float(grp_kn_w.attrs["xc"])
            yc_w = float(grp_kn_w.attrs["yc"])
        i_zz_x_w /= i_zz_x_w.max()
        i_zz_y_w /= i_zz_y_w.max()
        y1_w, y2_w = get_skymap_fwhm(grid_y_w, i_zz_y_w, yc_w)
        x1_w, x2_w = get_skymap_fwhm(grid_x_w, i_zz_x_w, xc_w)
        assert x2_w > x1_w
        if len(tmp["cm"].keys()) > 0: ax.plot(xc_w, yc_w, **tmp["cm"])
        if len(tmp["ysize"].keys()) > 0: ax.errorbar([xc_w, xc_w], [y1_w, y2_w],
                                                     xerr=[int_x_w.max() / 10, int_x_w.max() / 10], **tmp["ysize"])
        if len(tmp["xsize"].keys()) > 0: ax.errorbar([x1_w, x2_w], [yc_w, yc_w],
                                                     yerr=[int_y_w.max() / 10, int_y_w.max() / 10], **tmp["xsize"])
        # --------------------
        if len(tmp["pcolormesh"].keys()) > 0:
            norm = Normalize(int_zz_w[np.isfinite(int_zz_w)].min(), int_zz_w[np.isfinite(int_zz_w)].max())
            # levels = np.geomspace(int_ration[~np.isnan(int_ration)].min(), int_ration[~np.isnan(int_ration)].max(), 50)
            cmap = cm.get_cmap(tmp["pcolormesh"]["cmap"])
            int_zz_w = int_zz_w.T
            int_zz_w[~np.isfinite(int_zz_w)] = 0.
            im = ax.pcolormesh(int_x_w, int_y_w, int_zz_w, norm=norm, cmap=cmap, alpha=1.0)
            im.set_rasterized(True)
            ax.set_facecolor(cmap(0.0))

    ''' --- plot relation --- '''
    if len(settings["kn_skymap_ratio"].keys()) > 0:
        assert (len(settings["kn_skymap"].keys()) > 0)
        assert (len(settings["kn_w_skymap"].keys()) > 0)
        tmp = copy.deepcopy(settings["kn_skymap_ratio"])
        if settings["precompute"]:
            int_x_w, int_y_w, int_zz, int_zz_w = \
                pb.get_combained_ej_skymaps_adjusted_to_other(time=time, freq=freq, other_pb_instance=pb_w,
                                                              nx=tmp["hist_nx"], ny=tmp["hist_ny"])
            kn_adjusted = grp.create_group("kn_adjusted nx={} ny={}".format(tmp["hist_nx"], tmp["hist_ny"]))
            kn_adjusted.create_dataset("int_x", data=int_x_w)
            kn_adjusted.create_dataset("int_y", data=int_y_w)
            kn_adjusted.create_dataset("int_zz", data=int_zz)
            kn_adjusted.create_dataset("int_zz_w", data=int_zz_w)
        else:
            kn_adjusted = grp["kn_adjusted nx={} ny={}".format(tmp["hist_nx"], tmp["hist_ny"])]
            int_x_w = np.array(np.array(kn_adjusted["int_x"]))
            int_y_w = np.array(np.array(kn_adjusted["int_y"]))
            int_zz = np.array(np.array(kn_adjusted["int_zz"]))
            int_zz_w = np.array(np.array(kn_adjusted["int_zz_w"]))
        int_ration = int_zz_w / int_zz
        int_ration = int_ration.T
        int_ration[~np.isfinite(int_ration)] = 1.

        # norm = Normalize(int_ration[np.isfinite(int_ration)].min(), int_ration[np.isfinite(int_ration)].max())
        # # levels = np.geomspace(int_ration[~np.isnan(int_ration)].min(), int_ration[~np.isnan(int_ration)].max(), 50)
        # cmap = cm.get_cmap('viridis')
        # int_ration = int_ration.T
        # int_ration[~np.isfinite(int_ration)] = 0.
        # im = ax.pcolormesh(int_x_w, int_y_w, int_ration, norm=norm, cmap=cmap, alpha=1.0)
        # im.set_rasterized(True)
        # ax.set_facecolor(cmap(0.0))

        vmin = tmp["pcolormesh"]["vmin"]
        vmax = tmp["pcolormesh"]["vmax"]
        vcenter = tmp["pcolormesh"]["vcenter"]
        levels = MaxNLocator(nbins=40).tick_values(int_ration.min(), int_ration.max())
        # levels = MaxNLocator(nbins=40).tick_values(-5, 1)
        cmap = plt.get_cmap(tmp["pcolormesh"]["cmap"])
        norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
        norm = colors.TwoSlopeNorm(vmin=vmin, vcenter=vcenter, vmax=vmax)
        # im = aux_ax2.contourf(theta, r, values, levels=np.logspace(np.log10(1e-4),np.log10(1e0),20), cmap=cmap, norm=norm)
        # im = aux_ax2.contourf(theta, r, values, levels=np.linspace(val.min(),val.max(),20), cmap=cmap, norm=norm)
        im = ax.pcolormesh(int_x_w, int_y_w, int_ration, cmap=cmap, norm=norm, alpha=1.0)
        im.set_rasterized(True)

    tmp_tile.close()

    ax_histx.set_yscale("linear")
    ax_histx.set_xscale("linear")
    ax_histx.minorticks_on()
    ax_histx.set_yscale("log")
    # ax_histx.set_ylabel("$\sum_{z}{I}$ [$\mu$Jy/mas$^2$]", fontsize=12)  # , color='blue')
    ax_histx.set_ylabel(r"$I_{\nu;\rm m}(x)$", fontsize=12)  # , color='blue')
    ax_histx.tick_params(axis='both', which='both', labelleft=True, labelbottom=False,
                         labelright=False, tick1On=True, tick2On=True,
                         labelsize=12, color="white",
                         direction='in',
                         bottom=False, top=True, left=True, right=True)
    ax_histx.set_facecolor(settings["histx_backgound_color"])
    ax_histy.set_yscale("linear")
    ax_histy.set_xscale("linear")
    ax_histy.minorticks_on()
    ax_histy.set_xscale("log")
    ax_histy.set_xlim(ax_histx.get_ylim())
    # ax_histy.set_xlabel("$\sum_{x}{I}$ [$\mu$Jy/mas$^2$]", fontsize=12)  # , color='blue')
    ax_histy.set_xlabel(r"$I_{\nu;\rm m}(z)$", fontsize=12)  # , color='blue')
    ax_histy.tick_params(axis='both', which='both', labelleft=False,
                         labelright=False, tick1On=True, tick2On=True,
                         labelsize=12, color="white",
                         direction='in',
                         bottom=True, top=True, left=True, right=False)
    ax_histy.set_facecolor(settings["histy_backgound_color"])
    ax_main.set_yscale("linear")
    ax_main.set_xscale("linear")
    ax_main.minorticks_on()
    ax_main.set_xlabel(settings["xlabel"], fontsize=12)
    ax_main.set_ylabel(settings["ylabel"], fontsize=12)
    # ax_main.set_xlim(int_x.min() * 1.1, int_x.max() * 1.1)
    # ax_main.set_ylim(int_y.min() * 1.1, int_y.max() * 1.1)
    if "xlim" in settings.keys() and len(settings["xlim"]) == 2: ax_main.set_xlim(settings["xlim"])
    else: ax_main.set_xlim(int_x.min() * 1.1, int_x.max() * 1.1)
    if "ylim" in settings.keys() and len(settings["ylim"]) == 2: ax_main.set_ylim(settings["ylim"])
    else: ax_main.set_ylim(int_y.min() * 1.1, int_y.max() * 1.1)
    ax_main.tick_params(axis='both', which='both', labelleft=True,
                        labelright=False, tick1On=True, tick2On=True,
                        labelsize=12,
                        direction='in',
                        bottom=True, top=True, left=True, right=True)
    ax_main.minorticks_on()
    ax_main.axhline(y=0, linestyle=':', linewidth=0.4)
    ax_main.axvline(x=0, linestyle=':', linewidth=0.4)
    ax_main.tick_params(axis='both', which='both', labelleft=True,
                        labelright=False, tick1On=True, tick2On=True,
                        labelsize=12, color="white",
                        direction='in',
                        bottom=True, top=True, left=True, right=True)

    if "text" in settings.keys() and settings["text"] == "full":
        tot_flux = float(pb1.KN.get_skymap_totfluxes(freq) * 1e3)
        ax_main.text(1.05, 0.6, r"$\\I_{\rm max}$="
                     + r"${}$".format(latex_float(int_zz.max() * 1e3)) + r"\,[$\mu$Jy/mas$^2$]\\"
                     + r"$F_{\nu}=$" + r"${}$".format(latex_float(tot_flux)) + r"\,[$\mu$Jy]\\"
                     + r" $\nu=$" + r"${}$".format(latex_float(freq)) + r"\,[Hz] \\"
                     + r"$\theta_{\rm obs}=" + r"{:.1f}".format(theta_obs) + r"$ [deg]\\"
                     + r"$t_{\rm obs}=" + r"{:.1f}".format(time) + "$ [day]"
                     , color='black', transform=ax.transAxes, fontsize=fontsize - 2)

    if settings["plot_grids"]:
        ax_main.grid()
        ax_histx.grid()
        ax_histy.grid()
    ax_histx.set_ylim(*settings["histx_lim"])
    ax_histy.set_xlim(*settings["histy_lim"])

    # ax.tick_params(  # axis='both',
    #     which='both',  # labelleft=True,
    #     # labelright=False,
    #     tick1On=True, tick2On=True,
    #     labelsize=12,
    #     direction='in',
    #     color="black",
    #     #               bottom=True, top=True, left=True, right=True
    # )
    # ax.set_xlim(settings["xlim"])
    # ax.set_ylim(settings["ylim"])
    # if settings["grid"]:
    #     ax.set_xlim(min(settings["xlim"]) * 1.1, max(settings["xlim"]) * 1.1)
    #     ax.set_ylim(min(settings["ylim"]) * 1.1, max(settings["ylim"]) * 1.1)
    #
    # # ax.tick_params(axis='both', which='both', labelleft=True,
    # #                labelright=False, tick1On=True, tick2On=True,
    # #                labelsize=11, color="white",
    # #                direction='in',
    # #                bottom=True, top=True, left=True, right=True)
    # # ax.set_title(r"$t_{\rm obs}=$" + r"{:d} day".format(time))
    # # ax.set_title(line["label"])
    #
    # ax.axhline(y=0, linestyle=':', linewidth=0.4, color="black")
    # ax.axvline(x=0, linestyle=':', linewidth=0.4, color="black")
    #
    # if (len(settings["title"].keys()) > 0):
    #     if settings["title"]["title"] == "time_fluxratio":
    #         times = pb.get_ej_skymap_times()
    #         fluxes = pb.get_ej_skymap_totfluxes(freq=freq)
    #         fluxes_w = pb_w.get_ej_skymap_totfluxes(freq=freq)
    #         fluxes_ration = fluxes_w / fluxes
    #         idx = find_nearest_index(times, time * cgs.day)
    #         # ax.set_title("t={:.0f} d. F/F={}".format(time, fluxes_ration[idx]))
    #         ax.set_title("t={:.0f} d. ".format(time) +
    #                      r"$F_{\nu}^{\rm w}/F_{\nu}^{\rm w/o}=$" +
    #                      "{:.2f}".format(fluxes_ration[idx]))
    #     elif settings["title"]["title"] == "time":
    #         ax.set_title("t={:.0f} d. ".format(time))
    #
    # # ------------------------------------------------------------------------------------------------------------------
    #
    #
    # ''' ------------------------------------- '''
    #
    # fig = plt.figure(figsize=(5, 5))
    #
    # # Add a gridspec with two rows and two columns and a ratio of 2 to 7 between
    # # the size of the marginal axes and the main axes in both directions.
    # # Also adjust the subplot parameters for a square plot.
    # gs = fig.add_gridspec(2, 2, width_ratios=(4, 2), height_ratios=(2, 4),
    #                       left=0.13, right=0.9, bottom=0.1, top=0.9, wspace=0.05, hspace=0.05)
    #
    # ax_main = fig.add_subplot(gs[1, 0])
    # ax_histx = fig.add_subplot(gs[0, 0], sharex=ax_main)
    # ax_histy = fig.add_subplot(gs[1, 1], sharey=ax_main)
    #
    # ''' -------------------------------------- '''
    #
    # # nx = 200
    # # ny = 100
    # # nx = np.complex(0, nx + 1)
    # # ny = np.complex(0, ny + 1)
    # # edges_x = np.mgrid[all_x.min() * 1.2:all_x.max() * 1.2:nx]
    # # edges_y = np.mgrid[all_y.min() * 1.2:all_y.max() * 1.2:ny]
    # # # grid_x, grid_y = np.mgrid[xrs_i.min() * 1.2:xrs_i.max() * 1.2:nx,
    # # #                  yrs_i.min() * 1.2:yrs_i.max() * 1.2:ny]
    # # i_zz, _ = np.histogramdd(tuple([all_x, all_y]), bins=tuple([edges_x, edges_y]), weights=all_fluxes)
    # # grid_x = 0.5 * (edges_x[1:] + edges_x[:-1])
    # # grid_y = 0.5 * (edges_y[1:] + edges_y[:-1])
    # # tmp = []
    # # for i in range(len(grid_x)):
    # #     tmp.append(np.sum(i_zz[i,:] * np.diff(edges_y))/np.sum(np.diff(edges_y)))
    # # tmp = np.array(tmp)
    # #
    # # y_sum_ii = np.sum(i_zz * 1e3, axis=1)
    # # max_ii = np.max(y_sum_ii)
    # # x1 = grid_x[np.argmin(y_sum_ii < max_ii * 0.5)]
    # # x2 = grid_x[::-1][np.argmin(y_sum_ii[::-1] < max_ii * 0.5)]
    #
    # grid_x, i_zz_x, _ = pb.get_skymap_lat_dist(all_x, all_y, all_fluxes, collapse_axis="y")
    # xc, yc = pb.get_ej_skymap_cm(all_x, all_y, all_fluxes)
    # x2, x1 = pb.get_skymap_fwhm(grid_x, i_zz_x, xc)
    #
    # ax = ax_histx
    # ax.plot(grid_x, i_zz_x * 1e3, lw=1., ls='-', drawstyle='steps', color="black")  # , color='blue')
    # # ax.plot(grid_x, np.sum(i_zz * 1e3, axis=1) / np.sum(np.diff(edges_y)), lw=1., ls='--', drawstyle='steps')#, color='blue')
    # # ax.plot(int_x, int_zz.max(axis=1), lw=1., ls='-')
    # # ax.plot(int_x, np.sum(int_zz * int_y, axis=1)/np.sum(int_y,axis=1), lw=1., ls='--')
    # tot_flux = float(pb.get_ej_skymap_totfluxes(freq) * 1e3)
    # ax.text(1.05, 0.6, r"$\\I_{\rm max}$="
    #         + r"${}$".format(latex_float(int_zz.max() * 1e3)) + r"\,[$\mu$Jy/mas$^2$]\\"
    #         + r"$F_{\nu}=$" + r"${}$".format(latex_float(tot_flux)) + r"\,[$\mu$Jy]\\"
    #         + r" $\nu=$" + r"${}$".format(latex_float(freq)) + r"\,[Hz] \\"
    #         + r"$\theta_{\rm obs}=" + r"{:.1f}".format(theta_obs) + r"$ [deg]\\"
    #         + r"$t_{\rm obs}=" + r"{:.1f}".format(time) + "$ [day]"
    #         , color='black', transform=ax.transAxes, fontsize=fontsize - 2)
    # ax.set_yscale("log")
    # ax.set_ylabel("$\sum_{z}{I}$ [$\mu$Jy/mas$^2$]", fontsize=12)  # , color='blue')
    # ax.tick_params(axis='both', which='both', labelleft=True, labelbottom=False,
    #                labelright=False, tick1On=True, tick2On=True,
    #                labelsize=fontsize,
    #                direction='in',
    #                bottom=False, top=True, left=True, right=True)
    # # ax1=ax.twinx()
    # # ax1.plot(grid_x, tmp, lw=1., ls='-', drawstyle='steps', color='red')
    # # ax1.set_yscale("log")
    # # ax1.set_ylabel(r"$\sum_{y} I \Delta y / \sum_{y} \Delta y$ [$\mu$Jy/mas$^2$]", fontsize=12, color='red')
    # ax.axvline(x=x1, color='gray', linestyle='dotted')
    # ax.axvline(x=x2, color='gray', linestyle='dotted')
    # ax.axvline(x=xcs_m, color='gray', linestyle='solid')
    #
    # ''' --------------------------------------- '''
    #
    # # nx = 200
    # # ny = 100
    # # nx = np.complex(0, nx + 1)
    # # ny = np.complex(0, ny + 1)
    # # edges_x = np.mgrid[all_x.min() * 1.2:all_x.max() * 1.2:nx]
    # # edges_y = np.mgrid[all_y.min() * 1.2:all_y.max() * 1.2:ny]
    # # # grid_x, grid_y = np.mgrid[xrs_i.min() * 1.2:xrs_i.max() * 1.2:nx,
    # # #                  yrs_i.min() * 1.2:yrs_i.max() * 1.2:ny]
    # # i_zz, _ = np.histogramdd(tuple([all_x, all_y]), bins=tuple([edges_x, edges_y]), weights=all_fluxes)
    # # grid_x = 0.5 * (edges_x[1:] + edges_x[:-1])
    # # grid_y = 0.5 * (edges_y[1:] + edges_y[:-1])
    # # tmp = []
    # # for i in range(len(grid_y)):
    # #     tmp.append(np.sum(i_zz[:, i] * np.diff(edges_x)) / np.sum(np.diff(edges_x)))
    # # tmp = np.array(tmp)
    # #
    # # x_sum_ii = np.sum(i_zz * 1e3, axis=0)
    # # max_ii = np.max(x_sum_ii)
    # # y1 = grid_y[np.argmin(x_sum_ii < max_ii * 0.5)]
    # # y2 = grid_y[::-1][np.argmin(x_sum_ii[::-1] < max_ii * 0.5)]
    # # y1, y2, grid_y, i_zz_y = pb.get_ej_skymap_fwhm(all_x, all_y, all_fluxes,axis="z")
    # grid_y, i_zz_y, _ = pb.get_skymap_lat_dist(all_x, all_y, all_fluxes, collapse_axis="x")
    # xc, yc = pb.get_ej_skymap_cm(all_x, all_y, all_fluxes)
    # y2, y1 = pb.get_skymap_fwhm(grid_y, i_zz_y, yc)
    #
    # ax = ax_histy
    # ax.plot(i_zz_y * 1e3, grid_y, lw=1., ls='-', drawstyle='steps', color='black')
    # # ax.plot(int_x, int_zz.max(axis=1), lw=1., ls='-')
    # # ax.plot(int_x, np.sum(int_zz * int_y, axis=1)/np.sum(int_y,axis=1), lw=1., ls='--')
    # ax.set_xscale("log")
    # ax.set_xlabel("$\sum_{x}{I}$ [$\mu$Jy/mas$^2$]", fontsize=12)  # , color='blue')
    # ax.tick_params(axis='both', which='both', labelleft=False,
    #                labelright=False, tick1On=True, tick2On=True,
    #                labelsize=fontsize,
    #                direction='in',
    #                bottom=True, top=True, left=True, right=False)
    # # ax1 = ax.twiny()
    # # ax1.plot(tmp, grid_y, lw=1., ls='-', drawstyle='steps', color='red')
    # # ax1.set_xscale("log")
    # # ax1.set_xlabel(r"$\sum_{y} I \Delta y / \sum_{y} \Delta y$ [$\mu$Jy/mas$^2$]", fontsize=12, color='red')
    # ax.axhline(y=y1, color='gray', linestyle='dotted')
    # ax.axhline(y=y2, color='gray', linestyle='dotted')
    # ax.axhline(y=ycs_m, color='gray', linestyle='solid')
    #
    # # ax.yaxis.set_ticklabels([])
    #
    # ''' --------------------------------------- '''
    #
    # # my_norm = LogNorm(1e-12, 1e-8)
    # # my_norm = LogNorm(1e-2, 1e0)
    # # my_norm = Normalize(int_zz.min(),int_zz.max())
    # # my_norm = LogNorm(int_zz.max()*1e-2,int_zz.max())
    # cmap = cm.get_cmap('viridis')
    #
    # ax = ax_main
    #
    # # if not (title is None): ax.set_title(title)
    # im = _plot_image(ax, int_x, int_y, int_zz, LogNorm(int_zz.max() * 1e-2, int_zz.max()),
    #                  levels=np.geomspace(int_zz.max() * 1e-2, int_zz.max(), 50))
    #
    # # center of mass
    # ax.plot(xcs_m, ycs_m, "x", color='white')
    #
    # # _i = np.argmax(x_latAvDist)
    # # xm = (int_x)[np.argmax(x_latAvDist), np.argmax(y_latAvDist)]
    # # ym = (int_y)[np.argmax(x_latAvDist), np.argmax(y_latAvDist)]
    # # ax.plot(_x, fwhm, ls='none', marker='.')
    # # ax.scatter([xm], [ym], s=80, marker='o', color="white", facecolors='none')
    # # _fwhm = fwhm
    #
    # # im.cmap.set_over("black")
    # # im.cmap.set_over("gray")
    # #
    # ax.set_yscale("linear")
    # ax.set_xscale("linear")
    # #
    # ax.set_xlabel("X [mas]", fontsize=12)
    # ax.set_ylabel("Z [mas]", fontsize=12)
    # #
    # # ax.set_xlim(grid_x.min(), grid_x.max())
    # # ax.set_ylim(grid_y.min(), grid_y.max())
    # #
    # ax.set_xlim(int_x.min() * 1.1, int_x.max() * 1.1)
    # ax.set_ylim(int_y.min() * 1.1, int_y.max() * 1.1)
    # #
    # ax.tick_params(axis='both', which='both', labelleft=True,
    #                labelright=False, tick1On=True, tick2On=True,
    #                labelsize=12,
    #                direction='in',
    #                bottom=True, top=True, left=True, right=True)
    # ax.minorticks_on()
    # ax.axhline(y=0, linestyle=':', linewidth=0.4)
    # ax.axvline(x=0, linestyle=':', linewidth=0.4)
    # ax.tick_params(axis='both', which='both', labelleft=True,
    #                labelright=False, tick1On=True, tick2On=True,
    #                labelsize=fontsize, color="white",
    #                direction='in',
    #                bottom=True, top=True, left=True, right=True)
    # --------- CBAR ----------------
    divider = make_axes_locatable(ax_cbar)
    cax = divider.append_axes('left', size='99%', pad=0.9)
    plt.delaxes(ax_cbar)
    cbar = plt.colorbar(im, cax=cax,
                        # format='%.0e',ticks=ticks
                        orientation='vertical',
                        # label=r"$I_{\nu}$ [mJy/mas$^2$]"
                        )
    cbar.set_label(label=r"$I_{\nu}$ [mJy/mas$^2$]", fontsize=12)
    cax.tick_params(axis='both', which='both', labelleft=False,
                    labelright=True, tick1On=True, tick2On=True,
                    labelsize=12,
                    direction='in',
                    bottom=True, top=True, left=True, right=True)
    # print("plotted: \n")
    # plt.tight_layout()
    figname = settings["figfpath"]  # make_prefix_for_pars(pars)
    if settings["save_figs"]:
        print(r"Saving:\n{}".format(figname + ".png"))
        plt.savefig(".png", dpi=256)
    if settings["save_pdf"]:
        print(r"Saving:\n{}".format(figname + ".pdf"))
        plt.savefig(figname + ".pdf")
    if settings["show_figs"]: plt.show()