import h5py
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy.spatial import cKDTree
from glob import glob
import pandas as pd
import os

from PyBlastAfterglowMag.interface import PyBlastAfterglow
from PyBlastAfterglowMag.utils import cgs
from sklearn.metrics import mean_squared_error

data_dir = "/media/vsevolod/T7/work/prj_grb_afterglow/GRB190114C/working_dirs/"

class EBL:
    def __init__(self):
        pass
    def fill_nan_by_extrapolation(self):
        # Create meshgrid for interpolation
        X, Y = np.meshgrid(self.ebl_z, self.ebl_en)

        # Flatten the grids and data for processing with griddata
        points = np.column_stack((X.ravel(), Y.ravel()))
        values = self.ebl_table.ravel()

        # Mask for valid (non-NaN) and invalid (NaN) points
        valid_mask = ~np.isnan(values)
        missing_mask = np.isnan(values)

        # Only keep valid points and values
        valid_points = points[valid_mask]
        valid_values = values[valid_mask]

        # Initial interpolation (linear or cubic) - change method if needed
        filled_values = interpolate.griddata(valid_points, valid_values, points, method='linear', fill_value=np.nan, rescale=True)

        # Check if there are still nans after griddata
        if np.isnan(filled_values).any():
            # raise ValueError()
            # Find nearest non-NaN points using a KDTree and fill NaNs
            kdtree = cKDTree(valid_points)
            distances, indices = kdtree.query(points[missing_mask])
            filled_values[missing_mask] = valid_values[indices]

        if np.isnan(filled_values).any():
            raise ValueError()

        # Reshape back to the original shape of the input data array
        res = filled_values.reshape(self.ebl_table.shape)
        return res

    def load(self):
        # load EBL table and compute absorbtion
        ebl_table = np.loadtxt("../../../data/EBL/Franceschini18/table.txt")
        self.ebl_z = ebl_table[0, 1:]
        self.ebl_en = ebl_table[1:, 0] * 1e12 # TeV -> eV
        self.ebl_table = ebl_table[1:,1:]
    def interpolate_fill_nans(self):
        self.ebl_table = self.fill_nan_by_extrapolation()
    def set_interpolator(self, z:np.ndarray, en:np.ndarray):
        X, Y = np.meshgrid(self.ebl_z, self.ebl_en)
        # Flatten the grids and data for processing with griddata
        points = np.column_stack((X.ravel(), Y.ravel()))
        values = self.ebl_table.ravel()
        #
        x, y = np.meshgrid(z, en)
        new_points = np.column_stack((x.ravel(), y.ravel()))
        # Initial interpolation (linear or cubic) - change method if needed
        res = interpolate.griddata(points, values, new_points, method='linear', fill_value=np.nan, rescale=True)
        return np.reshape(res,newshape=(len(z),len(en)))

def integrated_spectrum(t1,t2,working_dir):
    pba = PyBlastAfterglow(workingdir=working_dir+'/',readparfileforpaths=True)
    freqs = pba.GRB.get_lc_freqs(unique=True,spec=False)
    times = pba.GRB.get_lc_times(unique=True,spec=False)
    mask = ((times >= t1) & (times <= t2))
    spec = [pba.GRB.get_lc(time=t) for t in times if t >= t1 and t <= t2]
    spec = np.reshape(spec, newshape=(len(spec),len(spec[0])))
    integ_spec = np.trapz(spec,x=times[mask],axis=0)
    return freqs, integ_spec

def compute_rms_scores_for_runs(data_freq, data_flux, t1, t2):
    # load lightcurves
    working_dirs = glob(data_dir + '*')
    scores = []
    for working_dir in working_dirs:
        # pba = PyBlastAfterglow(workingdir=working_dir+'/',readparfileforpaths=True)
        # freqs = pba.GRB.get_lc_freqs(unique=True,spec=False)
        # times = pba.GRB.get_lc_times(unique=True,spec=False)
        # mask = ((times >= 68) & (times <= 110))
        # spec = [pba.GRB.get_lc(time=t) for t in times if t >= 68 and t <= 110]
        # spec = np.reshape(spec, newshape=(len(spec),len(spec[0])))
        freqs, integ_spec = integrated_spectrum(t1=t1,t2=t2,working_dir=working_dir)
        interp_spec = interpolate.interp1d(freqs, integ_spec*freqs*1.e-26,kind='linear')(data_freq)
        rms = mean_squared_error(np.log10(data_flux), np.log10(interp_spec))
        scores.append(rms)

        # ax.plot(data[:,0],interp_spec,marker='x')
        # ax.plot(freqs,integ_spec*freqs*1e-20)
    df = pd.DataFrame.from_dict({"workingdir":working_dirs,"rms":scores})
    df.to_csv("./runs.csv")

def get_rms_scores_for_all_runs(data_freq:list[np.ndarray],
                                data_flux:list[np.ndarray],
                                t1s:list[float], t2s:list[float])->tuple[list[str],list[float]]:
    # load lightcurves
    working_dirs = glob(data_dir + '*')
    scores = []
    min_score = np.inf
    for i, working_dir in enumerate(working_dirs):
        all_freqs, all_fluxes, all_model_fluxes = [], [], []
        for (obs_freq, obs_flux, t1, t2) in zip(data_freq,data_flux,t1s,t2s):
            # compute spectrum for this time interval (integrated over time)
            freqs, model_flux_dens = integrated_spectrum(t1=t1,t2=t2,working_dir=working_dir)
            flux = model_flux_dens*freqs*1.e-26 # flux in cgs units (eg / cm^2 / s)
            # interpolate spectrum for observer freqs
            model_fluxes = interpolate.interp1d(freqs, flux,kind='linear')(obs_freq)
            # store result for later error calculation
            all_freqs.append(freqs)
            all_fluxes.append(obs_flux)
            all_model_fluxes.append(model_fluxes)
        # compute error for all data for this simulation
        rms = mean_squared_error(np.log10(np.concatenate(all_fluxes).flatten()),
                                 np.log10(np.concatenate(all_model_fluxes).flatten()))
        rms /= len(np.concatenate(all_fluxes).flatten())
        scores.append(rms)
        min_score = min(min_score, rms)
        if i % 100 == 0:
            print(f"Processing {i}/{len(working_dirs)} min_score={min_score}")
    return (working_dirs, scores)

def load_observations(fname:str)->tuple[np.ndarray,np.ndarray,np.ndarray,np.ndarray]:
    # load VHE objservations
    with open(fname) as f:
        lines = [line.rstrip('\n') for line in f]
    new_lines = []
    for line in lines:
        if not ((line == '') or (line[0] == '|') or (line[0] == "\\")):
            new_lines.append([float(val) for val in line.split()])
    data = np.reshape(np.array(new_lines),newshape=(len(new_lines),len(new_lines[0])))
    data[:,0] *= 2.417989242e14 # ev -> Hz
    # data[:,1] *= 1.e26
    # data[:,2] *= 1.e26
    # data[:,3] *= 1.e26
    return (data[:,0], data[:,1], data[:,2], data[:,3]) # Hz, erg/cm^2/s/Hz, up_err, low_err

def apply_ebl_absorption_to_data(freqs:np.ndarray, fluxes:np.ndarray, z:float)->np.ndarray:
    # load EBL table and compute absorbtion
    ebl = EBL()
    ebl.load()
    ebl.interpolate_fill_nans()
    tau = ebl.set_interpolator(z=np.full_like(freqs,z),en=freqs/2.417989242e14) # Hz->ev in input
    atten = np.exp(-tau)
    atten[~np.isfinite(atten)] = 1.
    atten = atten[:,0] # remove repeated dimension from interpolation

    return fluxes*atten

def plot():

    z = 0.4245 # https://www.aanda.org/articles/aa/pdf/2022/03/aa41788-21.pdf

    # VHE intervals
    data = [
        dict(name="int1",t1=68, t2=110, file='./magic_int1_points.txt',color='blue',label=r"$68-110$ [s]"),
        dict(name="int2",t1=110, t2=180, file='./magic_int2_points.txt',color='orange',label=r"$180-360$ [s]"),
        # dict(name="int3",t1=180, t2=360, file='./magic_int3_points.txt',color='red'),
        # dict(name="int4",t1=360, t2=625, file='./magic_int4_points.txt',color='green'),
        # dict(name="int5",t1=625, t2=2400, file='./magic_int5_points.txt',color='purple'),
    ]

    # process data

    # freqs, obs_data, model_data = [], [], []
    # # process data for each time interval
    # for data_i in data:
    #     # load observations for this interval
    #     freqs, fluxes, _, _ = load_observations(data_i["file"])
    #     # apply EBL absorption to data (same as in the synthetic light curves)
    #     fluxes = apply_ebl_absorption_to_data(freqs=freqs,fluxes=fluxes, z=z)
    #     data_i["freqs"] = freqs
    #     data_i["fluxes"] = fluxes
    #
    # # compute error metric
    # working_dirs, scores = get_rms_scores_for_all_runs(
    #     data_freq=[data_i["freqs"] for data_i in data],
    #     data_flux=[data_i["fluxes"] for data_i in data],
    #     t1s=[data_i["t1"] for data_i in data],
    #     t2s=[data_i["t2"] for data_i in data],
    # )
    # df = pd.DataFrame.from_dict({"workingdir":working_dirs,"rms":scores})
    # df.to_csv("./runs.csv")

    # plot results

    fig, ax = plt.subplots(ncols=1,nrows=1,figsize=(5,3),layout='constrained')
    data_to_save = []
    for data_i in data:
        # load observations for this interval
        freqs, fluxes, err_u, err_l = load_observations(data_i["file"])
        # apply EBL absorption to data (same as in the synthetic light curves)
        err_l = apply_ebl_absorption_to_data(freqs=freqs,fluxes=fluxes-err_l, z=z)
        err_u = apply_ebl_absorption_to_data(freqs=freqs,fluxes=err_u-fluxes, z=z)
        fluxes = apply_ebl_absorption_to_data(freqs=freqs,fluxes=fluxes, z=z)

        # plot
        # ax.plot(freqs,fluxes,color=data_i['color'],marker='o',ls='none',fillstyle='none')
        yerr = np.array(list(zip(err_l,err_u))).T
        ax.errorbar(freqs,fluxes,yerr=yerr,
                    fmt='o',mfc='white',marker='o',capsize=2,color=data_i['color'],mec= data_i['color'],
                    ecolor = data_i['color'],fillstyle='full',label=data_i['label'])
        data_to_save.append(
            {"freqs":freqs, "fluxes":fluxes, "err_l":err_l, "err_u":err_u} | data_i
        )
    # print(data_to_save)

    # plot best model (load, sort by score, plot)
    df = pd.read_csv("./runs.csv",index_col=0)
    df.sort_values(by='rms',inplace=True)
    df.drop_duplicates('rms',inplace=True)
    for i, ls in zip((0,),['-','--',':']):
        working_dir, rms = df.iloc[i] # best model (lowest score)
        print(working_dir, rms)
        for data_i in data:
            freqs, model_flux_dens = integrated_spectrum(t1=data_i['t1'],t2=data_i['t2'], working_dir=working_dir)
            ax.plot(freqs,model_flux_dens*freqs*1.e-26,color=data_i['color'],ls=ls)


    ax.set(xscale='log',yscale='log',ylim=(1e-10,1e-7),xlim=(1e18,1e27))
    ax.set_ylabel(r"Flux [erg cm$^{-2}$ s$^{-1}$]", fontsize=12)
    ax.set_xlabel(r"$\nu$ [Hz]", fontsize=12)
    ax.grid(ls=':')
    ax.legend(fancybox=False, loc='lower left', columnspacing=0.8,
              # bbox_to_anchor=(0.5, 0.5),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
              shadow=False, ncol=1, fontsize=12, labelcolor='black',
              framealpha=0.0, borderaxespad=0.)
    # ax.set_yscale('log')
    ax.minorticks_on()
    ax.tick_params(axis='both', which='both', direction='in', labelsize=11)
    ax.set_title("GRB 190114C",fontsize=14)
    name="tophat_fs_ssc_grid_GRB_190114C"
    plt.savefig(os.getcwd() + '/figs/' + name + '.png', dpi=256)
    plt.savefig(os.getcwd() + '/figs/' + name + '.pdf')
    if plot: plt.show()
    plt.close(fig)
    plt.show()

    # load VHE objservations
    # with open('./magic_int1_points.txt') as f:
    #     lines = [line.rstrip('\n') for line in f]
    # new_lines = []
    # for line in lines:
    #     if not ((line == '') or (line[0] == '|') or (line[0] == "\\")):
    #         new_lines.append([float(val) for val in line.split()])
    # data = np.reshape(np.array(new_lines),newshape=(len(new_lines),len(new_lines[0])))
    # freqs, fluxes = load_observations('./magic_int1_points.txt')
    # fluxes = apply_ebl_absorption_to_data(freqs=freqs,fluxes=fluxes, z=z)

    # load EBL table and compute absorbtion
    # ebl = EBL()
    # ebl.load()
    # ebl.interpolate_fill_nans()
    # tau = ebl.set_interpolator(z=np.full_like(data[:,0],z),en=data[:,0])
    # atten = np.exp(-tau)
    # atten[~np.isfinite(atten)] = 1.
    # data[:,1]*=atten[:,0]
    # data[:,2]*=atten[:,0]
    # data[:,3]*=atten[:,0]

    # convert units
    # data[:,0] *= 2.417989242e14 # ev -> Hz
    # data[:,1:5] *= 1e-23 * 1e3# erg s^-1 cm^-2 Hz^-1 -> Jansky -> mJy

    # compute_rms_scores_for_runs(data_freq=freqs, data_flux=fluxes, t1=68, t2=110)

    df = pd.read_csv("./runs.csv",index_col=0)
    df.sort_values(by='rms',inplace=True)

    fig, ax = plt.subplots(ncols=1,nrows=1)
    ax.plot(freqs,fluxes,marker='o',color='gray',ls='none')

    for i in range(2):
        working_dir, rms = df.iloc[i]
        freqs, integ_spec = integrated_spectrum(t1=68, t2=110,working_dir=working_dir)
        ax.plot(freqs,integ_spec*freqs*1.e-26)

    # scores = np.sort(scores)
    # for i in range(5):
    #     print(f'Best index = {i} working_dir={working_dirs[i]} rmse={scores[i]}')
    # freqs, integ_spec = integrated_spectrum(t1=68,t2=110,working_dir=working_dirs[scores[0]])
    # ax.plot(freqs,integ_spec*freqs*1.e-26)

    # ax.errorbar(data[:,0],data[:,1],yerr=np.array(list(zip(data[:,2]-data[:,1],data[:,1]-data[:,3]))).T,fmt='.', ecolor = 'red')

    ax.set(xscale='log',yscale='log',ylim=(1e-10,1e-7),xlim=(1e18,1e27))
    ax.set_ylabel(r"Flux [erg cm$^{-2}$ s$^{-1}$]")
    ax.set_xlabel(r"$\nu$ [Hz]")
    plt.show()

# def preprocess(fname="./model_lcs.h5"):
#     working_dirs = glob(data_dir + '*')
#     dfile = h5py.File(fname,'w')
#     for i, working_dir in enumerate(working_dirs):
#         pba = PyBlastAfterglow(workingdir=working_dir+'/',readparfileforpaths=True)
#         freqs = pba.GRB.get_lc_freqs(unique=True,spec=False)
#         times = pba.GRB.get_lc_times(unique=True,spec=False)
#         spec = pba.GRB.get_lc()
#         dfile.create_dataset()
#
#         if i % 100 == 0:
#             print(f"Processing {i}/{len(working_dirs)} min_score={min_score}")
#     return (working_dirs, scores)

if __name__ == '__main__':
    plot()