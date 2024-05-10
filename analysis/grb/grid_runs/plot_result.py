import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy.spatial import cKDTree
from glob import glob

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

def plot():

    z = 0.4245 # https://www.aanda.org/articles/aa/pdf/2022/03/aa41788-21.pdf

    # load VHE objservations
    with open('./magic_int1_points.txt') as f:
        lines = [line.rstrip('\n') for line in f]
    new_lines = []
    for line in lines:
        if not ((line == '') or (line[0] == '|') or (line[0] == "\\")):
            new_lines.append([float(val) for val in line.split()])
    data = np.reshape(np.array(new_lines),newshape=(len(new_lines),len(new_lines[0])))

    # load EBL table and compute absorbtion
    ebl = EBL()
    ebl.load()
    ebl.interpolate_fill_nans()
    tau = ebl.set_interpolator(z=np.full_like(data[:,0],z),en=data[:,0])
    atten = np.exp(-tau)
    atten[~np.isfinite(atten)] = 1.
    data[:,1]*=atten[:,0]
    data[:,2]*=atten[:,0]
    data[:,3]*=atten[:,0]

    # convert units
    data[:,0] *= 2.417989242e14 # ev -> Hz
    # data[:,1:5] *= 1e-23 * 1e3# erg s^-1 cm^-2 Hz^-1 -> Jansky -> mJy


    fig, ax = plt.subplots(ncols=1,nrows=1)
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
        freqs, integ_spec = integrated_spectrum(t1=68,t2=110,working_dir=working_dir)
        interp_spec = interpolate.interp1d(freqs, integ_spec*freqs*1e-20,kind='linear')(data[:,0])
        rms = mean_squared_error(data[:,1], interp_spec)
        scores.append(rms)

        # ax.plot(data[:,0],interp_spec,marker='x')
        # ax.plot(freqs,integ_spec*freqs*1e-20)

    index = np.argmin(np.array(scores))
    print(f'Best index = {index} working_dir={working_dirs[index]} rmse={scores[index]}')
    freqs, integ_spec = integrated_spectrum(t1=68,t2=110,working_dir=working_dirs[index])
    ax.plot(freqs,integ_spec*freqs*1e-20)

    ax.errorbar(data[:,0],data[:,1],yerr=np.array(list(zip(data[:,2]-data[:,1],data[:,1]-data[:,3]))).T,fmt='.', ecolor = 'red')
    ax.set(xscale='log',yscale='log')
    ax.set_ylabel(r"$F_{\nu}$ [Jy]")
    ax.set_xlabel(r"$\nu$ [Hz]")
    plt.show()

if __name__ == '__main__':
    plot()