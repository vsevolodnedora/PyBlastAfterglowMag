import copy

import numpy as np
from sympy.polys.benchmarks.bench_solvers import k2

np.set_printoptions(precision=2)

import h5py
import pandas as pd
from glob import glob
import os
import re
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.colors import LogNorm, Normalize
from matplotlib.colors import BoundaryNorm
from matplotlib.patches import Patch
from matplotlib.ticker import MaxNLocator
import json

from scipy.stats import linregress
from scipy.interpolate import interp1d

from scipy.optimize import curve_fit
from scipy.optimize import minimize
from scipy.optimize import differential_evolution

import pysr
import sympy
from pysr import PySRRegressor
from sklearn.model_selection import train_test_split

# settings
#plt.style.use("fivethirtyeight")

# Supress runtime warning
import warnings
warnings.filterwarnings("ignore")

try:
    import package.src.PyBlastAfterglowMag as PBA
except ImportError:
    try:
        import PyBlastAfterglowMag as PBA
    except:
        raise ImportError("Cannot import PyBlastAfterglowMag")

DATA_PATH = str(__file__).split("analysis/kn/")[0] + "analysis/kn/kenta_data/"

# load the metadata
with open(DATA_PATH+"metadata.json") as json_file:
    json_data = json.load(json_file)
    print(json_data)
SIMS = pd.DataFrame.from_dict(json_data).T
SIMS.set_index("name")
# select only new simulations
df = SIMS[SIMS["given_time"] == "new"]

get_ej_data = lambda name : DATA_PATH+name+'/'+"ej_collated.h5"

# -------------------------------------------------------------------------------
EJ_TEXT_PATH = str(__file__).split("analysis/kn/")[0] + "analysis/kn/ejecta/output/"
df_text = pd.read_csv(EJ_TEXT_PATH+"ejecta_fasttail_vals_at_massmax.csv",index_col=0)

def find_nearest_index(array, value):
    ''' Finds index of the value in the array that is the closest to the provided one '''
    idx = (np.abs(array - value)).argmin()
    return idx

def get_2d_ek_(ej:PBA.id_kenta.EjectaData, text:float, new_vinf_len:int=80)\
        ->tuple[np.ndarray,np.ndarray,np.ndarray, np.ndarray]:
    thetas = ej.get_theta()
    vinf = ej.get_vinf()
    masses = ej.get(v_n="mass",text=text)

    vinf = vinf[:-1]
    vinf_c = 0.5 * (vinf[1:] + vinf[:-1])
    mom = PBA.utils.MomFromBeta(vinf_c)
    ctheta = 0.5 * (thetas[1:] + thetas[:-1])
    # ctheta = 0.5 * (thetas[1:] + thetas[:-1])
    mass = 0.5 * (masses[:,1:] + masses[:,:-1]) [:,:-1]
    # mass = 0.5 * (mass[1:,:] + mass[:-1,:])
    ek = mass * PBA.cgs.solar_m * (vinf_c[np.newaxis,:] * PBA.cgs.c)**2

    # mask = ctheta < np.pi/2.
    return (mom, ctheta, ek.T, mass.T)
def get_2d_ek(ej:PBA.id_kenta.EjectaData, text:float, new_vinf_len:int or None=80) \
        ->tuple[np.ndarray,np.ndarray,np.ndarray, np.ndarray]:
    thetas = ej.get_theta()
    vinf = ej.get_vinf()
    masses = ej.get(v_n="mass",text=text)

    vinf, thetas, masses = PBA.id_kenta.reinterpolate_hist2(
        vinf, thetas, masses, new_theta_len=None, new_vinf_len=new_vinf_len, mass_conserving=True
    )

    vinf_c = 0.5 * (vinf[1:] + vinf[:-1])
    mom = PBA.utils.MomFromBeta(vinf_c)
    ctheta = 0.5 * (thetas[1:] + thetas[:-1])
    ek = masses * PBA.cgs.solar_m * (vinf_c[np.newaxis,:] * PBA.cgs.c)**2
    return (mom, ctheta, ek.T, masses.T)

    # # masses = masses.T
    # masses=masses[:-1,:]  # remove vinf > 1
    #
    # vinf = vinf[:-1]
    # vinf_c = 0.5 * (vinf[1:] + vinf[:-1])
    # mom = PBA.utils.MomFromBeta(vinf_c)
    # ctheta = 0.5 * (thetas[1:] + thetas[:-1])
    # # ctheta = 0.5 * (thetas[1:] + thetas[:-1])
    # mass = 0.5 * (masses[:,1:] + masses[:,:-1]) [:,:-1]
    # # mass = 0.5 * (mass[1:,:] + mass[:-1,:])
    # ek = mass * PBA.cgs.solar_m * (vinf_c[np.newaxis,:] * PBA.cgs.c)**2
    #
    # # mask = ctheta < np.pi/2.
    # return (mom, ctheta, ek.T, mass.T)

class FitBase:

    name = 'Base'

    def __init__(self):
        pass


    def compute_chi2(self, x, y_values, y_fit, n):
        # x = 10**x
        # y_values = 10**y
        # y_fit = 10**y_pred

        # The standard deviation is the square root of the average of the squared
        # deviations from the mean, i.e., ``std = sqrt(mean(x))``, where
        # ``x = abs(a - a.mean())**2``.
        ss_res = np.sum(((y_values - y_fit)) ** 2)
        reduced_chi_squared = ss_res / (len(y_values) - n)
        return reduced_chi_squared

    def compute_r2(self,x,y_values,y_fit,n):
        # x=10**x
        # y_values=10**y
        # y_fit =10**y_pred

        ss_res = np.sum((y_values - y_fit) ** 2)
        ss_tot = np.sum((y_values - np.mean(y_values)) ** 2)
        r_squared = 1 - (ss_res / ss_tot)
        adjusted_r_squared = 1 - (1 - r_squared) * (len(y_values) - 1) / (len(y_values) - n - 1)
        return adjusted_r_squared

    def compute_sse_res(self,x,y_values,y_fit):
        return np.sum((y_values - y_fit) ** 2)

    def print_table(self,infos:pd.DataFrame):
        info = copy.deepcopy(infos)
        info.to_csv(os.getcwd()+f'/output/piecewise_line_{self.name}.csv',index=True)
        info["y0"] = ["${}$".format(PBA.utils.latex_float(y0)) for y0 in info["y0"]]
        info["r2"] = [f"${r2:.3f}$" for r2 in info["r2"]]
        info["chi2"] = [f"${chi2:.3f}$" for chi2 in info["chi2"]]
        # info["chi2"] = [np.log10(chi2) for chi2 in info["chi2"]]
        # info["chi2"] = ["${}$".format(PBA.utils.latex_float(chi2)) for chi2 in info["chi2"]]
        # info["sse"] = ["${}$".format(PBA.utils.latex_float(sse)) for sse in info["sse"]]
        # del info["chi2"]

        keys = ["label"]
        keys += copy.deepcopy( self.keys )
        keys += ['sse']
        # keys += ['chi2','r2', 'sse']


        print(info[keys].to_latex(float_format="%.2f",index=False))

class Fit1D_2seg(FitBase):
    name = "2segFit"
    keys = ["x0","y0","k1","k2"]
    labels = [r"$\mathcal{M}_0$" , r"$\mathcal{E}_0$" , r"$k_1$" , r"$k_2$" ]
    def __init__(self):
        super().__init__()
        pass

    @staticmethod
    def piecewise_linear(x, x0, y0, k1, k2):
        condlist = [x < x0,
                    x >= x0]
        funclist = [lambda x: k1 * x + y0 - k1 * x0,
                    lambda x: k2 * x + y0 - k2 * x0 # k2 * x + (y0 - k1 * x0) + (k1 - k2) * x0, #
                    # lambda x: k3 * x + y0 - k1*x0 + k1*x0 - k2*x0 + k2*x1 - k3*x1 #k3 * x + (y0 - k1 * x0) + (k1 - k2) * x0 + (k2 - k3) * x1
                    ]
        return np.piecewise(x, condlist, funclist)
    @staticmethod
    def piecewise_power(x, x0, y0, k1, k2):
        condlist = [x < x0,
                    x >= x0]
        funclist = [lambda x: y0*(x/x0)**k1,# k1 * x + y0 - k1 * x0,
                    lambda x: y0*(x/x0)**k2#k2 * x + y0 - k2 * x0, # k2 * x + (y0 - k1 * x0) + (k1 - k2) * x0, #
                    # lambda x: k3 * x + y0 - k1*x0 + k1*x0 - k2*x0 + k2*x1 - k3*x1 #k3 * x + (y0 - k1 * x0) + (k1 - k2) * x0 + (k2 - k3) * x1
                    ]
        return np.piecewise(x, condlist, funclist)

    def fit(self, x_values:np.ndarray, y_values:np.ndarray, undo_x, undo_y):
        # Find the index of the maximum y-value
        max_index = np.argmax(y_values)
        # Corresponding x-value for the maximum point
        x_max = x_values[max_index]
        # Finding the closest index to x=0
        zero_index = np.argmin(np.abs(x_values))
        # The x-value closest to zero
        x_zero = x_values[zero_index]
        # Initial guesses for the parameters
        initial_guesses = [0.5*x_max, 0.9*max(y_values), -7, -10]

        # Fit the model to the data
        # popt, _ = curve_fit(piecewise_linear, x_values, y_values, p0=initial_guesses,absolute_sigma=True)
        def residuals(x, *args):
            # res = np.zeros_like(x)
            # mask1 = x<x[np.argmax(y_values)]
            # res[mask1] = (piecewise_linear(x[mask1], *args) - y_values[mask1])
            # mask2 = x>=x[np.argmax(y_values)]
            # res[mask2] = (piecewise_linear(x[mask2], *args) - y_values[mask2])
            res = (self.piecewise_linear(x, *args) - y_values)

            res = (self.piecewise_linear(x, *args) - y_values) ** 2 * y_values ** 3
            mask = np.where((x>args[1])&(x<args[2]))
            res[mask] *= 1.2

            # res[mask] = (10**piecewise_linear(x, *args) - 10**y_values)[mask]
            # res[np.where(x<args[0])]**4
            # res[np.where((x>args[1])&(x<args[2]))]*=2
            # res[np.where((x>args[2]))]*=2
            # res[np.argmax(y_values)-4:np.argmax(y_values)+4] *= 10
            return res
        popt, _ = curve_fit(residuals, x_values, np.zeros_like(y_values), p0=initial_guesses,absolute_sigma=False)
        # print(popt)
        # Extract the optimized parameters
        x0_opt, y0_opt, k1_opt, k2_opt = popt

        # Generate x-values for plotting the fitted function
        x_fit = x_values#np.linspace(np.min(x_values), np.max(x_values), 1000)
        y_fit = self.piecewise_linear(x_fit, *popt)

        # reduced_chi_squared = self.compute_chi2(undo_x(x_values),undo_y(y_values),undo_y(y_fit),len(popt))
        # adjusted_r_squared = self.compute_r2(undo_x(x_values),undo_y(y_values),undo_y(y_fit),len(popt))
        # sse_res = self.compute_sse_res(undo_x(x_values),undo_y(y_values),undo_y(y_fit))
        reduced_chi_squared = self.compute_chi2(x_values,y_values,y_fit,len(popt))
        adjusted_r_squared = self.compute_r2(x_values,y_values,y_fit,len(popt))
        sse_res = self.compute_sse_res(x_values,y_values,y_fit)


        print(f"Chi2={reduced_chi_squared} R2={adjusted_r_squared} SSE={sse_res}")

        # Plot the original data and the fitted curve
        # plt.figure(figsize=(10, 6))
        # plt.scatter(x_values, y_values, color='blue', label='Data Points')
        # plt.plot(x_fit, y_fit, 'r-', label='Fitted Piece-wise Linear Function')
        # plt.title('Fit of Piece-wise Linear Function to Data')
        # plt.xlabel('x')
        # plt.ylabel('y')
        # plt.axvline(x=x0_opt, color='green', linestyle='--', label='Breakpoint at x0')
        # plt.axvline(x=x1_opt, color='purple', linestyle='--', label='Breakpoint at x_zero')
        # plt.legend()
        # plt.grid(True)
        # plt.show()

        coeffs = [undo_x(x0_opt),undo_y(y0_opt),k1_opt,k2_opt]
        return x_fit, y_fit, coeffs, dict(chi2=reduced_chi_squared, r2=adjusted_r_squared, sse=sse_res,
                                          x0=coeffs[0],y0=coeffs[1],
                                          k1=coeffs[2],k2=coeffs[3])

class Fit1D_3seg(FitBase):
    name = "3segFit"
    keys = ["x0","x1","y0","k1","k2","k3"]
    labels = [r"$\mathcal{M}_0$" , r"$\mathcal{M}_1$" , r"$\mathcal{E}_0$" , r"$k_1$" , r"$k_2$" , r"$k_3$"]
    def __init__(self):
        super().__init__()
        pass

    @staticmethod
    def piecewise_linear(x, x0, x1, y0, k1, k2, k3):
        condlist = [x < x0,
                    (x >= x0) & (x < x1),
                    x >= x1]
        funclist = [lambda x: k1 * x + y0 - k1 * x0,
                    lambda x: k2 * x + y0 - k2 * x0, # k2 * x + (y0 - k1 * x0) + (k1 - k2) * x0, #
                    # lambda x: k3 * x + y0 - k1*x0 + k1*x0 - k2*x0 + k2*x1 - k3*x1 #k3 * x + (y0 - k1 * x0) + (k1 - k2) * x0 + (k2 - k3) * x1
                    lambda x: k3 * x + y0 - k2*x0 + k2*x1 - k3*x1 #k3 * x + (y0 - k1 * x0) + (k1 - k2) * x0 + (k2 - k3) * x1
                    ]
        return np.piecewise(x, condlist, funclist)
    @staticmethod
    def piecewise_power(x, x0, x1, y0, k1, k2, k3):
        condlist = [x < x0,
                    (x >= x0) & (x < x1),
                    x >= x1]
        funclist = [lambda x: y0*(x/x0)**k1,# k1 * x + y0 - k1 * x0,
                    lambda x: y0*(x/x0)**k2,#k2 * x + y0 - k2 * x0, # k2 * x + (y0 - k1 * x0) + (k1 - k2) * x0, #
                    # lambda x: k3 * x + y0 - k1*x0 + k1*x0 - k2*x0 + k2*x1 - k3*x1 #k3 * x + (y0 - k1 * x0) + (k1 - k2) * x0 + (k2 - k3) * x1
                    lambda x: y0*(x/x1)**k3 * (x1/x0)**k2 #k3 * x + y0 - k2*x0 + k2*x1 - k3*x1 #k3 * x + (y0 - k1 * x0) + (k1 - k2) * x0 + (k2 - k3) * x1
                    ]
        return np.piecewise(x, condlist, funclist)

    def fit(self, x_values:np.ndarray, y_values:np.ndarray, undo_x, undo_y):
        # Find the index of the maximum y-value
        max_index = np.argmax(y_values)
        # Corresponding x-value for the maximum point
        x_max = x_values[max_index]
        # Finding the closest index to x=0
        zero_index = np.argmin(np.abs(x_values))
        # The x-value closest to zero
        x_zero = x_values[zero_index]
        # Initial guesses for the parameters
        initial_guesses = [x_max, max(x_values)*0.75, max(y_values), -7, -7, -20]

        # Fit the model to the data
        # popt, _ = curve_fit(piecewise_linear, x_values, y_values, p0=initial_guesses,absolute_sigma=True)
        def residuals(x, *args):
            # res = np.zeros_like(x)
            # mask1 = x<x[np.argmax(y_values)]
            # res[mask1] = (piecewise_linear(x[mask1], *args) - y_values[mask1])
            # mask2 = x>=x[np.argmax(y_values)]
            # res[mask2] = (piecewise_linear(x[mask2], *args) - y_values[mask2])
            res = (self.piecewise_linear(x, *args) - y_values)

            res = (self.piecewise_linear(x, *args) - y_values) ** 2 * y_values ** 3
            mask = np.where((x>args[1])&(x<args[2]))
            res[mask] *= 1.2

            # res[mask] = (10**piecewise_linear(x, *args) - 10**y_values)[mask]
            # res[np.where(x<args[0])]**4
            # res[np.where((x>args[1])&(x<args[2]))]*=2
            # res[np.where((x>args[2]))]*=2
            # res[np.argmax(y_values)-4:np.argmax(y_values)+4] *= 10
            return res
        popt, _ = curve_fit(residuals, x_values, np.zeros_like(y_values), p0=initial_guesses,absolute_sigma=False)
        # print(popt)
        # Extract the optimized parameters
        x0_opt, x1_opt, y0_opt, k1_opt, k2_opt, k3_opt = popt

        # Generate x-values for plotting the fitted function
        x_fit = x_values#np.linspace(np.min(x_values), np.max(x_values), 1000)
        y_fit = self.piecewise_linear(x_fit, *popt)

        # reduced_chi_squared = self.compute_chi2(undo_x(x_values),undo_y(y_values),undo_y(y_fit),len(popt))
        # adjusted_r_squared = self.compute_r2(undo_x(x_values),undo_y(y_values),undo_y(y_fit),len(popt))
        # sse_res = self.compute_sse_res(undo_x(x_values),undo_y(y_values),undo_y(y_fit))
        reduced_chi_squared = self.compute_chi2(x_values,y_values,y_fit,len(popt))
        adjusted_r_squared = self.compute_r2(x_values,y_values,y_fit,len(popt))
        sse_res = self.compute_sse_res(x_values,y_values,y_fit)

        print(f"Chi2={reduced_chi_squared} R2={adjusted_r_squared} SSE={sse_res}")

        # Plot the original data and the fitted curve
        # plt.figure(figsize=(10, 6))
        # plt.scatter(x_values, y_values, color='blue', label='Data Points')
        # plt.plot(x_fit, y_fit, 'r-', label='Fitted Piece-wise Linear Function')
        # plt.title('Fit of Piece-wise Linear Function to Data')
        # plt.xlabel('x')
        # plt.ylabel('y')
        # plt.axvline(x=x0_opt, color='green', linestyle='--', label='Breakpoint at x0')
        # plt.axvline(x=x1_opt, color='purple', linestyle='--', label='Breakpoint at x_zero')
        # plt.legend()
        # plt.grid(True)
        # plt.show()

        coeffs = [undo_x(x0_opt),undo_x(x1_opt),undo_y(y0_opt),k1_opt,k2_opt,k3_opt]
        return x_fit, y_fit, coeffs, dict(chi2=reduced_chi_squared, r2=adjusted_r_squared,sse=sse_res,
                                          x0=coeffs[0],x1=coeffs[1],y0=coeffs[2],
                                          k1=coeffs[3],k2=coeffs[4],k3=coeffs[5])

class Fit1D_4seg(FitBase):
    name = "4segFit"
    keys = ["x0","x1","x2","y0","k1","k2","k3"]
    labels = [r"$\mathcal{M}_0$" , r"$\mathcal{M}_1$" , r"$\mathcal{M}_2$" , r"$\mathcal{E}_0$" , r"$k_1$" , r"$k_2$" , r"$k_3$"]
    def __init__(self):
        super().__init__()
        pass

    @staticmethod
    def piecewise_linear(x, x0, x1, x2, y0, k1, k2, k3):
        condlist = [x < x0,
                    (x >= x0) & (x < x1),
                    (x >= x1) & (x < x2),
                    x >= x2]
        funclist = [lambda x: k1 * x + y0 - k1 * x0,
                    lambda x: y0,
                    lambda x: k2 * x + y0 - k2 * x1,
                    # k2 * x + (y0 - k1 * x0) + (k1 - k2) * x0, #
                    # lambda x: k3 * x + y0 - k1*x0 + k1*x0 - k2*x0 + k2*x1 - k3*x1 #k3 * x + (y0 - k1 * x0) + (k1 - k2) * x0 + (k2 - k3) * x1
                    lambda x: k3 * x + y0 - k2*x1 + k2*x2 - k3*x2 #k3 * x + (y0 - k1 * x0) + (k1 - k2) * x0 + (k2 - k3) * x1
                    ]
        return np.piecewise(x, condlist, funclist)
    @staticmethod
    def piecewise_power(x, x0, x1, x2, y0, k1, k2, k3):
        condlist = [x < x0,
                    (x >= x0) & (x < x1),
                    (x >= x1) & (x < x2),
                    x >= x2]
        funclist = [lambda x: y0*(x/x0)**k1, # k1 * x + y0 - k1 * x0,
                    lambda x: y0,
                    lambda x: y0*(x/x1)**k2, #k2 * x + y0 - k2 * x0, # k2 * x + (y0 - k1 * x0) + (k1 - k2) * x0, #
                    # lambda x: k3 * x + y0 - k1*x0 + k1*x0 - k2*x0 + k2*x1 - k3*x1 #k3 * x + (y0 - k1 * x0) + (k1 - k2) * x0 + (k2 - k3) * x1
                    lambda x: y0*(x/x2)**k3 * (x2/x1)**k2 #k3 * x + y0 - k2*x0 + k2*x1 - k3*x1 #k3 * x + (y0 - k1 * x0) + (k1 - k2) * x0 + (k2 - k3) * x1
                    ]
        return np.piecewise(x, condlist, funclist)

    def fit(self, x_values:np.ndarray, y_values:np.ndarray, undo_x, undo_y):
        # Find the index of the maximum y-value
        max_index = np.argmax(y_values)
        # Corresponding x-value for the maximum point
        x_max = x_values[max_index]
        # Finding the closest index to x=0
        zero_index = np.argmin(np.abs(x_values))
        # The x-value closest to zero
        x_zero = x_values[zero_index]
        # Initial guesses for the parameters
        initial_guesses = [x_max, max(x_values)*0.05, max(x_values)*0.75, max(y_values), -7, -7, -20]

        # Fit the model to the data
        # popt, _ = curve_fit(piecewise_linear, x_values, y_values, p0=initial_guesses,absolute_sigma=True)
        def residuals(x, *args):
            # res = np.zeros_like(x)
            # mask1 = x<x[np.argmax(y_values)]
            # res[mask1] = (piecewise_linear(x[mask1], *args) - y_values[mask1])
            # mask2 = x>=x[np.argmax(y_values)]
            # res[mask2] = (piecewise_linear(x[mask2], *args) - y_values[mask2])
            res = (Fit1D_4seg.piecewise_linear(x, *args) - y_values) ** 2

            res = (Fit1D_4seg.piecewise_linear(x, *args) - y_values) ** 2 * y_values ** 3
            mask = np.where((x>args[1])&(x<args[2]))
            res[mask] *= 1.1

            # res[mask] = (10**piecewise_linear(x, *args) - 10**y_values)[mask]
            # res[np.where(x<args[0])]**4
            # res[np.where((x>args[1])&(x<args[2]))]*=2
            # res[np.where((x>args[2]))]*=2
            # res[np.argmax(y_values)-4:np.argmax(y_values)+4] *= 10
            return res
        popt, _ = curve_fit(residuals, x_values, np.zeros_like(y_values), p0=initial_guesses,absolute_sigma=False)
        # print(popt)
        # Extract the optimized parameters
        x0_opt, x1_opt, x2_opt, y0_opt, k1_opt, k2_opt, k3_opt = popt

        # Generate x-values for plotting the fitted function
        x_fit = x_values#np.linspace(np.min(x_values), np.max(x_values), 1000)
        y_fit = Fit1D_4seg.piecewise_linear(x_fit, *popt)

        # reduced_chi_squared = self.compute_chi2(undo_x(x_values),undo_y(y_values),undo_y(y_fit),len(popt))
        # adjusted_r_squared = self.compute_r2(undo_x(x_values),undo_y(y_values),undo_y(y_fit),len(popt))
        # sse_res = self.compute_sse_res(undo_x(x_values),undo_y(y_values),undo_y(y_fit))
        reduced_chi_squared = self.compute_chi2(x_values,y_values,y_fit,len(popt))
        adjusted_r_squared = self.compute_r2(x_values,y_values,y_fit,len(popt))
        sse_res = self.compute_sse_res(x_values,y_values,y_fit)

        print(f"Chi2={reduced_chi_squared} R2={adjusted_r_squared} SSE={sse_res}")
        # Plot the original data and the fitted curve
        # plt.figure(figsize=(10, 6))
        # plt.scatter(x_values, y_values, color='blue', label='Data Points')
        # plt.plot(x_fit, y_fit, 'r-', label='Fitted Piece-wise Linear Function')
        # plt.title('Fit of Piece-wise Linear Function to Data')
        # plt.xlabel('x')
        # plt.ylabel('y')
        # plt.axvline(x=x0_opt, color='green', linestyle='--', label='Breakpoint at x0')
        # plt.axvline(x=x1_opt, color='purple', linestyle='--', label='Breakpoint at x_zero')
        # plt.legend()
        # plt.grid(True)
        # plt.show()

        coeffs = [undo_x(x0_opt),undo_x(x1_opt),undo_x(x2_opt),undo_y(y0_opt),k1_opt,k2_opt,k3_opt]
        return x_fit, y_fit, coeffs, dict(chi2=reduced_chi_squared, r2=adjusted_r_squared, sse=sse_res,
                                          x0=coeffs[0],x1=coeffs[1],x2=coeffs[2],y0=coeffs[3],
                                          k1=coeffs[4],k2=coeffs[5],k3=coeffs[6])


def fit_data(x, y, undo_x, undo_y, name:str, fitting_obj, save_res:bool=True):
    # np.savetxt(os.getcwd()+f'/output/+name+"log_mom_log_ek.txt",X=np.column_stack((x,y)),fmt="%.3f")

    # o_fit = Fit1D()
    x_fit, y_fit, coeffs, fit_dict = fitting_obj.fit(x, y, undo_x, undo_y)

    # using fitting function
    # mom = np.linspace(*mom_lim, n_shells[sim])
    vinf_ = np.linspace(0.005,1, len(y)+1)
    mom = PBA.MomFromBeta(vinf_)[:-1]
    ek = fitting_obj.piecewise_power(mom, *coeffs)
    # mass = ek / (PBA.cgs.solar_m * PBA.cgs.c**2 * PBA.BetaFromMom(mom)**2)#piecewise_power(mom, *coeffs) / PBA.cgs.solar_m
    # plt.close()
    # plt.loglog(x_fit,y_fit,marker='.',ls='None')
    # plt.loglog(np.log10(mom), np.log10(ek),marker='x',ls='None')
    # plt.show()

    if save_res:
        np.savetxt(os.getcwd()+f'/output/'+name+f"_log_mom_log_ek_sph_and_fit_{fitting_obj.name}.txt",
                   X=np.column_stack((x,y, np.log10(mom), np.log10(ek))),
                   fmt="%.3f")

    # return (np.log10(mom), np.log10(ek), fit_dict)
    return (x_fit, y_fit, fit_dict)


def plot_all_sim_ejecta_mass_evol(crit=None,yscale="linear",ylim=(0,0.04),title="Volume integrated ejecta mass",
                             figname="figname"):
    fig, ax = plt.subplots(ncols=1,nrows=1,figsize=(7,4))
    #ax2 = ax.twinx()
    for idx, sim_dic in enumerate(df.iterrows()):
        sim_dic = sim_dic[1]
        ej = PBA.id_kenta.EjectaData(get_ej_data(sim_dic["name"]),verbose=True)
        if (crit is None):
            mass = ej.total_mass()
        else:
            mass = ej.total_mass_vs_text(crit=crit)
        tmerg = sim_dic["tmerg"]
        if sim_dic["given_time"] == "new":
            ax.plot(ej.getText()-tmerg, mass,
                    color=sim_dic["color"], ls=sim_dic["ls"], label=sim_dic["label"])
        else:
            ax.plot(ej.getText()-tmerg, mass, color=sim_dic["color"], ls=sim_dic["ls"],
                    label=sim_dic["label"])
        # data = PBA.id_kenta.Data(fpath_rhomax=sim_dic["datadir"]+sim_dic["rhomax"],
        #                          fpath_mdot=sim_dic["datadir"]+sim_dic["mdot_extract"])
        # ax.plot(ej.getText(), ej.total_mass_fasttail(),color=color_pal[sim_dic["idx"]],
        #         label=r"$\Gamma\beta>1$",ls="--")
        # ax2.plot(*data.get_rhomax(),color=color_pal[sim_dic["idx"]],label=r"$\rho_{\rm max}$",ls=":")
        # ax2.plot(*data.get_mdot(),color=color_pal[sim_dic["idx"]],label=r"$\dot{M}$",ls="-.")
    ax.tick_params(labelsize=12)
    ax.set_ylim(*ylim)
    ax.set_yscale(yscale)
    ax.legend(fontsize=12,ncol=2)
    # ax.set_yscale("log")
    ax.set_xlabel(r"$t_{\rm ext} - t_{\rm merg}$ [ms]",fontsize=14)
    ax.set_ylabel(r"Ejecta mass $[M_{\odot}]$",fontsize=14)
    ax.set_title(title,fontsize=14)
    # ax.grid(which="both",axis="both")
    plt.tight_layout()
    plt.savefig(os.getcwd()+f'/figs/{figname}.png',dpi=256)
    plt.savefig(os.getcwd()+f'/figs/{figname}.pdf')
    plt.show()


def plot_all_sim_ejecta_mass(
        o_fit, xlim=(1e-3, 4),ylim0=(1e43, 1e51),ylim1=(1e43, 1e51),ylim2=(-.5,.5),
        figname="figname",figname_coeffs="figname", sim: str or None=None, plot_fit_coeffs:bool=True):

    do_cumulative = True
    log_type = 10
    log_type_y = 10

    get_cumulative = lambda val : np.cumsum(val[::-1])[::-1] if do_cumulative else val
    do_log = lambda val : (val if not log_type else (np.log2(val) if log_type==2 else np.log10(val)))
    un_log = lambda val : (val if not log_type else (2**(val) if log_type==2 else 10**(val)))

    do_log_y = lambda val : (val if not log_type_y else (np.log2(val) if log_type_y==2 else np.log10(val)))
    un_log_y = lambda val : (val if not log_type_y else (2**(val) if log_type_y==2 else 10**(val)))

    fig, axes = plt.subplots(ncols=1,nrows=3,figsize=(4.6,2*3.2),
                             layout='constrained',sharex='col',#sharex='col',sharey='row',
                             gridspec_kw={'height_ratios': [2,2,1]})

    # df_fit = pd.read_csv(os.getcwd()+'/'+'piecewise_line_fits.csv',index_col=0)

    infos = dict()
    idx = 0

    x1, x2, y1, y2 = -1.5, -0.9, -2.5, -1.9  # subregion of the original image
    axins = axes[1].inset_axes(
        [0.05, 0.1, 0.5, 0.5],
        xlim=(x1, x2), ylim=(y1, y2), xticklabels=[], yticklabels=[])

    # sub region of the original image
    x1, x2, y1, y2 = 8e-2, 5e-1, 2e48, 4e49
    axins.set_xlim(x1, x2)
    axins.set_ylim(y1, y2)
    axins.set_xscale('log')
    axins.set_yscale('log')
    axins.grid(ls=':')
    axins.set_xticklabels([])
    axins.set_yticklabels([])
    axins.set_xlabel('')
    axins.xaxis.set_visible(False)
    axins.set_xticks([])
    axins.minorticks_on()
    axins.tick_params(axis='both', which='both', labelleft=False,
                      labelright=True, tick1On=True, tick2On=True,
                      labelsize=12,
                      direction='in',
                      bottom=True, top=True, left=True, right=True)
    axes[1].indicate_inset_zoom(axins, edgecolor="black")

    # for idx, sim_dic in enumerate(df.iterrows()):
    # for (sim_dic,fit_dic) in zip(df.iterrows(), df_fit.iterrows()):
    # o_fit = Fit1D_4seg()
    # o_fit = Fit1D_3seg()

    for (name, sim_dic) in df.iterrows():
        text_dict = df_text.loc[name]
        # sim_dic = sim_dic[1]
        if sim and name != sim:
            continue

        # name = sim_dic[0]
        ej = PBA.id_kenta.EjectaData(get_ej_data(sim_dic['name']),verbose=True)
        mom, _, ek, mass = get_2d_ek(ej=ej,
                               # text=get_text_max_mass(ej=ej,crit='fast')
                               text=float(sim_dic['tmerg'])+float(text_dict['text'])
                               )
        # ek = np.sum(mass,axis=1)
        ek = np.sum(ek,axis=1)
        # ek = np.sum(mass,axis=1)*PBA.cgs.solar_m#*(PBA.cgs.c**2*PBA.BetaFromMom(mom)**2)

        # ej_mass = ej.get(v_n="mass",text=get_text_max_mass(ej=ej,crit='fast'))
        # ek = np.cumsum(np.sum(ek,axis=0)[::-1] * PBA.cgs.solar_m * vinf_c**2 * PBA.cgs.c**2)[::-1]/np.sum(ek * PBA.cgs.solar_m * vinf_c**2 * PBA.cgs.c**2)
        # ek = np.cumsum(ek[::-1])[::-1]#/np.sum(ek)


        l_mom, l_ek = do_log(mom), do_log_y(ek)
        mask = np.isfinite(l_ek)
        l_mom = l_mom[mask]
        l_ek = l_ek[mask]

        # N = 5
        # l_ek = np.convolve(l_ek, np.ones(N)/N, mode='valid')
        # l_mom = np.convolve(l_mom, np.ones(N)/N, mode='valid')


        axes[0].plot(un_log(l_mom), get_cumulative(un_log_y(l_ek)), color=sim_dic["color"], ls=sim_dic["ls"], label=sim_dic["label"],lw=1.2)#, lw=0.7, drawstyle='steps')
        axes[1].plot(un_log(l_mom), un_log_y(l_ek), color=sim_dic["color"], ls=sim_dic["ls"], label=sim_dic["label"],lw=1.2)#, lw=0.7, drawstyle='steps')


        l_mom_pred, l_ek_pred, info = fit_data(l_mom, l_ek, un_log, un_log_y,
                                               fitting_obj=o_fit, name=sim_dic["name"])

        info["label"] = sim_dic["label"]
        infos[name] = info
        axes[0].plot(un_log(l_mom_pred), get_cumulative(un_log_y(l_ek_pred)), color=sim_dic["color"], ls=sim_dic["ls"],lw=.6)#, lw=0.7, drawstyle='steps')
        axes[1].plot(un_log(l_mom_pred), un_log_y(l_ek_pred), color=sim_dic["color"], ls=sim_dic["ls"],lw=.6)#, lw=0.7, drawstyle='steps')
        axes[2].plot(un_log(l_mom_pred), (l_ek-l_ek_pred), color=sim_dic["color"], ls=sim_dic["ls"], lw=0.7)

        axins.plot(un_log(l_mom), un_log_y(l_ek), color=sim_dic["color"], ls=sim_dic["ls"], label=sim_dic["label"],lw=1.2)
        axins.plot(un_log(l_mom_pred), un_log_y(l_ek_pred), color=sim_dic["color"], ls=sim_dic["ls"],lw=.6)#, lw=0.7, drawstyle='steps')

        '''
        fit_dic = df_fit.loc[name]
        y_pred = piecewise_power(10**l_mom,x0=fit_dic['x0'],x1=fit_dic['x1'],y0=fit_dic['y0'],
                                 k1=fit_dic['k1'],k2=fit_dic['k2'],k3=fit_dic['k3'])
        # y_pred = 10**piecewise_linear(l_mom,x0=fit_dic['x0'],x1=fit_dic['x1'],y0=fit_dic['y0'],
        #                          k1=fit_dic['k1'],k2=fit_dic['k2'],k3=fit_dic['k3'])
        axes[0].plot(10**l_mom, get_cumulative(y_pred), color='gray', ls=sim_dic["ls"],lw=1.2)#, lw=0.7, drawstyle='steps')
        axes[1].plot(10**l_mom, y_pred, color='gray', ls=sim_dic["ls"],lw=1.2)#, lw=0.7, drawstyle='steps')
        '''

        # axes[0].plot(mom, get_cumulative(ek), color=sim_dic["color"], ls=sim_dic["ls"], label=sim_dic["label"],lw=1.2)#, lw=0.7, drawstyle='steps')

        # mom_pred, ek_pred, info = fit_data(mom, ek,name=sim_dic["name"])
        # infos[sim_dic["name"]] = info
        # axes[0].plot(mom_pred, get_cumulative(ek_pred), color=sim_dic["color"], ls=sim_dic["ls"],lw=.6)#, lw=0.7, drawstyle='steps')
        # axes[1].plot(mom_pred, (ek-ek_pred) / ek, color=sim_dic["color"], ls=sim_dic["ls"], lw=0.7)

        # break

    axes[0].plot([1e-4,1e-2], [1e39,1e41], color='gray', ls='-', label='Simulation',lw=1.2)#, lw=0.7, drawstyle='steps')
    axes[0].plot([1e-4,1e-2], [1e39,1e41], color='gray', ls='-',label='Fit',lw=.6)#, lw=0.7, drawstyle='steps')

    infos = pd.DataFrame.from_dict(infos).T
    o_fit.print_table(infos=infos)

    for ax in axes:
        ax.tick_params(labelsize=12)
        # ax.set_ylim(*ylim)

        ax.tick_params(labelsize=12,axis='both', which='both',direction='in',tick1On=True, tick2On=True)
        ax.minorticks_on()
        ax.grid(ls=':',lw=0.8)
        ax.set_xscale('log')

    #for idx in range(len(df)):
    axes[0].set_yscale('log')
    axes[1].set_yscale('log')
    # axes[2].set_yscale('log')
    axes[0].set_ylim(*ylim0)
    axes[1].set_ylim(*ylim1)
    axes[2].set_ylim(*ylim2)
    axes[-1].set_xlim(*xlim)

    # ax.set_yscale(yscale)
    # ax.set_xscale(xscale)

    n = 2
    ax = axes[0]
    han, lab = ax.get_legend_handles_labels()
    ax.add_artist(ax.legend(han[:-1 * n], lab[:-1 * n],
                            **dict(fancybox=False,loc= 'lower left',columnspacing=0.4,
                            #"bbox_to_anchor": (0.5, 1.2),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
                            shadow=False, ncol= 1, fontsize= 12,framealpha=0., borderaxespad= 0., frameon=False)))
    ax.add_artist(ax.legend(han[len(han) - n:], lab[len(lab) - n:],
                            **dict(fancybox=False,loc= 'upper right',columnspacing=0.4,
                            #"bbox_to_anchor": (0.5, 1.2),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
                            shadow=False, ncol= 1, fontsize= 12,framealpha=0., borderaxespad= 0., frameon=False)))
    # ax.set_yscale("log")
    axes[-1].set_xlabel(r"$\Gamma\beta$",fontsize=12)
    axes[0].set_ylabel(r"$E_{\rm k} (> \Gamma \beta)$ [erg]",fontsize=12)
    axes[1].set_ylabel(r"$E_{\rm k}$ [erg]",fontsize=12)
    # ax.set_title(title,fontsize=12)
    # axes[-1].set_ylabel(r"$\Delta E_{\rm k}$ [erg]",fontsize=12)
    axes[-1].set_ylabel(r"$\log_{10}(E_{\rm k;\,sph})-\log_{10}(E_{\rm k;\,fit})$",fontsize=12)

    plt.tight_layout()
    figname_ = os.getcwd()+f'/figs/{figname}_{o_fit.name}'
    plt.savefig(figname_+'.png',dpi=256)
    plt.savefig(figname_+'.pdf')
    plt.show()


    if (True & plot_fit_coeffs):
        fig,axes = plt.subplots(ncols=len(o_fit.keys),nrows=1,figsize=(12,3),layout='constrained',sharex='all')
        for i, key in enumerate(o_fit.keys):
            for j, (name, sim_dict) in enumerate(df.iterrows()):
                # df_i/ = result_1seg[name]
                # axes[i].axhline(y=df_i.iloc[0][key],color=sim_dict['color'],linestyle=sim_dict['ls'],linewidth=0.6)

                # df_i = result_9seg[name]
                if key == 'y0':
                    val = np.log10(infos.loc[name][key])
                    # lbl = r"$\log_{10}($"+sim_dict["label"]+"$)$"
                else:
                    val = infos.loc[name][key]
                    # lbl = sim_dict["label"]
                axes[i].plot([sim_dict["label"]], val, marker=sim_dict['marker'],color=sim_dict['color'],fillstyle='none')

        # # plot fit to y0 coefficient (get average coefficients for each angular segment (assume that average is good))
        # for (name, sim_dict) in df.iterrows():
        #
        #     angles = np.array( result_9seg[name].index, dtype=np.float64 )
        #     slope = np.float64( result_1seg[name]["slope"] )
        #     intercept = np.float64( result_1seg[name]["intercept"] )
        #     y0_fit = slope * np.sin( angles*np.pi/180 ) + intercept
        #
        #     axes[o_fit.keys.index('y0')].plot(angles, un_log_y( y0_fit ),ls=sim_dict['ls'],color=sim_dict['color'],fillstyle='none')

        for i, key in enumerate(o_fit.keys):
            axes[i].set_title(o_fit.labels[i],fontsize=12)

        # for ax in axes:
        # axes[o_fit.keys.index('y0')].set_yscale("log")
        # axes[o_fit.keys.index('y0')].set_ylabel("log")
        axes[o_fit.keys.index('y0')].set_title(r"$\log_{10}($"+o_fit.labels[o_fit.keys.index('y0')]+"$)$",fontsize=12)

        for ax in axes:
            # for ax in ax:
            # ax.set_xlim(.1,90.1)
            ax.tick_params(labelsize=12,which='both',direction='in',tick1On=True, tick2On=True)
            ax.tick_params(axis='x', which='minor', bottom=False)
            # ax.minorticks_on()
            ax.grid(lw=0.6)
            # ax.set_xlabel("Polar angle [deg]",fontsize=12)
            for tick in ax.get_xticklabels():
                tick.set_rotation(75)

        figname_ = os.getcwd()+f'/figs/{figname_coeffs}_{o_fit.name}'
        plt.savefig(figname_+'.png',dpi=256)
        plt.savefig(figname_+'.pdf')

        plt.show()

def plot_all_sim_ejecta_mass_row(
        xlim=(1e-3, 4),ylim0=(1e43, 1e51),ylim1=(1e43, 1e51),ylim2=(-.5,.5),
       figname="figname",figname_coeffs="figname", sim: str or None=None, plot_fit_coeffs:bool=True):
    o_fits = [Fit1D_2seg(),Fit1D_3seg(),Fit1D_4seg()]
    fit_colors=['blue','green','red']
    fit_lcs = [':','-.','--']
    do_cumulative = True
    log_type = 10
    log_type_y = 10

    get_cumulative = lambda val : np.cumsum(val[::-1])[::-1] if do_cumulative else val
    do_log = lambda val : (val if not log_type else (np.log2(val) if log_type==2 else np.log10(val)))
    un_log = lambda val : (val if not log_type else (2**(val) if log_type==2 else 10**(val)))

    do_log_y = lambda val : (val if not log_type_y else (np.log2(val) if log_type_y==2 else np.log10(val)))
    un_log_y = lambda val : (val if not log_type_y else (2**(val) if log_type_y==2 else 10**(val)))

    fig, axes = plt.subplots(ncols=len(df),nrows=3,figsize=(12.,5.),
                             layout='constrained',sharex='col',sharey='row',#sharey='row',
                             gridspec_kw={'height_ratios': [2,2,1]})

    # df_fit = pd.read_csv(os.getcwd()+'/'+'piecewise_line_fits.csv',index_col=0)

    infos_fits = dict()
    infos_fits2 = dict()
    for o_fit, color, lc in zip(o_fits,fit_colors,fit_lcs):
        infos_fits[o_fit.name] = dict()
    for i_s, (name, sim_dic) in enumerate(df.iterrows()):
        infos_fits2[name] = dict(label=sim_dic['label'])
    # for o_fit, color, lc in zip(o_fits,fit_colors,fit_lcs):
    #     for q in ["chi2","r2","sse"]:
    #         infos_fits[f"{o_fit.name} {q}"] = 0
    idx = 0

    # x1, x2, y1, y2 = -1.5, -0.9, -2.5, -1.9  # subregion of the original image
    # axins = axes[1].inset_axes(
    #     [0.05, 0.1, 0.5, 0.5],
    #     xlim=(x1, x2), ylim=(y1, y2), xticklabels=[], yticklabels=[])

    # sub region of the original image
    # x1, x2, y1, y2 = 8e-2, 5e-1, 2e48, 4e49
    # axins.set_xlim(x1, x2)
    # axins.set_ylim(y1, y2)
    # axins.set_xscale('log')
    # axins.set_yscale('log')
    # axins.grid(ls=':')
    # axins.set_xticklabels([])
    # axins.set_yticklabels([])
    # axins.set_xlabel('')
    # axins.xaxis.set_visible(False)
    # axins.set_xticks([])
    # axins.minorticks_on()
    # axins.tick_params(axis='both', which='both', labelleft=False,
    #                   labelright=True, tick1On=True, tick2On=True,
    #                   labelsize=12,
    #                   direction='in',
    #                   bottom=True, top=True, left=True, right=True)
    # axes[1].indicate_inset_zoom(axins, edgecolor="black")

    # for idx, sim_dic in enumerate(df.iterrows()):
    # for (sim_dic,fit_dic) in zip(df.iterrows(), df_fit.iterrows()):
    # o_fit = Fit1D_4seg()
    # o_fit = Fit1D_3seg()

    for i_s, (name, sim_dic) in enumerate(df.iterrows()):
        text_dict = df_text.loc[name]
        # sim_dic = sim_dic[1]
        if sim and name != sim:
            continue

        # name = sim_dic[0]
        ej = PBA.id_kenta.EjectaData(get_ej_data(sim_dic['name']),verbose=True)
        mom, _, ek, mass = get_2d_ek(ej=ej,
                                     # text=get_text_max_mass(ej=ej,crit='fast')
                                     text=float(sim_dic['tmerg'])+float(text_dict['text']))
        # ek = np.sum(mass,axis=1)
        ek = np.sum(ek,axis=1)
        # ek = np.sum(mass,axis=1)*PBA.cgs.solar_m#*(PBA.cgs.c**2*PBA.BetaFromMom(mom)**2)

        # ej_mass = ej.get(v_n="mass",text=get_text_max_mass(ej=ej,crit='fast'))
        # ek = np.cumsum(np.sum(ek,axis=0)[::-1] * PBA.cgs.solar_m * vinf_c**2 * PBA.cgs.c**2)[::-1]/np.sum(ek * PBA.cgs.solar_m * vinf_c**2 * PBA.cgs.c**2)
        # ek = np.cumsum(ek[::-1])[::-1]#/np.sum(ek)


        l_mom, l_ek = do_log(mom), do_log_y(ek)
        mask = np.isfinite(l_ek)
        l_mom = l_mom[mask]
        l_ek = l_ek[mask]

        # N = 5
        # l_ek = np.convolve(l_ek, np.ones(N)/N, mode='valid')
        # l_mom = np.convolve(l_mom, np.ones(N)/N, mode='valid')


        axes[0][i_s].plot(un_log(l_mom), get_cumulative(un_log_y(l_ek)), color='black', ls='-',lw=1)#, lw=0.7, drawstyle='steps')
        axes[1][i_s].plot(un_log(l_mom), un_log_y(l_ek), color='black', ls='-',lw=1)#, lw=0.7, drawstyle='steps')

        for o_fit, color, ls in zip(o_fits,fit_colors,fit_lcs):
            l_mom_pred, l_ek_pred, info = fit_data(l_mom, l_ek, un_log, un_log_y,
                                                   fitting_obj=o_fit, name=sim_dic["name"])

            info["label"] = sim_dic["label"]
            infos_fits[o_fit.name][name] = info
            axes[0][i_s].plot(un_log(l_mom_pred), get_cumulative(un_log_y(l_ek_pred)), color=color, ls=ls,lw=1)#, lw=0.7, drawstyle='steps')
            axes[1][i_s].plot(un_log(l_mom_pred), un_log_y(l_ek_pred), color=color, ls=ls,lw=1)#, lw=0.7, drawstyle='steps')
            axes[2][i_s].plot(un_log(l_mom_pred), (l_ek-l_ek_pred), color=color, ls=ls, lw=1)


            for q in ["chi2","r2","sse"]:
                infos_fits2[name][f"{o_fit.name} {q}"] = info[q]

        # axins.plot(un_log(l_mom), un_log_y(l_ek), color=sim_dic["color"], ls=sim_dic["ls"], label=sim_dic["label"],lw=1.2)
        # axins.plot(un_log(l_mom_pred), un_log_y(l_ek_pred), color=sim_dic["color"], ls=sim_dic["ls"],lw=.6)#, lw=0.7, drawstyle='steps')

        '''
        fit_dic = df_fit.loc[name]
        y_pred = piecewise_power(10**l_mom,x0=fit_dic['x0'],x1=fit_dic['x1'],y0=fit_dic['y0'],
                                 k1=fit_dic['k1'],k2=fit_dic['k2'],k3=fit_dic['k3'])
        # y_pred = 10**piecewise_linear(l_mom,x0=fit_dic['x0'],x1=fit_dic['x1'],y0=fit_dic['y0'],
        #                          k1=fit_dic['k1'],k2=fit_dic['k2'],k3=fit_dic['k3'])
        axes[0].plot(10**l_mom, get_cumulative(y_pred), color='gray', ls=sim_dic["ls"],lw=1.2)#, lw=0.7, drawstyle='steps')
        axes[1].plot(10**l_mom, y_pred, color='gray', ls=sim_dic["ls"],lw=1.2)#, lw=0.7, drawstyle='steps')
        '''

        # axes[0].plot(mom, get_cumulative(ek), color=sim_dic["color"], ls=sim_dic["ls"], label=sim_dic["label"],lw=1.2)#, lw=0.7, drawstyle='steps')

        # mom_pred, ek_pred, info = fit_data(mom, ek,name=sim_dic["name"])
        # infos[sim_dic["name"]] = info
        # axes[0].plot(mom_pred, get_cumulative(ek_pred), color=sim_dic["color"], ls=sim_dic["ls"],lw=.6)#, lw=0.7, drawstyle='steps')
        # axes[1].plot(mom_pred, (ek-ek_pred) / ek, color=sim_dic["color"], ls=sim_dic["ls"], lw=0.7)

        # break


    axes[0][0].plot([1e-4,1e-2], [1e39,1e41], color='black', ls='-', label='NR',lw=1)#, lw=0.7, drawstyle='steps')
    for o_fit, color, ls in zip(o_fits,fit_colors,fit_lcs):
        axes[0][0].plot([1e-4,1e-2], [1e39,1e41], color=color, ls=ls,label=o_fit.name,lw=1)#, lw=0.7, drawstyle='steps')

    for o_fit, color in zip(o_fits,fit_colors):
        print('-'*20+o_fit.name+'-'*20)
        infos = pd.DataFrame.from_dict(infos_fits[o_fit.name]).T
        o_fit.print_table(infos=infos)
        print('-'*21+'-'*25)
    print('-'*20+"Total"+'-'*20)

    df_sum = pd.DataFrame.from_dict(infos_fits2).T
    print(df_sum.to_latex(index=False,float_format="%.3f"))
    print('-'*21+'-'*25)

    for ax in axes:
        for ax in ax:
            ax.tick_params(labelsize=12)
            # ax.set_ylim(*ylim)

            ax.tick_params(labelsize=12,axis='both', which='both',direction='in',tick1On=True, tick2On=True)
            ax.minorticks_on()
            # ax.grid(ls=':',lw=0.8)
            ax.set_xscale('log')

    #for idx in range(len(df)):
    for ax in axes[0]: ax.set_yscale('log')
    for ax in axes[1]: ax.set_yscale('log')
    # axes[2].set_yscale('log')
    for ax in axes[0]: ax.set_ylim(*ylim0)
    for ax in axes[1]: ax.set_ylim(*ylim1)
    for ax in axes[2]: ax.set_ylim(*ylim2)
    for ax in axes[-1]:
        ax.set_xlim(*xlim)
        ax.grid(ls='-',lw=0.6)
    for ax, (name, sim_dic) in zip(axes[0],df.iterrows()):
        ax.set_title(sim_dic['label'],fontsize=12)

    # ax.set_yscale(yscale)
    # ax.set_xscale(xscale)


    n = 2
    ax = axes[0][0]
    ax.legend(**dict(fancybox=False,loc= 'lower left',columnspacing=0.4,
                     #"bbox_to_anchor": (0.5, 1.2),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
                     shadow=False, ncol= 1, fontsize= 12,framealpha=0., borderaxespad= 0., frameon=False))
    # han, lab = ax.get_legend_handles_labels()
    # ax.add_artist(ax.legend(han[:-1 * n], lab[:-1 * n],
    #                         **dict(fancybox=False,loc= 'lower left',columnspacing=0.4,
    #                                #"bbox_to_anchor": (0.5, 1.2),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
    #                                shadow=False, ncol= 1, fontsize= 12,framealpha=0., borderaxespad= 0., frameon=False)))
    # ax.add_artist(ax.legend(han[len(han) - n:], lab[len(lab) - n:],
    #                         **dict(fancybox=False,loc= 'upper right',columnspacing=0.4,
    #                                #"bbox_to_anchor": (0.5, 1.2),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
    #                                shadow=False, ncol= 1, fontsize= 12,framealpha=0., borderaxespad= 0., frameon=False)))
    # ax.set_yscale("log")
    for ax in axes[-1]:
        ax.set_xlabel(r"$\Gamma\beta$",fontsize=12)
    axes[0][0].set_ylabel(r"$E_{\rm k} (> \Gamma \beta)$ [erg]",fontsize=12)
    axes[1][0].set_ylabel(r"$E_{\rm k}$ [erg]",fontsize=12)
    # ax.set_title(title,fontsize=12)
    # axes[-1].set_ylabel(r"$\Delta E_{\rm k}$ [erg]",fontsize=12)
    # axes[-1][0].set_ylabel(r"$\log_{10}(E_{\rm k;\,sph})-\log_{10}(E_{\rm k;\,fit})$",fontsize=12)
    axes[-1][0].set_ylabel(r"$\Delta\log_{10}(E_{\rm k})$",fontsize=12)

    plt.tight_layout()
    figname_ = os.getcwd()+f'/figs/{figname}_{"all_fits"}'
    plt.savefig(figname_+'.png',dpi=256)
    plt.savefig(figname_+'.pdf')
    plt.show()


    if (True & plot_fit_coeffs):
        for o_fit, color in zip(o_fits,fit_colors):


            fig,axes = plt.subplots(ncols=len(o_fit.keys),nrows=1,figsize=(12,3),layout='constrained',sharex='all')
            for i, key in enumerate(o_fit.keys):
                for j, (name, sim_dict) in enumerate(df.iterrows()):
                    # df_i/ = result_1seg[name]
                    # axes[i].axhline(y=df_i.iloc[0][key],color=sim_dict['color'],linestyle=sim_dict['ls'],linewidth=0.6)

                    # df_i = result_9seg[name]
                    if key == 'y0':
                        val = np.log10(infos_fits[o_fit.name][name][key])
                        # lbl = r"$\log_{10}($"+sim_dict["label"]+"$)$"
                    else:
                        val = infos_fits[o_fit.name][name][key]
                        # lbl = sim_dict["label"]
                    axes[i].plot([sim_dict["label"]], val, marker=sim_dict['marker'],color=sim_dict['color'],fillstyle='none')

            # keys to fit
            print("-"*10+o_fit.name+10*'-')
            keys_to_fit_linear = ['x0','x1','k1', 'k2','k3']
            keys_to_fit_loglinear = ['y0']
            keys_to_fit_parabolic = ['k2','k3']
            for i, key in enumerate(o_fit.keys):
                if key in keys_to_fit_linear:
                    q_ = []
                    y_ = []
                    lbls = []
                    for j, (name, sim_dict) in enumerate(df.iterrows()):
                        if name.__contains__("SFHo_"):
                            q_.append(float(sim_dict["q"]))
                            y_.append(float(infos_fits[o_fit.name][name][key]))
                            lbls.append(sim_dict['label'])
                    slope, intercept, r_value, p_value, std_err = linregress(np.array(q_), np.array(y_))

                    q_range = np.array(q_)#np.mgrid[np.min(q_):np.max(q_):10]
                    y_fit = slope * q_range + intercept

                    print(f"\tLinear for {key}: "
                          f"{o_fit.labels[o_fit.keys.index(key)][:-1]} = "
                          f"{slope:.2f} \, q + {intercept:.2f}$")

                    axes[i].plot(lbls, y_fit, color='gray',ls='-',lw=1.0)
                if key in keys_to_fit_loglinear:
                    q_ = []
                    y_ = []
                    lbls = []
                    for j, (name, sim_dict) in enumerate(df.iterrows()):
                        if name.__contains__("SFHo_"):
                            q_.append(float(sim_dict["q"]))
                            y_.append(float(infos_fits[o_fit.name][name][key]))
                            lbls.append(sim_dict['label'])
                    slope, intercept, r_value, p_value, std_err = linregress(np.array(q_), np.array(y_))

                    q_range = np.array(q_)#np.mgrid[np.min(q_):np.max(q_):10]
                    y_fit = slope * q_range + intercept

                    print(f"\tLoglinear for {key}: "
                          f"{o_fit.labels[o_fit.keys.index(key)][:-1]} = "
                          f"{PBA.latex_float(slope)} \, q + {PBA.latex_float(intercept)}$")

                    axes[i].plot(lbls, np.log10(y_fit) , color='gray',ls='-',lw=1.0)
                if key in keys_to_fit_parabolic:
                    points = []
                    lbls = []
                    for j, (name, sim_dict) in enumerate(df.iterrows()):
                        if name.__contains__("SFHo_"):
                            points.append((float(sim_dict["q"]),
                                           float(infos_fits[o_fit.name][name][key])))
                            lbls.append(sim_dict['label'])
                    # print(points)
                    # Create matrices for the linear system
                    A = np.array([[x**2, x, 1] for (x, y) in points])
                    b = np.array([y for x, y in points])

                    # Solve the linear system
                    coefficients = np.linalg.solve(A, b)

                    a, b, c = coefficients
                    print(f"\tParabolic for {key}: "
                          f"{o_fit.labels[o_fit.keys.index(key)][:-1]} = "
                          f"{a:.2f} \, q^2 + {b:.2f} \, q + {c:.2f}$")
                    func = lambda q : a*q**2 + b*q + c
                    y_fit = [func(q) for (q,y) in points]
                    axes[i].plot(lbls, y_fit , color='gray',ls='--',lw=1.0)
            # # plot fit to y0 coefficient (get average coefficients for each angular segment (assume that average is good))
            # for (name, sim_dict) in df.iterrows():
            #
            #     angles = np.array( result_9seg[name].index, dtype=np.float64 )
            #     slope = np.float64( result_1seg[name]["slope"] )
            #     intercept = np.float64( result_1seg[name]["intercept"] )
            #     y0_fit = slope * np.sin( angles*np.pi/180 ) + intercept
            #
            #     axes[o_fit.keys.index('y0')].plot(angles, un_log_y( y0_fit ),ls=sim_dict['ls'],color=sim_dict['color'],fillstyle='none')

            for i, key in enumerate(o_fit.keys):
                axes[i].set_title(o_fit.labels[i],fontsize=12)

            # for ax in axes:
            # axes[o_fit.keys.index('y0')].set_yscale("log")
            # axes[o_fit.keys.index('y0')].set_ylabel("log")
            axes[o_fit.keys.index('y0')].set_title(r"$\log_{10}($"+o_fit.labels[o_fit.keys.index('y0')]+"$)$",fontsize=12)

            for ax in axes:
                # for ax in ax:
                # ax.set_xlim(.1,90.1)
                ax.tick_params(labelsize=12,which='both',direction='in',tick1On=True, tick2On=True)
                ax.tick_params(axis='x', which='minor', bottom=False)
                # ax.minorticks_on()
                ax.grid(lw=0.6)
                # ax.set_xlabel("Polar angle [deg]",fontsize=12)
                for tick in ax.get_xticklabels():
                    tick.set_rotation(75)

            figname_ = os.getcwd()+f'/figs/{figname_coeffs}_{o_fit.name}'
            plt.savefig(figname_+'.png',dpi=256)
            plt.savefig(figname_+'.pdf')

            plt.show()

def plot_ej_mom_theta_dist_for_text():
    fig, axes = plt.subplots(ncols=len(df),nrows=1,figsize=(3*4.6,3.2),layout='constrained',sharey='all')

    i = 0
    for sim, sim_dict in df.iterrows():
        text_dict = df_text.loc[sim]
        ej = PBA.id_kenta.EjectaData(get_ej_data(sim_dict['name']),verbose=True)
        mom,ctheta,ek = get_2d_ek(ej=ej,text=sim_dict['tmerg']+text_dict['text'])

        cmap = plt.get_cmap('viridis')
        norm = LogNorm(vmin=ek.max()*1e-7, vmax=ek.max())

        ax = axes[i]
        im = ax.pcolormesh(ctheta*180./np.pi,mom, ek, cmap=cmap, norm=norm)
        # ax0.set_yscale('log')
        # fig.colorbar(im, ax=ax)
        ax.set_yscale('log')
        ax.set_title(sim_dict['label'])
        # ax.set_ylim()
        i+=1
    axes[0].set_ylabel(r"$\Gamma\beta$",fontsize=12)
    for ax in axes:
        ax.tick_params(labelsize=12,axis='both', which='both',direction='in',tick1On=True, tick2On=True)
        ax.minorticks_on()
        ax.set_xlabel(r"Polar Angle [deg]",fontsize=12)
    plt.show()

def plot_sim_ekecta_mass_OLD(sim_:str or None):

    # fig, axes = plt.subplots(nrows=2,ncols=len(df) if not sim_ else 1,sharex='col',sharey='row',
    #                          layout='constrained',figsize=(12,4.5))

    do_cumulative = True
    log_type = 10
    log_type_y = 10

    get_cumulative = lambda val : np.cumsum(val[::-1])[::-1] if do_cumulative else val
    do_log = lambda val : (val if not log_type else (np.log2(val) if log_type==2 else np.log10(val)))
    un_log = lambda val : (val if not log_type else (2**(val) if log_type==2 else 10**(val)))

    do_log_y = lambda val : (val if not log_type_y else (np.log2(val) if log_type_y==2 else np.log10(val)))
    un_log_y = lambda val : (val if not log_type_y else (2**(val) if log_type_y==2 else 10**(val)))
    o_fit = Fit1D_3seg()
    keys = ["x0","x1","y0","k1","k2","k3"]

    segs = [(-.1,91)]
    csegs = [0.5*(seg[0]+seg[1]) for seg in segs]
    # colors = ["blue","lime","green","orange","red"]
    result_1seg = dict()
    for (name, sim_dict) in df.iterrows():
        text_dict = df_text.loc[name]
        ej = PBA.id_kenta.EjectaData(get_ej_data(sim_dict['name']),verbose=False)
        mom, ctheta, ek, mass = get_2d_ek(ej=ej,
                                          # text=get_text_max_mass(ej=ej,crit='fast')
                                          text=float(sim_dict['tmerg'])+float(text_dict['text']))

        result_ = dict()

        for segment in segs:
            if sim_ and name != sim_:
                continue

            c_seg = int( (segment[1]+segment[0])*0.5 )

            low = segment[0] * np.pi / 180.
            up  = segment[1] * np.pi / 180.
            ek_ = np.sum(ek[:,(ctheta > low) & (ctheta <= up)],axis=1)

            l_mom, l_ek = do_log(mom), do_log_y(ek_)
            mask = np.isfinite(l_ek)
            l_mom = l_mom[mask]
            l_ek = l_ek[mask]

            x_fit, y_fit, fit_dict = fit_data(x=l_mom, y=l_ek, undo_x=un_log, undo_y=un_log_y, name=name,
                                              fitting_obj=o_fit,save_res=False)

            result_[c_seg] = fit_dict
        df_i = pd.DataFrame.from_dict(result_).T
        result_1seg[name] = df_i

    segs = [(-.1,30.),(30,60.),(60.,91.)]
    csegs = [0.5*(seg[0]+seg[1]) for seg in segs]
    # colors = ["blue","lime","green","orange","red"]
    result_3seg = dict()
    for (name, sim_dict) in df.iterrows():
        text_dict = df_text.loc[name]
        ej = PBA.id_kenta.EjectaData(get_ej_data(sim_dict['name']),verbose=False)
        mom, ctheta, ek, mass = get_2d_ek(ej=ej,
                                          # text=get_text_max_mass(ej=ej,crit='fast')
                                          text=float(sim_dict['tmerg'])+float(text_dict['text']))

        result_ = dict()

        for segment in segs:
            if sim_ and name != sim_:
                continue

            c_seg = int( (segment[1]+segment[0])*0.5 )

            low = segment[0] * np.pi / 180.
            up  = segment[1] * np.pi / 180.
            ek_ = np.sum(ek[:,(ctheta > low) & (ctheta <= up)],axis=1)

            l_mom, l_ek = do_log(mom), do_log_y(ek_)
            mask = np.isfinite(l_ek)
            l_mom = l_mom[mask]
            l_ek = l_ek[mask]

            x_fit, y_fit, fit_dict = fit_data(x=l_mom, y=l_ek, undo_x=un_log, undo_y=un_log_y, name=name,
                                              fitting_obj=o_fit,save_res=False)

            result_[c_seg] = fit_dict
        df_i = pd.DataFrame.from_dict(result_).T
        result_3seg[name] = df_i

    segs = [(-.1,18.),(18,36.),(36.,54.),(54.,72,),(72.,91.)]
    csegs = [0.5*(seg[0]+seg[1]) for seg in segs]
    # colors = ["blue","lime","green","orange","red"]
    result_5seg = dict()
    for (name, sim_dict) in df.iterrows():
        text_dict = df_text.loc[name]
        ej = PBA.id_kenta.EjectaData(get_ej_data(sim_dict['name']),verbose=False)
        mom, ctheta, ek, mass = get_2d_ek(ej=ej,
                                          # text=get_text_max_mass(ej=ej,crit='fast')
                                          text=float(sim_dict['tmerg'])+float(text_dict['text']))

        result_ = dict()

        for segment in segs:
            if sim_ and name != sim_:
                continue

            c_seg = int( (segment[1]+segment[0])*0.5 )

            low = segment[0] * np.pi / 180.
            up  = segment[1] * np.pi / 180.
            ek_ = np.sum(ek[:,(ctheta > low) & (ctheta <= up)],axis=1)

            l_mom, l_ek = do_log(mom), do_log_y(ek_)
            mask = np.isfinite(l_ek)
            l_mom = l_mom[mask]
            l_ek = l_ek[mask]

            x_fit, y_fit, fit_dict = fit_data(x=l_mom, y=l_ek, undo_x=un_log, undo_y=un_log_y, name=name,
                                              fitting_obj=o_fit,save_res=False)

            result_[c_seg] = fit_dict
        df_i = pd.DataFrame.from_dict(result_).T
        result_5seg[name] = df_i


    segs = [(-.1,10.),(10,20.),(20.,30.),(30.,40,),(40.,50.),(50.,60.),(60.,70.),(70.,80.),(80.,91.)]
    csegs = [0.5*(seg[0]+seg[1]) for seg in segs]
    # colors = ["blue","lime","green","orange","red"]
    result_9seg = dict()
    for (name, sim_dict) in df.iterrows():
        text_dict = df_text.loc[name]
        ej = PBA.id_kenta.EjectaData(get_ej_data(sim_dict['name']),verbose=False)
        mom, ctheta, ek, mass = get_2d_ek(ej=ej,
                                          # text=get_text_max_mass(ej=ej,crit='fast')
                                          text=float(sim_dict['tmerg'])+float(text_dict['text']))

        result_ = dict()

        for segment in segs:
            if sim_ and name != sim_:
                continue

            c_seg = int( (segment[1]+segment[0])*0.5 )

            low = segment[0] * np.pi / 180.
            up  = segment[1] * np.pi / 180.
            ek_ = np.sum(ek[:,(ctheta > low) & (ctheta <= up)],axis=1)

            l_mom, l_ek = do_log(mom), do_log_y(ek_)
            mask = np.isfinite(l_ek)
            l_mom = l_mom[mask]
            l_ek = l_ek[mask]

            x_fit, y_fit, fit_dict = fit_data(x=l_mom, y=l_ek, undo_x=un_log, undo_y=un_log_y, name=name,
                                              fitting_obj=o_fit,save_res=False)

            result_[c_seg] = fit_dict
        df_i = pd.DataFrame.from_dict(result_).T
        result_9seg[name] = df_i

    fig,axes = plt.subplots(ncols=6,nrows=3,figsize=(12,9),layout='constrained',sharex='all')
    for i, key in enumerate(keys):
        for (name, sim_dict) in df.iterrows():
            df_i = result_1seg[name]
            axes[0][i].axhline(y=df_i.iloc[0][key],color=sim_dict['color'],linestyle=sim_dict['ls'],linewidth=0.6)
            axes[1][i].axhline(y=df_i.iloc[0][key],color=sim_dict['color'],linestyle=sim_dict['ls'],linewidth=0.6)
            axes[2][i].axhline(y=df_i.iloc[0][key],color=sim_dict['color'],linestyle=sim_dict['ls'],linewidth=0.6)
            df_i = result_3seg[name]
            for (cth, fit_dict_cth) in df_i.iterrows():
                axes[0][i].plot(cth, fit_dict_cth[key],marker=sim_dict['marker'],color=sim_dict['color'],fillstyle='none')
            df_i = result_5seg[name]
            for (cth, fit_dict_cth) in df_i.iterrows():
                axes[1][i].plot(cth, fit_dict_cth[key],marker=sim_dict['marker'],color=sim_dict['color'],fillstyle='none')
            df_i = result_9seg[name]
            for (cth, fit_dict_cth) in df_i.iterrows():
                axes[2][i].plot(cth, fit_dict_cth[key],marker=sim_dict['marker'],color=sim_dict['color'],fillstyle='none')

    for ax in axes:
        ax[2].set_yscale("log")
    # for ax in axes:
    #     for ax in ax[3]:
    #         # ax.set_xscale("log")
    #         ax.set_yscale("log")
    for ax in axes:
        for ax in ax:
            ax.set_xlim(0,90)
    plt.show()






















    fit_data_segs = dict()
    for cseg in csegs:
        for (name, sim_dic) in df.iterrows():
            fit_data_segs[f"{int(cseg)} {name}"] = dict()


    i = 0
    for (name, sim_dic) in df.iterrows():
        text_dict = df_text.loc[name]
        # sim_dic = sim_dic[1]
        if sim_ and name != sim_:
            continue

        sim_dict = df.loc[name]
        text_dict = df_text.loc[name]

        ej = PBA.id_kenta.EjectaData(get_ej_data(sim_dict['name']),verbose=True)
        mom, ctheta, ek, mass = get_2d_ek(ej=ej,
                                     # text=get_text_max_mass(ej=ej,crit='fast')
                                     text=float(sim_dict['tmerg'])+float(text_dict['text'])
                                     )
        # ek = np.sum(mass,axis=1)
        # ek = np.sum(ek,axis=1)

        low = 30. * np.pi / 180.
        mid = 60. * np.pi / 180.

        ek_tot = np.sum(ek, axis=1)
        ek_equatorial = np.sum(ek[:,(ctheta >= mid)],axis=1) #np.sum(ek[:,(ctheta >= mid)], axis=1)
        ek_mid = np.sum(ek[:,(ctheta > low) & (ctheta <= mid)], axis=1)
        ek_polar = np.sum(ek[:,(ctheta <= low)], axis=1)


        ax = axes[0][i]
        ax.plot(mom,ek_equatorial, color='red')
        ax.plot(mom,ek_mid, color='green')
        ax.plot(mom,ek_polar, color='blue')

        ax = axes[1][i]
        # data = np.insert(np.log10(mom), 0, -1)
        # for segment, color in zip([(0.,15.),(15,30.),(30.,45.),(45.,60,),(60.,75.),(75.,90.)],
        #                           ["blue","cyan","lime","green","orange","red"]):


        for segment, color in zip([(0.,18.),(18,36.),(36.,54.),(54.,72,),(72.,90.)],
                                  ["blue","lime","green","orange","red"]):
            cseg = int( (segment[1]+segment[0])*0.5 )
            low = segment[0] * np.pi / 180.
            up  = segment[1] * np.pi / 180.
            ek_ = np.sum(ek[:,(ctheta > low) & (ctheta <= up)],axis=1)

            l_mom, l_ek = do_log(mom), do_log_y(ek_)
            mask = np.isfinite(l_ek)
            l_mom = l_mom[mask]
            l_ek = l_ek[mask]

            x_fit, y_fit, fit_dict = fit_data(x=l_mom, y=l_ek, undo_x=un_log, undo_y=un_log_y, name=name,
                                              fitting_obj=o_fit,save_res=False)

            ax.plot(mom,ek_,color=color)
            ax.plot(un_log(x_fit),un_log_y(y_fit),color=color,ls=':')

            fit_data_segs[f"{int(cseg)} {name}"] = fit_dict

        #     data = np.column_stack((data, np.insert(np.log10(ek_),0, 0.5*(up+low) )))
        #
        # np.savetxt(os.getcwd()+f'/output/+name+"_2D_6seg_log_ek_log_mom.txt",X=data,fmt="%.3f")

        # low = 30. * np.pi / 180.
        # mid = 60. * np.pi / 180.
        #
        # ek_tot = np.sum(ek, axis=1)
        # ek_equatorial = np.sum(ek[:,(ctheta >= mid)],axis=1) #np.sum(ek[:,(ctheta >= mid)], axis=1)
        # ek_mid = np.sum(ek[:,(ctheta > low) & (ctheta <= mid)], axis=1)
        # ek_polar = np.sum(ek[:,(ctheta <= low)], axis=1)
        i+=1
    for ax in axes:
        for ax in ax:
            ax.set_xscale("log")
            ax.set_yscale("log")
    plt.show()

    print(f"ctheta:{ctheta.shape} mom:{mom.shape} ek:{ek.shape}")


    # plt.figure(figsize=(6, 6))
    # plt.title("Smoothed Raster Image")
    # plt.imshow(image, cmap='viridis')
    # plt.colorbar()
    # plt.show()
    #
    # fig,ax = plt.subplots(ncols=1,nrows=1)
    # ax.set_title("Smoothed Raster Image")
    # im = ax.imshow(image, cmap='viridis')
    # ax.set_rasterized()
    # fig.colorbar(im, ax=ax)
    # plt.savefig(outputfile,dpi=600)
    # plt.show()



    ''' -------------------------------------------------------------------------- '''

    default_pysr_params = dict(
        procs=4,
        populations=8,
        population_size=50,
        ncycles_per_iteration=500,
        early_stop_condition=(
            "stop_if(loss, complexity) = loss < 1e-6 && complexity < 10"
            # Stop early if we find a good and simple equation
        ),
        model_selection="best",
        weight_randomize=0.1, # ^ Randomize the tree much more frequently
    )
    model = PySRRegressor(
        niterations=100,
        binary_operators=["+", "*", "-", "/"],
        unary_operators=["cos", "exp", "sin", "square", "cube", "cos2(x)=cos(x)^2", "sin2(x)=sin(x)^2"],
        extra_sympy_mappings={"cos2": lambda x: sympy.cos(x)**2,
                              "sin2": lambda x: sympy.sin(x)**2},
        **default_pysr_params,
    )

    log_ek = np.log10(ek_tot)
    mom_ = mom[np.isfinite(log_ek)]
    log_ek = log_ek[np.isfinite(log_ek)]
    log_ek = log_ek[mom_>0.1]
    mom_ = mom_[mom_>0.1]
    model.fit(mom_.reshape(-1, 1),log_ek.reshape(-1, 1))#(X, y)
    print("Best:", model.sympy())
    print(model)
    # for i in range(4):
    #     print(f"Eq{i}: ", model.sympy(i))
    ypredict = model.predict(mom_.reshape(-1, 1))
    print("Default selection MSE:", np.power(ypredict - log_ek, 2).mean())
    ypredict = 10**ypredict.flatten()

    # exit(1)

    fig, axes = plt.subplots(nrows=2,ncols=1,sharex='all',layout='constrained')

    ax = axes[0]
    do_cumulative = True
    get_cumulative = lambda val : np.cumsum(val[::-1])[::-1] if do_cumulative else val
    ax.plot(mom,get_cumulative(ek_tot),color='black',label=r'Total')
    ax.plot(mom_,get_cumulative(ypredict),color='gray')
    ax.plot(mom,get_cumulative(ek_equatorial),label=r'$\theta > 60$')
    ax.plot(mom,get_cumulative(ek_mid),label=r'$30 < \theta < 60$')
    ax.plot(mom,get_cumulative(ek_polar),label=r'$\theta < 30$')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(1e-2,4)
    # ax.set_ylim(1e21,1e31)
    ax.legend(fancybox=False,loc= 'lower left',columnspacing=0.4,
              #"bbox_to_anchor": (0.5, 1.2),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
              shadow=False, ncol= 1, fontsize= 12,
              framealpha=0., borderaxespad= 0., frameon=False)
    ax.set_ylabel(r"$E_k$")

    ax = axes[1]
    ax.pcolormesh(mom, ctheta*180/np.pi, ek.T, cmap=plt.get_cmap('jet'),
                  norm=LogNorm(vmin=ek.max()*1e-7, vmax=ek.max()))
    ax.set_xscale('log')
    ax.set_xlabel(r"$\Gamma\beta$")
    ax.set_ylabel(r"Polar angle [deg]")

    plt.show()
def plot_sim_ekecta_mass(
        o_fit,sim_:str or None, plot_fit_coeffs:bool, plot_fit_ek:bool,plot_box_plots:bool,
        save_fit_result:bool,figname:str):

    # fig, axes = plt.subplots(nrows=2,ncols=len(df) if not sim_ else 1,sharex='col',sharey='row',
    #                          layout='constrained',figsize=(12,4.5))

    do_cumulative = True
    log_type = 10
    log_type_y = 10

    get_cumulative = lambda val : np.cumsum(val[::-1])[::-1] if do_cumulative else val
    do_log = lambda val : (val if not log_type else (np.log2(val) if log_type==2 else np.log10(val)))
    un_log = lambda val : (val if not log_type else (2**(val) if log_type==2 else 10**(val)))

    do_log_y = lambda val : (val if not log_type_y else (np.log2(val) if log_type_y==2 else np.log10(val)))
    un_log_y = lambda val : (val if not log_type_y else (2**(val) if log_type_y==2 else 10**(val)))



    ''' ----------- DEFINE ANGULAR SEGMENTS TO ANALYZE -------------  '''

    # segs = [(-.1,10.),(10,20.),(20.,30.),(30.,40,),(40.,50.),(50.,60.),(60.,70.),(70.,80.),(80.,91.)]
    # segs = np.arange(start=0,stop=90+18,step=18)
    # segs = [(low,up) for (low,up) in zip(segs[:-1],segs[1:])]
    segs = [(0, 36.87),
            (36.87, 53.13),
            (53.13, 66.42),
            (66.42, 78.46),
            (78.46, 90.00)]
    csegs = [0.5*(seg[0]+seg[1]) for seg in segs]

    n_v = 80
    vinf_fit = np.linspace(0.005,1, n_v+1)
    mom_fit = PBA.MomFromBeta(vinf_fit)[:-1]


    result_1seg = dict()
    result_9seg = dict()
    result_9seg_xy = dict()


    # colors = ["blue",'cyan',"lime","green","orange","red"]
    colors = ["blue","green","magenta","orange","red"]
    # cmap = plt.get_cmap('tab10')
    # norm = Normalize(vmin=0,vmax=6)
    # colors = [cmap(norm(i)) for i in range(len(csegs))]

    ''' ----------- COLLECT DATA FOR EACH SIMULATION FOR EACH SEGMENT -------------  '''

    for (name, sim_dict) in df.iterrows():
        if sim_ and name != sim_:
            continue

        text_dict = df_text.loc[name]
        ej = PBA.id_kenta.EjectaData(get_ej_data(sim_dict['name']),verbose=False)
        mom, ctheta, ek, mass = get_2d_ek(ej=ej,
                                          # text=get_text_max_mass(ej=ej,crit='fast')
                                          text=float(sim_dict['tmerg'])+float(text_dict['text']),
                                          new_vinf_len=n_v+1)


        # --------------------- N Segments (independent fit of each segment)  ------------------

        result_ = dict()
        for segment in segs:

            c_seg = float( (segment[1] + segment[0])*0.5 )

            low = segment[0] * np.pi / 180.
            up  = segment[1] * np.pi / 180.
            ek_ = np.sum(ek[:,(ctheta > low) & (ctheta <= up)],axis=1)

            l_mom, l_ek = do_log(mom), do_log_y(ek_)
            mask = np.isfinite(l_ek)
            l_mom = l_mom[mask]
            l_ek = l_ek[mask]

            x_fit, y_fit, fit_dict = fit_data(x=l_mom, y=l_ek, undo_x=un_log, undo_y=un_log_y, name=name,
                                              fitting_obj=o_fit,save_res=False)

            ek_fit = o_fit.piecewise_power(mom_fit, *[np.float64(fit_dict[key]) for key in o_fit.keys])


            result_[f"{c_seg:.3f}"] = fit_dict
            result_9seg_xy[f"{name} {c_seg:.3f}"] = dict(
                nr=np.column_stack((l_mom, l_ek)),
                fit_nr=np.column_stack((x_fit, y_fit)),
                fit=np.column_stack((do_log(mom_fit), do_log_y(ek_fit))),
            )


        df_i = pd.DataFrame.from_dict(result_).T
        result_9seg[name] = df_i

        # -------------------- ONE SEGMENT (total quantities)  -----------------

        result_ = dict()
        for segment in [(-.1, 91)]:

            c_seg = float( (segment[1]+segment[0])*0.5 )

            low = segment[0] * np.pi / 180.
            up  = segment[1] * np.pi / 180.
            ek_ = np.sum(ek[:,(ctheta > low) & (ctheta <= up)],axis=1)

            l_mom, l_ek = do_log(mom), do_log_y(ek_)
            mask = np.isfinite(l_ek)
            l_mom = l_mom[mask]
            l_ek = l_ek[mask]

            x_fit, y_fit, fit_dict = fit_data(x=l_mom, y=l_ek, undo_x=un_log, undo_y=un_log_y, name=name,
                                              fitting_obj=o_fit,save_res=False)

            # -------- fit y0 from N SEGMENT fit ------
            angles = np.array( result_9seg[name].index, dtype=np.float64 )
            y0_s = do_log_y( np.array( result_9seg[name]["y0"] ) )
            # slope, intercept, r_value, p_value, std_err = linregress(np.sin(angles*np.pi/180), y0_s)
            slope, intercept, r_value, p_value, std_err = linregress(angles*np.pi/180, y0_s)
            fit_dict["slope"] = slope
            fit_dict["intercept"] = intercept

            result_[f"{c_seg:.3f}"] = fit_dict

        df_i = pd.DataFrame.from_dict(result_).T
        result_1seg[name] = df_i

        # -------------------- ADD NEW FIT VALUES TO THE N SEG FIT  -----------------
        slope = np.float64( result_1seg[name]["slope"] )
        intercept = np.float64( result_1seg[name]["intercept"] )
        coeffs_tot = [np.average(np.array(result_9seg[name][key])) for key in o_fit.keys]
        for segment in segs:
            c_seg = float( (segment[1] + segment[0])*0.5 )
            # y0_fit = slope * np.sin(c_seg*np.pi/180) + intercept
            y0_fit = slope * c_seg*np.pi/180 + intercept
            y0_fit = un_log_y(y0_fit)
            coeffs_tot[o_fit.keys.index("y0")] = y0_fit
            ek_fit = o_fit.piecewise_power(mom_fit, *np.array(coeffs_tot).flatten().tolist())

            result_9seg_xy[f"{name} {c_seg:.3f}"]["fit_u"] = np.column_stack((
                do_log(mom_fit), do_log_y(ek_fit) # fit using general func and general grids
            ))

        # -------------------- Save Result -----------------
        if save_fit_result:

            # fig, ax = plt.subplots(ncols=1,nrows=1,layout='constrained')
            # for i_seg, segment in enumerate(segs):
            #     c_seg = int( (segment[1]+segment[0])*0.5 )
            #     tbl = result_9seg_xy[f"{name} {c_seg}"]["nr"]
            #     ax.plot(un_log(tbl[:,0]),un_log_y(tbl[:,1]),color=colors[i_seg],ls='-')

            table, table_fit, table_ufit = [], [], []
            # l_mom_fit = np.zeros(0,)
            # vinf_ = np.linspace(0.005,1, len(l_ek)+1)
            # mom = PBA.MomFromBeta(vinf_)[:-1]
            angles = []

            for i_seg, segment in enumerate(segs):
                c_seg = float( (segment[1]+segment[0])*0.5 )
                angles.append(c_seg*np.pi/180) # radians
                tbl_ = result_9seg_xy[f"{name} {c_seg:.3f}"]["nr"]
                table.append(un_log_y(tbl_[:,1]))

            # save result whn eah segment was fitted separately
            for i_seg, segment in enumerate(segs):
                c_seg = float( (segment[1]+segment[0])*0.5 )

                # coeffs = [np.float64(result_9seg[name][key][c_seg]) for key in keys]
                # ek = o_fit.piecewise_power(mom, *np.array(coeffs).tolist())
                # table_fit.append(ek)

                tbl_ = result_9seg_xy[f"{name} {c_seg:.3f}"]["fit"]
                table_fit.append(un_log_y(tbl_[:,1]))

                # ax.plot(un_log(tbl_[:,0]),un_log_y(tbl_[:,1]),color=colors[i_seg],ls='--')

            # save result when all segments are fitted with the same (averaged) pars except y0
            # slope = np.float64( result_1seg[name]["slope"] )
            # intercept = np.float64( result_1seg[name]["intercept"] )
            # coeffs_tot = [np.average(np.array(result_9seg[name][key])) for key in keys]
            for i_seg, segment in enumerate(segs):
                c_seg = float( (segment[1]+segment[0])*0.5 )
                # y0_fit = slope * np.sin(c_seg*np.pi/180) + intercept
                # y0_fit = 10**y0_fit
                # coeffs_tot[keys.index("y0")] = y0_fit
                # ek_fit = o_fit.piecewise_power(mom, *np.array(coeffs_tot).flatten().tolist())
                # table_ufit.append(ek_fit)

                tbl_ = result_9seg_xy[f"{name} {c_seg:.3f}"]["fit_u"]
                table_ufit.append(un_log_y(tbl_[:,1]))

                # ax.plot(un_log(tbl_[:,0]),un_log_y(tbl_[:,1]),color=colors[i_seg],ls=':')

            # table = np.reshape(np.array(table),newshape=(len(table),len(table[0]))) # [ith, imom]
            table_fit = np.reshape(np.array(table_fit),newshape=(len(table_fit),len(table_fit[0]))) # [ith, imom]
            table_ufit = np.reshape(np.array(table_ufit),newshape=(len(table_ufit),len(table_ufit[0]))) # [ith, imom]
            # mom_fit = un_log(l_mom_fit)
            angles = np.array(angles)

            # ax.set(xscale='log',yscale='log',xlim=(1e-2,3),ylim=(1e44,1e49))
            # plt.show()


            fname = os.getcwd()+f'/output/{name}_log_mom_log_ek_asym_and_fit_{o_fit.name}.h5'
            print(f"Saving: {fname}")
            with h5py.File(fname,'w') as f:
                f.create_dataset(name="mom",data=mom_fit) # 4-momentum (dimensionless)
                f.create_dataset(name="ctheta",data=angles) # polar angles segments centers [rad]
                f.create_dataset(name="ek_fit",data=table_fit) # kinetic energy per momentum bin and per angular bin [erg]
                f.create_dataset(name="ek_ufit",data=table_ufit) # kinetic energy per momentum bin and per angular bin [erg] (universal fit)

    if (save_fit_result):
        for sim_name, series in result_1seg.items():
            # df_fit_1 = pd.DataFrame(series)
            series['label'] = df.loc[sim_name]['label']
            series = series[['label','slope','intercept']]
            print(series.to_latex(index=False,float_format="%.2f"))



    # plot coefficients for each simulation and each angular segment
    if (True & plot_fit_coeffs):

        lims = {"x1":[0.5,6.5],"k3":[0,-80]}

        fig,axes = plt.subplots(ncols=len(o_fit.keys),nrows=1,figsize=(12,4.0),
                                # layout='constrained',
                                sharex='all')
        for i, key in enumerate(o_fit.keys):
            for (name, sim_dict) in df.iterrows():
                df_i = result_1seg[name]
                if (key != 'y0'):
                    axes[i].axhline(y=df_i.iloc[0][key],color=sim_dict['color'],linestyle=sim_dict['ls'],linewidth=0.6)

                df_i = result_9seg[name]
                for j, (cth, fit_dict_cth) in enumerate(df_i.iterrows()):
                    label =sim_dict['label'] if (i == 0 and j == 0) else None
                    axes[i].plot(float(cth), fit_dict_cth[key],marker=sim_dict['marker'],color=sim_dict['color'],fillstyle='none',
                                 label=label,ls=sim_dict['ls'])

        # plot fit to y0 coefficient (get average coefficients for each angular segment (assume that average is good))
        for (name, sim_dict) in df.iterrows():

            angles = np.array( result_9seg[name].index, dtype=np.float64 )
            slope = np.float64( result_1seg[name]["slope"] )
            intercept = np.float64( result_1seg[name]["intercept"] )
            # y0_fit = slope * np.sin( angles*np.pi/180 ) + intercept
            y0_fit = slope * angles*np.pi/180 + intercept

            axes[o_fit.keys.index('y0')].plot(
                angles, un_log_y( y0_fit ),ls=sim_dict['ls'],color=sim_dict['color'],fillstyle='none'
            )

        for i, key in enumerate(o_fit.keys):
            axes[i].set_title(o_fit.labels[i],fontsize=12)
            if key in lims.keys():
                axes[i].set_ylim(*lims[key])

        # for ax in axes:
        axes[o_fit.keys.index('y0')].set_yscale("log")
        for ax in axes:
            # for ax in ax:
            ax.set_xlim(.1,90.1)
            ax.tick_params(labelsize=12,which='both',direction='in',tick1On=True, tick2On=True)
            ax.minorticks_on()
            ax.set_xlabel("Polar angle [deg]",fontsize=12)
        # plt.tight_layout()
        pass
        fig.legend(fancybox=False, loc='center', columnspacing=0.8,
                   bbox_to_anchor=(0.5,0.97),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
                   shadow=False, ncol=len(df), fontsize=12, labelcolor='black',
                   framealpha=0.0, borderaxespad=0.)
        plt.subplots_adjust(left=0.03,right=0.99,wspace=0.3)
        figname_ = os.getcwd()+f'/figs/{figname}_{o_fit.name}_coeffs'
        plt.savefig(figname_+'.png',dpi=256)
        plt.savefig(figname_+'.pdf')
        plt.show()

    if (True & plot_fit_ek):
        fig,axes = plt.subplots(ncols=len(df),nrows=1,figsize=(12,4),sharex='col',sharey='row') # layout='constrained',

        for i, (name, sim_dict) in enumerate(df.iterrows()):
            if sim_ and name != sim_:
                continue
            for (segment, color) in zip(segs, colors):
                c_seg = float( (segment[1]+segment[0])*0.5 )
                tbl = result_9seg_xy[f"{name} {c_seg:.3f}"]
                l_mom, l_ek = tbl["nr"][:,0], tbl["nr"][:,1]
                l_mom_fit, l_ek_fit = tbl["fit"][:,0], tbl["fit"][:,1]
                label = r"$\theta_{\rm c}="+f"{c_seg:.1f}$ deg." if i == 0 else None
                axes[i].plot(un_log(l_mom),un_log_y(l_ek),color=color,ls='-',lw=1,label=label)
                axes[i].plot(un_log(l_mom_fit),un_log_y(l_ek_fit),color=color,ls='--',lw=1)

                l_ek_fit = interp1d(l_mom_fit,l_ek_fit,kind='linear')(l_mom)
                # difference = l_ek-l_ek_fit
                # axes[1,i].plot(un_log(l_mom),difference,color=color,ls='-',lw=0.7)
            # pass
            # df_fit_total = result_1seg[name]
            # angles = np.array( result_9seg[name].index, dtype=np.float64 )
            # slope, intercept = y0_fit_dict["slope"], y0_fit_dict["intercept"]
            # slope = np.float64( result_1seg[name]["slope"] )
            # intercept = np.float64( result_1seg[name]["intercept"] )
            # coeffs_tot = [np.average(np.array(result_9seg[name][key])) for key in keys]

            # vinf_ = np.linspace(0.005,1, len(l_ek)+1)
            # mom_fit = PBA.MomFromBeta(vinf_)[:-1]
            # coeffs_tot = np.array(coeffs_tot).flatten().tolist()
            # ek_fit_tot = o_fit.piecewise_power(mom_fit, *coeffs_tot)

            # axes[0,i].plot(mom_fit,ek_fit_tot,color='gray',ls='-',lw=0.7)



            # y0_tot = float( df_fit_total["y0"] )
            # differences = []
            for j, (segment, color) in enumerate( zip(segs, colors) ):
                c_seg = float( (segment[1]+segment[0])*0.5 )
                tbl = result_9seg_xy[f"{name} {c_seg:.3f}"]
                l_mom, l_ek = tbl["nr"][:,0], tbl["nr"][:,1]
                # l_mom_fit, l_ek_fit = tbl["fit"][:,0], tbl["fit"][:,1]
                l_mom_ufit, l_ek_ufit = tbl["fit_u"][:,0], tbl["fit_u"][:,1]


                # y0_fit = slope * np.sin(c_seg*np.pi/180) + intercept
                # y0_fit = 10**y0_fit
                #
                # coeffs_tot[keys.index("y0")] = y0_fit
                # ek_fit = o_fit.piecewise_power(mom_fit, *np.array(coeffs_tot).flatten().tolist())

                axes[i].plot(un_log(l_mom_ufit), un_log_y(l_ek_ufit), color=color,ls=':',lw=1.)

                # tbl = result_9seg_xy[f"{name} {c_seg}"]
                # l_mom, l_ek = tbl["nr"][:,0], tbl["nr"][:,1]
                # l_mom_fit, l_ek_fit = tbl["fit"][:,0], tbl["fit"][:,1]
                # l_ek_fit = interp1d(l_mom_ufit,l_ek_ufit,kind='linear')(l_mom)
                # difference = l_ek-l_ek_fit
                # axes[1,i].plot(un_log(l_mom), difference,color=color,ls=':',lw=0.7)
                # differences.append(difference)
                # mean = np.mean(difference)
            # axes[1,i].axhline(y=mean,color=color,linestyle='-')
            # bplot = axes[1,i].boxplot(differences, patch_artist=True, positions=np.logspace(-2,0,len(segs)))
            # for patch, color in zip(bplot['boxes'], colors):
            #     patch.set_facecolor(color)

        for ax in axes:
            ax.tick_params(labelsize=12,which='both',direction='in',tick1On=True, tick2On=True)
            ax.minorticks_on()
        for ax in axes:
            ax.set_xscale("log")
            ax.set_yscale("log")
        for ax in axes:
            ax.set_ylim(2e44,2e49)
        for ax in axes:
            ax.set_xlim(2e-2,4)
        axes[0].set_ylabel(r"$E_{\rm k}$ [erg]",fontsize=12)
        # axes[1,0].set_ylabel(r"$\Delta \log_{10}(E_{\rm k})$",fontsize=12)
        for ax, (sim,sim_dict) in zip(axes,df.iterrows()):
            ax.set_title(sim_dict['label'],fontsize=12)

        for ax in axes:
            ax.set_xlabel(r"$\Gamma\beta$",fontsize=12)
        # plt.colorbar()
        fig.legend(fancybox=False, loc='center', columnspacing=0.8,
                   bbox_to_anchor=(0.5,0.97),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
                   shadow=False, ncol=len(segs), fontsize=12, labelcolor='black',
                   framealpha=0.0, borderaxespad=0.)

        n = 3
        ax = axes[0]
        ax.plot([1e-4,1e-2], [1e39,1e41], color='gray', ls='-', label='NR',lw=1.)#, lw=0.7, drawstyle='steps')
        ax.plot([1e-4,1e-2], [1e39,1e41], color='gray', ls='--',label=o_fit.name,lw=1)#, lw=0.7, drawstyle='steps')
        ax.plot([1e-4,1e-2], [1e39,1e41], color='gray', ls=':',label=o_fit.name+r"$^*$",lw=1)#, lw=0.7, drawstyle='steps')
        han, lab = ax.get_legend_handles_labels()
        ax.add_artist(ax.legend(han[len(han) - n:], lab[len(lab) - n:],
                                **dict(fancybox=False,loc= 'lower center',columnspacing=0.4,
                                       #"bbox_to_anchor": (0.5, 1.2),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
                                       shadow=False, ncol= 1, fontsize= 12,framealpha=0., borderaxespad= 0., frameon=False)))
        #
        # ax.add_artist(ax.legend(han[:-1 * n], lab[:-1 * n],
        #                         **dict(fancybox=False,loc= 'lower center',columnspacing=0.4,
        #                                #"bbox_to_anchor": (0.5, 1.2),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
        #                                shadow=False, ncol= 1, fontsize= 12,framealpha=0., borderaxespad= 0., frameon=False)))
        # ax.add_artist(ax.legend(han[len(han) - n:], lab[len(lab) - n:],
        #                         **dict(fancybox=False,loc= 'upper left',columnspacing=0.4,
        #                                #"bbox_to_anchor": (0.5, 1.2),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
        #                                shadow=False, ncol= 1, fontsize= 12,framealpha=0., borderaxespad= 0., frameon=False)))

        # plt.sub()
        # for ax in axes:
        # ax[0].set_ylim(1e42,1e49)
        # ax[1].set_ylim(-2,2)
        # for ax
        # ax.set_xlabel("Polar angle [deg]",fontsize=12)
        plt.subplots_adjust(left=0.05,right=0.99,wspace=0.05)
        figname_ = os.getcwd()+f'/figs/{figname}_{o_fit.name}_ek_tall'
        plt.savefig(figname_+'.png',dpi=256)
        plt.savefig(figname_+'.pdf')
        plt.show()

    if (True & plot_fit_ek):
        fig,axes = plt.subplots(ncols=len(df),nrows=2,figsize=(12,5),layout='constrained',sharex='col',sharey='row')

        for i, (name, sim_dict) in enumerate(df.iterrows()):
            if sim_ and name != sim_:
                continue
            for (segment, color) in zip(segs, colors):
                c_seg = float( (segment[1]+segment[0])*0.5 )
                tbl = result_9seg_xy[f"{name} {c_seg:.3f}"]
                l_mom, l_ek = tbl["nr"][:,0], tbl["nr"][:,1]
                l_mom_fit, l_ek_fit = tbl["fit"][:,0], tbl["fit"][:,1]
                axes[0,i].plot(un_log(l_mom),un_log_y(l_ek),color=color,ls='-',lw=0.7)
                axes[0,i].plot(un_log(l_mom_fit),un_log_y(l_ek_fit),color=color,ls='--',lw=0.7)

                l_ek_fit = interp1d(l_mom_fit,l_ek_fit,kind='linear')(l_mom)
                difference = l_ek-l_ek_fit
                axes[1,i].plot(un_log(l_mom),difference,color=color,ls='-',lw=0.7)
            # pass
            # df_fit_total = result_1seg[name]
            # angles = np.array( result_9seg[name].index, dtype=np.float64 )
            # slope, intercept = y0_fit_dict["slope"], y0_fit_dict["intercept"]
            # slope = np.float64( result_1seg[name]["slope"] )
            # intercept = np.float64( result_1seg[name]["intercept"] )
            # coeffs_tot = [np.average(np.array(result_9seg[name][key])) for key in keys]

            # vinf_ = np.linspace(0.005,1, len(l_ek)+1)
            # mom_fit = PBA.MomFromBeta(vinf_)[:-1]
            # coeffs_tot = np.array(coeffs_tot).flatten().tolist()
            # ek_fit_tot = o_fit.piecewise_power(mom_fit, *coeffs_tot)

            # axes[0,i].plot(mom_fit,ek_fit_tot,color='gray',ls='-',lw=0.7)



            # y0_tot = float( df_fit_total["y0"] )
            # differences = []
            for j, (segment, color) in enumerate( zip(segs, colors) ):
                c_seg = float( (segment[1]+segment[0])*0.5 )
                tbl = result_9seg_xy[f"{name} {c_seg:.3f}"]
                l_mom, l_ek = tbl["nr"][:,0], tbl["nr"][:,1]
                # l_mom_fit, l_ek_fit = tbl["fit"][:,0], tbl["fit"][:,1]
                l_mom_ufit, l_ek_ufit = tbl["fit_u"][:,0], tbl["fit_u"][:,1]


                # y0_fit = slope * np.sin(c_seg*np.pi/180) + intercept
                # y0_fit = 10**y0_fit
                #
                # coeffs_tot[keys.index("y0")] = y0_fit
                # ek_fit = o_fit.piecewise_power(mom_fit, *np.array(coeffs_tot).flatten().tolist())

                axes[0,i].plot(un_log(l_mom_ufit), un_log_y(l_ek_ufit), color=color,ls=':',lw=0.7)

                # tbl = result_9seg_xy[f"{name} {c_seg}"]
                # l_mom, l_ek = tbl["nr"][:,0], tbl["nr"][:,1]
                # l_mom_fit, l_ek_fit = tbl["fit"][:,0], tbl["fit"][:,1]
                l_ek_fit = interp1d(l_mom_ufit,l_ek_ufit,kind='linear')(l_mom)
                difference = l_ek-l_ek_fit
                axes[1,i].plot(un_log(l_mom), difference,color=color,ls=':',lw=0.7)
                # differences.append(difference)
                # mean = np.mean(difference)
            # axes[1,i].axhline(y=mean,color=color,linestyle='-')
            # bplot = axes[1,i].boxplot(differences, patch_artist=True, positions=np.logspace(-2,0,len(segs)))
            # for patch, color in zip(bplot['boxes'], colors):
            #     patch.set_facecolor(color)

        for ax in axes:
            for ax in ax:
            # ax.set_xlim(.1,90.1)
                ax.tick_params(labelsize=12,which='both',direction='in',tick1On=True, tick2On=True)
                ax.minorticks_on()
        for ax in axes[0]:
            ax.set_xscale("log")
            ax.set_yscale("log")
        axes[0][0].set_xlim(1e-2,3)
        for ax in axes[0]:
            ax.set_ylim(2e45,2e49)
        for ax in axes[1]:
            ax.set_ylim(-.9,.9)
            ax.grid(lw=0.5)
        for ax in axes[1]:
            ax.set_xlim(2e-2,3)
        axes[0,0].set_ylabel(r"$E_{\rm k}$ [erg]",fontsize=12)
        axes[1,0].set_ylabel(r"$\Delta \log_{10}(E_{\rm k})$",fontsize=12)
        for ax, (sim,sim_dict) in zip(axes[0],df.iterrows()):
            ax.set_title(sim_dict['label'],fontsize=12)

        for ax in axes[1]:
            ax.set_xlabel(r"$\Gamma\beta$",fontsize=12)
        # for ax in axes:
            # ax[0].set_ylim(1e42,1e49)
            # ax[1].set_ylim(-2,2)
        # for ax
        # ax.set_xlabel("Polar angle [deg]",fontsize=12)

        figname_ = os.getcwd()+f'/figs/{figname}_{o_fit.name}_ek_with_residuals'
        plt.savefig(figname_+'.png',dpi=256)
        plt.savefig(figname_+'.pdf')
        plt.show()

    if (True & plot_box_plots):

        fig,axes = plt.subplots(ncols=len(df),nrows=1,figsize=(12,3),layout='constrained',sharex='all',sharey='all')

        for i, (name, sim_dict) in enumerate(df.iterrows()):
            if sim_ and name != sim_:
                continue

            # y0_tot = float( df_fit_total["y0"] )
            differences_i = []
            differences_u = []
            for j, (segment, color) in enumerate( zip(segs, colors) ):
                c_seg = float( (segment[1]+segment[0])*0.5 )
                tbl = result_9seg_xy[f"{name} {c_seg:.3f}"]
                l_mom, l_ek = tbl["nr"][:,0], tbl["nr"][:,1]
                # l_mom_fit, l_ek_fit = tbl["fit"][:,0], tbl["fit"][:,1]
                # l_mom_ufit, l_ek_ufit = tbl["fit_u"][:,0], tbl["fit_u"][:,1]

                l_ek_fit = interp1d(tbl["fit_u"][:,0], tbl["fit_u"][:,1], kind='linear')(l_mom)
                difference = l_ek-l_ek_fit
                # axes[1,i].plot(un_log(l_mom), difference,color=color,ls=':',lw=0.7)
                differences_u.append(difference)
                # mean = np.mean(difference)

                l_ek_fit = interp1d(tbl["fit"][:,0], tbl["fit"][:,1], kind='linear')(l_mom)
                difference = l_ek-l_ek_fit
                # axes[1,i].plot(un_log(l_mom), difference,color=color,ls=':',lw=0.7)
                differences_i.append(difference)
                # mean = np.mean(difference)

            def set_box(bplot, color, alpha, lw, markersize=0):
                box_color = color
                whisker_color =color
                cap_color =color
                median_color = color
                flier_color = color

                for patch in bplot['boxes']:
                    patch.set_facecolor(box_color)
                    patch.set_alpha(alpha)

                for whisker in bplot['whiskers']:
                    whisker.set_color(whisker_color)
                    whisker.set_linewidth(lw)  # You can set the line width if you want

                for cap in bplot['caps']:
                    cap.set_color(cap_color)
                    cap.set_linewidth(lw)  # You can set the line width if you want

                for median in bplot['medians']:
                    median.set_color(median_color)
                    median.set_linewidth(lw)  # You can set the line width if you want

                for flier in bplot['fliers']:
                    flier.set_markeredgecolor(flier_color)
                    flier.set_marker('o')  # You can change the marker style if you want
                    flier.set_markersize(markersize)  # You can change the size of the markers

            # axes[1,i].axhline(y=mean,color=color,linestyle='-')
            bplot = axes[i].boxplot(differences_u, patch_artist=True, positions=np.array(csegs,dtype=int),widths=4)
            set_box(bplot,color='blue',alpha=0.4,lw=1,markersize=0)

            bplot = axes[i].boxplot(differences_i, patch_artist=True, positions=np.array(csegs,dtype=int),widths=8)
            set_box(bplot,color='red',alpha=0.4,lw=2,markersize=0)

        legend_elements = [
            Patch(facecolor='red', alpha=0.4, label=o_fit.name),
            Patch(facecolor='blue', alpha=0.4, label=o_fit.name+"$^{*}$"),
            # Patch(color='blue', label='Whiskers (Range)'),
            # Patch(color='blue', label='Median'),
            # Patch(color='blue', label='Caps (Max/Min)'),
            # Patch(color='blue', label='Outliers', marker='o', markersize=5)  # Only if you have fliers
        ]

        # Adding the legend to the splot
        axes[0].legend(handles=legend_elements,
                       fancybox=False,loc= 'lower left',columnspacing=0.4,
                       #"bbox_to_anchor": (0.5, 1.2),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
                       shadow=False, ncol= 1, fontsize= 12,framealpha=0., borderaxespad= 0., frameon=False)

        for ax in axes:
            ax.tick_params(labelsize=12,which='both',direction='in',tick1On=True, tick2On=True)
            ax.minorticks_on()
            ax.set_ylim(-1.,1.)
            for tick in ax.get_xticklabels():
                tick.set_rotation(45)
            ax.set_xlabel(r"$\theta_{\rm c}$ [deg]",fontsize=12)
            # ax.set_xlim(-1,91)
        axes[0].set_ylabel(r"$\Delta\log_{10}(E_{\rm k})$",fontsize=12)

        for ax, (sim,sim_dict) in zip(axes,df.iterrows()):
            ax.set_title(sim_dict['label'],fontsize=12)
        # for ax in axes[0]:
        #     ax.set_xscale("log")
        #     ax.set_yscale("log")
        # axes[0][0].set_xlim(1e-2,3)
        # for ax in axes[0]:
        #     ax.set_ylim(1e44,2e48)
        # for ax in axes[1]:
        #     ax.set_ylim(-1.,1.)
        # for ax in axes[1]:
        #     ax.set_xlim(1e-2,3)
        # for ax in axes:
        # ax[0].set_ylim(1e42,1e49)
        # ax[1].set_ylim(-2,2)
        # for ax
        # ax.set_xlabel("Polar angle [deg]",fontsize=12)
        figname_ = os.getcwd()+f'/figs/{figname}_{o_fit.name}_residuals_barplot'
        plt.savefig(figname_+'.png',dpi=256)
        plt.savefig(figname_+'.pdf')
        plt.show()


        exit(0)



    exit(0)











    fit_data_segs = dict()
    for cseg in csegs:
        for (name, sim_dic) in df.iterrows():
            fit_data_segs[f"{int(cseg)} {name}"] = dict()


    i = 0
    for (name, sim_dic) in df.iterrows():
        text_dict = df_text.loc[name]
        # sim_dic = sim_dic[1]
        if sim_ and name != sim_:
            continue

        sim_dict = df.loc[name]
        text_dict = df_text.loc[name]

        ej = PBA.id_kenta.EjectaData(get_ej_data(sim_dict['name']),verbose=True)
        mom, ctheta, ek, mass = get_2d_ek(ej=ej,
                                          # text=get_text_max_mass(ej=ej,crit='fast')
                                          text=float(sim_dict['tmerg'])+float(text_dict['text'])
                                          )
        # ek = np.sum(mass,axis=1)
        # ek = np.sum(ek,axis=1)

        low = 30. * np.pi / 180.
        mid = 60. * np.pi / 180.

        ek_tot = np.sum(ek, axis=1)
        ek_equatorial = np.sum(ek[:,(ctheta >= mid)],axis=1) #np.sum(ek[:,(ctheta >= mid)], axis=1)
        ek_mid = np.sum(ek[:,(ctheta > low) & (ctheta <= mid)], axis=1)
        ek_polar = np.sum(ek[:,(ctheta <= low)], axis=1)


        ax = axes[0][i]
        ax.plot(mom,ek_equatorial, color='red')
        ax.plot(mom,ek_mid, color='green')
        ax.plot(mom,ek_polar, color='blue')

        ax = axes[1][i]
        # data = np.insert(np.log10(mom), 0, -1)
        # for segment, color in zip([(0.,15.),(15,30.),(30.,45.),(45.,60,),(60.,75.),(75.,90.)],
        #                           ["blue","cyan","lime","green","orange","red"]):


        for segment, color in zip([(0.,18.),(18,36.),(36.,54.),(54.,72,),(72.,90.)],
                                  ["blue","lime","green","orange","red"]):
            cseg = int( (segment[1]+segment[0])*0.5 )
            low = segment[0] * np.pi / 180.
            up  = segment[1] * np.pi / 180.
            ek_ = np.sum(ek[:,(ctheta > low) & (ctheta <= up)],axis=1)

            l_mom, l_ek = do_log(mom), do_log_y(ek_)
            mask = np.isfinite(l_ek)
            l_mom = l_mom[mask]
            l_ek = l_ek[mask]

            x_fit, y_fit, fit_dict = fit_data(x=l_mom, y=l_ek, undo_x=un_log, undo_y=un_log_y, name=name,
                                              fitting_obj=o_fit,save_res=False)

            ax.plot(mom,ek_,color=color)
            ax.plot(un_log(x_fit),un_log_y(y_fit),color=color,ls=':')

            fit_data_segs[f"{int(cseg)} {name}"] = fit_dict

        #     data = np.column_stack((data, np.insert(np.log10(ek_),0, 0.5*(up+low) )))
        #
        # np.savetxt(os.getcwd()+'/'+name+"_2D_6seg_log_ek_log_mom.txt",X=data,fmt="%.3f")

        # low = 30. * np.pi / 180.
        # mid = 60. * np.pi / 180.
        #
        # ek_tot = np.sum(ek, axis=1)
        # ek_equatorial = np.sum(ek[:,(ctheta >= mid)],axis=1) #np.sum(ek[:,(ctheta >= mid)], axis=1)
        # ek_mid = np.sum(ek[:,(ctheta > low) & (ctheta <= mid)], axis=1)
        # ek_polar = np.sum(ek[:,(ctheta <= low)], axis=1)
        i+=1
    for ax in axes:
        for ax in ax:
            ax.set_xscale("log")
            ax.set_yscale("log")
    plt.show()

    print(f"ctheta:{ctheta.shape} mom:{mom.shape} ek:{ek.shape}")





    ''' -------------------------------------------------------------------------- '''

    default_pysr_params = dict(
        procs=4,
        populations=8,
        population_size=50,
        ncycles_per_iteration=500,
        early_stop_condition=(
            "stop_if(loss, complexity) = loss < 1e-6 && complexity < 10"
            # Stop early if we find a good and simple equation
        ),
        model_selection="best",
        weight_randomize=0.1, # ^ Randomize the tree much more frequently
    )
    model = PySRRegressor(
        niterations=100,
        binary_operators=["+", "*", "-", "/"],
        unary_operators=["cos", "exp", "sin", "square", "cube", "cos2(x)=cos(x)^2", "sin2(x)=sin(x)^2"],
        extra_sympy_mappings={"cos2": lambda x: sympy.cos(x)**2,
                              "sin2": lambda x: sympy.sin(x)**2},
        **default_pysr_params,
    )

    log_ek = np.log10(ek_tot)
    mom_ = mom[np.isfinite(log_ek)]
    log_ek = log_ek[np.isfinite(log_ek)]
    log_ek = log_ek[mom_>0.1]
    mom_ = mom_[mom_>0.1]
    model.fit(mom_.reshape(-1, 1),log_ek.reshape(-1, 1))#(X, y)
    print("Best:", model.sympy())
    print(model)
    # for i in range(4):
    #     print(f"Eq{i}: ", model.sympy(i))
    ypredict = model.predict(mom_.reshape(-1, 1))
    print("Default selection MSE:", np.power(ypredict - log_ek, 2).mean())
    ypredict = 10**ypredict.flatten()

    # exit(1)

    fig, axes = plt.subplots(nrows=2,ncols=1,sharex='all',layout='constrained')

    ax = axes[0]
    do_cumulative = True
    get_cumulative = lambda val : np.cumsum(val[::-1])[::-1] if do_cumulative else val
    ax.plot(mom,get_cumulative(ek_tot),color='black',label=r'Total')
    ax.plot(mom_,get_cumulative(ypredict),color='gray')
    ax.plot(mom,get_cumulative(ek_equatorial),label=r'$\theta > 60$')
    ax.plot(mom,get_cumulative(ek_mid),label=r'$30 < \theta < 60$')
    ax.plot(mom,get_cumulative(ek_polar),label=r'$\theta < 30$')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(1e-2,4)
    # ax.set_ylim(1e21,1e31)
    ax.legend(fancybox=False,loc= 'lower left',columnspacing=0.4,
              #"bbox_to_anchor": (0.5, 1.2),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
              shadow=False, ncol= 1, fontsize= 12,
              framealpha=0., borderaxespad= 0., frameon=False)
    ax.set_ylabel(r"$E_k$")

    ax = axes[1]
    ax.pcolormesh(mom, ctheta*180/np.pi, ek.T, cmap=plt.get_cmap('jet'),
                  norm=LogNorm(vmin=ek.max()*1e-7, vmax=ek.max()))
    ax.set_xscale('log')
    ax.set_xlabel(r"$\Gamma\beta$")
    ax.set_ylabel(r"Polar angle [deg]")

    plt.show()

def plot_sim_ekecta_mass_text_dependency(name:str):

    # fig, axes = plt.subplots(nrows=2,ncols=1,sharex='col',sharey='row',layout='constrained')


    sim_dict = df.loc[name]
    #
    ej = PBA.id_kenta.EjectaData(get_ej_data(sim_dict["name"]),verbose=True)

    mom, ctheta, = np.zeros(0,), np.zeros(0,)
    ek, mass = [], []
    for text in ej.getText():
        mom, ctheta, ek_, mass_ = get_2d_ek(ej=ej,text=text,new_vinf_len=None)
        ek.append(ek_)
        mass.append(mass_)

    ek_vs_mom = np.column_stack([np.sum(ek_,axis=1) for ek_ in ek])
    ek_vs_theta = np.column_stack([np.sum(ek_,axis=0) for ek_ in ek])

    cmap = plt.get_cmap('jet')

    fig, axes = plt.subplots(nrows=2,ncols=1,sharex='col',sharey='row',layout='constrained')

    # ax = axes[0]
    im = axes[0].pcolormesh(ej.getText()-sim_dict["tmerg"], mom, ek_vs_mom, cmap=cmap,
                            norm=LogNorm(vmin=ek_vs_mom.max()*1e-5, vmax=ek_vs_mom.max()))

    im = axes[1].pcolormesh(ej.getText()-sim_dict["tmerg"], ctheta*180/np.pi, ek_vs_theta, cmap=cmap,
                            norm=LogNorm(vmin=ek_vs_theta.max()*1e-5, vmax=ek_vs_theta.max()))
    axes[0].set_yscale('log')
    for ax in axes:
        ax.set_xlim(0.,ej.getText()[-1]-sim_dict["tmerg"])
    # ax0.set_yscale('log')
    fig.colorbar(im, ax=axes[0])
    fig.colorbar(im, ax=axes[1])
    plt.show()


    #
    # ej = PBA.id_kenta.EjStruct(get_ej_data(sim_dict["name"]),verbose=True)
    # data = ej.get_2D_id(text=ej.getText()[0], method_r0="from_beta",t0=1e3)
    #
    # ctheta = data["ctheta"][0,:]
    # mom = data["mom"][:,0]
    # ek = np.column_stack([np.cumsum(ej.get_2D_id(text=text, method_r0="from_beta",t0=1e3),axis=0)] for text in ej.getText())
    #
    # cmap = plt.get_cmap('jet')
    #
    # fig, axes = plt.subplots(nrows=2,ncols=1,sharex='all')
    #
    # # ax = axes[0]
    # im = axes[0].pcolormesh(ej.getText()-sim_dict["tmerg"], ctheta, ek.T, cmap=plt.get_cmap('jet'),
    #                         norm=LogNorm(vmin=ek.max()*1e-7, vmax=ek.max()))
    # # ax0.set_yscale('log')
    # fig.colorbar(im, ax=axes[0])
    # plt.show()

    texts = ej.getText()
    vinf = ej.get_vinf()
    thetas = ej.get_theta()
    print(ej.get_vinf().shape, ej.get_theta().shape,ej.get(v_n="mass",text=ej.getText()[0]).shape)

    # mass_for_vin/f = lambda vinf, text : ej.get(v_n="mass",text=text)[:,find_nearest_index(ej.get_vinf(),vinf)]

    # mass_vs_vinf = lambda text : np.sum(ej.get(v_n="mass",text=text), axis=0)
    # mass_vs_theta = lambda text : np.sum(ej.get(v_n="mass",text=text), axis=1)

    # res = []
    # for text in texts:
    #     res.append(mass_vs_vinf(text))
    # res = np.reshape(np.array(res),newshape=(len(texts),(len(vinf))))
    vinf_c = 0.5 * (vinf[1:] + vinf[:-1])
    mom_c = PBA.MomFromBeta(vinf_c) #get_Gamma(vinf_c) * vinf_c
    mom_c = mom_c[:-1]
    theta_c = 0.5 * (thetas[1:] + thetas[:-1])



    def get_mass_vs_vinf(angle:str):
        if angle == "equatorial":
            theta_l = 60. * np.pi / 180
            theta_h = 90. * np.pi / 180
            np.column_stack([
                np.sum(ej.get(v_n="mass",text=text)[thetas ], axis=0) for text in texts
            ]).T

    mass_vs_vinf = np.column_stack([np.sum(ej.get(v_n="mass",text=text), axis=0) for text in texts]).T
    mass_vs_theta = np.column_stack([np.sum(ej.get(v_n="mass",text=text), axis=1) for text in texts]).T

    mass_vs_vinf_c = 0.5 * (mass_vs_vinf[:,1:] + mass_vs_vinf[:,1:]) [:,:-1]
    mass_vs_theta_c = 0.5 * (mass_vs_theta[:,1:] + mass_vs_theta[:,1:])


    cmap = plt.get_cmap('jet')

    fig, axes = plt.subplots(nrows=2,ncols=1,sharex='all')

    # ax = axes[0]
    im = axes[0].pcolormesh(texts-sim_dict["tmerg"], mom_c, mass_vs_vinf_c.T, cmap=cmap,
                            norm=LogNorm(vmin=mass_vs_vinf_c.max()*1e-7, vmax=mass_vs_vinf_c.max()))

    im = axes[1].pcolormesh(texts-sim_dict["tmerg"], theta_c, mass_vs_theta_c.T, cmap=cmap,
                            norm=LogNorm(vmin=mass_vs_theta_c.max()*1e-7, vmax=mass_vs_theta_c.max()))
    # ax0.set_yscale('log')
    fig.colorbar(im, ax=axes[0])
    fig.colorbar(im, ax=axes[1])
    plt.show()
def plot_sim_ekecta_mass_text_dependency_all():
    fig, axes = plt.subplots(nrows=2,ncols=len(df),sharex='col',sharey='row',layout='constrained',figsize=(12,4))

    for i_s, (name, sim_dict) in enumerate(df.iterrows()):

        sim_dict = df.loc[name]
        #
        ej = PBA.id_kenta.EjectaData(get_ej_data(sim_dict["name"]),verbose=True)

        mom, ctheta, = np.zeros(0,), np.zeros(0,)
        ek, mass = [], []
        for text in ej.getText():
            mom, ctheta, ek_, mass_ = get_2d_ek(ej=ej,text=text,new_vinf_len=None)
            ek.append(ek_)
            mass.append(mass_)

        ek_vs_mom = np.column_stack([np.sum(ek_,axis=1)/np.sum(ek_,axis=1).max() for ek_ in ek])
        ek_vs_theta = np.column_stack([np.sum(ek_,axis=0)/np.sum(ek_,axis=0).max() for ek_ in ek])

        cmap = plt.get_cmap('jet')


        # ax = axes[0]
        im1 = axes[0][i_s].pcolormesh(ej.getText()-sim_dict["tmerg"], mom, ek_vs_mom, cmap=cmap,
                                # norm=LogNorm(vmin=ek_vs_mom.max()*1e-5, vmax=ek_vs_mom.max()))
                                norm=LogNorm(vmin=1e-7, vmax=1))

        im2 = axes[1][i_s].pcolormesh(ej.getText()-sim_dict["tmerg"], ctheta*180/np.pi, ek_vs_theta, cmap=cmap,
                                # norm=LogNorm(vmin=ek_vs_theta.max()*1e-5, vmax=ek_vs_theta.max()))
                                norm=LogNorm(vmin=1e-7, vmax=1))
        axes[0][i_s].set_yscale('log')
        # for ax in axes[i_s]:
        axes[-1][i_s].set_xlim(0.,ej.getText()[-1]-sim_dict["tmerg"])



    # ax0.set_yscale('log')
    fig.colorbar(im1, ax=axes[0][-1])
    fig.colorbar(im2, ax=axes[1][-1])
    plt.show()


    #
    # ej = PBA.id_kenta.EjStruct(get_ej_data(sim_dict["name"]),verbose=True)
    # data = ej.get_2D_id(text=ej.getText()[0], method_r0="from_beta",t0=1e3)
    #
    # ctheta = data["ctheta"][0,:]
    # mom = data["mom"][:,0]
    # ek = np.column_stack([np.cumsum(ej.get_2D_id(text=text, method_r0="from_beta",t0=1e3),axis=0)] for text in ej.getText())
    #
    # cmap = plt.get_cmap('jet')
    #
    # fig, axes = plt.subplots(nrows=2,ncols=1,sharex='all')
    #
    # # ax = axes[0]
    # im = axes[0].pcolormesh(ej.getText()-sim_dict["tmerg"], ctheta, ek.T, cmap=plt.get_cmap('jet'),
    #                         norm=LogNorm(vmin=ek.max()*1e-7, vmax=ek.max()))
    # # ax0.set_yscale('log')
    # fig.colorbar(im, ax=axes[0])
    # plt.show()

    texts = ej.getText()
    vinf = ej.get_vinf()
    thetas = ej.get_theta()
    print(ej.get_vinf().shape, ej.get_theta().shape,ej.get(v_n="mass",text=ej.getText()[0]).shape)

    # mass_for_vin/f = lambda vinf, text : ej.get(v_n="mass",text=text)[:,find_nearest_index(ej.get_vinf(),vinf)]

    # mass_vs_vinf = lambda text : np.sum(ej.get(v_n="mass",text=text), axis=0)
    # mass_vs_theta = lambda text : np.sum(ej.get(v_n="mass",text=text), axis=1)

    # res = []
    # for text in texts:
    #     res.append(mass_vs_vinf(text))
    # res = np.reshape(np.array(res),newshape=(len(texts),(len(vinf))))
    vinf_c = 0.5 * (vinf[1:] + vinf[:-1])
    mom_c = PBA.MomFromBeta(vinf_c) #get_Gamma(vinf_c) * vinf_c
    mom_c = mom_c[:-1]
    theta_c = 0.5 * (thetas[1:] + thetas[:-1])



    def get_mass_vs_vinf(angle:str):
        if angle == "equatorial":
            theta_l = 60. * np.pi / 180
            theta_h = 90. * np.pi / 180
            np.column_stack([
                np.sum(ej.get(v_n="mass",text=text)[thetas ], axis=0) for text in texts
            ]).T

    mass_vs_vinf = np.column_stack([np.sum(ej.get(v_n="mass",text=text), axis=0) for text in texts]).T
    mass_vs_theta = np.column_stack([np.sum(ej.get(v_n="mass",text=text), axis=1) for text in texts]).T

    mass_vs_vinf_c = 0.5 * (mass_vs_vinf[:,1:] + mass_vs_vinf[:,1:]) [:,:-1]
    mass_vs_theta_c = 0.5 * (mass_vs_theta[:,1:] + mass_vs_theta[:,1:])


    cmap = plt.get_cmap('jet')

    fig, axes = plt.subplots(nrows=2,ncols=1,sharex='all')

    # ax = axes[0]
    im = axes[0].pcolormesh(texts-sim_dict["tmerg"], mom_c, mass_vs_vinf_c.T, cmap=cmap,
                            norm=LogNorm(vmin=mass_vs_vinf_c.max()*1e-7, vmax=mass_vs_vinf_c.max()))

    im = axes[1].pcolormesh(texts-sim_dict["tmerg"], theta_c, mass_vs_theta_c.T, cmap=cmap,
                            norm=LogNorm(vmin=mass_vs_theta_c.max()*1e-7, vmax=mass_vs_theta_c.max()))
    # ax0.set_yscale('log')
    fig.colorbar(im, ax=axes[0])
    fig.colorbar(im, ax=axes[1])
    plt.show()


def plot_sim_ekecta_mass_old(sim_name:str):

    text_dict = df_text.loc[sim_name]
    sim_dict = df.loc[sim_name]

    ej = PBA.id_kenta.EjectaData(get_ej_data(sim_dict["name"]),verbose=True)
    texts = ej.getText()
    vinf = ej.get_vinf()
    thetas = ej.get_theta()
    print(ej.get_vinf().shape,
          ej.get_theta().shape,
          ej.get(v_n="mass",text=text_dict['text']+sim_dict['tmerg']).shape)
    vinf_centers = 0.5 * (vinf[1:]+vinf[:-1])
    vinf_centers = vinf_centers[:-1]
    moms = PBA.MomFromBeta(vinf_centers)
    ej_data_vs_mom,ej_data_vs_theta = [],[]
    for text in texts:
        ej_mass = ej.get(v_n="mass",text=text)
        ej_ek = ej_mass * PBA.utils.cgs.solar_m * (vinf[np.newaxis, :] * PBA.utils.cgs.c) ** 2
        ej_ek_centers = 0.5*(ej_ek[:,1:] + ej_ek[:,:-1])[:,:-1]
        ej_cumulative_ek_vs_mom = np.cumsum(np.sum(ej_ek_centers,axis=0)[::-1])[::-1]/np.sum(ej_ek_centers)
        ej_data_vs_mom.append(ej_cumulative_ek_vs_mom)

        ej_ek_centers = 0.5*(ej_ek[1:,:] + ej_ek[:-1,:])
        ej_data_vs_theta.append(np.sum(ej_ek_centers,axis=0)/np.sum(ej_ek_centers))

    ej_data_vs_mom = np.reshape(np.array(ej_data_vs_mom),newshape=(len(texts),len(ej_data_vs_mom[0])))
    ej_data_vs_theta = np.reshape(np.array(ej_data_vs_theta),newshape=(len(texts),len(ej_data_vs_theta[0])))
    ej_data_vs_mom[~np.isfinite(ej_data_vs_mom)] = 0.
    ej_data_vs_theta[~np.isfinite(ej_data_vs_theta)] = 0.

    levels = MaxNLocator(nbins=15).tick_values(ej_data_vs_mom.min(), ej_data_vs_mom.max())
    # pick the desired colormap, sensible levels, and define a normalization
    # instance which takes data values and translates those into levels.
    cmap = plt.get_cmap('viridis')
    norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
    norm = LogNorm(vmin=ej_data_vs_mom.max()*1e-7, vmax=ej_data_vs_mom.max())

    fig, axes = plt.subplots(nrows=2,ncols=1)

    ax = axes[0]
    im = ax.pcolormesh(texts-sim_dict["tmerg"], moms, ej_data_vs_mom.T, cmap=cmap, norm=norm)
    im = ax.pcolormesh(texts-sim_dict["tmerg"], thetas, ej_data_vs_theta.T, cmap=cmap, norm=norm)
    # ax0.set_yscale('log')
    fig.colorbar(im, ax=ax)
    ax.set_title('pcolormesh with levels')
    plt.show()



''' --------------------------------------------------- '''

def plot_rho_max(df:pd.DataFrame):

    fig,ax = plt.subplots(ncols=1,nrows=1)

    for name, sim_dic in df.iterrows():
        print(sim_dic["name"], sim_dic["rhomax"], sim_dic["mdot_extract"])
        if (sim_dic["rhomax"]!=-1):
            data = PBA.id_kenta.Data(fpath_rhomax=DATA_PATH+sim_dic["name"]+'/'+sim_dic["rhomax"],
                                     fpath_mdot=DATA_PATH+sim_dic["name"]+'/'+sim_dic["mdot_extract"])
            # df = data.df_mdot
            tmerg = sim_dic["tmerg"]
            # times = df["time"]
            # times = times - tmerg

            ax.plot(data.df_rho["time"]-tmerg, data.df_rho[f"rho_max"], color=sim_dic['color'],lw=1,alpha=.8,label=sim_dic['label'],ls=sim_dic['ls'])
    ax.grid(color='gray',ls='--',alpha=.5)
    ax.legend()
    ax.set_yscale("log")
    plt.show()

# Definition of a retarded time for ejecta mass flux
def get_tret(r : float, vinf_ave0 : float =.7) -> float:
    vinf_ave0 = .7
    return (r*1000.*100./vinf_ave0/PBA.utils.cgs.c) * 1e3 # ms
def get_tret_arr(r : float, vinf_ave0 : np.ndarray) -> np.ndarray:
    # vinf_ave0 = .7
    return (r*1000.*100./vinf_ave0/PBA.utils.cgs.c) * 1e3 # ms
def plot_rho_mdot_rext(ylim=(1e-7, 1e-2), ylim2=(1, 6), xlim=(-2, 2.), crit="fast") -> None:

    # color_pal_ = sns.color_palette(n_colors=10)
    ncols = len([(name, dic) for name, dic in df.iterrows() if dic["rhomax"]!=-1])
    fig, axes = plt.subplots(ncols=ncols,
                             nrows=1, figsize=(12,2.5), sharey='all',layout='constrained')
    i = 0
    for name, sim_dic in df.iterrows():

        if sim_dic["rhomax"] == -1 or sim_dic["mdot_extract"] == -1:
            continue

        data = PBA.id_kenta.Data(fpath_rhomax=DATA_PATH+'/'+sim_dic["name"]+'/'+sim_dic["rhomax"],
                                 fpath_mdot=DATA_PATH+'/'+sim_dic["name"]+'/'+sim_dic["mdot_extract"])
        rext = sorted(data.get_rext())
        rext = rext[:-2][::-1]
        norm = Normalize(vmin=0,vmax=len(rext)-2)
        colors=['red','orange','green','blue']

        time = data.df_mdot["time"]
        #color_pal = sns.color_palette(n_colors=len(rext)+1)

        tmerg = sim_dic["tmerg"]

        ax = axes[i]
        for j, r in enumerate(rext):
            tret = get_tret(r=r) # ms
            tret = get_tret_arr(r=r, vinf_ave0=data.df_mdot[f"vave_{crit}"]) # ms
            # ax.plot(time - tret- tmerg, data.df_mdot[f"mdot_fast_r{r}"],
            #             color=color_pal_[j], label=r"$R_{\rm ext}$"+f"={r} [km]",lw=1,alpha=.9)
            vals = data.df_mdot[f"mdot_{crit}_r{r}"]
            if np.sum(vals) > 1e-10:
                ax.plot(time - tret - tmerg, vals, #color=color_pal_[j],
                        color=colors[j],#plt.cm.get_cmap('seismic')(norm(j)),
                        label=r"$R_{\rm ext}$"+f"={r}",lw=0.7,alpha=.8)

        # for i, r in enumerate(rext):
        # ax.plot(data.df_mdot["time"], data.df_mdot[f"mdot_fast_r{r}"],color=color_pal[i])
        # ax.plot(data.ret_time(r0=r,vave_key="vave_fast"), data.df_mdot[f"mdot_fast_r{r}"],color=color_pal[i])
        # df.plot("time",f"mdot_fast_r{r}",ax=ax,color=color_pal[i])


        ax2 = ax.twinx()
        ax2.plot(data.df_rho["time"]-tmerg, data.df_rho[f"rho_max"], color='black',lw=1.0,alpha=.8)

        # df.plot("time",f"rho_max",ax=ax2,color="black")

        # ax.set_xlabel("time [ms]", fontsize=12)
        #ax.set_ylabel("Mdot [Msun/s]", fontsize=12)
        #ax2.set_ylabel("rho max / rho nuc", fontsize=12)
        ax.set_yscale("log")
        # ax2.set_yscale("log")
        # ax.grid(color='black')
        # ax2.grid(color='gray',ls=':',alpha=.5)
        ax2.tick_params(labelsize=12,which='both',direction='in',labelright=True,tick2On=True,right=True,tick1On=True)
        # ax.legend(loc='lower right', fontsize=12)
        ax.tick_params(labelsize=12,which='both',direction='in')
        if len(ylim) > 0 : ax.set_ylim(*ylim)
        if len(ylim2) > 0 : ax2.set_ylim(*ylim2)
        if len(xlim) > 0 : ax.set_xlim(*xlim)

        ax.set_title(sim_dic['label'], fontsize=12)#(r"$\dot{M}$ vs $\rho_{\rm max}$" + f" for {sim_dic['label']}", fontsize=14)

        if i!=ncols-1:
            ax2.yaxis.set_ticklabels([])
        if i == 0:
            ax.legend(fancybox=False,loc= 'upper left',columnspacing=0.2,
                      #"bbox_to_anchor": (0.5, 1.2),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
                      shadow=False, ncol= 1, fontsize= 11,
                      framealpha=0., borderaxespad= 0., frameon=False)
        i+=1

    # axes[-1].legend(loc='lower right', fontsize=12)
    axes[0].set_ylabel(r"$\dot{M}_{\rm ft}$ [$M_{\odot}$/s]", fontsize=12)
    for ax in axes:
        ax.grid(ls=':',lw=0.7)
        ax.set_xlabel(r"$t-t_{\rm merg}$ [ms]", fontsize=12)
        ax.minorticks_on()
        ax2.minorticks_on()

    ax2.set_ylabel(r"$\rho_{\rm max}/\rho_{\rm nuc}$", color="black",fontsize=12)

    # fig.text(0.5, 0.001, "time [ms]", fontsize=14, ha='center')
    # # plt.tight_layout()
    # # fig.text(0.5, 0.001, "time [ms]", fontsize=14, ha='center')
    # plt.savefig(os.getcwd()+'/figs/'+"mdot_rhomax.png",dpi=256)
    # plt.savefig(os.getcwd()+'/figs/s'+"mdot_rext_rhomax.pdf")
    # plt.show()

    # data.df_mdot.plot("time","mdot_fast_r48")
    # data.df_mdot.plot("time","mdot_fast_r144")
    # data.df_mdot.plot("time","mdot_fast_r482")
    # data.df_mdot.plot("time","mdot_fast_r144")
    plt.savefig(os.getcwd()+'/figs/'+'mdot_rhomax.pdf')
    plt.show()

def plot_init_profile(mom : np.ndarray, ctheta : np.ndarray, data2d : np.ndarray,
                      xmin:float=0, xmax:float=90, ymin:float=1e-2, ymax:float=6, vmin:float=1e-12, vmax:float=1e-6,
                      norm_mode:str="log", cmap:plt.cm = plt.get_cmap('RdYlBu_r'),
                      xscale:str="linear", yscale:str="linear", cbar_label:str=r"$M_{\rm ej}$",
                      subplot_mode:str="sum",
                      hist_ylabel:str=r"$M_{\rm ej}$ [M$_{\odot}$]",
                      hist_ylabel2:str=r"$M_{\rm ej}(\Gamma\beta>1) [M$_{\odot}$]",
                      title:str=None, figpath:str=None, fontsize:int=12):
    """

    :param mom: np.ndarray; 1D;
        mom = Gamma * beta, where bete is dimensionless terminal velocity, Gamma is lorentz factor [dimensionless]
    :param ctheta: np.ndarray; 1D;
        polar angle (0-pole, pi/2-equator) [radians]
    :param data2d: np.ndarray 2D [shape=(len(mom),len(ctheta)]
        Mass or kinetic energy to plot; Kinetic energy is defined as $\E_k = mass * (\beta * c)^2$ where
        $mass$ is the mass in grams, $\beta$ is dimensionless velocity and $c$ is the speed of light
    """

    # from radians to dgrees
    ctheta *= 180 / np.pi

    fig = plt.figure(figsize=(5, 5))

    ax0 = fig.add_axes([0.16, 0.12, 0.81 - 0.15, 0.59 - 0.12])
    ax1 = fig.add_axes([0.16, 0.61, 0.81 - 0.15, 0.91 - 0.61])
    cax = fig.add_axes([0.83, 0.12, 0.86 - 0.83, 0.59 - 0.12])

    # set how to deal with 1D histogram (top subplot) Summ all quantities for all angles or compute mass-averaged
    if subplot_mode == 'sum':
        res = np.sum(data2d, axis=0)
        ax1.plot(ctheta, res, color='blue', ls='-', drawstyle='steps')
        ax2 = ax1.twinx()
        ax2.plot(ctheta, np.sum(data2d[mom>1.], axis=0), color='red', ls='-', drawstyle='steps')
        ax2.set_yscale('log')
        ax2.set_ylabel(hist_ylabel2, fontsize=fontsize, color='red')

    elif subplot_mode == 'seg_sum':
        segs = np.arange(start=0,stop=90+18,step=18)
        seg_colors = ["blue","green","magenta","orange","red"]
        segs = [(low,up) for (low,up) in zip(segs[:-1],segs[1:])]

        for (seg_color, segment) in zip(seg_colors, segs):
            c_seg = int( (segment[1] + segment[0])*0.5 )
            low = segment[0]
            up  = segment[1]
            ek_ = np.sum(data2d[:,(ctheta > low) & (ctheta <= up)],axis=0)
            ax1.plot(ctheta, ek_, color=seg_color, ls='-', drawstyle='steps',
                     label=r"$\theta_{\rm c}=$"+f"{int(c_seg)} deg.")
        ax1.legend(fancybox=False,loc= 'lower right',columnspacing=0.4,
                      #"bbox_to_anchor": (0.5, 1.2),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
                      shadow=False, ncol= 2, fontsize= 12,
                      framealpha=0., borderaxespad= 0., frameon=False)
    else:
        vinf = mom / np.sqrt(1.0+mom*mom)
        data2d[~np.isfinite(data2d)] = 0.
        vinf[~np.isfinite(vinf)] = 0.
        res = np.sum(data2d * vinf[:, np.newaxis], axis=0) / np.sum(data2d, axis=0)

        # plot the 1D histogram
        ax1.plot(ctheta, res, color='black', ls='-', drawstyle='steps')


    # adjust subplot settings
    if (not title is None): ax1.set_title(title,fontsize=fontsize)
    if norm_mode=="log": ax1.set_yscale("log")
    else: ax1.set_yscale("linear")
    ax1.set_ylabel(hist_ylabel, fontsize=fontsize,color='blue')
    ax1.get_yaxis().set_label_coords(-0.15, 0.5)
    ax1.set_xlim(xmin, xmax)
    ax1.xaxis.set_ticklabels([])
    for ax in [ax1,ax2]:
        if ax:
            ax.tick_params(axis='both', which='both', labelsize=12, direction='in')
            ax.minorticks_on()


    # prepare for plotting 2D histogram
    if vmin is None: vmin = data2d[(data2d > 0) & (np.isfinite(data2d))].min()
    if vmax is None: vmax = data2d[(data2d > 0) & (np.isfinite(data2d))].max()
    if (norm_mode=="log"): norm = LogNorm(vmin, vmax)
    elif (norm_mode=="linear"): norm = Normalize(vmin, vmax)
    elif (norm_mode=="levels"):
        levels = MaxNLocator(nbins=15).tick_values(vmin, vmax)
        norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
    else:
        raise KeyError(" norm_mode is not recognized ")

    # plot colormesh
    im = ax0.pcolor(ctheta, mom, data2d, cmap=cmap, norm=norm, shading='auto')
    # ax0.axhline(y=1, linestyle='--', color='gray')

    # adjust the bottom subplot
    ax0.set_ylim(ymin, ymax)
    ax0.set_xlim(xmin, xmax)
    ax0.set_xscale(xscale)
    ax0.set_yscale(yscale)
    ax0.set_ylabel(r"$\Gamma\beta$", fontsize=fontsize)
    ax0.set_xlabel(r"Polar angle", fontsize=fontsize)
    ax0.get_yaxis().set_label_coords(-0.15, 0.5)
    # ax0.text(0.75, 0.88, title, color='white', transform=ax0.transAxes, fontsize=fontsize)
    ax0.minorticks_on()
    ax0.tick_params(axis='both', which='both', labelleft=True,
                    labelright=False, tick1On=True, tick2On=True,
                    labelsize=12,
                    direction='in',
                    bottom=True, top=True, left=True, right=True)

    cbar = plt.colorbar(im, cax=cax, norm=norm)
    # cbar.set_label(r"$M_{\rm ej}$ [M$_{\odot}$]", fontsize=fontsize)
    cbar.set_label(cbar_label, fontsize=fontsize)
    cbar.ax.tick_params(axis='both', which='both',
                        labelsize=12, direction='in')
    cbar.ax.minorticks_on()

    # save figures
    if not figpath is None: plt.savefig(figpath+".png", dpi=256)
    if not figpath is None: plt.savefig(figpath+".pdf")
    plt.show()
    plt.close()
def plot_init_profile_new(mom : np.ndarray, ctheta : np.ndarray, data2d : np.ndarray,
                      xmin:float=0, xmax:float=90, ymin:float=1e-2, ymax:float=6, vmin:float=1e-12, vmax:float=1e-6,
                      norm_mode:str="log", cmap:plt.cm = plt.get_cmap('RdYlBu_r'),
                      xscale:str="linear", yscale:str="linear", cbar_label:str=r"$M_{\rm ej}$",
                      subplot_mode:str="sum", hist_ylabel:str=r"$M_{\rm ej}$ [M$_{\odot}$]",
                      title:str=None, figpath:str=None, fontsize:int=12):
    """

    :param mom: np.ndarray; 1D;
        mom = Gamma * beta, where bete is dimensionless terminal velocity, Gamma is lorentz factor [dimensionless]
    :param ctheta: np.ndarray; 1D;
        polar angle (0-pole, pi/2-equator) [radians]
    :param data2d: np.ndarray 2D [shape=(len(mom),len(ctheta)]
        Mass or kinetic energy to plot; Kinetic energy is defined as $\E_k = mass * (\beta * c)^2$ where
        $mass$ is the mass in grams, $\beta$ is dimensionless velocity and $c$ is the speed of light
    """

    # from radians to dgrees
    ctheta *= 180 / np.pi

    fig,axes = plt.subplots(figsize=(5, 5),ncols=1,nrows=2,layout='constrained',sharex='col',
                            gridspec_kw={'height_ratios': [1.5,2]})

    # ax0 = fig.add_axes([0.16, 0.12, 0.81 - 0.15, 0.59 - 0.12])
    # ax1 = fig.add_axes([0.16, 0.61, 0.81 - 0.15, 0.91 - 0.61])
    # cax = fig.add_axes([0.83, 0.12, 0.86 - 0.83, 0.59 - 0.12])

    # set how to deal with 1D histogram (top subplot) Summ all quantities for all angles or compute mass-averaged
    ax1 = axes[0]

    if subplot_mode == 'sum':
        # res = np.sum(data2d, axis=0)
        ax1.plot(ctheta, np.sum(data2d, axis=0), color='black', ls='-', drawstyle='steps')
        ax2 = ax1.twinx()
        ax2.plot(ctheta, np.sum(data2d[mom>1.], axis=0), color='gray', ls='-', drawstyle='steps')
        ax2.set_yscale('log')

    elif subplot_mode == 'seg_sum':
        segs = np.arange(start=0,stop=90+18,step=18)
        seg_colors = ["blue","green","magenta","orange","red"]
        segs = [(low,up) for (low,up) in zip(segs[:-1],segs[1:])]

        for (seg_color, segment) in zip(seg_colors, segs):
            c_seg = int( (segment[1] + segment[0])*0.5 )
            low = segment[0]
            up  = segment[1]
            ek_ = np.sum(data2d[:,(ctheta > low) & (ctheta <= up)],axis=0)
            ax1.plot(ctheta, ek_, color=seg_color, ls='-', drawstyle='steps',
                     label=r"$\theta_{\rm c}=$"+f"{int(c_seg)} deg.")
        ax1.legend(fancybox=False,loc= 'lower right',columnspacing=0.4,
                   #"bbox_to_anchor": (0.5, 1.2),  # loc=(0.0, 0.6),  # (1.0, 0.3), # <-> |
                   shadow=False, ncol= 2, fontsize= 12,
                   framealpha=0., borderaxespad= 0., frameon=False)
    else:
        vinf = mom / np.sqrt(1.0+mom*mom)
        data2d[~np.isfinite(data2d)] = 0.
        vinf[~np.isfinite(vinf)] = 0.
        res = np.sum(data2d * vinf[:, np.newaxis], axis=0) / np.sum(data2d, axis=0)

        # plot the 1D histogram
        ax1.plot(ctheta, res, color='black', ls='-', drawstyle='steps')


    # adjust subplot settings
    if (not title is None): ax1.set_title(title,fontsize=fontsize)
    if norm_mode=="log": ax1.set_yscale("log")
    else: ax1.set_yscale("linear")
    ax1.set_ylabel(hist_ylabel, fontsize=fontsize)
    ax1.get_yaxis().set_label_coords(-0.15, 0.5)
    ax1.set_xlim(xmin, xmax)
    ax1.xaxis.set_ticklabels([])
    ax1.tick_params(axis='both', which='both', labelleft=True,
                    labelright=False, tick1On=True, tick2On=True,
                    labelsize=12, direction='in',
                    bottom=True, top=True, left=True, right=True)
    ax1.minorticks_on()

    # prepare for plotting 2D histogram
    if vmin is None: vmin = data2d[(data2d > 0) & (np.isfinite(data2d))].min()
    if vmax is None: vmax = data2d[(data2d > 0) & (np.isfinite(data2d))].max()
    if (norm_mode=="log"): norm = LogNorm(vmin, vmax)
    elif (norm_mode=="linear"): norm = Normalize(vmin, vmax)
    elif (norm_mode=="levels"):
        levels = MaxNLocator(nbins=15).tick_values(vmin, vmax)
        norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
    else:
        raise KeyError(" norm_mode is not recognized ")

    ax0 = axes[1]
    # plot colormesh
    im = ax0.pcolor(ctheta, mom, data2d, cmap=cmap, norm=norm, shading='auto')
    # ax0.axhline(y=1, linestyle='--', color='gray')

    # adjust the bottom subplot
    ax0.set_ylim(ymin, ymax)
    ax0.set_xlim(xmin, xmax)
    ax0.set_xscale(xscale)
    ax0.set_yscale(yscale)
    ax0.set_ylabel(r"$\Gamma\beta$", fontsize=fontsize)
    ax0.set_xlabel(r"Polar angle", fontsize=fontsize)
    ax0.get_yaxis().set_label_coords(-0.15, 0.5)
    # ax0.text(0.75, 0.88, title, color='white', transform=ax0.transAxes, fontsize=fontsize)
    ax0.minorticks_on()
    ax0.tick_params(axis='both', which='both', labelleft=True,
                    labelright=False, tick1On=True, tick2On=True,
                    labelsize=12,
                    direction='in',
                    bottom=True, top=True, left=True, right=True)

    cbar = plt.colorbar(im, ax=ax0, norm=norm)
    # cbar.set_label(r"$M_{\rm ej}$ [M$_{\odot}$]", fontsize=fontsize)
    cbar.set_label(cbar_label, fontsize=fontsize)
    cbar.ax.tick_params(axis='both', which='both',
                        labelsize=12, direction='in')
    cbar.ax.minorticks_on()

    # save figures
    if not figpath is None: plt.savefig(figpath+".png", dpi=256)
    if not figpath is None: plt.savefig(figpath+".pdf")
    plt.show()
    plt.close()
def plot_ejecta_angular_distribution(name:str):

    # sim_dic = SIMULATIONS["DD2_135_135_res150"]
    sim_dic = df.loc[name]
    text_dict = df_text.loc[name]
    text = float(sim_dic['tmerg'])+float(text_dict['text'])
    ej = PBA.id_kenta.EjectaData(get_ej_data(sim_dic["name"]),verbose=True)
    # ej_mass = ej.get(v_n="mass",text=text)
    # # ej_mass = ej.get(v_n="ye",text=text)
    # theta = ej.get_theta()
    # mom = PBA.MomFromBeta(ej.get_vinf())

    mom, ctheta, ek, masses = get_2d_ek(ej=ej,text=text,new_vinf_len=80)

    plot_init_profile(mom=mom, ctheta=ctheta, data2d=ek, yscale='log', title=sim_dic["label"],
                      vmin=1e41, vmax=1e48, norm_mode='log', subplot_mode='sum',
                      hist_ylabel=r"$E_{\rm k}$ [erg]",
                      hist_ylabel2=r"$E_{\rm k}(\Gamma\beta>1)$ [erg]",
                      cbar_label=r"$E_{\rm k}$ [erg]",
                      # vmin=0.,vmax=0.6,norm_mode='linear', subplot_mode='ave',
                      #   hist_ylabel=r"$\langle Y_{e;\,\rm ej} \rangle$", cbar_label=r"$Y_{e;\,\rm ej}$",
                      figpath=os.getcwd()+"/figs/"+f"{name}_angular_profile")
    # with h5py.File(os.getcwd()+"/"+"SFHo_13_14_res150_text20_pm.h5",'w') as f:
    #     f.create_dataset(name="Gamma_beta",data=mom)
    #     f.create_dataset(name="polar_angle",data=theta)
    #     f.create_dataset(name="mass",data=ej_mass) # [imom,itheta]
    # print(mom.shape,theta.shape,ej_mass.shape)

    # sim_dic = SIMULATIONS["SFHo_13_14_res150"]
    # text = 13 + sim_dic["tmerg"]
    # ej = PBA.id_kenta.EjectaData(sim_dic["datadir"]+"ej_collated.h5",verbose=True)
    # ej_mass = ej.get(v_n="mass",text=text)
    # theta = ej.get_theta()
    # mom = ej.get_vinf() * PBA.utils.get_Gamma(ej.get_vinf())
    # plot_init_profile(mom=mom,ctheta=theta,mass=ej_mass.T, yscale='log',title=sim_dic["label"],
    #                   figpath=figs+sim_dic["name"]+"_initial")
    #
    # sim_dic = SIMULATIONS["DD2_135_135_res150_floor"]
    # text = 13 + sim_dic["tmerg"]
    # ej = PBA.id_kenta.EjectaData(sim_dic["datadir"]+"ej_collated.h5",verbose=True)
    # ej_mass = ej.get(v_n="mass",text=text)
    # theta = ej.get_theta()
    # mom = ej.get_vinf() * PBA.utils.get_Gamma(ej.get_vinf())
    # plot_init_profile(mom=mom,ctheta=theta,mass=ej_mass.T, yscale='log',title=sim_dic["label"],
    #                   figpath=figs+sim_dic["name"]+"_initial")


if __name__ == '__main__':

    # process new simulation data -> collated.h5
    # task_process()
    # plot_sim_ekecta_mass_text_dependency(name="SFHo_135_135_res150_new")
    # plot_sim_ekecta_mass_text_dependency_all()
    # plot_sim_ekecta_mass_old(sim_name="SFHo_13_14_res150")

    ''' ejecta properties at a given extraction time '''
    # plot_ej_mom_theta_dist_for_text()
    # plot_all_sim_ejecta_mass_evol(crit=None,yscale="linear",ylim=(0., 0.01),figname="ejecta_mass_evol",
    #                          title="Volume integrated ejecta mass")

    # plot_all_sim_ejecta_mass(o_fit=Fit1D_4seg(),figname="ej_mom_ek_nr_",figname_coeffs="ej_mom_ek_coeffs_nr_",
    #                          ylim0=(1e43, 1e51),ylim1=(1e43, 1e51),ylim2=(-0.5,0.5), xlim=(1e-2, 4))
    # plot_all_sim_ejecta_mass(o_fit=Fit1D_3seg(),figname="ej_mom_ek_nr_",figname_coeffs="ej_mom_ek_coeffs_nr_",
    #                          ylim0=(1e43, 1e51),ylim1=(1e43, 1e51),ylim2=(-0.5,0.5), xlim=(1e-2, 4))
    # plot_all_sim_ejecta_mass_row(figname="ej_mom_ek_nr_",figname_coeffs="ej_mom_ek_coeffs_nr_",
    #                              ylim0=(2e43, 1e51),ylim1=(2e43, 4e49),ylim2=(-0.75,0.75), xlim=(2e-2, 4))

    plot_sim_ekecta_mass(figname="angular_dist",o_fit = Fit1D_3seg(),sim_=None,
                         plot_fit_coeffs=True,plot_box_plots=True,plot_fit_ek=True,save_fit_result=True)
    # plot_sim_ekecta_mass(o_fit = Fit1D_4seg(),sim_=None,
    #                      plot_fit_coeffs=True,plot_box_plots=True,plot_fit_ek=True,save_fit_result=True)

    # plot_ejecta_angular_distribution(name="SFHo_13_14_res150")
    # plot_sim_ekecta_mass(sim_dic=df.loc["DD2_135_135_res150"])

    # plot_rho_max(df=df)

    # for (sim_name, _) in df.iterrows():
    #     plot_ejecta_angular_distribution(name=sim_name)
    # plot_ejecta_angular_distribution_all()


    ''' Mass flux & rho_max as a function of time '''
    # plot_rho_mdot_rext()