"""
    Uses PyBlastAfterglowMag code
    From iside `\package\' run:
    pip uninstall --no-cache-dir PyBlastAfterglowMag & pip install .
    to install the postprocessing/run utilities

"""
import copy

import numpy as np
import h5py
import shutil
from glob import glob
from multiprocessing import Pool
from scipy.integrate import ode
import matplotlib.pyplot as plt
from scipy import interpolate
# import pytest
from pathlib import Path
from matplotlib.colors import Normalize, LogNorm
from pathlib import Path
from matplotlib.colors import Normalize, LogNorm
from mpl_toolkits.axes_grid1 import ImageGrid
from mpl_toolkits.axisartist.grid_finder import MaxNLocator
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import Normalize, LogNorm
from matplotlib import ticker, cm, rc, rcParams
from matplotlib.patches import Rectangle
import matplotlib.colors as colors
from matplotlib import cm
from matplotlib import cm
import os
import re
from scipy.optimize import curve_fit
from scipy.optimize import leastsq

import PyBlastAfterglowMag as PBA

def model(params, v, theta):
    a, v0, sigma_v, theta0, sigma_theta, b, alpha = params
    gaussian_part = a * np.exp(-(v - v0)**2 / (2 * sigma_v**2)) * np.exp(-(theta - theta0)**2 / (2 * sigma_theta**2))
    exponential_part = b * np.exp(-alpha * v)
    return gaussian_part + exponential_part

def error_function(params, v, theta, mass):
    prediction = model(params, v, theta)
    residuals = prediction - mass
    penalized_residuals = residuals * (1 + v)
    return penalized_residuals.flatten()

def _1():
    dir = "/media/vsevolod/T7/work/KentaData/"
    simname = "SFHoTim276_12_15_0025_150mstg_B0_HLLC" + "/"
    files = glob(dir + simname + "ejecta_*.h5")
    if (len(files) == 0):
        raise FileNotFoundError("Files not found")
    id = PBA.id_kenta.ProcessRawFiles(files=files, verbose=True)
    collated_ej_data_fpath = dir + simname + "ej_collated.h5"
    # id.process_save(collated_ej_data_fpath)

    # extract data from collated file for a given extraction time and save it into working dir.
    text = 25
    id = PBA.id_kenta.EjStruct(fpath=collated_ej_data_fpath, verbose=True)
    id_dict = id.get_2D_id(text=text, method_r0="from_beta", t0=1e3, new_theta_len=None, new_vinf_len=None)


    V, Theta = np.meshgrid(id_dict["mom"], id_dict["ctheta"])
    mass = np.log10(id_dict["mass"].T)

    # Initial guesses for parameters
    params0 = [1e-5, 0.5, 0.1, 0.8, 0.1, 1e-5, 1.0]
    result, _ = leastsq(error_function, params0, args=(V, Theta, mass))

    fig, axes = plt.subplots(3, 1, figsize=(6, 12))

    # Plotting the original data
    cax1 = axes[0].pcolormesh(V, Theta, mass, shading='auto')
    axes[0].set_title('Original Data')
    axes[0].set_ylabel('Theta')
    fig.colorbar(cax1, ax=axes[0])

    # Plotting the predicted data
    cax2 = axes[1].pcolormesh(V, Theta, model(result, V, Theta), shading='auto')
    axes[1].set_title('Predicted Data')
    axes[1].set_ylabel('Theta')
    fig.colorbar(cax2, ax=axes[1])

    # Plotting the residuals
    residuals = mass - model(result, V, Theta)
    cax3 = axes[2].pcolormesh(V, Theta, residuals, shading='auto')
    axes[2].set_title('Residuals')
    axes[2].set_xlabel('Velocity')
    axes[2].set_ylabel('Theta')
    fig.colorbar(cax3, ax=axes[2])

    plt.tight_layout()
    plt.show()


def _2(V, Theta, mass_2D):
    def gaussian_2d(vtheta, a, v0, sigma_v, theta0, sigma_theta):
        v, theta = vtheta
        return a * np.exp(-(v - v0) ** 2 / (2 * sigma_v ** 2)) * np.exp(-(theta - theta0) ** 2 / (2 * sigma_theta ** 2))

    # Set initial parameter guesses
    params0 = [1, 0.1, 1, 0.0, 1]

    # Flattening the 2D data arrays for curve_fit
    V_flat = V.ravel()
    Theta_flat = Theta.ravel()
    mass_flat = mass_2D.ravel()

    # Performing the fit
    V_flat = V_flat[np.isfinite(mass_flat)]
    Theta_flat = Theta_flat[np.isfinite(mass_flat)]
    mass_flat = mass_flat[np.isfinite(mass_flat)]
    # params, covariance = curve_fit(gaussian_2d, (V_flat, Theta_flat), mass_flat, p0=params0)

    weights = V_flat**2
    params, covariance = curve_fit(gaussian_2d, (V_flat, Theta_flat), mass_flat, p0=params0,
                                   sigma=1.0/weights, absolute_sigma=True,maxfev=int(1e5))

    # Predicting using the fitted parameters
    predicted_data = gaussian_2d((V, Theta), *params).reshape(V.shape)

    residuals = mass_2D - predicted_data
    residuals[~np.isfinite(residuals)] = 0.
    ss_res = np.sum(residuals ** 2)
    ss_tot = np.sum((mass_flat - np.mean(mass_flat)) ** 2)
    r_squared = 1 - (ss_res / ss_tot)

    print(f"R-squared value: {r_squared}")
    print(f"Total_massbefore: {np.sum(mass_2D)} predicted: {np.sum(np.exp(predicted_data))}")
    print(f"Total mass fast_tail before: {np.sum(np.exp(mass_2D[V > 1]))} predicted: {np.sum(np.exp(predicted_data[V > 1]))}")

    fig, axes = plt.subplots(3, 1, figsize=(6, 12), sharex="all")

    # Original Data
    cax1 = axes[0].pcolormesh(V, Theta, mass_2D, shading='auto')
    axes[0].set_title('Original Data')
    axes[0].set_ylabel('Theta')
    fig.colorbar(cax1, ax=axes[0])

    # Predicted Data
    cax2 = axes[1].pcolormesh(V, Theta, predicted_data, shading='auto')
    axes[1].set_title('Predicted Data')
    axes[1].set_ylabel('Theta')
    fig.colorbar(cax2, ax=axes[1])

    # Residuals
    cax3 = axes[2].pcolormesh(V, Theta, residuals, shading='auto')
    axes[2].set_title('Residuals')
    axes[2].set_xlabel('Velocity')
    axes[2].set_ylabel('Theta')
    fig.colorbar(cax3, ax=axes[2])

    plt.tight_layout()
    plt.show()

def main():
    dir = "/media/vsevolod/T7/work/KentaData/"
    simname = "SFHoTim276_12_15_0025_150mstg_B0_HLLC" + "/"
    files = glob(dir + simname + "ejecta_*.h5")
    if (len(files) == 0):
        raise FileNotFoundError("Files not found")
    id = PBA.id_kenta.ProcessRawFiles(files=files, verbose=True)
    collated_ej_data_fpath = dir + simname + "ej_collated.h5"
    # id.process_save(collated_ej_data_fpath)

    # extract data from collated file for a given extraction time and save it into working dir.
    text = 25
    id = PBA.id_kenta.EjStruct(fpath=collated_ej_data_fpath, verbose=True)
    id_dict = id.get_2D_id(text=text, method_r0="from_beta", t0=1e3, new_theta_len=None, new_vinf_len=None)

    V, Theta = np.meshgrid(id_dict["mom"], id_dict["ctheta"])
    mass = np.log(id_dict["mass"].T)
    _2(V,Theta,mass)

if __name__ == '__main__':
    main()