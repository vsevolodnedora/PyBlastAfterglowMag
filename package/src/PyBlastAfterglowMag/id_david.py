'''

    This script reads output files from postprocessing outflowed.cc code that
    in turn postprocesses the output of the WhiskyTHC NR GRHD code.
    Input files can be: * hist_vel_inf.dat or "corr_vel_inf_theta.h5" *
    These files contain one- or two-dimensional histograms of the ejecta
    velocity distribution.
    The script removes one half of the sphere from the angular distribution,
    and parts of it where mass is 0 or velocity is negative, e.g., it clears the
    data so PyBlastAfterglow does not need to deal with unphysical initial data.

    Usage
    python3 ./id_maker_from_thc_outflow.py -i path/to/hist_or_corr.h5 -o path/to/output.h5 -m corr_or_hist -l 30 --factor 2.

'''
import matplotlib.pyplot as plt
import numpy as np
import h5py
from scipy import interpolate
import copy
import os
import sys
import argparse

# import package.src.PyBlastAfterglowMag.rebin
# from .rebin import rebin


from .id_tools import (reinterpolate_hist, compute_ek_corr)
from .utils import (cgs, MomFromBeta,GammaFromMom,BetaFromMom,BetaFromGamma,MomFromGamma,GammaFromBeta)


def load_vinf_hist(hist_fpath):
    hist1 = np.loadtxt(hist_fpath)
    vinf, mass = hist1[:, 0], hist1[:, 1]
    return (vinf, mass)



def load_corr_file2(corr_fpath,
                    reinterpolate_theta=True,
                    new_theta_len=None,
                    beta_min=0., beta_max=0.,
                    dist="pw",
                    ):
    # load the corr file
    if not os.path.isfile(corr_fpath):
        raise FileNotFoundError(f"file not found: {corr_fpath}")
    dfile = h5py.File(corr_fpath, mode="r")
    for key in dfile.keys():
        print(key, np.array(dfile[key]).shape)

    # extract arrays of data
    theta_edges = np.array(dfile["theta"], dtype=np.double)
    vinf_edges = np.array(dfile["vel_inf"], dtype=np.double)
    mass = np.array(dfile["mass"], dtype=np.double)  # [i_theta, i_vinf]
    mtot = np.sum(mass)

    print("Data: mass={} theta_edges={} vinf_edges={}".format(mass.shape, theta_edges.shape, vinf_edges.shape))
    if (len(theta_edges) - 1 == len(mass[0, :])) and (len(vinf_edges) - 1 == len(mass[:, 0])):
        print("Transposing data")
        mass = np.array(dfile["mass"], dtype=np.double).T

    # check data
    if (len(theta_edges) - 1 == len(mass[:, 0])) and (len(vinf_edges) - 1 == len(mass[0, :])):
        pass
    elif (len(theta_edges) == len(mass[:, 0])) and (len(vinf_edges) == len(mass[0, :])):
        print("WARNING in datafile. Expected mass[n_theta-1, n_vinf-1]"
              " Got: mass[{}, {}] with ntheta-1={} and n_vinf-1={}"
              " Averaging Mass to use grid as edges "
                      .format(len(mass[:, 0]), len(mass[0, :]), len(theta_edges) - 1, len(vinf_edges) - 1))

        mass = 0.5 * (mass[:,1:] + mass[:,:-1])
        mass = 0.5 * (mass[1:,:] + mass[:-1,:])
    else:
        raise IOError("Error in datafile. Expected mass[n_theta-1, n_vinf-1]"
                      " Got: mass[{}, {}] with ntheta-1={} and n_vinf-1={}"
                      .format(len(mass[:, 0]), len(mass[0, :]), len(theta_edges) - 1, len(vinf_edges) - 1))

    theta_centers = 0.5 * (theta_edges[1:] + theta_edges[:-1])
    velocity_centers = 0.5 * (vinf_edges[1:] + vinf_edges[:-1])

    # check again
    if not len(mass[:, 0]) == len(theta_centers):
        raise ValueError("Histogram mismatch :: mass[:, 0]={} theta={}".format(len(mass[:, 0]), len(theta_edges)))
    if not len(mass[0, :]) == len(velocity_centers):
        raise ValueError("Histogram mismatch :: mass[0, :]={} vinf={}".format(len(mass[0, :]), len(vinf_edges)))

    print("Initial data theta edges: [{}, {}] * pi".format(theta_edges[0] / np.pi, theta_edges[-1] / np.pi))
    print("Initial data vinf edges : [{}, {}] c".format(vinf_edges.min(), vinf_edges.max()))
    print("total mass = {} [before interpolation from {} ]".format(mtot, mass.shape))

    # remove the low beta data if needed
    if (beta_min > 0):
        print("Removing data with beta < {}".format(beta_min))
        mass = mass[:, ~(velocity_centers < beta_min)]
        vinf_edges = vinf_edges[~(vinf_edges < beta_min)]
    if (beta_max > 0):
        print("Removing data with beta > {}".format(beta_max))
        mass = mass[:, ~(velocity_centers > beta_max)]
        vinf_edges = vinf_edges[~(vinf_edges > beta_max)]

    # the data is usually stored for both hemispheres (but with 'z' symmetry it is the same)
    if (theta_edges[-1] > 3. * np.pi / 4.):
        print("Removing data with theta > pi/2.")
        nmax = len(theta_edges[theta_edges < np.pi / 2.]) + 1  # ...

        theta_one_hemisphere = theta_edges[: nmax]
        mass = mass[: nmax - 1]

        thetas_pol_edges = theta_one_hemisphere

    else:
        raise IOError("Error! Expected theta to go from ~0 to ~pi But the end is {}".format(theta_edges[-1]))

    if (reinterpolate_theta):
        print("Rebinning historgram")
        thetas_pol_edges, mass = reinterpolate_hist(thetas_pol_edges, mass, new_theta_len, dist=dist)

    thetas_pol_centers = 0.5 * (thetas_pol_edges[1:] + thetas_pol_edges[:-1])
    vinf_centers = 0.5 * (vinf_edges[1:] + vinf_edges[:-1])

    assert thetas_pol_centers[0] < thetas_pol_edges[-1]
    assert thetas_pol_centers[-1] <= np.pi / 2.

    # assert len(thetas_pol_edges)-1 == len(mass[:, 0])
    # assert len(vinf_centers) == len(mass[0, :])

    return (thetas_pol_edges, thetas_pol_centers, vinf_centers, mass)

def clean_data_hist(_vinf, _mass):
    tmp_vinf = []
    tmp_mass = []
    for i in range(len(_vinf)):
        if ((_vinf[i] > 0.) and (_mass[i] > 0.) and (_vinf[i] < 1.)):
            tmp_vinf.append(_vinf[i])
            tmp_mass.append(_mass[i])
    if (len(tmp_mass) < 2):
        raise ValueError("Failed clean the data")
    return (np.array(tmp_vinf), np.array(tmp_mass))
def clean_data_corr(thetas_pol, beta, mass, remove_pi_over_2=False):
    assert len(thetas_pol) == len(mass[:, 0])
    assert len(beta) == len(mass[0, :])
    tmp_beta, tmp_mass = [], []
    # remove all data where there is no mass
    for i in range(len(beta)):
        if (beta[i] <= 0.):
            continue
        if (np.sum(mass[:, i]) <= 0.):
            continue
        tmp_beta.append(beta[i])
        tmp_mass.append(mass[:, i])
    tmp_beta = np.array(tmp_beta)
    tmp_mass = np.reshape(np.array(tmp_mass), newshape=(len(tmp_beta), len(thetas_pol))).T
    tmp_theta_pol, tmp_tmp_mass = [], []
    # remove all data where there is no mass
    for i in range(len(thetas_pol)):
        if (thetas_pol[i] < 0):
            continue
        if (np.sum(tmp_mass[i, :]) <= 0.):
            continue
        tmp_theta_pol.append(thetas_pol[i])
        tmp_tmp_mass.append(tmp_mass[i, :])
    tmp_theta_pol = np.array(tmp_theta_pol)
    tmp_tmp_mass = np.reshape(np.array(tmp_tmp_mass), newshape=(len(tmp_theta_pol), len(tmp_beta)))
    # remove data at pi/4, replace it with interpolation. Reason -- WiskyTHC has some weird spike at pi/4.
    if (remove_pi_over_2):
        for i in range(len(tmp_beta)):
            idx = abs((tmp_theta_pol - thetas_pol[0]) / (np.pi / 4.) - 1.) < 1e-2
            tmp_tmp_mass[idx, i] = interpolate.interp1d(
                np.delete(tmp_theta_pol, idx),
                np.delete(tmp_tmp_mass[:, i], idx), kind="linear")(tmp_theta_pol[idx])
    return (tmp_theta_pol, tmp_beta, tmp_tmp_mass)
def compute_ek_hist(_vinf, _mass):
    return _mass * cgs.solar_m * (_vinf * _vinf * cgs.c * cgs.c)


def prepare_kn_ej_id_2d(nlayers, corr_fpath, outfpath,
                        r0type="fromrho", t0=100, r0frac=0.5,
                        dist="pw"):
    thetas, cthetas, betas, masses = load_corr_file2(
        corr_fpath=corr_fpath,
        reinterpolate_theta=True,
        new_theta_len=nlayers if ((not nlayers is None) and (nlayers > 0)) else None,
        dist=dist)

    # theta_corr, vinf_corr, mass_corr = clean_data_corr(thetas, betas, masses, remove_pi_over_2=True)
    ctheta_corr, vinf_corr, mass_corr = clean_data_corr(cthetas, betas, masses, remove_pi_over_2=True)
    ek_corr = compute_ek_corr(vinf_corr, mass_corr).T  # [n_beta, n_theta]
    print(ctheta_corr.shape, vinf_corr.shape, ek_corr.shape)
    # EjectaEk.plot_corr(theta_corr, vinf_corr, ek_corr)
    ek_corr2 = np.copy(ek_corr[(vinf_corr > 0.) & (vinf_corr < 1), :])
    mass_corr2 = np.copy(mass_corr.T[(vinf_corr > 0.) & (vinf_corr < 1), :])
    vinf_corr2 = np.copy(vinf_corr[(vinf_corr > 0.) & (vinf_corr < 1)])
    ctheta_corr2 = np.copy(ctheta_corr)
    # EjectaEk.plot_corr(theta_corr, vinf_corr, ek_corr)
    print(ctheta_corr2.shape, vinf_corr2.shape, ek_corr2.shape)

    # self.o_pba.setEjectaStructNumeric(theta_corr2, vinf_corr2, ek_corr2, fac, True, self.pars_kn["eats_method"])

    # vinf_corr2 = vinf_corr2[25:]
    # ek_corr2 = ek_corr2[25:,:]
    # mass_corr2 = mass_corr2[25:,:]

    ctheta_corr3 = np.zeros_like(ek_corr2)
    theta_corr3 = np.zeros_like(ek_corr2)
    for imom in range(len(vinf_corr2)):
        ctheta_corr3[imom,:] = ctheta_corr2
        theta_corr3[imom,:] = ctheta_corr2
    mom_corr3 = np.zeros_like(ek_corr2)
    for ith in range(len(ctheta_corr2)):
        mom_corr3[:,ith]=np.array( vinf_corr2*GammaFromBeta(vinf_corr2))

    r = np.zeros_like(mass_corr2)
    t = t0
    for ith in range(len(ctheta_corr3[0,:])):
        for ir in range(len(mom_corr3[:,0])):
            r[ir,ith] =  BetaFromMom(mom_corr3[ir,ith])*cgs.c * t

    # print(len(ctheta_corr2), ctheta_corr2)

    def plot_final_profile(ctheta, mom, data2d, cmap='jet',vmin=None,vmax=None):
        from matplotlib.colors import LogNorm
        if vmin is None: vmin = data2d[(data2d > 0) & (np.isfinite(data2d))].min()
        if vmax is None: vmax = data2d[(data2d > 0) & (np.isfinite(data2d))].max()
        norm = LogNorm(vmin, vmax)
        fig,ax = plt.subplots(ncols=1,nrows=1, figsize=(4.6,3.2))
        ctheta *= (180. / np.pi)
        im = ax.pcolor(mom, ctheta, data2d, cmap=cmap, norm=norm, shading='auto')
        # ax0.axhline(y=1, linestyle='--', color='gray')

        # adjust the bottom subplot
        ax.set_ylim(0, 90)
        ax.set_xlim(1e-2, 4)
        ax.set_xscale('log')
        ax.set_yscale('linear')
        ax.set_xlabel(r"$\Gamma\beta$", fontsize=12)
        ax.set_ylabel(r"Polar angle", fontsize=12)
        plt.show()

    # plot_final_profile(ctheta=ctheta_corr3,mom=mom_corr3,data2d=ek_corr2)
    # plot_final_profile(ctheta=ctheta_corr3,mom=mom_corr3,data2d=mass_corr2*cgs.solar_m)

    dfile = h5py.File(outfpath, "w")
    dfile.attrs.create("theta_wind", data=np.pi/2, dtype=np.float64)
    dfile.attrs.create("theta_core", data=np.pi/2, dtype=np.float64)
    dfile.create_dataset("r",data=r)# used for theta_w in PW method
    dfile.create_dataset("theta",data=ctheta_corr3)# used for theta_w in PW method
    dfile.create_dataset("ctheta",data=ctheta_corr3)
    dfile.create_dataset("mom",data=mom_corr3)
    dfile.create_dataset("ek",data=ek_corr2)
    dfile.create_dataset("mass",data=mass_corr2*cgs.solar_m)
    dfile.create_dataset("ye",data=np.zeros_like(ek_corr2))
    dfile.create_dataset("s",data=np.zeros_like(ek_corr2))
    dfile.close()
    print("file saved: {}".format(outfpath))
def prepare_kn_ej_id_1d(nlayers, hist_fpath, outfpath):
    betas, masses = load_vinf_hist(hist_fpath=hist_fpath)
    # thetas = np.full_like(betas, np.pi / 2.)

    # vinf_hist, mass_hist = o_data.load_vinf_hist()
    assert len(betas) == len(masses)
    vinf_hist, mass_hist = clean_data_hist(betas, masses)
    ek_hist = compute_ek_hist(vinf_hist, mass_hist)
    # ek_hist = np.cumsum(ek_hist)
    # vinf_hist = np.sqrt(ek_hist/mass_hist/cgs.solar_m)/cgs.c
    # theta_hist = np.zeros_like(ek_hist)
    theta_hist = np.full_like(vinf_hist, np.pi / 2.)
    gam_hist = GammaFromBeta(vinf_hist)
    if (len(gam_hist[~np.isfinite(gam_hist)]) > 0):
        raise ValueError("nan in gammas for beta={}".format(vinf_hist))

    # self.o_pba.setEjectaStructNumericUniformInTheta(theta_hist, ek_hist, gam_hist,
    #                                                 mass_hist * cgs.solar_m, nlayers, fac,
    #                                                 self.pars_kn["eats_method"])

    dfile = h5py.File(outfpath, "w")
    dfile.create_dataset("theta", data=theta_hist)
    dfile.create_dataset("mom", data=gam_hist*vinf_hist)
    dfile.create_dataset("ek", data=ek_hist)
    dfile.close()
    print("file saved: {}".format(outfpath))

# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type=str,
                        help="fpath to histogram or 2d histogram of the ejecta velocity distripution")
    parser.add_argument("-m", "--mode", type=str, choices=["hist", "corr"],
                        help="what data is given: histogram or a 2d histogram")
    parser.add_argument("-m", "--output", type=str,
                        help="fpath of the output file")
    parser.add_argument("-l", "--nlayers", type=int,
                        help="number of angular layers to make", default=30)
    parser.add_argument("--factor", type=float,
                        help="multiplier for ejecta mass", default=2.)
    args = parser.parse_args()
    mode = args.mode
    infpath = args.input
    outfpath = args.output
    nlayers = args.nlayers
    fac = args.factor
    # ---
    if (mode not in ["hist","corr"]):
        raise KeyError(" option '--mode' can be 'hist' or 'corr' only. Given:{}".format(mode))
    if (not os.path.isfile(infpath)):
        raise FileNotFoundError(" input file is not found {} ".format(infpath))
    if (os.path.isfile(outfpath)):
        raise FileExistsError(" output file already exists {}".format(outfpath))
    # ---
    if (mode == "hist"):
        prepare_kn_ej_id_1d(nlayers, hist_fpath=infpath, outfpath=outfpath)
    elif (mode == "corr"):
        prepare_kn_ej_id_2d(nlayers=nlayers, corr_fpath=infpath, outfpath=outfpath, dist="a")
    else:
        exit(1)
    # ---
    print("Ejecta initial data is prepared {}".format(outfpath))

