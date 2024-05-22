import numpy as np
import h5py
import os
from glob import glob

import matplotlib.pyplot as plt
from matplotlib import ticker, cm, rc, rcParams
from matplotlib.colors import Normalize, LogNorm, BoundaryNorm
from matplotlib.ticker import MaxNLocator

rc('text', usetex=True) # $ sudo apt-get install cm-super
rc('font', family='serif')
rcParams['font.size'] = 10

import package.src.PyBlastAfterglowMag as PBA

def ejecta_analysis():
    dir = "/media/vsevolod/T7/work/KentaData/"
    simname = "SFHoTim276_12_15_0025_150mstg_B0_HLLC"+"/"

    files = glob(dir+simname+"ejecta_*.h5")
    if (len(files) == 0):
        raise FileNotFoundError("Files not found")
    id = PBA.id_kenta.ProcessRawFiles(files=files,verbose=True)

    outfnmae = dir+simname+"ej_collated.h5"
    id.process_save(outfnmae)

    ej = PBA.id_kenta.EjectaData(fpath=outfnmae,verbose=True)

    id = PBA.id_kenta.EjStruct(fpath=outfnmae,verbose=True)
    id_dict = id.get_2D_id(text=25,method_r0="from_beta",t0=1e3,new_theta_len=None,new_vinf_len=None)
    mom = id_dict["mom"]
    ctheta = id_dict["ctheta"]
    ek = id_dict["ek"]
    mass = id_dict["mass"]
    id.plot_init_profile(mom=id_dict["mom"][:,0], ctheta=id_dict["ctheta"][0,:], mass=id_dict["ek"])

    outfnmae = dir+simname+"ej_id.h5"
    with h5py.File(outfnmae, "w") as dfile:
        for key, data in id_dict.items():
            dfile.create_dataset(name=key, data=data)

    id.plot_init_profile(mom=ej.vinf[ej.vinf<1]*PBA.utils.get_Gamma(ej.vinf[ej.vinf<1]), ctheta=ej.theta, mass=ej.get(v_n="mass",text=25)[:,ej.vinf<1].T)


    data = PBA.id_maker_from_kenta_bns.Data(fpath_rhomax=dir+simname+"rhomax_SFHo_12_15.txt",
                                            fpath_mdot=dir+simname+"Mdot_extraction_SFHo_12_15.txt")

    thetas, betas, masses = ej.theta, ej.vinf, ej.get(v_n="mass",text=25)


    xmin = 0
    xmax = 90
    ymin = 1e-2
    ymax = 6
    vmin = 1e-12
    vmax = 1e-6

    fontsize = 12

    gammas = PBA.utils.get_Gamma(betas)
    moms = gammas * betas
    thetas = thetas * 180 / PBA.utils.cgs.pi

    mask = moms > 1
    moms = moms[mask]
    masses = masses[:,mask]

    fig, ax = plt.subplots(figsize=(4.6, 3.0), ncols=1, nrows=1, constrained_layout=True)

    levels = MaxNLocator(nbins=15).tick_values(masses.min(), masses.max())
    cmap = plt.get_cmap('RdYlBu_r')
    norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
    norm = LogNorm(masses[(masses > 0) & (np.isfinite(masses))].min(), masses[(masses > 0) & (np.isfinite(masses))].max())

    im = ax.pcolor(thetas, moms, masses.T, cmap=cmap, norm=norm, shading='auto')
    # cbar = fig.colorbar(im, ax=ax, extend='both')
    # cbar.ax.set_title(r"Mass [$M_{\odot}$]")

    ax.axhline(y=1, linestyle='--', color='gray')

    ax.set_ylim(ymin, ymax)
    ax.set_xlim(xmin, xmax)

    ax.set_ylabel(r"$\Gamma\beta$", fontsize=fontsize)
    ax.set_xlabel(r"Polar angle [deg]", fontsize=fontsize)
    # ax.get_yaxis().set_label_coords(-0.15, 0.5)
    # ax0.text(0.75, 0.88, lbl, color='white', transform=ax0.transAxes, fontsize=fontsize)
    ax.minorticks_on()
    ax.tick_params(axis='both', which='both', labelleft=True,
                   labelright = False, tick1On = True, tick2On = True,
                   labelsize = fontsize-2,
                   direction = 'in',
                   bottom = True, top = True, left = True, right = True)

    # ax.text(0.05, 0.05, task["line"]["label"], color='black', transform=ax.transAxes)
    ax.text(0.05, 0.85, "Sim", color='black', transform=ax.transAxes)


    # ekdens = np.sum(eks, axis=0)
    # ax1.step(0.5 * (theta[1:] + theta[:-1]), mdens, where='mid', color='red')
    # hist_vinf, hist_mass = o_data.load_vinf_hist()
    # hist_mom = hist_vinf * get_Gamma(hist_vinf)
    # hist_eks = 0.5 * (hist_vinf * cgs.c) ** 2 * hist_mass * cgs.solar_m
    # ax1.step(mom, ekdens, where='mid', color='red')
    # ax1.step(hist_mom, hist_eks, where='mid', color='black')

    # cbar = plt.colorbar(im, cax=cax, norm=norm)
    cbar = fig.colorbar(im, ax=ax, norm=norm)
    cbar.ax.set_title(r"$M_{\rm ej}$ [M$_{\odot}$]", fontsize=fontsize)
    cbar.ax.minorticks_off()
    # plt.savefig(PAPERPATH + "kinetic_energy_struct_models.pdf")
    # plt.savefig(FIGPATH + "kinetic_energy_struct_models.png")
    # if save_figs: plt.savefig(FIGPATH + figname + ".png", dpi=256)
    # if save_figs and save_pdfs: plt.savefig(PAPERPATH + figname + ".pdf")
    # plt.savefig(sys.argv[0].replace(".py", "_") + name.lower() + ".pdf")
    # if show_figs:
    plt.show()
    plt.close()

    # ------------------------------------------------------------------------------------------------------------------

    times = ej.texts
    masses = np.column_stack([np.sum(ej.get(v_n="mass",text=t)[:,mask],axis=1) for t in ej.texts])

    fig, ax = plt.subplots(figsize=(4.6, 3.0), ncols=1, nrows=1, constrained_layout=True)

    levels = MaxNLocator(nbins=15).tick_values(masses.min(), masses.max())
    cmap = plt.get_cmap('RdYlBu_r')
    norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
    norm = LogNorm(masses[(masses > 0) & (np.isfinite(masses))].max()*1e-3, masses[(masses > 0) & (np.isfinite(masses))].max())

    im = ax.pcolor(times, thetas, masses, cmap=cmap, norm=norm, shading='auto')
    # cbar = fig.colorbar(im, ax=ax, extend='both')
    # cbar.ax.set_title(r"Mass [$M_{\odot}$]")

    ax.axhline(y=1, linestyle='--', color='gray')

    ax.set_ylim(xmin, xmax)
    # ax.set_xlim(times.min(), times.max())

    ax.set_ylabel(r"Polar angle [deg]", fontsize=fontsize)
    ax.set_xlabel(r"Extraction time [ms]", fontsize=fontsize)
    # ax.get_yaxis().set_label_coords(-0.15, 0.5)
    # ax0.text(0.75, 0.88, lbl, color='white', transform=ax0.transAxes, fontsize=fontsize)
    ax.minorticks_on()
    ax.tick_params(axis='both', which='both', labelleft=True,
                   labelright = False, tick1On = True, tick2On = True,
                   labelsize = fontsize-2,
                   direction = 'in',
                   bottom = True, top = True, left = True, right = True)

    # ax.text(0.05, 0.05, task["line"]["label"], color='black', transform=ax.transAxes)
    ax.text(0.05, 0.85, "Sim", color='black', transform=ax.transAxes)


    # ekdens = np.sum(eks, axis=0)
    # ax1.step(0.5 * (theta[1:] + theta[:-1]), mdens, where='mid', color='red')
    # hist_vinf, hist_mass = o_data.load_vinf_hist()
    # hist_mom = hist_vinf * get_Gamma(hist_vinf)
    # hist_eks = 0.5 * (hist_vinf * cgs.c) ** 2 * hist_mass * cgs.solar_m
    # ax1.step(mom, ekdens, where='mid', color='red')
    # ax1.step(hist_mom, hist_eks, where='mid', color='black')

    # cbar = plt.colorbar(im, cax=cax, norm=norm)
    cbar = fig.colorbar(im, ax=ax, norm=norm)
    cbar.ax.set_title(r"$M_{\rm ej}$ [M$_{\odot}$]", fontsize=fontsize)
    cbar.ax.minorticks_off()
    # plt.savefig(PAPERPATH + "kinetic_energy_struct_models.pdf")
    # plt.savefig(FIGPATH + "kinetic_energy_struct_models.png")
    # if save_figs: plt.savefig(FIGPATH + figname + ".png", dpi=256)
    # if save_figs and save_pdfs: plt.savefig(PAPERPATH + figname + ".pdf")
    # plt.savefig(sys.argv[0].replace(".py", "_") + name.lower() + ".pdf")
    # if show_figs:
    plt.show()
    plt.close()





    # plot ej.mass evolution
    fig, ax = plt.subplots(ncols=1,nrows=1)
    ax.plot(ej.texts, ej.total_mass(), color="gray",label="Total")
    ax.plot(ej.texts, ej.total_mass_fasttail(),color="black",label=r"$\Gamma\beta>1$")
    ax2 = ax.twinx()
    ax2.plot(*data.get_rhomax(),color="green",label=r"$\rho_{\rm max}$")
    ax2.plot(*data.get_mdot(),color="green",label=r"$\rho_{\rm max}$")
    ax2.set_yscale("log")
    plt.legend()
    ax.set_yscale("log")
    plt.show()



def main():
    ejecta_analysis()

if __name__ == '__main__':
    main()