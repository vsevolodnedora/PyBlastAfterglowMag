import numpy as np
import matplotlib.pyplot as plt
import h5py
import os

from mpl_toolkits.axes_grid1.inset_locator import (inset_axes, InsetPosition, mark_inset)
from matplotlib.colors import Normalize, LogNorm
from matplotlib import ticker, cm, rc, rcParams
import matplotlib.colors as colors
import matplotlib as mpl


import os
import math
import mpl_toolkits.axisartist.floating_axes as floating_axes
from matplotlib.projections import PolarAxes
from mpl_toolkits.axisartist.grid_finder import FixedLocator, MaxNLocator, DictFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable
import random

from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle

from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator

from matplotlib.patches import Wedge

rc('text', usetex=True) # $ sudo apt-get install cm-super
rc('font', family='serif')
rcParams['font.size'] = 10

curdir = os.getcwd() + '/'

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

def plot_ejecta_dyn_evol_for_movie():

    dfile_ej = h5py.File(curdir+"dyn_kn_gauss.h5")

    print(dfile_ej.keys())
    print(dfile_ej[list(dfile_ej.keys())[0]].keys())
    print(dfile_ej.attrs.keys())


    all_shells = int(dfile_ej.attrs["nshells"])
    all_layers = int(dfile_ej.attrs["nlayers"])
    print("all_shells={} all_layers={}".format(all_shells, all_shells))

    r_grid = np.array(dfile_ej["shell={} layer={}".format(0, 0)]["R"])

    cthetas = np.zeros((all_layers, all_shells))
    rs = np.zeros((all_layers, all_shells))
    val = np.zeros((all_layers, all_shells))

    ir = 100

    for ilayer in range(int(all_layers)):
        for ishell in range(int(all_shells)):
            key = "shell={} layer={}".format(ishell, ilayer)
            cthetas[ilayer, ishell] = np.array(dfile_ej[key]["ctheta"])[ir]
            rs[ilayer, ishell] = np.array(dfile_ej[key]["R"])[ir]
            i_val = np.array(dfile_ej[key]["mom"])
            val[ilayer, ishell] = i_val[ir]

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
    cthetas = np.array(cthetas)  # [::-1] # INVERT DATA ----------------------------------------------------------

    print("ctheta: {}".format(cthetas.shape))
    print("Rs: {}".format(rs.shape))


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

    angle_ticks = range(0, 100, 10)
    angle_ticks_rads = [a * np.pi / 180.0 for a in angle_ticks]  # [::-1] # INVERT TICKS -------------------
    angle_ticks_rads_plus_offset = [a for a in angle_ticks_rads]
    # angle_ticks_rads_plus_offset = angle_ticks_rads_plus_offset[::-1] # Polar Angle
    angle_ticks_for_plot = []
    for i in range(len(angle_ticks)):
        angle_ticks_for_plot.append((angle_ticks_rads_plus_offset[i], r"$" + str(angle_ticks[i]) + "$"))

    print("Angle ticks prepared")

    # lRs = np.log10(Rs/Rs0[0])
    # lRs0 = np.log10(Rs0/Rs0[0])

    # rs = rs/rs.min()
    lrs = np.log10(rs)#- np.log2(rs)[0, 0]
    val[~np.isfinite(lrs)] = np.nan
    lrs[~np.isfinite(lrs)] = 0.
    # lrs /= lrs.min()
    # radius_ticks = range(int(lrs[0, 0]), int(lrs[0, -1]), 1)
    radius_ticks = range(int(lrs[np.isfinite(lrs)].min()), int(lrs[np.isfinite(lrs)].max()), 1)
    radius_ticks_for_plot = []
    for i in range(len(radius_ticks)):
        radius_ticks_for_plot.append((radius_ticks[i], r"$" + str(radius_ticks[i]) + "$"))

    print("Radial ticks prepared")

    # ---------------------------------------

    scale = 1.5
    aspect = 1.20
    height = 3.0
    fig = plt.figure(1, figsize=(height * aspect * scale, height * scale))
    fig.subplots_adjust(wspace=0.1, left=0.05, right=0.95, top=0.94)
    fig.subplots_adjust()

    ax2, aux_ax2 = setup_arc_radial_axes(fig, 111, angle_ticks_for_plot, radius_ticks_for_plot, 1,
                                              radius_ticks[-1])

    # r, theta = np.meshgrid(lRs,cthetas)
    values = val

    # levels = ticker.LogLocator(base=10,numticks=1000).tick_values(1e-4, 0.8)
    levels = MaxNLocator(nbins=20).tick_values(val.min(), val.max())
    cmap = plt.get_cmap('viridis')  # seismic
    norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
    norm = colors.TwoSlopeNorm(vmin=val.min(), vcenter=0.1, vmax=val.max())

    # im = aux_ax2.pcolormesh(cthetas, lrs, values, cmap=cmap, norm=norm, alpha=0.7)
    im = aux_ax2.pcolormesh(cthetas, lrs, values, cmap=cmap, alpha=0.7, norm=LogNorm(vmin=1e-3, vmax=5e-1),
                            rasterized=True)

    cbar = plt.colorbar(im, orientation='vertical')
    cbar.ax.set_ylabel(r'$\log_{10}( \rho/\rho_{\rm ISM} )$', fontsize=12)

    # ax2.axis["bottom"].axes.set_xlabel('Angle [deg]', fontsize=20, weight="bold")

    plt.suptitle(' Jet layer dynamics ', fontsize=14, weight="bold")
    plt.legend(loc=3, prop={'size': 14})
    # plt.xlabel('Angle [deg]', fontsize=20, weight="bold", rotation=30)
    plt.ylabel(r'Radius $[R/R_0]$', fontsize=14, weight="bold")

    # plt.show()
    plt.savefig('plot.png', dpi=256)
    plt.show()
    plt.close()

    fig, ax = plt.subplots(figsize=(4, 3), subplot_kw={'projection': 'polar'}, ncols=1, nrows=1)
    # ax = plt.subplot(111, polar=True)
    # levels = MaxNLocator(nbins=25).tick_values(val.min(), val.max())
    levels = ticker.LogLocator(base=10, numticks=100).tick_values(val[np.isfinite(val)].min(), val[np.isfinite(val)].max())

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


def _load_data(ir = 10, v_n="mom"):
    dfile_ej = h5py.File(curdir+"magnetar_driven_ej.h5")

    # print(dfile_ej.keys())
    # print(dfile_ej[list(dfile_ej.keys())[0]].keys())
    # print(dfile_ej.attrs.keys())

    keys = dfile_ej.keys()
    ej_keys = [key for key in keys if key.__contains__("EJ")]
    print(ej_keys)


    all_shells = int(dfile_ej.attrs["nshells"])
    all_layers = 4;# int(dfile_ej.attrs["nlayers"])
    print("all_shells={} all_layers={}".format(all_shells, all_shells))

    r_grid = np.array(dfile_ej["shell={} layer={} key=R".format(0, 0)])
    ctheta_grid = np.array(dfile_ej["shell=0 layer=0 key=ctheta"])
    print(f"R size = {r_grid.shape} ctheta={ctheta_grid.shape}")

    cthetas = np.zeros((all_layers, all_shells))
    rs = np.zeros((all_layers, all_shells))
    ds = np.zeros((all_layers, all_shells))
    val = np.zeros((all_layers, all_shells))
    t = np.zeros((all_layers, all_shells))




    for ilayer in range(int(all_layers)):
        for ishell in range(int(all_shells)):
            key_t = "shell={} layer={} key={}".format(ishell, ilayer, "tburst")
            key_x = "shell={} layer={} key={}".format(ishell, ilayer, "EJr")
            key_d = "shell={} layer={} key={}".format(ishell, ilayer, "EJdelta")
            key_y = "shell={} layer={} key={}".format(ishell, ilayer, "ctheta")
            key_z = "shell={} layer={} key={}".format(ishell, ilayer, v_n)
            cthetas[ilayer, ishell] = np.array(dfile_ej[key_y])[ir]#np.array(dfile_ej[key]["ctheta"])[ir]
            rs[ilayer, ishell] = np.array(dfile_ej[key_x])[ir]
            ds[ilayer, ishell] = np.array(dfile_ej[key_d])[ir]
            t[ilayer, ishell] = np.array(dfile_ej[key_t])[ir]
            i_val = np.array(dfile_ej[key_z])
            val[ilayer, ishell] = i_val[ir]

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
    cthetas = np.array(cthetas)  # [::-1] # INVERT DATA ----------------------------------------------------------

    print("ctheta: {}".format(cthetas.shape))
    print("Rs: {}".format(rs.shape))

    return (t, rs, ds, cthetas, val)
def plot_polar_mpl(ax, theta, r):
    ax.plot(theta, r)
    ax.set_rlim(0)
    ax.set_title('polar matplotlib')
def plot_logpolar_mpl(ax, theta, r):
    ax.plot(theta, r)
    ax.set_rlim(0)
    ax.set_rscale('log')
    ax.set_title('log-polar matplotlib')
def plot_logpolar(ax, theta, r_, bullseye=None, **kwargs):
    min10 = np.log10(np.min(r_))
    max10 = np.log10(np.max(r_))
    if bullseye is None:
        bullseye = min10 - np.log10(0.5 * np.min(r_))
    r = np.log10(r_) - min10 + bullseye
    ax.plot(theta, r, **kwargs)
    l = np.arange(np.floor(min10), max10)
    ax.set_rticks(l - min10 + bullseye)
    ax.set_yticklabels(["1e%d" % x for x in l])
    ax.set_rlim(0, max10 - min10 + bullseye)
    ax.set_title('log-polar manual')
    return ax
def plot_logpolar_bar(ax, theta, theta_width, r_, height, bullseye=None, **kwargs):
    min10 = np.log10(np.min(r_))
    max10 = np.log10(np.max(r_))
    if bullseye is None:
        bullseye = min10 - np.log10(0.5 * np.min(r_))
    # r1 = np.log10(r_) - min10 + bullseye
    # ax.plot(theta, r, **kwargs)
    ax.bar(x=theta, height = height, width=theta_width, bottom=np.log10(r_) - min10 + bullseye)
    l = np.arange(np.floor(min10), max10)
    ax.set_rticks(l - min10 + bullseye)
    ax.set_yticklabels(["1e%d" % x for x in l])
    ax.set_rlim(0, max10 - min10 + bullseye)
    ax.set_title('log-polar manual')
    return ax
def plot_for_movie2():
    t, rs, ds, cthetas, val = _load_data(ir=1)
    nshells = rs.shape[1]
    nlayers = rs.shape[0]
    heights = np.log10(ds)#np.log10(np.diff(rs,axis=1))
    heights[~np.isfinite(heights)] = 0.#np.log10(np.diff(rs,axis=1))
    lrs = np.log10(rs)
    # heights = np.insert(heights,-1,0.5*heights[-1])
    # theta_widths = np.diff(heights,axis=0)
    # theta_widths = np.insert(theta_widths,-1,theta_widths[-1])

    # plt.plot(rs[0,:],val[0,:],'x')
    # plt.plot(rs[1,:],val[1,:],'x')
    # plt.plot(rs[2,:],val[2,:],'x')
    # plt.plot(rs[3,:],val[3,:],'x')

    # plt.plot(ds[0,:],val[0,:],'o')
    # plt.plot(ds[1,:],val[1,:],'x')
    # plt.plot(ds[2,:],val[2,:],'p')
    # plt.plot(ds[3,:],val[3,:],'d')

    plt.loglog(rs[0,:],ds[0,:],'o')
    plt.loglog(rs[1,:],ds[1,:],'x')
    plt.loglog(rs[2,:],ds[2,:],'p')
    plt.loglog(rs[3,:],ds[3,:],'d')

    # plt.plot(rs[:,0],val[:,0],'x')
    # plt.plot(rs[:,1],val[:,1],'o')
    # plt.plot(rs[:,2],val[:,2],'p')
    # plt.plot(rs[:,3],val[:,3],'v')

    plt.show()

    # plt.pcolormesh(range(nshells),range(nlayers),cthetas)
    # plt.pcolormesh(rs,cthetas,val)
    # plt.show()

    ctheta_grid = [ctheta_shell.max() for ctheta_shell in cthetas[:]]
    ctheta_grid = np.insert(ctheta_grid, len(ctheta_grid), np.pi/2.)

    # ctheta_width = np.diff(ctheta_grid)

    # plt.plot(ctheta_grid,np.zeros_like(ctheta_grid),marker = 'x',ls='none')
    # plt.plot(ctheta_grid[:-1]+np.diff(ctheta_grid)/2,np.zeros_like(ctheta_grid[:-1]),marker = 'o',ls='none')
    # plt.show()

    cthetas_centers = ctheta_grid[:-1]+np.diff(ctheta_grid)/2
    cthetas_widths = np.diff(ctheta_grid)*.999


    print(f"{ctheta_grid}, PI/2={np.pi/2.}")
    # ctheta_width = np.zeros(len(ctheta_grid))
    # for il in range(nlayers):
    #     if (il == 0):
    #         ctheta_width[0] = (ctheta_grid[1]-ctheta_grid[0])
    #     elif (il == nlayers-1):
    #         ctheta_width[il] = (ctheta_grid[il]-ctheta_grid[il-1])
    #     else:
    #         ctheta_width[il] = (ctheta_grid[il+1]-ctheta_grid[il-1])





    fig, ax = plt.subplots(nrows=1,ncols=1,figsize=[6,6],
                           subplot_kw={'projection': 'polar'})

    # val = np.log10(val)
    norm = LogNorm(vmin=val[np.isfinite(val) * val>0].min(), vmax=val[np.isfinite(val) * val>0].max())
    # norm = Normalize(val[np.isfinite(val)].min(), val[np.isfinite(val)].max())
    cmap = cm.get_cmap("inferno")
    cmap.set_under("lime")
    cmap.set_over("black")

    # ax = axes[0]
    for il in range(nlayers):
        for ish in range(nshells):
            bar_theta_mid = cthetas_centers[il]
            bar_theta_width = cthetas_widths[il]
            bar_x = rs[il,ish]
            bar_height = ds[il,ish]

            if (bar_x > 0 and bar_height > 0):
                bar_height = np.log10(bar_height)
                bar_x = np.log10(bar_x)
                rgba_color = cmap(norm(val[il,ish]))
                ax.bar(x=bar_theta_mid, height = bar_height, width=bar_theta_width, bottom=bar_x,# edgecolor='black',
                       color=rgba_color)
    ax.set_thetamax(np.rad2deg(np.pi/2.))
    ax.set_rlabel_position(90)
    # ax.grid(False)
    ax.set_xticks(ctheta_grid)
    tick = [ax.get_rmax(),ax.get_rmax()*0.97]
    for t in np.deg2rad(np.arange(0,90,5)):
        ax.plot([t,t], tick, lw=0.72, color="k")
    # ax.set_rticks(ctheta_grid)

    # fig.colorbar(pc)
    ax3 = fig.add_axes([0.9, 0.1, 0.03, 0.8])

    cb1 = mpl.colorbar.ColorbarBase(ax3, cmap=cmap, norm=norm, orientation='vertical', extend='both')

    ax.set_title('t={:.2e} s'.format(t.max()))
    plt.show()
    exit(1)




    theta1=15
    theta2=45
    r1 = 5
    r2 = 7

    fig, ax = plt.subplots(figsize=[6,6],
                           subplot_kw={'projection': 'polar'})
    theta_mid = np.deg2rad((theta1 + theta2)/2)
    theta_width = np.deg2rad(theta2 - theta1)
    height = r2 - r1
    ax.bar(x=theta_mid, height = height, width=theta_width, bottom=r1, edgecolor='black')
    ax.set_rlim(4, 8)
    # ax.set_rmin(1e4)
    # plt.rgrids([0.1, 0.2, 0.4, 0.6, 0.8])
    # ax.set_rscale("log")

    plt.show()













    r = 10 ** np.arange(-3, 1.0, 0.0001)
    theta = 2 * np.pi * np.log10(r)

    ax = plt.subplots(ncols=1, nrows=1, subplot_kw=dict(polar=True))[1]
    # plot_polar_mpl(ax[0], theta, r)
    # plot_logpolar_mpl(ax[1], theta, r)
    # plot_logpolar(ax, theta, r)

    rmin = 1e4
    rmax = 1e9

    min10 = np.log10(rmin)
    max10 = np.log10(rmax)
    # if bullseye is None:
    #     bullseye = min10 - np.log10(0.5 * np.min(r_))
    # r = np.log10(r_) - min10 + bullseye
    # ax.plot(theta, r)
    theta1=45
    theta2=80
    r1 = np.log10(1.e5)
    r2 = np.log10(1.e7)
    theta_mid = np.deg2rad((theta1 + theta2)/2)
    height = r2 - r1
    theta_width = np.deg2rad(theta2 - theta1)

    plot_logpolar_bar(ax,theta,theta_width,r1,r2-r1,None)

    l = np.arange(np.floor(min10), max10)
    ax.set_rticks(l - min10)
    ax.set_yticklabels(["1e%d" % x for x in l])
    ax.set_rlim(0, max10 - min10)
    ax.set_title('log-polar manual')
    # return ax

    plt.show()


def _plot_one_quarter(ax, cax, ir, v_n, theta_offset):

    t, rs, ds, cthetas, val = _load_data(ir=ir, v_n=v_n)

    nshells = rs.shape[1]
    nlayers = rs.shape[0]

    norm = LogNorm(vmin=val[np.isfinite(val) * val>0].min(), vmax=val[np.isfinite(val) * val>0].max())
    # norm = Normalize(val[np.isfinite(val)].min(), val[np.isfinite(val)].max())
    cmap = cm.get_cmap("inferno")
    cmap.set_under("lime")
    cmap.set_over("black")

    ctheta_grid = np.array([ctheta_shell.max() for ctheta_shell in cthetas[:]])+theta_offset
    ctheta_grid = np.insert(ctheta_grid, len(ctheta_grid), np.pi/2.+theta_offset)
    cthetas_centers = ctheta_grid[:-1]+np.diff(ctheta_grid)/2
    cthetas_widths = np.diff(ctheta_grid)*.999

    # ax = axes[0]
    for il in range(nlayers):
        for ish in range(nshells):
            bar_theta_mid = cthetas_centers[il]
            bar_theta_width = cthetas_widths[il]
            bar_x = rs[il,ish]
            bar_height = ds[il,ish]

            if (bar_x > 0 and bar_height > 0):
                bar_height = np.log10(bar_height)
                bar_x = np.log10(bar_x)
                rgba_color = cmap(norm(val[il,ish]))
                ax.bar(x=bar_theta_mid, height = bar_height, width=bar_theta_width, bottom=bar_x,# edgecolor='black',
                       color=rgba_color)
    # ax.set_thetamax(np.rad2deg(np.pi/2.))
    # ax.set_rlabel_position(90)
    # divider = make_axes_locatable(ax)
    # cax = divider.new_vertical(size='1%', pad=0.1)#, pack_start = True)
    # fig.add_axes(cax)
    return (cmap, norm)



def plot_for_movie3(ir):

    fig, axes = plt.subplots(nrows=2,ncols=2,figsize=[6,6],
                           subplot_kw={'projection': 'polar'})
    if not hasattr(axes,"__len__"): axes = np.array([axes])
    axes = axes.flatten()
    v_ns = ["mom","mom","mom","mom"]

    for i in range(len(axes)):
        ax = axes[i]
        theta_offset = 0
        if i == 0:
            theta_offset = 90
        elif i == 1:
            theta_offset = 0.
        elif i == 2:
            theta_offset = 180.
        elif i == 3:
            theta_offset = 270.
        else:
            pass

        ax.set_thetamin(theta_offset)
        ax.set_thetamax(theta_offset+90.)
        (cmap, norm) = _plot_one_quarter(ax=ax, cax=None, ir=ir, v_n=v_ns[i], theta_offset=np.deg2rad(theta_offset))

        if i == 0:
            theta_offset = 90 # (left, bottom, width, height)
            cax = fig.add_axes([0.01, 0.5, 0.03, 0.8/2.])
            # cax.axis('off')
            cb1 = mpl.colorbar.ColorbarBase(cax,cmap=cmap, norm=norm, orientation='vertical', extend='both')
            cax.yaxis.set_ticks_position('left')
        elif i == 1:
            theta_offset = 0.
            # cax = fig.add_axes([0.9, 0.1, 0.03, 0.8])
        elif i == 2:
            theta_offset = 180.
            # cax = fig.add_axes([0.9, 0.1, 0.03, 0.8])
        elif i == 3:
            theta_offset = 270.
            # cax = fig.add_axes([0.9, 0.1, 0.03, 0.8])
        else:
            pass

    plt.show()

    for icol in range(2):
        for irow in range(2):
            ax = axes[icol][irow]
            ax.set_thetamax(np.rad2deg(np.pi/2.))
            ax.set_thetamin(0)



# val = np.log10(val)

    # ax.grid(False)
    ax.set_xticks(ctheta_grid)
    tick = [ax.get_rmax(),ax.get_rmax()*0.97]
    for t in np.deg2rad(np.arange(0,90,5)):
        ax.plot([t,t], tick, lw=0.72, color="k")
    # ax.set_rticks(ctheta_grid)

    # fig.colorbar(pc)
    ax3 = fig.add_axes([0.9, 0.1, 0.03, 0.8])

    cb1 = mpl.colorbar.ColorbarBase(ax3, cmap=cmap, norm=norm, orientation='vertical', extend='both')

    ax.set_title('t={:.2e} s'.format(t.max()))
    plt.show()
    exit(1)




    theta1=15
    theta2=45
    r1 = 5
    r2 = 7

    fig, ax = plt.subplots(figsize=[6,6],
                           subplot_kw={'projection': 'polar'})
    theta_mid = np.deg2rad((theta1 + theta2)/2)
    theta_width = np.deg2rad(theta2 - theta1)
    height = r2 - r1
    ax.bar(x=theta_mid, height = height, width=theta_width, bottom=r1, edgecolor='black')
    ax.set_rlim(4, 8)
    # ax.set_rmin(1e4)
    # plt.rgrids([0.1, 0.2, 0.4, 0.6, 0.8])
    # ax.set_rscale("log")

    plt.show()













    r = 10 ** np.arange(-3, 1.0, 0.0001)
    theta = 2 * np.pi * np.log10(r)

    ax = plt.subplots(ncols=1, nrows=1, subplot_kw=dict(polar=True))[1]
    # plot_polar_mpl(ax[0], theta, r)
    # plot_logpolar_mpl(ax[1], theta, r)
    # plot_logpolar(ax, theta, r)

    rmin = 1e4
    rmax = 1e9

    min10 = np.log10(rmin)
    max10 = np.log10(rmax)
    # if bullseye is None:
    #     bullseye = min10 - np.log10(0.5 * np.min(r_))
    # r = np.log10(r_) - min10 + bullseye
    # ax.plot(theta, r)
    theta1=45
    theta2=80
    r1 = np.log10(1.e5)
    r2 = np.log10(1.e7)
    theta_mid = np.deg2rad((theta1 + theta2)/2)
    height = r2 - r1
    theta_width = np.deg2rad(theta2 - theta1)

    plot_logpolar_bar(ax,theta,theta_width,r1,r2-r1,None)

    l = np.arange(np.floor(min10), max10)
    ax.set_rticks(l - min10)
    ax.set_yticklabels(["1e%d" % x for x in l])
    ax.set_rlim(0, max10 - min10)
    ax.set_title('log-polar manual')
    # return ax

    plt.show()

if __name__ == '__main__':
    # plot_ejecta_dyn_evol_for_movie()
    # plot_for_movie2()
    plot_for_movie3(ir=100)