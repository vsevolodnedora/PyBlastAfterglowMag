import numpy as np
import os
import sys
import copy
from .utils import cgs

rebin_paths = [
    "../../../../../GIT/GitHub/rebin",
    "/home/vsevolod/Work/GIT/GitHub/rebin",
    "/home/enlil/vnedora/work/afterglow/rebin"
]

imported = False
for path in rebin_paths:
    if os.path.isdir(path):
        sys.path.insert(1, path)
        import rebin
        imported = True
        break
if not imported:
    raise ImportError("Faild to import rebin. Tried: these paths:{}".format(rebin_paths))

def compute_ek_corr(_vinf, _mass):
    tmp_ek = []
    _vinf = copy.deepcopy(np.asarray(_vinf))
    _mass = copy.deepcopy(np.asarray(_mass))
    for i in range(len(_mass[:, 0])):
        # tmp = _mass[i, :]
        # tmp *= cgs.solar_m
        # tmp *= _vinf
        # tmp *= _vinf
        # tmp *= cgs.c
        # tmp *= cgs.c
        # x = 1

        # arr = _mass[i, :] * cgs.solar_m
        # arr1 = (_vinf * _vinf * cgs.c * cgs.c)
        # arr *= arr1
        tmp_ek.append(_mass[i, :] * cgs.solar_m * _vinf * _vinf * cgs.c * cgs.c)
        # tmp_ek.append( copy.deepcopy( tmp ) )
    res = np.reshape(np.array(tmp_ek), newshape=(len(_mass[:, 0]), len(_vinf)))
    assert res.shape == _mass.shape
    return res

def _generate_grid_cthetas(nlayers, theta0):
    fac = np.arange(0, nlayers + 1) / float(nlayers)
    thetas = 2. * np.arcsin(fac * np.sin(theta0 / 2.))
    cthetas = 0.5 * (thetas[1:] + thetas[:-1])
    return (thetas, cthetas)
def _generate_grid_cthetas2(nlayers, theta0, theta1):
    dtheta = theta1 / (nlayers + 1)
    thetas_l, thetas_h, cthetas = [], [], []
    for i in range(nlayers+1):
        cthetas.append(i * dtheta + dtheta / 2.) # = (double)i * dtheta + dtheta / 2.;
        thetas_l.append(i * dtheta)# = (double) i * dtheta;
        thetas_h.append((i + 1) * dtheta) # double i_theta_c_h = (double) (i + 1) * dtheta;
    cthetas = np.array(cthetas)
    thetas_h = np.array(thetas_h)
    thetas_l = np.array(thetas_l)
    # fac = np.arange(0, nlayers + 1) / float(nlayers)
    # thetas = 2. * np.arcsin(fac * np.sin(theta0 / 2.))
    # cthetas = 0.5 * (thetas[1:] + thetas[:-1])
    return (thetas_h, cthetas)
def reinterpolate_hist(thetas_pol_edges, mass_2d_hist, new_theta_len=None, dist="pw"):
    print("Rebinning historgram")
    if (new_theta_len is None):
        new_theta_len = len(thetas_pol_edges)
    # if ()
    if dist == "pw":
        new_thetas_edges, new_theta_centers = _generate_grid_cthetas(new_theta_len - 1, theta0=np.pi / 2.)
    elif dist == "a":
        new_thetas_edges, new_theta_centers = _generate_grid_cthetas2(new_theta_len-1, theta0=0., theta1=np.pi / 2.)
    else:
        raise KeyError("Only 'pw' or 'a' are supported")
    # x = range(new_thetas_edges)
    # plt.close()
    # plt.plot(range(len(new_thetas_edges)), new_thetas_edges, ls='none', marker='x',  color="black")
    # plt.plot(range(len(new_thetas_edges2)), new_thetas_edges2, ls='none', marker='.', color="red")
    # plt.axhline(y=np.pi/2.)
    # plt.show()
    if (len(thetas_pol_edges) - 1 != len(mass_2d_hist[:, 0])):
        raise ValueError("something is wrong")

    if (len(new_thetas_edges) != len(thetas_pol_edges)):
        print("Change theta_grid {}->{}".format(len(thetas_pol_edges), len(new_thetas_edges)))

    new_mass = np.zeros((len(new_thetas_edges) - 1, len(mass_2d_hist[0, :])))
    for ibeta in range(len(mass_2d_hist[0, :])):
        tmp = rebin.rebin(thetas_pol_edges, mass_2d_hist[:, ibeta], new_thetas_edges,
                          interp_kind='piecewise_constant')
        new_mass[:, ibeta] = tmp

    thetas_pol_edges = new_thetas_edges
    mass = new_mass
    return (thetas_pol_edges, mass)