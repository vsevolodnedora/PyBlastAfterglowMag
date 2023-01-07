import matplotlib.pyplot as plt
import numpy as np
import h5py
from scipy import interpolate
import copy
import os
import sys
import argparse

def gauss_eneregy_dist(E_iso_c, theta, theta_c):
    # E_iso_c * np.exp(-0.5 * ( thetas_c[i] * thetas_c[i] / theta_c / theta_c ) )
    # E_iso_c * np.exp( -1. * cthetas0[i] * cthetas0[i] / (theta_c * theta_c) )
    return E_iso_c * np.exp(-0.5 * ( theta * theta / theta_c / theta_c ) )

def make_gaussian_dist_a(E_iso_c, Gamma0c, theta_c, theta_w, M0c, n_layers_a, gflat=False):
    c = 2.9979e10
    # set grid
    thetas_c_l = np.zeros( n_layers_a )
    thetas_c_h = np.zeros( n_layers_a )
    thetas_c = np.zeros( n_layers_a )
    dtheta = theta_w / n_layers_a
    for i in range(n_layers_a):
        theta_c_i = i * dtheta + dtheta / 2.
        i_theta_c_l = i * dtheta
        i_theta_c_h = (i + 1) * dtheta
        # i_theta_h = i_theta_c_h
        thetas_c[i] = theta_c_i
        thetas_c_l[i] = i_theta_c_l
        thetas_c_h[i] = i_theta_c_h
    # ---
    dist_E0_a = np.zeros ( n_layers_a )
    dist_G0_a = np.zeros ( n_layers_a )
    dist_M0_a = np.zeros ( n_layers_a )
    for i in range(n_layers_a):
        frac_of_solid_ang = 2 * np.sin(0.5 * thetas_c_h[i]) * np.sin(0.5 * thetas_c_h[i])
        dist_E0_a[i] = gauss_eneregy_dist(E_iso_c, thetas_c[i], theta_c)#E_iso_c * np.exp(-0.5 * ( thetas_c[i] * thetas_c[i] / theta_c / theta_c ) )
        if (gflat):
            dist_G0_a[i] = Gamma0c
        else:
            dist_G0_a[i] = 1.0 + (Gamma0c - 1) * dist_E0_a[i] / E_iso_c
        dist_M0_a[i] = dist_E0_a[i] / (( dist_G0_a[i] - 1.0) * c * c )
        # dist_E0_a[i] *= ( frac_of_solid_ang / 2. )
        # dist_M0_a[i] *= ( frac_of_solid_ang / 2. )

    return (thetas_c, dist_G0_a, dist_M0_a, dist_E0_a)

def make_gaussian_dist_pw(E_iso_c, Gamma0c, theta_c, theta_w, M0c, n_layers_pw, gflat=False):

    c = 2.9979e10

    theta_pw = np.zeros( n_layers_pw + 1 )
    cthetas0 = np.zeros( n_layers_pw )
    for i in range(n_layers_pw + 1):
        fac = i / n_layers_pw
        theta_pw[i] = 2.0 * np.arcsin( fac * np.sin(theta_w / 2.0 ) )
    thetas_h0_pw = np.zeros( n_layers_pw )
    for i in range(n_layers_pw):
        cthetas0[i] = 0.5 * ( theta_pw[i+1] + theta_pw[i] )
        thetas_h0_pw[i] = theta_pw[i + 1]

    def CellsInLayer(i_layer):
        return 2 * i_layer + 1

    cil = np.zeros( n_layers_pw )
    for i in range(n_layers_pw):
        cil[i] = CellsInLayer(i)
    ncells = cil.sum() # total number of cells

    dist_E0_pw = np.zeros( n_layers_pw )
    dist_G0_pw = np.zeros( n_layers_pw )
    dist_M0_pw = np.zeros( n_layers_pw )
    ang_size_layer = 2.0 * np.pi * ( 2.0 * np.sin(0.5 * theta_w) * np.sin(0.5 * theta_w) ) / (4.0 * np.pi)
    for i in range(n_layers_pw):
        dist_E0_pw[i] = gauss_eneregy_dist(E_iso_c, cthetas0[i], theta_c)#E_iso_c * np.exp( -1. * cthetas0[i] * cthetas0[i] / (theta_c * theta_c) )
        if (gflat):
            dist_G0_pw[i] = Gamma0c
        else:
            dist_G0_pw[i] = 1. + (Gamma0c - 1.) * np.exp( -1. * cthetas0[i] * cthetas0[i] / (2. * theta_c * theta_c) )
        # dist_E0_pw[i] *= ang_size_layer
        dist_M0_pw[i] = dist_E0_pw[i] / (dist_G0_pw[i] * c * c)
        # dist_E0_pw[i] /= ncells
        # dist_M0_pw[i] /= ncells
    return (cthetas0, dist_G0_pw, dist_M0_pw, dist_E0_pw)

def main():

    thetas_c, dist_G0_a, dist_M0_a, dist_E0_a = make_gaussian_dist_a(E_iso_c = 1.e52,
        Gamma0c = 300.,
        M0c = -1.,
        theta_c = 0.085,
        theta_w = 0.2618,
        n_layers_a = 20)

    cthetas0, dist_G0_pw, dist_M0_pw, dist_E0_pw = make_gaussian_dist_pw(E_iso_c = 1.e52,
        Gamma0c = 300.,
        M0c = -1.,
        theta_c = 0.085,
        theta_w = 0.2618,
        n_layers_pw = 20)

    plt.semilogy(thetas_c, dist_E0_a, marker='x', label='a')
    plt.semilogy(cthetas0, dist_E0_pw, marker='.', label='pw')

    # plt.semilogy(thetas_c, dist_G0_a, marker='x', label='a')
    # plt.semilogy(cthetas0, dist_G0_pw, marker='.', label='pw')

    plt.legend()
    plt.show()

if __name__ == '__main__':
    main()