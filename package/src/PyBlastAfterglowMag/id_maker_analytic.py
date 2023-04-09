import matplotlib.pyplot as plt
import numpy as np
import h5py
from scipy import interpolate
import copy
import os
import sys
import argparse

''' UNUSED NOT IMPLEMENTED NOT SUPPORTED '''

def get_Beta(Gamma):
    return (1. / Gamma) * np.sqrt((Gamma - 1.) * (Gamma + 1.))

def MomFromGam(gam):
    return np.sqrt(gam*gam - 1.)

def GamFromMom(mom):
    return np.sqrt(1.0+mom*mom)

def BetFromMom(mom):
    return mom / GamFromMom(mom)

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

    i_theta_c_h = np.hstack((np.array([0.]),thetas_c_h)).flatten()
    return (i_theta_c_h, thetas_c, dist_G0_a, dist_M0_a, dist_E0_a)

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
        dist_E0_pw[i] =  ang_size_layer * gauss_eneregy_dist(E_iso_c, cthetas0[i], theta_c)#E_iso_c * np.exp( -1. * cthetas0[i] * cthetas0[i] / (theta_c * theta_c) )
        if (gflat):
            dist_G0_pw[i] = Gamma0c
        else:
            dist_G0_pw[i] = 1. + (Gamma0c - 1.) * np.exp( -1. * cthetas0[i] * cthetas0[i] / (2. * theta_c * theta_c) )
        # dist_E0_pw[i] *= ang_size_layer
        dist_M0_pw[i] = dist_E0_pw[i] / (dist_G0_pw[i] * c * c)
        dist_E0_pw[i] /= ncells * CellsInLayer(i_layer=i) # TODO workaround (inside the code i devide by CellsInLyaer()
        dist_M0_pw[i] /= ncells * CellsInLayer(i_layer=i)

        '''
        
( 6.80844e+46, 6.75726e+46, 6.65604e+46, 6.50705e+46, 6.31357e+46, 6.07977e+46, 5.8106e+46, 5.51158e+46, 5.18861e+46, 4.84782e+46, 4.49531e+46, 4.13706e+46, 3.77868e+46, 3.42535e+46, 3.08165e+46, 2.75154e+46, 2.43827e+46, 2.14435e+46, 1.87163e+46, 1.62126e+46, 1.39377e+46, 1.18914e+46, 1.00688e+46, 8.46098e+45, 7.05607e+45, 5.83984e+45, 4.79659e+45, 3.90981e+45, 3.16277e+45, 2.53902e+45, 2.02277e+45, 1.59923e+45, 1.25473e+45, 9.76938e+44, 7.54841e+44, 5.78779e+44, 4.40389e+44, 3.32523e+44, 2.49154e+44, 1.85254e+44, 1.36685e+44, 1.00074e+44, 7.27058e+43, 5.24151e+43, 3.74956e+43, 2.66156e+43, 1.87466e+43, 1.31018e+43, 9.08577e+42, 6.25183e+42, )
( 299.857, 298.732, 296.494, 293.168, 288.791, 283.412, 277.09, 269.892, 261.895, 253.181, 243.839, 233.962, 223.643, 212.978, 202.062, 190.987, 179.845, 168.719, 157.691, 146.834, 136.215, 125.895, 115.925, 106.35, 97.2057, 88.5215, 80.3185, 72.6108, 65.4057, 58.7047, 52.5034, 46.7927, 41.5592, 36.7858, 32.4527, 28.5374, 25.0159, 21.8631, 19.0528, 16.5591, 14.3559, 12.4177, 10.7197, 9.23825, 7.95083, 6.83621, 5.87461, 5.0477, 4.33867, 3.73217, )
( 2.52632e+23, 2.51678e+23, 2.4978e+23, 2.46959e+23, 2.43247e+23, 2.38685e+23, 2.33322e+23, 2.27218e+23, 2.20435e+23, 2.13045e+23, 2.05122e+23, 1.96744e+23, 1.87992e+23, 1.78947e+23, 1.69689e+23, 1.60297e+23, 1.50846e+23, 1.41411e+23, 1.32058e+23, 1.2285e+23, 1.13844e+23, 1.05092e+23, 9.66366e+22, 8.85164e+22, 8.07619e+22, 7.33978e+22, 6.64418e+22, 5.99062e+22, 5.37971e+22, 4.81159e+22, 4.28588e+22, 3.80182e+22, 3.35828e+22, 2.95382e+22, 2.58677e+22, 2.25523e+22, 1.95718e+22, 1.6905e+22, 1.45301e+22, 1.24251e+22, 1.05681e+22, 8.93794e+21, 7.51387e+21, 6.27618e+21, 5.20617e+21, 4.2863e+21, 3.50025e+21, 2.83294e+21, 2.27052e+21, 1.80031e+21, )
( 0.00261053, 0.00783162, 0.0130528, 0.018274, 0.0234953, 0.0287168, 0.0339385, 0.0391605, 0.0443827, 0.0496052, 0.054828, 0.0600513, 0.0652749, 0.070499, 0.0757235, 0.0809486, 0.0861742, 0.0914004, 0.0966273, 0.101855, 0.107083, 0.112312, 0.117542, 0.122772, 0.128003, 0.133236, 0.138469, 0.143703, 0.148938, 0.154174, 0.159411, 0.164649, 0.169889, 0.175129, 0.180371, 0.185614, 0.190858, 0.196104, 0.201351, 0.206599, 0.211849, 0.2171, 0.222353, 0.227607, 0.232863, 0.238121, 0.24338, 0.248641, 0.253903, 0.259167, )

        '''

    return (theta_pw, cthetas0, dist_G0_pw, dist_M0_pw, dist_E0_pw)



class JetStruct:
    def __init__(self, n_layers_pw,n_layers_a):
        self.m_theta_w = 0
        self.m_theta_c = 0
        self.dist_E0_pw = np.zeros(n_layers_pw)
        self.dist_Mom0_pw = np.zeros(n_layers_pw)
        self.dist_M0_pw = np.zeros(n_layers_pw)
        self.dist_Ye_pw = np.zeros(n_layers_pw)
        self.dist_s_pw = np.zeros(n_layers_pw)
        self.cthetas0 = np.zeros(n_layers_pw)
        self.theta_pw = np.zeros(n_layers_pw)
        self.cil = np.zeros(n_layers_pw)
        self.nlayers_pw = n_layers_pw
        self.ncells = 0

        self.nlayers_a = n_layers_a
        self.thetas_c_l = np.zeros(n_layers_a)
        self.thetas_c_h = np.zeros(n_layers_a)
        self.thetas_c = np.zeros(n_layers_a)
        # self.thetas_c_l = np.zeros(n_layers_a)


    @staticmethod
    def CellsInLayer(i_layer):
        return 2 * i_layer + 1

    def setThetaGridA(self):

        self.thetas_c_l=np.zeros( self.nlayers_a )
        self.thetas_c_h=np.zeros( self.nlayers_a )
        self.thetas_c=np.zeros( self.nlayers_a )

        dtheta = self.m_theta_w / self.nlayers_a
        #        double _tmp=0;
        for i in range(self.nlayers_a):

            # account for geometry
            theta_c_i = i * dtheta + dtheta / 2.
            i_theta_c_l = i * dtheta
            i_theta_c_h = (i + 1) * dtheta
            # double i_theta_h = i_theta_c_h;

            self.thetas_c[i] = theta_c_i
            self.thetas_c_l[i] = i_theta_c_l
            self.thetas_c_h[i] = i_theta_c_h
            #std::cout << "ilayer=" << i << "theta" \ _tmp+=i_theta_c_h;


    def setThetaGridPW(self):
        self.theta_pw = np.zeros( self.nlayers_pw + 1 )
        self.cthetas0 = np.zeros( self.nlayers_pw )
        for i in range(self.nlayers_pw + 1):
            fac = i / self.nlayers_pw
            self.theta_pw[i] = 2.0 * np.arcsin( fac * np.sin(self.m_theta_w / 2.0 ) )


        thetas_h0_pw = np.zeros( self.nlayers_pw )
        for i in range(self.nlayers_pw):
            self.cthetas0[i] = 0.5 * ( self.theta_pw[i+1] + self.theta_pw[i] )
            thetas_h0_pw[i] = self.theta_pw[i + 1]

        #        ncells_pw = 0;

        # compute the number of phi cells in each 'theta' layer
        cil = np.zeros( self.nlayers_pw )
        for i in range(self.nlayers_pw):
            cil[i] = JetStruct.CellsInLayer(i)
        self.ncells = cil.sum() # total number of cells

    def gaussian(self, E_iso_c, Gamma0c, theta_c, theta_w, M0c, gflat=False):

        self.m_theta_w = theta_w
        self.m_theta_c = theta_c

        c = 2.99792458e10

        # set piece-wise

        self.setThetaGridPW()
        ang_size_layer = 2.0 * np.pi * ( 2.0 * np.sin(0.5 * theta_w) * np.sin(0.5 * theta_w) );
        for i  in range(len(self.cthetas0)):
            self.dist_E0_pw[i] = E_iso_c * ang_size_layer / (4.0 * np.pi) * np.exp( -1. * self.cthetas0[i] * self.cthetas0[i] / (theta_c * theta_c) )
            Gamma = 0
            if (gflat):
                Gamma = MomFromGam( Gamma0c )
            else:
                Gamma = 1. + (Gamma0c - 1.) * np.exp(-1. * self.cthetas0[i] * self.cthetas0[i] / (2. * theta_c * theta_c))

            self.dist_Mom0_pw[i] = MomFromGam(Gamma)
            self.dist_M0_pw[i] = self.dist_E0_pw[i] / (Gamma * c * c)
            self.dist_E0_pw[i] /= self.ncells
            self.dist_E0_pw[i] *= JetStruct.CellsInLayer(i)
            self.dist_M0_pw[i] /= self.ncells
            self.dist_M0_pw[i] *= JetStruct.CellsInLayer(i)


        # set adaptive
        # self.nlayers_a = n_layers_a
        self.setThetaGridA()
        self.dist_E0_a = np.zeros( self.nlayers_a )
        self.dist_Mom0_a = np.zeros( self.nlayers_a )
        self.dist_M0_a = np.zeros( self.nlayers_a )
        self.dist_Ye_a = np.zeros( self.nlayers_a )
        self.dist_s_a = np.zeros( self.nlayers_a )
        for i in range(self.nlayers_a):
            frac_of_solid_ang = 2 * np.sin(0.5 * self.thetas_c_h[i]) * np.sin(0.5 * self.thetas_c_h[i])
            self.dist_E0_a[i] = E_iso_c * np.exp(-0.5 * ( self.thetas_c[i] * self.thetas_c[i] / theta_c / theta_c ) )
            Gamma = 0
            if (gflat):
                Gamma = Gamma0c
            else:
                Gamma = 1.0 + (Gamma0c - 1) * self.dist_E0_a[i] / E_iso_c
            self.dist_Mom0_a[i] = MomFromGam( Gamma )
            self.dist_M0_a[i] = self.dist_E0_a[i] / (( Gamma - 1.0) * c*c )

            # self.dist_E0_a[i] /= ( frac_of_solid_ang / 2. )
            # self.dist_M0_a[i] /= ( frac_of_solid_ang / 2. )

        i = 0

    def tophat(self, E_iso, Gamma0, theta_h, M0, Ye, s):

        self.m_theta_w = theta_h
        self.m_theta_c = theta_h

        c = 2.99792458e10


        # set piece-wise
        # self.nlayers_pw = n_layers;
        self.dist_E0_pw = np.zeros( self.nlayers_pw )
        self.dist_Mom0_pw = np.zeros( self.nlayers_pw )
        self.dist_M0_pw = np.zeros( self.nlayers_pw )
        self.dist_Ye_pw = np.zeros( self.nlayers_pw )
        self.dist_s_pw = np.zeros( self.nlayers_pw )
        self.setThetaGridPW()

        one_min_cos = 2 * np.sin(0.5 * theta_h) * np.sin(0.5 * theta_h)
        ang_size_layer = 2.0 * np.pi * one_min_cos
        for i in range(len(self.cthetas0)):
            self.dist_E0_pw[i] = E_iso * ang_size_layer / (4.0 * np.pi)
            self.dist_Mom0_pw[i] = MomFromGam(Gamma0)
            # Gamma = GamFromMom(self.dist_Mom0_pw[i])
            # self.dist_M0_pw[i] = M0 < 0 ? dist_E0_pw[i] / (Gamma * CGS::c * CGS::c) : M0 * ang_size_layer / (4 * np.pi);
            self.dist_M0_pw[i] = self.dist_E0_pw[i] / (Gamma0 * c * c) if M0 < 0 else M0 * ang_size_layer / (4 * np.pi)

            self.dist_E0_pw[i] /= self.ncells
            self.dist_E0_pw[i] *= JetStruct.CellsInLayer(i)
            self.dist_M0_pw[i] /= self.ncells
            self.dist_M0_pw[i] *= JetStruct.CellsInLayer(i)
            self.dist_Ye_pw[i] = Ye
            self.dist_s_pw[i] = s



        # set adaptive
        # self.nlayers_a = n_layers;
        self.setThetaGridA()
        self.dist_E0_a = np.zeros( self.nlayers_a )
        self.dist_Mom0_a= np.zeros( self.nlayers_a )
        self.dist_M0_a=np.zeros( self.nlayers_a )
        self.dist_Ye_a=np.zeros( self.nlayers_a )
        self.dist_s_a=np.zeros( self.nlayers_a )
        frac_of_solid_ang = 2 * np.sin(0.5 * theta_h) * np.sin(0.5 * theta_h) # for pi/2 -> 1.
        self.dist_E0_a[0] = E_iso * frac_of_solid_ang / 2.
        self.dist_Mom0_a[0] = MomFromGam(Gamma0)
        Gamma = GamFromMom(self.dist_Mom0_a[0])
        # self.dist_M0_a[0] = M0 < 0 ? E_iso / ((Gamma - 1.0) * CGS::c*CGS::c) * frac_of_solid_ang / 2. : M0 * frac_of_solid_ang / 2.;
        self.dist_M0_a[0] = E_iso / ((Gamma - 1.0) * c*c) * frac_of_solid_ang / 2. if M0 < 0 else M0 * frac_of_solid_ang / 2.
        self.dist_Ye_a[0] = 0.
        self.dist_s_a[0] = 0.

        # with as file:
        #     string = readline().split(" ")

        # os.path.isfile()
        self.dist_E0_a /= (frac_of_solid_ang / 2.) # TODO Get rid of shis (used in code)
        self.dist_M0_a /= (frac_of_solid_ang / 2.)

    def saveCurrentStructure(self, outfpath, type="pw"):
        if type == "pw":
            dfile = h5py.File(outfpath, "w")
            dfile.create_dataset("r", data=np.zeros_like(self.theta_pw[1:]))
            dfile.create_dataset("theta", data=self.theta_pw[1:])
            dfile.create_dataset("ctheta", data=self.cthetas0)
            dfile.create_dataset("mom", data=self.dist_Mom0_pw)
            dfile.create_dataset("mass", data=self.dist_M0_pw)
            dfile.create_dataset("ek", data=self.dist_E0_pw)
            dfile.create_dataset("ye", data=self.dist_Ye_pw)
            dfile.create_dataset("s", data=self.dist_s_pw)
            dfile.close()
        elif(type=="a"):
            dfile = h5py.File(outfpath, "w")
            dfile.create_dataset("r", data=np.zeros_like(self.thetas_c_h))
            dfile.create_dataset("theta", data=self.thetas_c_h)
            dfile.create_dataset("ctheta", data=self.thetas_c)
            dfile.create_dataset("mom", data=self.dist_Mom0_a)
            dfile.create_dataset("mass", data=self.dist_M0_a)
            dfile.create_dataset("ek", data=self.dist_E0_a)
            dfile.create_dataset("ye", data=self.dist_Ye_a)
            dfile.create_dataset("s", data=self.dist_s_a)
            dfile.close()
        else:
            raise KeyError("not implemented")



def prepare_grb_ej_id_1d(pars, outfpath, type="pw"):

    o_jet = JetStruct(pars["nlayers_pw"],pars["nlayers_a"])
    if (pars["struct"] == "gaussian"):
        o_jet.gaussian(E_iso_c=pars["Eiso_c"],Gamma0c=pars["Gamma0c"],
                       theta_c=pars["theta_c"],theta_w=pars["theta_w"],M0c=pars["M0c"])
    elif (pars["struct"] == "tophat"):
        o_jet.tophat(E_iso=pars["Eiso_c"],Gamma0=pars["Gamma0c"],theta_h=pars["theta_c"],M0=pars["M0c"],Ye=0,s=0)
    else:
        raise KeyError("Not implemented")

    o_jet.saveCurrentStructure(outfpath=outfpath,type=type)

    #
    # if (type == "a"):
    #     thetas0, cthetas0, dist_G0, dist_M0, dist_E0 = make_gaussian_dist_a(
    #         pars["Eiso_c"], pars["Gamma0c"], pars["theta_c"], pars["theta_w"], pars["M0c"], pars["nlayers_a"])
    # elif (type == "pw"):
    #     thetas0, cthetas0, dist_G0, dist_M0, dist_E0 = make_gaussian_dist_pw(
    #         pars["Eiso_c"], pars["Gamma0c"], pars["theta_c"], pars["theta_w"], pars["M0c"], pars["nlayers_pw"])
    # else:
    #     raise KeyError(" use 'a' or 'pw' ")
    #
    # # ek = np.zeros((len(dist_G0),len(cthetas0)))
    # # for i in range(len(dist_G0)):
    # #     for j in range(len(cthetas0)):
    # #         ek[i,j] = dist_G0
    # ek = np.array(dist_E0)
    # mom = np.sqrt(dist_G0*dist_G0 - 1.)
    #
    # dfile = h5py.File(outfpath, "w")
    # dfile.create_dataset("theta", data=thetas0)
    # dfile.create_dataset("ctheta", data=cthetas0)
    # dfile.create_dataset("mom", data=mom)
    # dfile.create_dataset("ek", data=ek)
    # dfile.create_dataset("ye", data=np.zeros_like(ek))
    # dfile.create_dataset("s", data=np.zeros_like(ek))
    # dfile.close()
    # print("file saved: {}".format(outfpath))

def prepare_grb_ej_id_2d(pars, outfpath, type="pw"):

    if (type == "a"):
        thetas0,cthetas0, dist_G0, dist_M0, dist_E0 = make_gaussian_dist_a(
            pars["Eiso_c"], pars["Gamma0c"], pars["theta_c"], pars["theta_w"], pars["M0c"], pars["nlayers_a"])
    elif (type == "pw"):
        thetas0,cthetas0, dist_G0, dist_M0, dist_E0 = make_gaussian_dist_pw(
            pars["Eiso_c"], pars["Gamma0c"], pars["theta_c"], pars["theta_w"], pars["M0c"], pars["nlayers_pw"])
    else:
        raise KeyError(" use 'a' or 'pw' ")

    # ek = np.zeros((len(dist_G0),len(cthetas0)))
    # for i in range(len(dist_G0)):
    #     for j in range(len(cthetas0)):
    #         ek[i,j] = dist_G0
    ek = np.diag(dist_E0)


    dfile = h5py.File(outfpath, "w")
    dfile.create_dataset("theta", data=thetas0)
    dfile.create_dataset("ctheta", data=cthetas0)
    dfile.create_dataset("mom", data=np.sqrt(dist_G0*dist_G0 - 1.))
    dfile.create_dataset("ek", data=ek)
    dfile.create_dataset("ye", data=np.zeros_like(ek))
    dfile.create_dataset("s", data=np.zeros_like(ek))
    dfile.close()
    print("file saved: {}".format(outfpath))

def main():

    o_struct = JetStruct(Gamma0c = 300.,
                         M0c = -1.,
                         theta_c = 0.085,
                         theta_w = 0.2618,
                         n_layers_pw = 50,
                         n_layers_a = 20)

    thetas0, thetas_c, dist_G0_a, dist_M0_a, dist_E0_a = make_gaussian_dist_a(E_iso_c = 1.e52,
        Gamma0c = 300.,
        M0c = -1.,
        theta_c = 0.085,
        theta_w = 0.2618,
        n_layers_a = 20)

    thetas0, cthetas0, dist_G0_pw, dist_M0_pw, dist_E0_pw = make_gaussian_dist_pw(E_iso_c = 1.e52,
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