import matplotlib.pyplot as plt
import numpy as np
import h5py
from scipy import interpolate
import copy
import os
import sys
import argparse


from .utils import *

''' UNUSED NOT IMPLEMENTED NOT SUPPORTED '''


class JetStruct:
    def __init__(self, n_layers_pw, n_layers_a):
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

        self.beta0_min = 1.e-5 # if structure gives layer with too small velocity -- code takes too long


    @staticmethod
    def _CellsInLayer(i_layer):
        return 2 * i_layer + 1

    def _setThetaGridA(self):

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


    def _setThetaGridPW(self):
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
            cil[i] = JetStruct._CellsInLayer(i)
        self.ncells = cil.sum() # total number of cells

    def _gaussian(self, E_iso_c, Gamma0c, theta_c, theta_w, M0c, gflat=False):

        self.m_theta_w = theta_w
        self.m_theta_c = theta_c

        c = 2.99792458e10

        # set piece-wise

        self._setThetaGridPW()
        ang_size_layer = 2.0 * np.pi * ( 2.0 * np.sin(0.5 * theta_w) * np.sin(0.5 * theta_w) )
        for i  in range(len(self.cthetas0)):
            self.dist_E0_pw[i] = E_iso_c * ang_size_layer / (4.0 * np.pi) * np.exp( -1. * self.cthetas0[i] * self.cthetas0[i] / (theta_c * theta_c) )

            Gamma = 0
            if (gflat):
                Gamma = MomFromGamma( Gamma0c )
            else:
                Gamma = 1. + (Gamma0c - 1.) * np.exp(-1. * self.cthetas0[i] * self.cthetas0[i] / (2. * theta_c * theta_c))

            self.dist_Mom0_pw[i] = MomFromGamma(Gamma)
            self.dist_M0_pw[i] = self.dist_E0_pw[i] / (Gamma * c * c)
            self.dist_E0_pw[i] /= self.ncells
            self.dist_E0_pw[i] *= JetStruct._CellsInLayer(i)
            self.dist_M0_pw[i] /= self.ncells
            self.dist_M0_pw[i] *= JetStruct._CellsInLayer(i)
        # plt.semilogy(self.cthetas0,self.dist_E0_pw)
        # plt.show()

        # set adaptive
        self._setThetaGridA()
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
            self.dist_Mom0_a[i] = MomFromGamma( Gamma )
            self.dist_M0_a[i] = self.dist_E0_a[i] / (( Gamma - 1.0) * c*c )

            self.dist_E0_a[i] *= ( frac_of_solid_ang / 2. )
            self.dist_M0_a[i] *= ( frac_of_solid_ang / 2. )

        # remove data that might cause issue in C++ code (too slow blastwaves)
        mask = BetaFromMom(self.dist_Mom0_a) > self.beta0_min
        self.dist_E0_a = self.dist_E0_a[mask]
        self.dist_Mom0_a = self.dist_Mom0_a[mask]
        self.dist_M0_a = self.dist_M0_a[mask]
        self.dist_Ye_a = self.dist_Ye_a[mask]
        self.dist_s_a = self.dist_s_a[mask]
        self.thetas_c_l = self.thetas_c_l[mask]
        self.thetas_c_h = self.thetas_c_h[mask]
        self.thetas_c = self.thetas_c[mask]

        self.nlayers_a = len(self.dist_Mom0_a)

        i = 0

    def _tophat(self, E_iso, Gamma0, theta_h, M0, Ye, s):

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
        self._setThetaGridPW()

        one_min_cos = 2. * np.sin(0.5 * theta_h) * np.sin(0.5 * theta_h)
        ang_size_layer = 2.0 * np.pi * one_min_cos / (4.0 * np.pi)
        for i in range(len(self.cthetas0)):
            self.dist_E0_pw[i] = E_iso * ang_size_layer
            self.dist_Mom0_pw[i] = MomFromGamma(Gamma0)
            # Gamma = GamFromMom(self.dist_Mom0_pw[i])
            # self.dist_M0_pw[i] = M0 < 0 ? dist_E0_pw[i] / (Gamma * CGS::c * CGS::c) : M0 * ang_size_layer / (4 * np.pi);
            self.dist_M0_pw[i] = self.dist_E0_pw[i] / (Gamma0 * c * c) if M0 < 0 else M0 * ang_size_layer / (4 * np.pi)

            self.dist_E0_pw[i] /= self.ncells
            self.dist_E0_pw[i] *= JetStruct._CellsInLayer(i)
            self.dist_M0_pw[i] /= self.ncells
            self.dist_M0_pw[i] *= JetStruct._CellsInLayer(i)
            self.dist_Ye_pw[i] = Ye
            self.dist_s_pw[i] = s

        # set adaptive
        # self.nlayers_a = n_layers;
        self._setThetaGridA()
        self.dist_E0_a = np.zeros( self.nlayers_a )
        self.dist_Mom0_a = np.zeros( self.nlayers_a )
        self.dist_M0_a = np.zeros( self.nlayers_a )
        self.dist_Ye_a = np.zeros( self.nlayers_a )
        self.dist_s_a = np.zeros( self.nlayers_a )
        frac_of_solid_ang = 2 * np.sin(0.5 * theta_h) * np.sin(0.5 * theta_h) # for pi/2 -> 1.
        self.dist_E0_a[0] = E_iso * frac_of_solid_ang / 2.
        self.dist_Mom0_a[0] = MomFromGamma(Gamma0)
        Gamma = GammaFromMom(self.dist_Mom0_a[0])
        # self.dist_M0_a[0] = M0 < 0 ? E_iso / ((Gamma - 1.0) * CGS::c*CGS::c) * frac_of_solid_ang / 2. : M0 * frac_of_solid_ang / 2.;
        self.dist_M0_a[0] = E_iso / ((Gamma - 1.0) * c*c) * frac_of_solid_ang / 2. if M0 < 0 else M0 * frac_of_solid_ang / 2.
        self.dist_Ye_a[0] = 0.
        self.dist_s_a[0] = 0.

        # with as file:
        #     string = readline().split(" ")

        # os.path.isfile()
        # self.dist_E0_a /= (frac_of_solid_ang / 2.) # TODO Get rid of shis (used in code)
        # self.dist_M0_a /= (frac_of_solid_ang / 2.)

    def get_1D_id(self, pars : dict, type : str) -> tuple[dict,dict]:

        if (pars["struct"] == "gaussian"):
            self._gaussian(E_iso_c=pars["Eiso_c"], Gamma0c=pars["Gamma0c"], theta_c=pars["theta_c"], theta_w=pars["theta_w"], M0c=pars["M0c"])
        elif (pars["struct"] == "tophat"):
            self._tophat(E_iso=pars["Eiso_c"], Gamma0=pars["Gamma0c"], theta_h=pars["theta_c"], M0=pars["M0c"], Ye=0, s=0)
        # elif (pars["struct"] == "numeric"):
        #     self._numeric(mom=pars["mom"],ek=pars["ek"],theta_max=np.pi/2.)
        else:
            raise KeyError("Not implemented")

        if (type == "pw" or type=="piece-wise"):
            res = {
                "r":np.zeros_like(self.theta_pw[1:]),
                "theta":self.theta_pw[1:],
                "ctheta":self.cthetas0,
                "mom":self.dist_Mom0_pw,
                "mass":self.dist_M0_pw,
                "ek": self.dist_E0_pw,
                "ye":self.dist_Ye_pw,
                "s":self.dist_s_pw
            }
        elif(type=="a" or type=="adaptive"):
            res = {
                "r":np.zeros_like(self.thetas_c_h),
                "theta":self.thetas_c_h,
                "ctheta":self.thetas_c,
                "mom":self.dist_Mom0_a,
                "mass":self.dist_M0_a,
                "ek": self.dist_E0_a,
                "ye":self.dist_Ye_a,
                "s":self.dist_s_a
            }
        else: raise KeyError(f"type={type} is not recognized. Supported: piece-wise or adaptive")

        pars["eats_type"] = type
        pars["theta_wing"] = self.m_theta_w
        pars["theta_core"] = self.m_theta_c
        return (res, pars)

    def save_1d_id(self, id_dict : dict, id_pars : dict, outfpath : str):
        if os.path.isfile(outfpath):
            os.remove(outfpath)
        with h5py.File(outfpath, "w") as dfile:
            for key, data in id_dict.items():
                dfile.create_dataset(name=key, data=np.array(data))
            for key, data in id_pars.items():
                dfile.attrs.create(key, data=id_pars[key])

    def OLD_saveCurrentStructure(self, outfpath, type="pw"):
        if type == "pw":
            dfile = h5py.File(outfpath, "w")
            dfile.attrs.create(name="theta_wing",data=self.m_theta_w)
            dfile.attrs.create(name="theta_core",data=self.m_theta_c)
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

            # dfile.create_dataset("r", data=np.column_stack((self.thetas_c_h, np.zeros_like(self.thetas_c_h))))
            # dfile.create_dataset("theta", data=np.column_stack((self.thetas_c_h, np.zeros_like(self.thetas_c_h))))
            # dfile.create_dataset("ctheta", data=np.column_stack((self.thetas_c, np.zeros_like(self.thetas_c_h))))
            # dfile.create_dataset("mom", data=np.column_stack((self.dist_Mom0_a, np.zeros_like(self.thetas_c_h))))
            # dfile.create_dataset("mass", data=np.column_stack((self.dist_M0_a, np.zeros_like(self.thetas_c_h))))
            # dfile.create_dataset("ek", data=np.column_stack((self.dist_E0_a, np.zeros_like(self.thetas_c_h))))
            # dfile.create_dataset("ye", data=np.column_stack((self.dist_Ye_a, np.zeros_like(self.thetas_c_h))))
            # dfile.create_dataset("s", data=np.column_stack((self.dist_s_a, np.zeros_like(self.thetas_c_h))))
            dfile.attrs.create(name="theta_wing",data=self.m_theta_w)
            dfile.attrs.create(name="theta_core",data=self.m_theta_c)
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


from .id_tools import _generate_grid_cthetas, reinterpolate_hist2
class EjectaStruct:
    def __init__(self, data:dict, verbose:bool=False):
        self.mom = data["mom"]
        self.Gamma = GammaFromMom(self.mom)
        self.beta = BetaFromMom(self.mom)
        ek = data["ek"]
        self.masses = ek / (self.beta*cgs.c)**2 / cgs.solar_m # Msun
        if "mass" in data.keys(): self.masses = data["mass"]
        self.theta_max = data["theta_max"] if "theta_max" in data.keys() else np.pi/2.
        self.thetas = data["theta"] if "theta" in data.keys() else np.zeros(0,)
        # self.theta = np.linspace(0., np.pi/2., 30)
        self.verbose = verbose

        vinf_ = np.array([0.00505   , 0.01495   , 0.02485   , 0.03475   , 0.04465   ,
               0.05455   , 0.06445   , 0.07435   , 0.08425   , 0.09415   ,
               0.10405   , 0.11395   , 0.12385   , 0.13375001, 0.14365   ,
               0.15355   , 0.16345   , 0.17335001, 0.18325   , 0.19315   ,
               0.20305   , 0.21295001, 0.22284999, 0.23275   , 0.24265   ,
               0.25255001, 0.26245001, 0.27235001, 0.28224999, 0.29214999,
               0.30204999, 0.31195   , 0.32185   , 0.33175001, 0.34165001,
               0.35155001, 0.36144999, 0.37134999, 0.38124999, 0.39115   ,
               0.40105   , 0.41095001, 0.42085001, 0.43075001, 0.44064999,
               0.45054999, 0.46044999, 0.47035   , 0.48025   , 0.49015   ,
               0.50005001, 0.50994998, 0.51985002, 0.52974999, 0.53965002,
               0.54955   , 0.55944997, 0.56935   , 0.57924998, 0.58915001,
               0.59904999, 0.60895002, 0.61884999, 0.62875003, 0.63865   ,
               0.64854997, 0.65845001, 0.66834998, 0.67825001, 0.68814999,
               0.69805002, 0.70795   , 0.71785003, 0.72775   , 0.73764998,
               0.74755001, 0.75744998, 0.76735002, 0.77724999, 0.78715003,
               0.79705   , 0.80694997, 0.81685001, 0.82674998, 0.83665001,
               0.84654999, 0.85645002, 0.86635   , 0.87625003, 0.88615   ,
               0.89604998, 0.90595001, 0.91584998, 0.92575002, 0.93564999,
               0.94555002, 0.95545   , 0.96534997, 0.97525001, 0.98514998,
               0.99505001, 1.00495005])
        mass_ = np.array([1.42934993e-05, 2.81558612e-05, 4.18902864e-05, 5.49883686e-05,
               6.71728753e-05, 7.93477796e-05, 9.03650074e-05, 9.99638161e-05,
               1.06393039e-04, 1.14078934e-04, 1.20442004e-04, 1.25621365e-04,
               1.31380212e-04, 1.37401962e-04, 1.35773252e-04, 1.31803339e-04,
               1.28578022e-04, 1.21958766e-04, 1.12421781e-04, 1.05865845e-04,
               1.04305444e-04, 1.00847123e-04, 9.77009174e-05, 9.28183251e-05,
               8.64190741e-05, 8.04303496e-05, 7.34467442e-05, 6.45315520e-05,
               5.84138289e-05, 5.41905847e-05, 5.18406007e-05, 5.00567093e-05,
               4.79984981e-05, 4.52172200e-05, 4.36353079e-05, 4.12852484e-05,
               3.86902403e-05, 3.55739052e-05, 3.19736234e-05, 2.94090203e-05,
               2.65221517e-05, 2.32450411e-05, 2.06010387e-05, 1.82979120e-05,
               1.60633013e-05, 1.41810794e-05, 1.26227698e-05, 1.10825812e-05,
               9.98310092e-06, 8.59366857e-06, 7.68986024e-06, 6.72754388e-06,
               5.80179122e-06, 5.02050985e-06, 4.26874428e-06, 3.65377231e-06,
               3.09202510e-06, 2.60035665e-06, 2.17417334e-06, 1.88691697e-06,
               1.57081152e-06, 1.26042424e-06, 1.04541498e-06, 8.82787611e-07,
               7.88398391e-07, 6.76416664e-07, 6.05813300e-07, 5.20005786e-07,
               4.58172834e-07, 3.93091524e-07, 3.70763685e-07, 3.30565199e-07,
               2.69349669e-07, 2.35668985e-07, 2.01745396e-07, 1.63239616e-07,
               1.53146018e-07, 1.39251430e-07, 8.78833921e-08, 6.90512113e-08,
               5.36619981e-08, 3.69990454e-08, 2.74241617e-08, 2.09313324e-08,
               1.62076965e-08, 1.24396036e-08, 9.17189324e-09, 6.36312659e-09,
               4.96148054e-09, 3.84213972e-09, 2.76035849e-09, 1.89559595e-09,
               1.31431044e-09, 5.77015496e-10, 4.49675220e-11, 1.92655784e-13,
               7.57131716e-14, 1.55272744e-14, 0.00000000e+00, 0.00000000e+00,
               0.00000000e+00, 0.00000000e+00])

        # plt.semilogy(vinf_,mass_)
        # plt.semilogy(self.beta,self.masses)
        # plt.show()
        #
        # self.beta=vinf_
        # self.mom=get_Gamma(self.beta)*self.beta
        # self.masses=mass_

    def get_2D_spherical_id(self,dist:str, t0:float, new_theta_len: int = 30) -> dict:


        theta = np.linspace(0.,self.theta_max,new_theta_len)
        masses2d = np.row_stack(
            [self.masses/(new_theta_len-1) for _ in range(new_theta_len-1)]
        )
        v_inf, thetas, masses = reinterpolate_hist2(self.beta, theta, masses2d,
                                                    new_theta_len=None,
                                                    new_vinf_len=None,
                                                    mass_conserving=True)
        masses = masses.T # -> [i_vinf, i_theta]
        thetas = 0.5 * (thetas[1:] + thetas[:-1])
        v_inf  = 0.5 * (v_inf[1:] + v_inf[:-1])
        mask = ((v_inf > 0) & (v_inf < 1) & (np.sum(masses, axis=1) > 0))

        res = {}
        res["mass"] = masses[mask,:] * cgs.solar_m
        res["r"] = np.zeros_like(res["mass"])
        for ith in range(len(thetas)):
            for i_vinf in range(len(v_inf[mask])):
                res["r"][i_vinf,ith] = v_inf[mask][i_vinf] * cgs.c * t0
        res["ek"] = res["mass"] * (v_inf[mask,np.newaxis]*cgs.c)**2 * cgs.solar_m #np.column_stack([masses[mask, i] * cgs.solar_m * (v_inf[mask]*cgs.c)**2 for i in range(len(thetas))])
        mom = MomFromBeta(v_inf[mask]) #v_inf[mask] * get_Gamma(v_inf[mask])
        thetas,mom = np.meshgrid(thetas,mom)
        if not (mom.shape == res["ek"].shape):
            raise ValueError("{} != {}".format(mom.shape,res["ek"].shape))
        res["mom"] = mom
        res["theta"] = thetas
        res["ctheta"] = thetas
        return res


        # masses = np.zeros((len(new_theta_centers), len(self.mom)))
        #
        # for imom in range(len(masses[0,:])):
        #     masses[:,imom] = self.ek[imom] / (cgs.c * cgs.c * self.beta[imom]) / len(masses[:,imom]) / cgs.solar_m # Msun
        #
        # mask = ((self.beta > 0) & (self.beta < 1) & (np.sum(masses, axis=0) > 0))
        #
        # masses = masses[:,mask].T
        # mom = self.mom[mask]
        # vinf = self.beta[mask]
        #
        #
        # # masses = masses.T
        # thetas, mom = np.meshgrid(new_theta_centers, mom)
        #
        # res["ek"] = np.column_stack([masses[:,i] * cgs.solar_m * vinf**2 * cgs.c * cgs.c
        #                              for i in range(len(new_theta_centers))])
        # res["mom"] = mom
        # res["theta"] = thetas
        # res["ctheta"] = thetas
        # res["mass"] = masses*cgs.solar_m
        #
        # res["r"] = np.zeros_like(masses)
        # for ith in range(len(new_theta_centers)):
        #     for i_vinf in range(len(vinf)):
        #         res["r"][i_vinf,ith] = BetFromMom(vinf[i_vinf]) * cgs.c * t0 # TODO THis is theta independent!
        #
        # if not (res["mom"].shape == res["ek"].shape ):
        #     raise ValueError("{} != {}".format(res["mom"].shape,res["ek"].shape))
        # if not (res["ek"].shape == res["mass"].shape ):
        #     raise ValueError("{} != {}".format(res["ek"].shape,res["mass"].shape))
        #
        # return res

    def get_2D_id(self,dist:str, t0:float) -> dict:
        masses, thetas, v_inf = self.masses.T,self.thetas,self.beta
        res = {}
        mask = ((v_inf > 0) & (v_inf < 1) & (np.sum(masses, axis=1) > 0))
        res["mass"] = masses[mask,:] * cgs.solar_m
        res["r"] = np.zeros_like(res["mass"])
        for ith in range(len(thetas)):
            for i_vinf in range(len(v_inf[mask])):
                res["r"][i_vinf,ith] = v_inf[mask][i_vinf] * cgs.c * t0
        res["ek"] = res["mass"] * (v_inf[mask,np.newaxis]*cgs.c)**2 * cgs.solar_m #np.column_stack([masses[mask, i] * cgs.solar_m * (v_inf[mask]*cgs.c)**2 for i in range(len(thetas))])
        mom = MomFromBeta(v_inf[mask]) #v_inf[mask] * get_Gamma(v_inf[mask])
        thetas,mom = np.meshgrid(thetas,mom)
        if not (mom.shape == res["ek"].shape):
            raise ValueError("{} != {}".format(mom.shape,res["ek"].shape))
        res["mom"] = mom
        res["theta"] = thetas
        res["ctheta"] = thetas
        return res


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
def OLD_prepare_grb_ej_id_1d(pars, outfpath, type="pw"):

    o_jet = JetStruct(pars["nlayers_pw"],pars["nlayers_a"])
    if (pars["struct"] == "gaussian"):

        o_jet._gaussian(E_iso_c=pars["Eiso_c"],
                        Gamma0c=pars["Gamma0c"],
                        theta_c=pars["theta_c"],
                        theta_w=pars["theta_w"],
                        M0c=pars["M0c"])

    elif (pars["struct"] == "tophat"):
        o_jet._tophat(E_iso=pars["Eiso_c"],
                      Gamma0=pars["Gamma0c"],
                      theta_h=pars["theta_c"],
                      M0=pars["M0c"],
                      Ye=0,
                      s=0)
    else:
        raise KeyError("Not implemented")

    o_jet.OLD_saveCurrentStructure(outfpath=outfpath,type=type)

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

def OLD_prepare_grb_ej_id_2d(pars, outfpath, type="pw"):

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