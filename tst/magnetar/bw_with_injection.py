from magnetar_evol import Magnetar

import numpy as np
from scipy.integrate import ode
from tqdm import tqdm
from scipy import interpolate
from scipy.optimize import fsolve,leastsq

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, LogNorm
from matplotlib import ticker, cm, rc
rc('text', usetex=True) # $ sudo apt-get install cm-super
rc('font', family='serif')

get_beta = lambda Gamma: np.sqrt(1. - np.power(Gamma, -2.))
get_Gamma = lambda beta: np.sqrt(1. / (1. - np.power(beta, 2.) ))

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

def find_nearest_index(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx
def find_nearest(array, value):
    return array[find_nearest_index(array, value)]

def rho_dlnrho1dR(R, nn, A0, s, R_EJ, R_ISM):
    if not nn is None:
        rho = nn * cgs.mppme
        dlnrho1dR = 0.
    else:
        if R < R_EJ:
            rho = A0 * R_EJ ** (-s) * cgs.mppme
            dlnrho1dR = 0.
        elif R >= R_EJ and R < R_ISM:
            rho = A0 * R ** (-s) * cgs.mppme
            dlnrho1dR = -s / R
        else:
            rho = A0 * R_ISM ** (-s) * cgs.mppme
            # rho = pars.A0 / pars.M0 * pars.R_ISM ** (-pars.s) * mppme
            dlnrho1dR = 0.
    # if scalebyM0: rho = rho / M0
    return (rho, dlnrho1dR)

class EqOpts:

    @staticmethod
    def dthetadr_None(*args):
        return 0.

    @staticmethod
    def dthetadr_Adi(gammaAdi, Gamma, beta, R, theta, aa, useSpread, Rd, thetamax=np.pi / 2.):
        """ basic method """
        vperp = np.sqrt(gammaAdi * (gammaAdi - 1) * (Gamma - 1) / (1 + gammaAdi * (Gamma - 1))) * cgs.c
        return vperp / R / Gamma / beta / cgs.c if (theta < thetamax) & useSpread else 0.

    @staticmethod
    def dthetadr_Adi_Rd(gammaAdi, Gamma, beta, R, theta, aa, useSpread, Rd, thetamax=np.pi / 2.):
        """ basic method that starts working after Rd """
        return EqOpts.dthetadr_Adi(gammaAdi, Gamma, beta, R, theta, aa, useSpread, Rd, thetamax) if (R > Rd) else 0.

    @staticmethod
    def dthetadr_AA(gammaAdi, Gamma, beta, R, theta, aa, useSpread, Rd, thetamax=np.pi / 2.):
        """ Method of Granot & Piran 2012 with 'a' parameter """
        return 1. / (R * Gamma ** (1. + aa) * theta ** (aa)) if (theta < thetamax) & useSpread else 0.

    @staticmethod
    def dthetadr_AA_Rd(gammaAdi, Gamma, beta, R, theta, aa, useSpread, Rd, thetamax=np.pi / 2.):
        """ Method of Granot & Piran 2012 with 'a' parameter that starts working after Rd """
        return EqOpts.dthetadr_AA(gammaAdi, Gamma, beta, R, theta, aa, useSpread, Rd, thetamax) if (R > Rd) else 0.
    #
    # @staticmethod
    # def dthetadr(gammaAdi, Gamma, R, theta, aa=np.nan):
    #     """
    #     Source:
    #     1. / (R * Gamma ** (1. + aa) * theta ** (aa))
    #     is from https://academic.oup.com/mnras/article/421/1/570/990386
    #     vperp / R / Gamma * one_over_beta / cgs.c
    #     is from ???
    #     :param gammaAdi:
    #     :param Gamma:
    #     :param R:
    #     :param one_over_beta:
    #     :param theta:
    #     :param aa:
    #     :return:
    #     """
    #
    #     if theta < np.pi:
    #         if np.isfinite(aa):
    #             return 1. / (R * Gamma ** (1. + aa) * theta ** (aa))
    #         else:
    #             vperp = np.sqrt(gammaAdi * (gammaAdi - 1) * (Gamma - 1) / (1 + gammaAdi * (Gamma - 1))) * cgs.c
    #             one_over_beta = 1. / np.power(1 - np.power(Gamma, -2.), 0.5)
    #             return vperp / R / Gamma * one_over_beta / cgs.c
    #     else:
    #         return 0.

    @staticmethod
    def dmdr(Gamma, RR, thetaE, theta, rho, aa=-1.):
        """
        https://arxiv.org/pdf/1203.5797.pdf
        https://academic.oup.com/mnras/article/421/1/570/990386 (for introduction of 'a' parameter)

        basic form was dmdr = 2 * cgs.pi * R ** 2. * rho * one_minus_costheta / pars["m_scale"]
        """
        # First term: change in swept-up mass due to the change in solid angle
        t1 = 0. if (aa < 0) else (1. / 3.) * np.sin(theta) / (Gamma ** (1. + aa) * theta ** (aa))
        # Second term: change in swept-up mass due to radial expansion
        t2 = (np.cos(thetaE) - np.cos(theta))  # -> always (1 - cos(theta))
        return 2. * np.pi * rho * (t1 + t2) * RR ** 2.

    # --- Adiabatic index (EOS) ---

    @staticmethod
    def gamma_adi_nava(Gamma, beta):
        """ From Nava 2013 paper """
        return (4 + 1 / Gamma) / 3.

    @staticmethod
    def gamma_adi_peer(Gamma, beta):
        """ From Peer 2012 arxiv:1203.5797 """
        mom = Gamma * beta
        theta = mom / 3. * (mom + 1.07 * mom ** 2.) / (1 + mom + 1.07 * mom ** 2.)
        zz = theta / (0.24 + theta)
        gamma_adi = (5. - 1.21937 * zz + 0.18203 * zz ** 2. - 0.96583 * zz ** 3. + 2.32513 * zz ** 4. -
                     2.39332 * zz ** 5. + 1.07136 * zz ** 6.) / 3.
        return gamma_adi

    # --- rho2 (comoving density) ---

    @staticmethod
    def rho2_rel(Gamma, beta, rho, gammaAdi):
        """ Eq. 36 in Rumar + 2014 arXiv:1410.0679 if Gamma >> 1"""
        return 4. * rho * Gamma

    @staticmethod
    def rho2_transrel(Gamma, beta, rho, gammaAdi):
        """ Eq. 36 in Rumar + 2014 arXiv:1410.0679 """
        return (gammaAdi * Gamma + 1.) / (gammaAdi - 1.) * rho

    # --- thickness ---

    @staticmethod
    def shock_thickness(mass, rhoprime, theta, Gamma, R, ncells):
        """ eq. from Johannesson 2006 : astro-ph/0605299 """
        one_min_costheta = (1 - np.cos(theta)/ncells)
        # rhoprime = 4. * rho * Gamma
        delta = mass / (2 * np.pi * one_min_costheta * rhoprime * Gamma * R ** 2)
        return delta

class Nava_fs_rhs:

    def __init__(self):
        pass

    def __call__(self, R, params, pars, Eqs):
        """
        [0] 1 / beta / cgs.c,
        [1] 1 / beta / Gamma / cgs.c,
        [2] dGammadR,
        [3] dEint2dR,
        [4] dthetadR,
        [5] dErad2dR,
        [6] dEsh2dR,
        [7] dEad2dR,
        [8] dM2dR
        :param R:
        :param params:
        :param pars:
        :return:
        """

        tburst = params[0]
        tcomoving = params[1]
        Gamma = params[2]
        Eint2 = params[3]
        theta = params[4]
        M2 = params[8]
        #
        M0 = pars["M0"]
        Gamma0 = pars["Gamma0"]
        theta0 = pars["theta0"]
        rho = pars["rho"] / M0
        dlnrho1dR = pars["dlnrho1dR"]
        #
        # if Eint2 < 0:
        #     beta = np.sqrt(1. - Gamma ** -2)
        #     print("ERROR! Eint2 < 0 Gamma < 1 [Gamma:{} beta:{} Eint2:{} M2:{} ] Resetting to 0".format(
        #         Gamma, beta, Eint2, M2,
        #     ))
        #     Eint2 = 0
        # if Gamma < 1:
        #     Gamma = 1.0001
        #     print("ERROR Gamma < 1 [Gamma:{} Eint2:{} M2:{} ] Resetting to 1.0001".format(
        #         Gamma, Eint2,  M2,
        #     ))
        #     # raise ValueError("Gamma < 1")
        # if Gamma > Gamma0:
        #     print("ERROR! Gamma({}) > Gamma0({}) -- resetting to Gamma0"
        #           .format(Gamma, Gamma0))
        #     Gamma = Gamma0


        #  one_over_beta = 1. / np.power(1 - np.power(Gamma, -2.), 0.5)
        beta = np.sqrt(1. - np.power(Gamma,-2))
        # beta0 = np.sqrt(1. - Gamma0 ** -2)

        gammaAdi = Eqs["eq_gammaAdi"](Gamma, beta)
        # gammaAdi = self.gammaAdi(Gamma, beta)  # (4 + 1 / Gamma) / 3.

        # one_minus_costheta = EquationsNava.one_minus_costheta(theta)
        # one_minus_costheta0 = 1. - np.cos(theta0)

        # Spreading
        # dthetadR = self.dthetadr(gammaAdi, Gamma, R, theta, pars["aa"]) * int(pars["useSpread"])
        dthetadR = Eqs["eq_dthetadr"](gammaAdi, Gamma, beta, R, theta, pars["aa"],
                                      pars["useSpread"], pars["Rd"], pars["thetaMax"])
        # Densities and volumes

        # rho2 = 4. * Gamma * rho
        # V2 = M2 / rho2
        dM2dR = EqOpts.dmdr(Gamma, R, pars["thetaE"], theta, rho, aa=pars["aa"]) / pars["ncells"] # 2 * cgs.pi * R ** 2. * rho * one_minus_costheta / (pars["m_scale"])

        # # # # dGammadR
        dGammadR = self.dGammadR_fs(Gamma, gammaAdi, dlnrho1dR, M2, dM2dR, Eint2)

        if dGammadR > 0:
            if Gamma > 0.95 * Gamma0:
                # raise ValueError("Gamma > 0.95 Gamma0 after RSISING")
                dGammadR = 0.

        # # # # Energies # # # #
        dEsh2dR = (Gamma - 1.) * dM2dR  # Shocked energy

        # --- Expansion energy
        dlnV2dR = dM2dR / M2 - dlnrho1dR - dGammadR / Gamma
        if pars["adiabLoss"]:
            dEad2dR = -(gammaAdi - 1.) * Eint2 * dlnV2dR
        else:
            dEad2dR = 0.

        # -- Radiative losses
        dErad2dR = pars["epsilon_e_rad"] * dEsh2dR

        ### -- Energy equation
        dEint2dR = dEsh2dR + dEad2dR - dErad2dR  # / (pars.M0 * c ** 2)


        return np.array([1 / beta / cgs.c,
                         1 / beta / Gamma / cgs.c,
                         dGammadR,
                         dEint2dR,
                         dthetadR,
                         dErad2dR,
                         dEsh2dR,
                         dEad2dR,
                         dM2dR
                         ])

    # @staticmethod
    # def gammaAdi(Gamma, beta):
    #     """ From Peer 2012 arxiv:1203.5797 """
    #     mom = Gamma * beta
    #     theta = mom / 3. * (mom + 1.07 * mom ** 2.) / (1 + mom + 1.07 * mom ** 2.)
    #     zz = theta / (0.24 + theta)
    #     gamma_adi = (5. - 1.21937 * zz + 0.18203 * zz ** 2. - 0.96583 * zz ** 3. + 2.32513 * zz ** 4. -
    #             2.39332 * zz ** 5. + 1.07136 * zz ** 6.) / 3.
    #     return gamma_adi

    # @staticmethod
    # def gammaAdi(Gamma):
    #     return (4 + 1 / Gamma) / 3.

    # @staticmethod
    # def dthetadr(gammaAdi, Gamma, R, theta, aa=np.nan):
    #     """
    #     Source:
    #     1. / (R * Gamma ** (1. + aa) * theta ** (aa))
    #     is from https://academic.oup.com/mnras/article/421/1/570/990386
    #     vperp / R / Gamma * one_over_beta / cgs.c
    #     is from ???
    #     :param gammaAdi:
    #     :param Gamma:
    #     :param R:
    #     :param one_over_beta:
    #     :param theta:
    #     :param aa:
    #     :return:
    #     """
    #
    #     if theta < np.pi:
    #         if np.isfinite(aa):
    #             return 1. / (R * Gamma ** (1. + aa) * theta ** (aa))
    #         else:
    #             vperp = np.sqrt(gammaAdi * (gammaAdi - 1) * (Gamma - 1) / (1 + gammaAdi * (Gamma - 1))) * cgs.c
    #             one_over_beta = 1. / np.power(1 - np.power(Gamma, -2.), 0.5)
    #             return vperp / R / Gamma * one_over_beta / cgs.c
    #     else:
    #         return 0.

    @staticmethod
    def GammaEff(Gamma, gammaAdi):
        return (gammaAdi * Gamma ** 2. - gammaAdi + 1.) / Gamma

    @staticmethod
    def dGammaEffdGamma(Gamma, gammaAdi):
        return (gammaAdi * Gamma ** 2. + gammaAdi - 1.) / Gamma ** 2
        #return 4. / 3. + 1. / Gamma ** 2. / 3. + 2. / Gamma ** 3. / 3.

    @staticmethod
    def dGammadR_fs(Gamma, gammaAdi, dlnrho1dR, M2, dM2dR, Eint2):

        GammaEff = Nava_fs_rhs.GammaEff(Gamma, gammaAdi)#, gammaAdi) #(gammaAdi * Gamma ** 2. - gammaAdi + 1.) / Gamma
        dGammaEffdGamma = Nava_fs_rhs.dGammaEffdGamma(Gamma, gammaAdi) #4. / 3. + 1. / Gamma ** 2. / 3. + 2. / Gamma ** 3. / 3.

        f_2 = GammaEff * (gammaAdi - 1.) * Eint2 / Gamma
        h_2 = GammaEff * (gammaAdi - 1.) * Eint2 * (dM2dR / M2 - dlnrho1dR)

        dGammadR = -((Gamma - 1.) * (GammaEff + 1.) * dM2dR - h_2) / ((1. + M2) + Eint2 * dGammaEffdGamma + f_2)

        return dGammadR

class Driver_t:

    def __init__(self, tstart, ode_rtol=1e-5, ode_nsteps=1000):

        self.tstart = tstart
        self.odeinstance = ode(self.rhs)
        self.odeinstance.set_integrator("dop853", rtol=ode_rtol, nsteps=ode_nsteps, first_step=tstart)
        self.odeinstance.set_f_params(self.pars_ode_rhs, self.eqs_ode_rhs)
        self.odeinstance.set_initial_value(self.initial_data, tstart * 0.9999)

    @classmethod
    def from_obj(cls, shell, **kwargs):
        return cls(
            shell.E0,
            shell.Gamma0,
            shell.thetaE,
            shell.theta0,
            shell.M0,
            shell.Rstart,
            shell.rho0,
            shell.dlnrho1dR_0,
            **kwargs
        )

    def update_ode_pars(self, **kwargs):

        if len(kwargs.keys()) > 0:
            for key in kwargs.keys():
                # if key.__contains__("eq_"):
                #     self.eqs_ode_rhs[key] = kwargs[key]
                # else:
                self.pars_ode_rhs[key] = kwargs[key]

        self.odeinstance.set_f_params(self.pars_ode_rhs, self.eqs_ode_rhs)

    def evolove(self, x, rho, dlnrho1dR, deinjdt):
        # Integrate ODEs
        self.update_ode_pars(rho = rho, dlnrho1dR = dlnrho1dR, deinjdt=deinjdt)
        self.dynamics = np.vstack((self.dynamics, np.zeros(len(self.all_v_ns))))
        i_res = np.copy(self.odeinstance.integrate(x))
        # print(self.odeinstance.set_solout((x, np.array(i_res))))
        if (x > self.tstart):
            assert np.isfinite(self.odeinstance.t)
            assert self.odeinstance.t > self.tstart
        self.dynamics[-1, :len(self.initial_data)] = self.apply_units(i_res)
        self._additional_quantities()
    # "eq_gammaAdi": EqOpts.gamma_adi_peer, "eq_rhoprime": EqOpts.rho2_transrel,
    def _additional_quantities(self):
        self.set_last("beta", get_beta(self.get_last("Gamma")))
        self.set_last("rho", self.pars_ode_rhs["rho"]) # the one used for evolution
        self.set_last("gammaAdi", self.kwargs["eq_gammaAdi"](self.get_last("Gamma"), self.get_last("beta")))
        self.set_last("rho2", self.kwargs["eq_rhoprime"](self.get_last("Gamma"),
                                                         self.get_last("beta"),
                                                         self.get_last("rho"),
                                                         self.get_last("gammaAdi")))
        #self.get_rhoprime(self.get_last("rho"), self.get_last("Gamma")))
        # self.set_last("R", self.odeinstance.t)
        self.set_last("tburst", self.odeinstance.t)
        self.set_last("tt", self.get_elapsed_time(self.get("R"),
                                                  self.get("Gamma"),
                                                  self.get("theta")))
        # mass, rhoprime, theta, Gamma, R
        self.set_last("thickness", EqOpts.shock_thickness(self.get_last("M2"),
                                                          self.get_last("rho2"),
                                                          self.get_last("theta"),
                                                          self.get_last("Gamma"),
                                                          self.get_last("R"),
                                                          self.kwargs["ncells"]))
        self.set_last("U_e", self.get_U_e())
        if "beta_ej" in self.all_v_ns:
            self.set_last("beta_ej", get_beta(self.get_last("Gamma_ej")))
            self.pars_ode_rhs["Rj_arr"].append(self.get_last("R"))
            self.pars_ode_rhs["tburst_arr"].append(self.get_last("tburst"))
            self.pars_ode_rhs["ctheta_j_arr"].append( self.kwargs["ctheta_j"] + 0.5 * (2.0 * self.get_last("theta") - 2.0 * self.kwargs["theta_w_j"]) )
            self.pars_ode_rhs["rho2_j_arr"].append( self.get_last("rho2") )
            if(self.pars_ode_rhs["t_ent"] > 0): self.pars_ode_rhs["is_entered"] = True

        #
        # self.set("rho", beta(self.get("Gamma")[-1]), -1)
        # self.set("rho", self.pars_ode_rhs["rho"], -1)
        # self.set("R", self.odeinstance.t, -1)
        # self.set("tt", self.get_elapsed_time(self.get("R"), self.get("Gamma"), self.get("theta")), -1)
        # delta = self.get_shock_front_thickness(self.get("R")[-1], self.get("rho")[-1], self.get("Gamma")[-1],
        #                                        self.get("theta")[-1], self.get("M2")[-1])
        # self.set("delta", delta, -1)
        # self.set("U_e", self.get_U_e(), -1)

    def _set_ode_inits(self, v_ns_dict):
        vals = [0. for v_n in self.v_ns_init_vals]
        for key in v_ns_dict.keys():
            if key in v_ns_dict:
                i = self.v_ns_init_vals.index(key)
            else:
                raise NameError("Key : {} : is not in init.data dict:"
                                "{}".format(key, v_ns_dict.keys()))
            vals[i] = v_ns_dict[key]
        assert len(vals) == len(self.v_ns_init_vals)
        return (vals)

    def i_nv(self, v_n):
        return self.all_v_ns.index(v_n)

    def set(self, v_n, array):
        self.dynamics[:, self.i_nv(v_n)] = array

    def set_last(self, v_n, value):
        self.dynamics[-1, self.i_nv(v_n)] = value

    def set_first(self, v_n, value):
        self.dynamics[0, self.i_nv(v_n)] = value

    def get(self, v_n):
        return self.dynamics[:, self.i_nv(v_n)]

    def get_last(self, v_n):
        return self.dynamics[-1, self.i_nv(v_n)]

    def get_first(self, v_n):
        return self.dynamics[0, self.i_nv(v_n)]

    def get_init_val(self, v_n):
        if self.dynamics.ndim == 2:
            return self.dynamics[0, self.i_nv(v_n)]
        else:
            return self.dynamics[self.i_nv(v_n)]

    def init_elapsed_time(self, Rstart, beta0, Gamma0, **kwargs):
        if kwargs["useSpread"]:
            dThetadr0 = 0
            tt = Rstart / (cgs.c * beta0) * (np.sqrt(1. + Rstart ** 2. * dThetadr0 ** 2.)) - Rstart / cgs.c
        else:
            tt = Rstart / (cgs.c * Gamma0 ** 2. * beta0 * (1. + beta0))
        return tt

    def get_elapsed_time(self, Rs, Gammas, thetas):
        beta = np.sqrt(1. - np.power(Gammas, -2))
        if self.kwargs["useSpread"]:
            dThetadr = np.concatenate([np.zeros(1), np.diff(thetas) / np.diff(Rs)])
            integrand = 1. / (cgs.c * beta) * np.sqrt(1. + Rs ** 2. * dThetadr ** 2.) - 1. / (cgs.c)
        else:
            integrand = 1. / (cgs.c * Gammas ** 2. * beta * (1. + beta))
        tti = np.trapz(integrand, Rs) + self.get_init_val("tt")
        return tti

    # def get_shock_front_thickness(self, R, rhoprime, Gamma, theta, m):
    #     """
    #     assuming compresion ration '4' by default. as rho2 = 4 * Gamma * rho is embedded here
    #
    #     delta = M2 / (8 * pi * (1-cos(theta)) * Gamma^2 * rho * r^2)
    #     where rho2 = 4 * rho * Gamma
    #     which obays that M2 / (4 * pi) = rho2 * Gamma * r^2 * delta
    #
    #     """
    #
    #     one_min_costheta = 2.#1. - np.cos(theta)/self.kwargs["ncells"]
    #     # delta = m / (8. * np.pi * one_min_costheta * Gamma ** 2. * rho * R ** 2.)
    #     delta = m / (2 * np.pi * one_min_costheta * rhoprime * Gamma * R ** 2)
    #     return delta

    # def get_rhoprime(self, rho, Gamma):
    #     """ 4. is the compression ratio """
    #     return 4. * rho * Gamma
    #     return

class SedovDimensionLess():
    def __init__(self,
                 gamma=4/3.,
                 nu=3,
                 w=0.,
                 epsilon=1e-50):
        self._nDim = nu
        self._w = w

        # Constants for the parametric equations:
        self.w1 = (3 * nu - 2 + gamma * (2 - nu)) / (gamma + 1.)
        self.w2 = (2. * (gamma - 1) + nu) / gamma
        self.w3 = nu * (2. - gamma)

        self.b0 = 1. / (nu * gamma - nu + 2)
        self.b2 = (gamma - 1.) / (gamma * (self.w2 - w))
        self.b3 = (nu - w) / (float(gamma) * (self.w2 - w))
        self.b5 = (2. * nu - w * (gamma + 1)) / (self.w3 - w)
        self.b6 = 2. / (nu + 2 - w)
        self.b1 = self.b2 + (gamma + 1.) * self.b0 - self.b6
        self.b4 = self.b1 * (nu - w) * (nu + 2. - w) / (self.w3 - w)
        self.b7 = w * self.b6
        self.b8 = nu * self.b6

        self.c0 = 2 * (nu - 1) * np.pi + (nu - 2) * (
                nu - 3)  # simple interpolation of correct function (only for nu=1,2,3)
        self.c5 = 2. / (gamma - 1)
        self.c6 = (gamma + 1) / 2.
        self.c1 = self.c5 * gamma
        self.c2 = self.c6 / gamma
        self.c3 = (nu * gamma - nu + 2.) / ((self.w1 - w) * self.c6)
        self.c4 = (nu + 2. - w) * self.b0 * self.c6

        # Characterize the solution
        f_min = self.c2 if self.w1 > w else self.c6

        f = np.logspace(np.log10(f_min), 0, np.int64(1e5))

        # Sort the etas for our interpolation function
        eta = self.parametrized_eta(f)
        f = f[eta.argsort()]
        eta.sort()

        d = self.parametrized_d(f)
        p = self.parametrized_p(f)
        v = self.parametrized_v(f)

        # If min(eta) != 0 then all values for eta < min(eta) = 0
        if eta[0] > 0:
            e01 = [0., eta[0] * (1 - 1e-10)]
            d01 = [0., 0]
            p01 = [0., 0]
            v01 = [0., 0]

            eta = np.concatenate([np.array(e01), eta])
            d = np.concatenate([np.array(d01), d])
            p = np.concatenate([np.array(p01), p])
            v = np.concatenate([np.array(v01), v])

        # Set up our interpolation functions
        self._d = interpolate.interp1d(eta, d, bounds_error=False, fill_value=np.nan, copy=False)  # 1. / self._rho1)
        self._p = interpolate.interp1d(eta, p, bounds_error=False, fill_value=np.nan, copy=False)  # 0.)
        self._v = interpolate.interp1d(eta, v, bounds_error=False, fill_value=np.nan, copy=False)  # 0.)

        # Finally Calculate the normalization of R_s:
        integral = eta ** (nu - 1) * (d * v ** 2 + p)
        integral = 0.5 * (integral[1:] + integral[:-1])
        d_eta = (eta[1:] - eta[:-1])

        # calculate integral and multiply by factor
        alpha = (integral * d_eta).sum() * (8 * self.c0) / ((gamma ** 2 - 1.) * (nu + 2. - w) ** 2)
        self._c = (1. / alpha) ** (1. / (nu + 2 - w))

    def parametrized_eta(self, var):
        return (var ** -self.b6) * ((self.c1 * (var - self.c2)) ** self.b2) * (
                (self.c3 * (self.c4 - var)) ** (-self.b1))

    def parametrized_d(self, var):
        return (var ** -self.b7) * \
            ((self.c1 * (var - self.c2)) ** (self.b3 - self._w * self.b2)) * \
            ((self.c3 * (self.c4 - var)) ** (self.b4 + self._w * self.b1)) * \
            ((self.c5 * (self.c6 - var)) ** -self.b5)

    def parametrized_p(self, var):
        return (var ** self.b8) * ((self.c3 * (self.c4 - var)) ** (self.b4 + (self._w - 2) * self.b1)) * \
            ((self.c5 * (self.c6 - var)) ** (1 - self.b5))

    def parametrized_v(self, var):
        return self.parametrized_eta(var) * var


    def rho_profile(self, r, r_shock, rho2):
        # density at radius r and time t
        eta = r / r_shock
        return rho2 * self._d(eta)

class SedovSolution():
    """
    see: [The Sedov self-similar point blast solutions in nonuniform media](
    https://link.springer.com/content/pdf/10.1007/BF01414626.pdf)
    rho0 = A*r**(-w)
    R_s = ((e * t**2)/(alpha * A))**(1/(nu + 2 - w))
    """



    def __init__(self,
                 e,
                 rho,
                 gamma=4/3.,
                 nu=3,
                 w=0.,
                 epsilon=1e-50):

        # w = 0 --> uniform background

        if not any(nu == np.array([1, 2, 3])):
            raise ValueError("nu (dimension of problem) need to be 1, 2 or 3!")

        self._epsilon = epsilon
        self._e = e
        self._gamma = gamma

        self._rho0 = rho # background density
        self._rho1 = ((gamma + 1.) / (gamma - 1.)) * rho # shock dens (Eq.3)

        self._nDim = nu
        self._w = w

        # Constants for the parametric equations:
        self.w1 = (3 * nu - 2 + gamma * (2 - nu)) / (gamma + 1.)
        self.w2 = (2. * (gamma - 1) + nu) / gamma
        self.w3 = nu * (2. - gamma)

        self.b0 = 1. / (nu * gamma - nu + 2)
        self.b2 = (gamma - 1.) / (gamma * (self.w2 - w))
        self.b3 = (nu - w) / (float(gamma) * (self.w2 - w))
        self.b5 = (2. * nu - w * (gamma + 1)) / (self.w3 - w)
        self.b6 = 2. / (nu + 2 - w)
        self.b1 = self.b2 + (gamma + 1.) * self.b0 - self.b6
        self.b4 = self.b1 * (nu - w) * (nu + 2. - w) / (self.w3 - w)
        self.b7 = w * self.b6
        self.b8 = nu * self.b6

        self.c0 = 2 * (nu - 1) * np.pi + (nu - 2) * ( nu - 3)  # simple interpolation of correct function (only for nu=1,2,3)
        self.c5 = 2. / (gamma - 1)
        self.c6 = (gamma + 1) / 2.
        self.c1 = self.c5 * gamma
        self.c2 = self.c6 / gamma
        self.c3 = (nu * gamma - nu + 2.) / ((self.w1 - w) * self.c6)
        self.c4 = (nu + 2. - w) * self.b0 * self.c6

        # Characterize the solution
        f_min = self.c2 if self.w1 > w else self.c6

        f = np.logspace(np.log10(f_min), 0, np.int64(1e5))

        # Sort the etas for our interpolation function
        eta = self.parametrized_eta(f)
        f = f[eta.argsort()]
        eta.sort()

        d = self.parametrized_d(f)
        p = self.parametrized_p(f)
        v = self.parametrized_v(f)

        # If min(eta) != 0 then all values for eta < min(eta) = 0
        if eta[0] > 0:
            e01 = [0., eta[0] * (1 - 1e-10)]
            d01 = [0., 0]
            p01 = [0., 0]
            v01 = [0., 0]

            eta = np.concatenate([np.array(e01), eta])
            d = np.concatenate([np.array(d01), d])
            p = np.concatenate([np.array(p01), p])
            v = np.concatenate([np.array(v01), v])

        # Set up our interpolation functions
        self._d = interpolate.interp1d(eta, d, bounds_error=False, fill_value=np.nan)#1. / self._rho1)
        self._p = interpolate.interp1d(eta, p, bounds_error=False, fill_value=np.nan)#0.)
        self._v = interpolate.interp1d(eta, v, bounds_error=False, fill_value=np.nan)#0.)

        # Finally Calculate the normalization of R_s:
        integral = eta ** (nu - 1) * (d * v ** 2 + p)
        integral = 0.5 * (integral[1:] + integral[:-1])
        d_eta = (eta[1:] - eta[:-1])

        # calculate integral and multiply by factor
        alpha = (integral * d_eta).sum() * (8 * self.c0) / ((gamma ** 2 - 1.) * (nu + 2. - w) ** 2)
        self._c = (1. / alpha) ** (1. / (nu + 2 - w))

    def parametrized_eta(self, var):
        return (var ** -self.b6) * ((self.c1 * (var - self.c2)) ** self.b2) * (
                (self.c3 * (self.c4 - var)) ** (-self.b1))

    def parametrized_d(self, var):
        return (var ** -self.b7) * ((self.c1 * (var - self.c2)) ** (self.b3 - self._w * self.b2)) * \
            ((self.c3 * (self.c4 - var)) ** (self.b4 + self._w * self.b1)) * (
                    (self.c5 * (self.c6 - var)) ** -self.b5)

    def parametrized_p(self, var):
        return (var ** self.b8) * ((self.c3 * (self.c4 - var)) ** (self.b4 + (self._w - 2) * self.b1)) * \
            ((self.c5 * (self.c6 - var)) ** (1 - self.b5))

    def parametrized_v(self, var):
        return self.parametrized_eta(var) * var

    # Shock properties
    def shock_radius(self, t):
        # outer radius at time t
        t = np.maximum(t, self._epsilon)
        return self._c * (self.e * t ** 2 / self.rho0) ** (1. / (self._nDim + 2 - self._w))

    def shock_velocity(self, t):
        # velocity of the shock wave
        t = np.maximum(t, self._epsilon)
        return (2. / (self._nDim + 2 - self._w)) * self.shock_radius(t) / t

    def post_shock_pressure(self, t):
        # post shock pressure
        return (2. / (self.gamma + 1)) * self.rho0 * self.shock_velocity(t) ** 2

    @property
    def post_shock_density(self, t=0):
        # post shock density
        return self._rho1

    def rho(self, r, t):
        # density at radius r and time t
        eta = r / self.shock_radius(t)
        return self.post_shock_density * self._d(eta)

    def pressure(self, r, t):
        # pressure at radius r and time t
        eta = r / self.shock_radius(t)
        return self.post_shock_pressure(t) * self._p(eta)

    def velocity(self, r, t):
        # velocity at radius r, and time t
        eta = r / self.shock_radius(t)
        return self._v(eta) * (2 / (self.gamma + 1)) * self.shock_velocity(t)

    def internal_energy(self, r, t):
        # internal energy at radius r and time t
        return self.pressure(r, t) / (self.rho(r, t) * (self.gamma - 1))

    def entropy(self, r, t):
        # entropy at radius, r, and time, t
        return self.pressure(r, t) / self.rho(r, t) ** self.gamma

    # Other properties
    @property
    def e(self):
        # total energy
        return self._e

    @property
    def gamma(self):
        # ratio of specific heats
        return self._gamma

    @property
    def rho0(self):
        # background density
        return self._rho0


class Nava_fs_rhs_t2:

    def __init__(self):
        pass

    def ej_rhs(self, pars, params, ii, ieq, einj):
        R_ej_0 = params[0 + ieq]
        tcomoving_ej_0 = params[1 + ieq]
        Gamma_ej_0 = params[2 + ieq]
        if (Gamma_ej_0 < 1):
            Gamma_ej_0 = 1.
        tmp_0 = np.power(Gamma_ej_0, -2.)
        beta_ej_0 = np.sqrt(1. - tmp_0)
        Eint2_ej_0 = params[3 + ieq]
        theta_ej_0 = params[4 + ieq]
        M2_ej_0 = params[8 + ieq]
        #
        M0_ej_0 = pars[f"M0_ej_{ii}"]
        Gamma0_ej_0 = pars[f"Gamma0_ej_{ii}"]
        theta0_ej_0 = pars[f"theta0_ej_{ii}"]

        # if (Eint2_ej_0 < Eint2):
        #     print("problem")

        ctheta_ej_0 = pars[f"ctheta_ej_{ii}"] + 0.5 * (2.0 * theta_ej_0 - 2.0 * pars[f"theta_w_ej_{ii}"])
        gammaAdi_ej_0 = EqOpts.gamma_adi_peer(Gamma_ej_0, beta_ej_0)

        # rho_ej_default_0 = pars["rho_ej_0"] / M0_ej_0
        # dlnrho1dR_ej_default_0 = pars["dlnrho1dR_ej_0"]
        # if (pars["use_ctheta_lim"]):
        #     if(pars["use_st_rho"]):
        #         if ((ctheta > ctheta_ej_0)&(R > R_ej_0)):
        #             if (not pars["is_entered"]) :
        #                 pars["t_ent"] = t
        #                 pars["r_ej_ent"] = R_ej_0
        #             rho_ej = 1e-90#rho_ej_default_0 * np.exp( - (R_ej_0 - pars["r_ej_ent"]) / R_ej_0 )
        #             dlnrho1dR_ej = 0.#rho_ej_default_0 * ( -1. * R_ej_0**-3 * pars["r_ej_ent"] * np.exp(pars["r_ej_ent"]/R_ej_0 - 1) )
        #
        #             if ((~np.isfinite(rho_ej))&(rho_ej < 0)):
        #                 raise ValueError()
        #             if (~np.isfinite(dlnrho1dR_ej)):
        #                 raise ValueError()
        #
        #             # o_sedov = SedovDimensionLess(gamma=gammaAdi_ej_0)
        #             # inside the jet cone & behind the head
        #             # From the expanding jet head (decaying)
        #
        #             # arc_jet = ctheta * R
        #             # arc_ej = ctheta_ej_0 * R_ej_0
        #             # assert arc_ej < arc_jet
        #             # rho_from_arc = pars["o_sedov"].rho_profile(arc_ej, arc_jet, rho2) / M0_ej_0
        #             #
        #             #
        #             #
        #             # rho_from_arc_m1 = pars["o_sedov"].rho_profile(arc_ej - 0.01 * arc_ej, arc_jet, rho2) / M0_ej_0
        #             # dlnrho1dR_ej_from_arc = (1. / arc_ej) * (rho_from_arc - rho_from_arc_m1) / (arc_ej - arc_ej * 0.01)
        #             # assert dlnrho1dR_ej_from_arc < 1
        #             #
        #             # rho_ej = rho_from_arc if (rho_from_arc < rho_ej_default_0) else rho_ej_default_0
        #             # dlnrho1dR_ej = dlnrho1dR_ej_from_arc if (rho_from_arc < rho_ej_default_0) else dlnrho1dR_ej_default_0
        #
        #             #
        #             # # From jet head itself (rising)
        #             rho_ej_behind = pars["o_sedov"].rho_profile(R_ej_0, R, rho2) / M0_ej_0
        #             rho_ej_behind_im1 = pars["o_sedov"].rho_profile(R_ej_0 - 0.01 * R_ej_0, R, rho2) / M0_ej_0
        #             dlnrho1dR_ej_behind = (1. / R_ej_0) * (rho_ej_behind - rho_ej_behind_im1) / (R_ej_0 - 0.01 * R_ej_0)
        #             rho_ej = rho_ej_behind
        #             dlnrho1dR_ej = dlnrho1dR_ej_behind
        #             # # chose the largest
        #             # if ((rho_from_arc < rho_ej_behind)&(rho_ej_default_0 > rho_from_arc)):
        #             #     rho_ej_final = rho_from_arc
        #             #     dlnrho1dR_ej_final = dlnrho1dR_ej_from_arc
        #             # # elif((rho_from_arc > rho_ej_behind)&(pars["rho_ej"] / M0_ej_0 <= rho_from_arc)):
        #             # #     rho_ej_final = rho_from_arc
        #             # #     dlnrho1dR_ej_final = dlnrho1dR_ej_from_arc
        #             # else:
        #             #     rho_ej_final = rho_ej_behind
        #             #     dlnrho1dR_ej_final = dlnrho1dR_ej_behind
        #             #
        #             # if (rho_ej_behind < rho_ej_default_0):
        #             #     rho_ej = rho_ej_default_0
        #             #     dlnrho1dR_ej = dlnrho1dR_ej_default_0
        #             # else:
        #             #     rho_ej = rho_ej_final
        #             #     dlnrho1dR_ej = dlnrho1dR_ej_final
        #         elif ((ctheta > ctheta_ej_0) & (R < R_ej_0)):
        #             # inside & above
        #             rho_ej = rho_ej_default_0
        #             dlnrho1dR_ej = dlnrho1dR_ej_default_0
        #         elif((ctheta < ctheta_ej_0)):
        #             # outside the jet cone
        #             rho_ej = rho_ej_default_0
        #             dlnrho1dR_ej = dlnrho1dR_ej_default_0
        #         else:
        #             rho_ej = rho_ej_default_0
        #             dlnrho1dR_ej = dlnrho1dR_ej_default_0
        #
        #
        #         # if (ctheta > ctheta_ej_0):
        #         #     o_sedov = SedovDimensionLess(gamma=gammaAdi_ej_0)
        #         #     # assert R > R_ej_0
        #         #     # print(pars["ctheta_j_arr"])
        #         #     # idx = np.argmax( pars["ctheta_j_arr"] > ctheta_ej_0 ) + 1
        #         #     assert pars["Rj_arr"][-1] > R_ej_0
        #         #     idx = find_nearest_index(np.asarray(pars["Rj_arr"]),R_ej_0)+1#np.argmax( pars["Rj_arr"] < R_ej_0 ) - 1
        #         #     # lendth of the 'arc' of the jet when jet had the same radius as currently ejecta has
        #         #     arc = pars["ctheta_j_arr"][idx] * pars["Rj_arr"][idx]
        #         #     # current arc of the ejecta
        #         #     arc_ej = ctheta_ej_0 * R_ej_0
        #         #     if not ( arc > arc_ej ):
        #         #         raise ValueError("this is essencially ctheta > ctheta_ej_0 requirement not fulfilled "
        #         #                          "R[idx]={:.2e} Rej={:.2e} Arc={:.2e} Arc_ej={:.2e}".format( pars["Rj_arr"][idx], R_ej_0, arc, arc_ej ))
        #         #     # rho2 of the jet when it had radius that ejecta now has
        #         #     rho2_ = pars["rho2_j_arr"][idx]
        #         #     # Sedov profile for laterally moving piece... yeah, I know
        #         #     rho_from_arc = o_sedov.rho_profile(arc_ej, arc, rho2_) / M0_ej_0
        #         #     rho_from_arc_m1 = o_sedov.rho_profile(arc_ej-0.1*arc_ej, arc, rho2_) / M0_ej_0
        #         #     dlnrho1dR_ej_from_arc= (1. / arc_ej) * (rho_from_arc - rho_from_arc_m1) / (arc_ej - arc_ej-0.1)
        #         #     assert dlnrho1dR_ej_from_arc < 1
        #         #
        #         #     rho_ej_behind = o_sedov.rho_profile(R_ej_0, R, rho2) / M0_ej_0
        #         #     rho_ej_behind_im1 = o_sedov.rho_profile(R_ej_0 - 0.01 * R_ej_0, R, rho2) / M0_ej_0
        #         #     dlnrho1dR_ej_behind = (1. / R_ej_0) * (rho_ej_behind - rho_ej_behind_im1) / (R_ej_0 - 0.01 * R_ej_0)
        #         #
        #         #     if ((rho_from_arc > rho_ej_behind)&(rho_ej_default_0 > rho_from_arc)):
        #         #         rho_ej_final = rho_from_arc
        #         #         dlnrho1dR_ej_final = dlnrho1dR_ej_from_arc
        #         #     # elif((rho_from_arc > rho_ej_behind)&(pars["rho_ej"] / M0_ej_0 <= rho_from_arc)):
        #         #     #     rho_ej_final = rho_from_arc
        #         #     #     dlnrho1dR_ej_final = dlnrho1dR_ej_from_arc
        #         #     else:
        #         #         rho_ej_final = rho_ej_behind
        #         #         dlnrho1dR_ej_final = dlnrho1dR_ej_behind
        #         #
        #         #     if (rho_ej_behind < rho_ej_default_0):
        #         #         rho_ej = rho_ej_default_0
        #         #         dlnrho1dR_ej = dlnrho1dR_ej_default_0
        #         #     else:
        #         #         rho_ej = rho_ej_final
        #         #         dlnrho1dR_ej = dlnrho1dR_ej_final
        #         # else:
        #         #     rho_ej = rho_ej_default_0
        #         #     dlnrho1dR_ej = dlnrho1dR_ej_default_0
        #     else:
        #
        #         if ((ctheta > ctheta_ej_0)&(R > R_ej_0)):
        #             rho_ej = 0.
        #             dlnrho1dR_ej = 0.
        #         else:
        #             rho_ej = pars["rho_ej"] / M0_ej_0
        #             dlnrho1dR_ej = pars["dlnrho1dR_ej"]
        # else:
        #     if(pars["use_st_rho"]):
        #         o_sedov = SedovDimensionLess(gamma=gammaAdi_ej_0)
        #         rho_ej = o_sedov.rho_profile(R_ej_0, R, rho2) / M0_ej_0
        #         rho_ej_ = o_sedov.rho_profile(R_ej_0 - 0.01 * R_ej_0, R, rho2) / M0_ej_0
        #         dlnrho1dR_ej = (1. / R_ej_0) * (rho_ej - rho_ej_) / (R_ej_0 - 0.01 * R_ej_0)
        #         if (rho_ej < pars["rho_ej"] / M0_ej_0):
        #             rho_ej = pars["rho_ej"] / M0_ej_0
        #             dlnrho1dR_ej = pars["dlnrho1dR_ej"]
        #         else:
        #             pass
        #     else:
        #         rho_ej = pars["rho_ej"] / M0_ej_0
        #         dlnrho1dR_ej = pars["dlnrho1dR_ej"]
        #
        rho_ej_0 = pars[f"rho_ej_{ii}"] / M0_ej_0
        dlnrho1dR_ej_0 = pars[f"dlnrho1dR_ej_{ii}"]


        if (not np.isfinite(Gamma_ej_0)):
            raise ValueError("nan in Gamma")
        if ((not np.isfinite(M2_ej_0)) | (M2_ej_0 < 0.)):
            raise ValueError("nan in M2_ej_0")
        if (not np.isfinite(beta_ej_0)):
            raise ValueError("nan in beta")

        gammaAdi_ej_0 = EqOpts.gamma_adi_peer(Gamma_ej_0, beta_ej_0)
        dthetadR_ej_0 = EqOpts.dthetadr_None(gammaAdi_ej_0, Gamma_ej_0, beta_ej_0, R_ej_0, theta_ej_0,
                                             pars["aa_ej"], pars["useSpread_ej"], pars["Rd_ej"], pars["thetaMax_ej"])
        dM2dR_ej_0 = EqOpts.dmdr(Gamma_ej_0, R_ej_0, pars[f"thetaE_ej_{ii}"], theta_ej_0, rho_ej_0, aa=pars["aa_ej"]) / pars[ "ncells_ej"]  # 2 * cgs.pi * R ** 2. * rho * one_minus_costheta / (pars["m_scale"])
        # # # # dGammadR
        dRdt_ej_0 = beta_ej_0 * cgs.c
        dEinj_0 = einj#pars["deinjdt"]
        dEinj_0 = dEinj_0 / (M0_ej_0 * cgs.c * cgs.c) / pars["ncells"]
        dEingdR_0 = dEinj_0 / dRdt_ej_0
        Doppler_0 = Gamma_ej_0 / (1. - beta_ej_0 * np.cos(theta_ej_0))
        dEingdR_abs_0 = dEingdR_0 # TODO add optical depth
        # dGammadR_ej_0 = self.dGammadR_fs(Gamma_ej_0, gammaAdi_ej_0, dlnrho1dR_ej, M2_ej_0, dM2dR_ej_0, Eint2_ej_0, dEingdR_abs)
        dGammadR_ej_0 = self.dGammadR_fs_withInj(
            Gamma_ej_0,gammaAdi_ej_0,0.,M2_ej_0,dM2dR_ej_0,Eint2_ej_0,dEingdR_abs_0,1.,Gamma_ej_0,0.,1
        )
        # # # # Energies # # # #
        dEsh2dR_ej_0 = (Gamma_ej_0 - 1.) * dM2dR_ej_0  # Shocked energy

        # --- Expansion energy
        dlnV2dR_ej_0 = dM2dR_ej_0 / M2_ej_0 - dlnrho1dR_ej_0 - dGammadR_ej_0 / Gamma_ej_0
        if pars["adiabLoss_ej"]: dEad2dR_ej_0 = -(gammaAdi_ej_0 - 1.) * Eint2_ej_0 * dlnV2dR_ej_0
        else: dEad2dR_ej_0 = 0.

        # -- Radiative losses
        dErad2dR_ej_0 = pars["epsilon_e_rad_ej"] * dEsh2dR_ej_0

        ### -- Energy equation
        dEingdR_abs_dop = dEingdR_abs_0 / Doppler_0 / Doppler_0
        dEint2dR_ej_0 = dEsh2dR_ej_0 + dEad2dR_ej_0 - dErad2dR_ej_0 + dEingdR_abs_dop# / (pars.M0 * c ** 2)

        if pars["useSpread"]:
            dttdr_ej_0 = 1. / (cgs.c * beta_ej_0) * np.sqrt(1. + R_ej_0 ** 2. * dthetadR_ej_0 ** 2.) - (1. / cgs.c)
        else:
            dttdr_ej_0 = 1. / (cgs.c * Gamma_ej_0 ** 2. * beta_ej_0 * (1. + beta_ej_0))

        assert np.isfinite(dGammadR_ej_0)
        assert dM2dR_ej_0 > 0.
        assert dthetadR_ej_0 >= 0.
        assert np.isfinite(dEint2dR_ej_0)

        out = np.array([
            dRdt_ej_0,  # 1 / beta / cgs.c,
            1 / beta_ej_0 / Gamma_ej_0 / cgs.c,
            dGammadR_ej_0 * dRdt_ej_0,
            dEint2dR_ej_0 * dRdt_ej_0,
            dthetadR_ej_0 * dRdt_ej_0,
            dErad2dR_ej_0 * dRdt_ej_0,
            dEsh2dR_ej_0 * dRdt_ej_0,
            dEad2dR_ej_0 * dRdt_ej_0,
            dM2dR_ej_0 * dRdt_ej_0,
            dttdr_ej_0 * dRdt_ej_0
        ])

        state = {"R":R_ej_0, "Gamma":Gamma_ej_0, "M0":M0_ej_0, "Eint2":Eint2_ej_0*M0_ej_0*cgs.c*cgs.c, "Adi":gammaAdi_ej_0}

        return (out,state)

    def __call__(self, t, params, pars, Eqs):
        """
        [0] 1 / beta / cgs.c,
        [1] 1 / beta / Gamma / cgs.c,
        [2] dGammadR,
        [3] dEint2dR,
        [4] dthetadR,
        [5] dErad2dR,
        [6] dEsh2dR,
        [7] dEad2dR,
        [8] dM2dR
        :param R:
        :param params:
        :param pars:
        :return:
        """

        # tburst = params[0]
        R = params[0]
        tcomoving = params[1]
        Gamma = params[2]
        beta = get_beta(Gamma)
        Eint2 = params[3]
        theta = params[4]
        M2 = params[8]
        #
        M0 = pars["M0"]
        Gamma0 = pars["Gamma0"]
        theta0 = pars["theta0"]
        rho = pars["rho"] / M0
        dlnrho1dR = pars["dlnrho1dR"]

        # R = t * beta * cgs.c

        #
        # if Eint2 < 0:
        #     beta = np.sqrt(1. - Gamma ** -2)
        #     print("ERROR! Eint2 < 0 Gamma < 1 [Gamma:{} beta:{} Eint2:{} M2:{} ] Resetting to 0".format(
        #         Gamma, beta, Eint2, M2,
        #     ))
        #     Eint2 = 0
        # if Gamma < 1:
        #     Gamma = 1.0001
        #     print("ERROR Gamma < 1 [Gamma:{} Eint2:{} M2:{} ] Resetting to 1.0001".format(
        #         Gamma, Eint2,  M2,
        #     ))
        #     # raise ValueError("Gamma < 1")
        # if Gamma > Gamma0:
        #     print("ERROR! Gamma({}) > Gamma0({}) -- resetting to Gamma0"
        #           .format(Gamma, Gamma0))
        #     Gamma = Gamma0

        #  one_over_beta = 1. / np.power(1 - np.power(Gamma, -2.), 0.5)
        # beta = np.sqrt(1. - np.power(Gamma,-2))
        # beta0 = np.sqrt(1. - Gamma0 ** -2)

        gammaAdi = EqOpts.gamma_adi_peer(Gamma, beta)
        # gammaAdi = self.gammaAdi(Gamma, beta)  # (4 + 1 / Gamma) / 3.

        # one_minus_costheta = EquationsNava.one_minus_costheta(theta)
        # one_minus_costheta0 = 1. - np.cos(theta0)

        # Spreading
        # dthetadR = self.dthetadr(gammaAdi, Gamma, R, theta, pars["aa"]) * int(pars["useSpread"])
        dthetadR = 0.#EqOpts.dthetadr_Adi(gammaAdi, Gamma, beta, R, theta, pars["aa"], pars["useSpread"], pars["Rd"], pars["thetaMax"])
        # Densities and volumes

        # rho2 = 4. * Gamma * rho
        # V2 = M2 / rho2
        dM2dR = EqOpts.dmdr(Gamma, R, pars["thetaE"], theta, rho, aa=pars["aa"]) / pars["ncells"]  # 2 * cgs.pi * R ** 2. * rho * one_minus_costheta / (pars["m_scale"])

        # # # # dGammadR
        dGammadR = self.dGammadR_fs(Gamma, gammaAdi, dlnrho1dR, M2, dM2dR, Eint2)

        if dGammadR > 0:
            if Gamma > 0.95 * Gamma0:
                raise ValueError("Gamma > 0.95 Gamma0 after RSISING")
                # dGammadR = 0.

        # # # # Energies # # # #
        dEsh2dR = (Gamma - 1.) * dM2dR  # Shocked energy

        # --- Expansion energy
        dlnV2dR = dM2dR / M2 - dlnrho1dR - dGammadR / Gamma
        if pars["adiabLoss"]:
            dEad2dR = -(gammaAdi - 1.) * Eint2 * dlnV2dR
        else:
            dEad2dR = 0.
        # dEad2dR = -(gammaAdi - 1.) * (3 / R - (1/Gamma)*dGammadR ) * Eint2 // from Gang


        # -- Radiative losses
        dErad2dR = pars["epsilon_e_rad"] * dEsh2dR

        ### -- Energy equation
        dEint2dR = dEsh2dR + dEad2dR - dErad2dR  # / (pars.M0 * c ** 2)

        dRdt = beta * cgs.c

        if pars["useSpread"]:
            dttdr = 1. / (cgs.c * beta) * np.sqrt(1. + R ** 2. * dthetadR ** 2.) - (1. / cgs.c)
        else:
            dttdr = 1. / (cgs.c * Gamma ** 2. * beta * (1. + beta))

        rho2 = EqOpts.rho2_transrel(Gamma, beta, pars["rho"], gammaAdi)
        dr_shock = EqOpts.shock_thickness(M2*M0,rho2,theta,Gamma,R, pars["ncells"])
        ctheta = pars["ctheta_j"] + 0.5 * (2.0 * theta - 2.0 * pars["theta_w_j"])

        out_grb = np.array([dRdt,  # 1 / beta / cgs.c,
                        1 / beta / Gamma / cgs.c,
                        dGammadR * dRdt,
                        dEint2dR * dRdt,
                        dthetadR * dRdt,
                        dErad2dR * dRdt,
                        dEsh2dR * dRdt,
                        dEad2dR * dRdt,
                        dM2dR * dRdt,
                        dttdr * dRdt
                       ])

        ''' ------------------------------------------------------------------------------------------------------- '''

        # i1 = pars["neq"]
        #
        # R_ej_0 = params[0 + i1]
        # tcomoving_ej_0 = params[1 + i1]
        # Gamma_ej_0 = params[2 + i1]
        # if (Gamma_ej_0 < 1):
        #     Gamma_ej_0 = 1.
        # tmp_0 = np.power(Gamma_ej_0, -2.)
        # beta_ej_0 = np.sqrt(1. - tmp_0)
        # Eint2_ej_0 = params[3 + i1]
        # theta_ej_0 = params[4 + i1]
        # M2_ej_0 = params[8 + i1]
        # #
        # M0_ej_0 = pars["M0_ej_0"]
        # Gamma0_ej_0 = pars["Gamma0_ej_0"]
        # theta0_ej_0 = pars["theta0_ej_0"]
        #
        # # if (Eint2_ej_0 < Eint2):
        # #     print("problem")
        #
        # ctheta_ej_0 = pars["ctheta_ej_0"] + 0.5 * (2.0 * theta_ej_0 - 2.0 * pars["theta_w_ej_0"])
        # gammaAdi_ej_0 = EqOpts.gamma_adi_peer(Gamma_ej_0, beta_ej_0)
        #
        # # rho_ej_default_0 = pars["rho_ej_0"] / M0_ej_0
        # # dlnrho1dR_ej_default_0 = pars["dlnrho1dR_ej_0"]
        # # if (pars["use_ctheta_lim"]):
        # #     if(pars["use_st_rho"]):
        # #         if ((ctheta > ctheta_ej_0)&(R > R_ej_0)):
        # #             if (not pars["is_entered"]) :
        # #                 pars["t_ent"] = t
        # #                 pars["r_ej_ent"] = R_ej_0
        # #             rho_ej = 1e-90#rho_ej_default_0 * np.exp( - (R_ej_0 - pars["r_ej_ent"]) / R_ej_0 )
        # #             dlnrho1dR_ej = 0.#rho_ej_default_0 * ( -1. * R_ej_0**-3 * pars["r_ej_ent"] * np.exp(pars["r_ej_ent"]/R_ej_0 - 1) )
        # #
        # #             if ((~np.isfinite(rho_ej))&(rho_ej < 0)):
        # #                 raise ValueError()
        # #             if (~np.isfinite(dlnrho1dR_ej)):
        # #                 raise ValueError()
        # #
        # #             # o_sedov = SedovDimensionLess(gamma=gammaAdi_ej_0)
        # #             # inside the jet cone & behind the head
        # #             # From the expanding jet head (decaying)
        # #
        # #             # arc_jet = ctheta * R
        # #             # arc_ej = ctheta_ej_0 * R_ej_0
        # #             # assert arc_ej < arc_jet
        # #             # rho_from_arc = pars["o_sedov"].rho_profile(arc_ej, arc_jet, rho2) / M0_ej_0
        # #             #
        # #             #
        # #             #
        # #             # rho_from_arc_m1 = pars["o_sedov"].rho_profile(arc_ej - 0.01 * arc_ej, arc_jet, rho2) / M0_ej_0
        # #             # dlnrho1dR_ej_from_arc = (1. / arc_ej) * (rho_from_arc - rho_from_arc_m1) / (arc_ej - arc_ej * 0.01)
        # #             # assert dlnrho1dR_ej_from_arc < 1
        # #             #
        # #             # rho_ej = rho_from_arc if (rho_from_arc < rho_ej_default_0) else rho_ej_default_0
        # #             # dlnrho1dR_ej = dlnrho1dR_ej_from_arc if (rho_from_arc < rho_ej_default_0) else dlnrho1dR_ej_default_0
        # #
        # #             #
        # #             # # From jet head itself (rising)
        # #             rho_ej_behind = pars["o_sedov"].rho_profile(R_ej_0, R, rho2) / M0_ej_0
        # #             rho_ej_behind_im1 = pars["o_sedov"].rho_profile(R_ej_0 - 0.01 * R_ej_0, R, rho2) / M0_ej_0
        # #             dlnrho1dR_ej_behind = (1. / R_ej_0) * (rho_ej_behind - rho_ej_behind_im1) / (R_ej_0 - 0.01 * R_ej_0)
        # #             rho_ej = rho_ej_behind
        # #             dlnrho1dR_ej = dlnrho1dR_ej_behind
        # #             # # chose the largest
        # #             # if ((rho_from_arc < rho_ej_behind)&(rho_ej_default_0 > rho_from_arc)):
        # #             #     rho_ej_final = rho_from_arc
        # #             #     dlnrho1dR_ej_final = dlnrho1dR_ej_from_arc
        # #             # # elif((rho_from_arc > rho_ej_behind)&(pars["rho_ej"] / M0_ej_0 <= rho_from_arc)):
        # #             # #     rho_ej_final = rho_from_arc
        # #             # #     dlnrho1dR_ej_final = dlnrho1dR_ej_from_arc
        # #             # else:
        # #             #     rho_ej_final = rho_ej_behind
        # #             #     dlnrho1dR_ej_final = dlnrho1dR_ej_behind
        # #             #
        # #             # if (rho_ej_behind < rho_ej_default_0):
        # #             #     rho_ej = rho_ej_default_0
        # #             #     dlnrho1dR_ej = dlnrho1dR_ej_default_0
        # #             # else:
        # #             #     rho_ej = rho_ej_final
        # #             #     dlnrho1dR_ej = dlnrho1dR_ej_final
        # #         elif ((ctheta > ctheta_ej_0) & (R < R_ej_0)):
        # #             # inside & above
        # #             rho_ej = rho_ej_default_0
        # #             dlnrho1dR_ej = dlnrho1dR_ej_default_0
        # #         elif((ctheta < ctheta_ej_0)):
        # #             # outside the jet cone
        # #             rho_ej = rho_ej_default_0
        # #             dlnrho1dR_ej = dlnrho1dR_ej_default_0
        # #         else:
        # #             rho_ej = rho_ej_default_0
        # #             dlnrho1dR_ej = dlnrho1dR_ej_default_0
        # #
        # #
        # #         # if (ctheta > ctheta_ej_0):
        # #         #     o_sedov = SedovDimensionLess(gamma=gammaAdi_ej_0)
        # #         #     # assert R > R_ej_0
        # #         #     # print(pars["ctheta_j_arr"])
        # #         #     # idx = np.argmax( pars["ctheta_j_arr"] > ctheta_ej_0 ) + 1
        # #         #     assert pars["Rj_arr"][-1] > R_ej_0
        # #         #     idx = find_nearest_index(np.asarray(pars["Rj_arr"]),R_ej_0)+1#np.argmax( pars["Rj_arr"] < R_ej_0 ) - 1
        # #         #     # lendth of the 'arc' of the jet when jet had the same radius as currently ejecta has
        # #         #     arc = pars["ctheta_j_arr"][idx] * pars["Rj_arr"][idx]
        # #         #     # current arc of the ejecta
        # #         #     arc_ej = ctheta_ej_0 * R_ej_0
        # #         #     if not ( arc > arc_ej ):
        # #         #         raise ValueError("this is essencially ctheta > ctheta_ej_0 requirement not fulfilled "
        # #         #                          "R[idx]={:.2e} Rej={:.2e} Arc={:.2e} Arc_ej={:.2e}".format( pars["Rj_arr"][idx], R_ej_0, arc, arc_ej ))
        # #         #     # rho2 of the jet when it had radius that ejecta now has
        # #         #     rho2_ = pars["rho2_j_arr"][idx]
        # #         #     # Sedov profile for laterally moving piece... yeah, I know
        # #         #     rho_from_arc = o_sedov.rho_profile(arc_ej, arc, rho2_) / M0_ej_0
        # #         #     rho_from_arc_m1 = o_sedov.rho_profile(arc_ej-0.1*arc_ej, arc, rho2_) / M0_ej_0
        # #         #     dlnrho1dR_ej_from_arc= (1. / arc_ej) * (rho_from_arc - rho_from_arc_m1) / (arc_ej - arc_ej-0.1)
        # #         #     assert dlnrho1dR_ej_from_arc < 1
        # #         #
        # #         #     rho_ej_behind = o_sedov.rho_profile(R_ej_0, R, rho2) / M0_ej_0
        # #         #     rho_ej_behind_im1 = o_sedov.rho_profile(R_ej_0 - 0.01 * R_ej_0, R, rho2) / M0_ej_0
        # #         #     dlnrho1dR_ej_behind = (1. / R_ej_0) * (rho_ej_behind - rho_ej_behind_im1) / (R_ej_0 - 0.01 * R_ej_0)
        # #         #
        # #         #     if ((rho_from_arc > rho_ej_behind)&(rho_ej_default_0 > rho_from_arc)):
        # #         #         rho_ej_final = rho_from_arc
        # #         #         dlnrho1dR_ej_final = dlnrho1dR_ej_from_arc
        # #         #     # elif((rho_from_arc > rho_ej_behind)&(pars["rho_ej"] / M0_ej_0 <= rho_from_arc)):
        # #         #     #     rho_ej_final = rho_from_arc
        # #         #     #     dlnrho1dR_ej_final = dlnrho1dR_ej_from_arc
        # #         #     else:
        # #         #         rho_ej_final = rho_ej_behind
        # #         #         dlnrho1dR_ej_final = dlnrho1dR_ej_behind
        # #         #
        # #         #     if (rho_ej_behind < rho_ej_default_0):
        # #         #         rho_ej = rho_ej_default_0
        # #         #         dlnrho1dR_ej = dlnrho1dR_ej_default_0
        # #         #     else:
        # #         #         rho_ej = rho_ej_final
        # #         #         dlnrho1dR_ej = dlnrho1dR_ej_final
        # #         # else:
        # #         #     rho_ej = rho_ej_default_0
        # #         #     dlnrho1dR_ej = dlnrho1dR_ej_default_0
        # #     else:
        # #
        # #         if ((ctheta > ctheta_ej_0)&(R > R_ej_0)):
        # #             rho_ej = 0.
        # #             dlnrho1dR_ej = 0.
        # #         else:
        # #             rho_ej = pars["rho_ej"] / M0_ej_0
        # #             dlnrho1dR_ej = pars["dlnrho1dR_ej"]
        # # else:
        # #     if(pars["use_st_rho"]):
        # #         o_sedov = SedovDimensionLess(gamma=gammaAdi_ej_0)
        # #         rho_ej = o_sedov.rho_profile(R_ej_0, R, rho2) / M0_ej_0
        # #         rho_ej_ = o_sedov.rho_profile(R_ej_0 - 0.01 * R_ej_0, R, rho2) / M0_ej_0
        # #         dlnrho1dR_ej = (1. / R_ej_0) * (rho_ej - rho_ej_) / (R_ej_0 - 0.01 * R_ej_0)
        # #         if (rho_ej < pars["rho_ej"] / M0_ej_0):
        # #             rho_ej = pars["rho_ej"] / M0_ej_0
        # #             dlnrho1dR_ej = pars["dlnrho1dR_ej"]
        # #         else:
        # #             pass
        # #     else:
        # #         rho_ej = pars["rho_ej"] / M0_ej_0
        # #         dlnrho1dR_ej = pars["dlnrho1dR_ej"]
        # #
        # rho_ej_0 = pars["rho_ej_0"] / M0_ej_0
        # dlnrho1dR_ej_0 = pars["dlnrho1dR_ej_0"]
        #
        #
        # if (not np.isfinite(Gamma_ej_0)):
        #     raise ValueError("nan in Gamma")
        # if ((not np.isfinite(M2_ej_0)) | (M2_ej_0 < 0.)):
        #     raise ValueError("nan in M2_ej_0")
        # if (not np.isfinite(beta_ej_0)):
        #     raise ValueError("nan in beta")
        #
        # gammaAdi_ej_0 = EqOpts.gamma_adi_peer(Gamma_ej_0, beta_ej_0)
        # dthetadR_ej_0 = EqOpts.dthetadr_None(gammaAdi_ej_0, Gamma_ej_0, beta_ej_0, R_ej_0, theta_ej_0, pars["aa_ej"], pars["useSpread_ej"], pars["Rd_ej"], pars["thetaMax_ej"])
        # dM2dR_ej_0 = EqOpts.dmdr(Gamma_ej_0, R_ej_0, pars["thetaE_ej_0"], theta_ej_0, rho_ej_0, aa=pars["aa_ej"]) / pars[ "ncells_ej"]  # 2 * cgs.pi * R ** 2. * rho * one_minus_costheta / (pars["m_scale"])
        # # # # # dGammadR
        # dEinj_0 = pars["deinjdt"]
        # dEinj_0 = dEinj_0 / (M0_ej_0 * cgs.c * cgs.c) / pars["ncells"]
        # dEingdR_0 = dEinj_0 / dRdt
        # Doppler_0 = Gamma_ej_0 / (1. - beta * np.cos(theta_ej_0))
        # dEingdR_abs_0 = dEingdR_0 # TODO add optical depth
        # # dGammadR_ej_0 = self.dGammadR_fs(Gamma_ej_0, gammaAdi_ej_0, dlnrho1dR_ej, M2_ej_0, dM2dR_ej_0, Eint2_ej_0, dEingdR_abs)
        # dGammadR_ej_0 = self.dGammadR_fs_withInj(
        #     Gamma_ej_0,gammaAdi_ej_0,0.,M2_ej_0,dM2dR_ej_0,Eint2_ej_0,dEingdR_abs_0,1.,Gamma_ej_0,0.,1
        # )
        # # # # # Energies # # # #
        # dEsh2dR_ej_0 = (Gamma_ej_0 - 1.) * dM2dR_ej_0  # Shocked energy
        #
        # # --- Expansion energy
        # dlnV2dR_ej_0 = dM2dR_ej_0 / M2_ej_0 - dlnrho1dR_ej_0 - dGammadR_ej_0 / Gamma_ej_0
        # if pars["adiabLoss_ej"]: dEad2dR_ej_0 = -(gammaAdi_ej_0 - 1.) * Eint2_ej_0 * dlnV2dR_ej_0
        # else: dEad2dR_ej_0 = 0.
        #
        # # -- Radiative losses
        # dErad2dR_ej_0 = pars["epsilon_e_rad_ej"] * dEsh2dR_ej_0
        #
        # ### -- Energy equation
        # dEingdR_abs_dop = dEingdR_abs_0 / Doppler_0 / Doppler_0
        # dEint2dR_ej_0 = dEsh2dR_ej_0 + dEad2dR_ej_0 - dErad2dR_ej_0 + dEingdR_abs_dop# / (pars.M0 * c ** 2)
        #
        # dRdt_ej_0 = beta_ej_0 * cgs.c
        #
        # assert np.isfinite(dGammadR)
        # assert np.isfinite(dGammadR_ej_0)
        # assert dM2dR_ej_0 > 0.
        # assert dthetadR_ej_0 >= 0.
        # assert np.isfinite(dEint2dR_ej_0)

        ''' ------------------------------------------------------------------------------------------------------- '''
        # i1 = pars["neq"] * (ii+1)

        out0,state0 = self.ej_rhs(pars,params,0,pars["neq"],pars["deinjdt"])
        out1,state1 = self.ej_rhs(pars,params,1,pars["neq"]*2,0.)
        out2,state2 = self.ej_rhs(pars,params,2,pars["neq"]*3,0.)



        if (False and state0["R"] > state1["R"]):
            # compute the stuff
            # def eqs(p,
            #         gJ, MJ, mJ, eJ, gAdiJ,
            #         gEJ, MEJ, mEJ, eEJ, gAdiEJ):
            #     gM, eM = p
            #     # gAdiJ = EqOpts.gamma_adi_peer(gJ,get_beta(gJ))
            #     # gAdiEJ = EqOpts.gamma_adi_peer(gEJ,get_beta(gEJ))
            #     gAdiM = EqOpts.gamma_adi_nava(gM,get_beta(gM))
            #     mM = (MJ + mJ) + (MEJ + mEJ)
            #     return (
            #             # total energy consercation (Ek + Eint)
            #             (gJ * cgs.c * cgs.c * (MJ + mJ) + Nava_fs_rhs.GammaEff(gJ, gAdiJ) * eJ)
            #             + (gEJ * cgs.c * cgs.c * (MEJ + mEJ) + Nava_fs_rhs.GammaEff(gEJ, gAdiEJ) * eEJ)
            #             - ( gM * cgs.c * cgs.c * mM + Nava_fs_rhs.GammaEff(gM,gAdiM) * eM ),
            #             # total momentum conservation
            #             cgs.c * np.sqrt(gJ * gJ - 1) * (MJ + mJ + gAdiJ * eJ / cgs.c ** 2) +
            #             cgs.c * np.sqrt(gEJ * gEJ - 1) * (MEJ + mEJ + gAdiEJ * eEJ / cgs.c ** 2) -
            #             cgs.c * np.sqrt(gM * gM - 1) * (mM + gAdiM * eM / cgs.c ** 2),
            #             # total mass conservation
            #             # (MJ + mJ) + (MEJ + mEJ) - mM
            #
            #     )
            def eqs(p,
                    gJ, mJ, eJ, gAdiJ,
                    gEJ, mEJ, eEJ, gAdiEJ):
                gM, eM = p
                # gAdiJ = EqOpts.gamma_adi_peer(gJ,get_beta(gJ))
                # gAdiEJ = EqOpts.gamma_adi_peer(gEJ,get_beta(gEJ))
                gAdiM = EqOpts.gamma_adi_peer(gM,get_beta(gM))
                mM = mJ + mEJ
                return (
                    # total energy consercation (Ek + Eint)
                    (gJ * mJ + Nava_fs_rhs.GammaEff(gJ, gAdiJ) * eJ)
                    + (gEJ * mEJ + Nava_fs_rhs.GammaEff(gEJ, gAdiEJ) * eEJ)
                    - ( gM * mM + Nava_fs_rhs.GammaEff(gM,gAdiM) * eM ),
                    # total momentum conservation
                    np.sqrt(gJ * gJ - 1) * (mJ + gAdiJ * eJ ) +
                    np.sqrt(gEJ * gEJ - 1) * (mEJ + gAdiEJ * eEJ ) -
                    np.sqrt(gM * gM - 1) * (mM + gAdiM * eM ),
                    # total mass conservation
                    # (MJ + mJ) + (MEJ + mEJ) - mM
                )

            gJ=state0["Gamma"]; mJ=state0["M0"]; eJ=state0["Eint2"]; gAdiJ=state0["Adi"]
            gEJ=state1["Gamma"]; mEJ=state1["M0"]; eEJ=state1["Eint2"]; gAdiEJ=state1["Adi"]
            i_gM=state1["Gamma"]; i_mM=state0["M0"]+state1["M0"]; i_eM=state0["Eint2"]+state1["Eint2"]
            #
            arr_gM = np.linspace(1.16, 1.17, 250)
            arr_gM = np.linspace(1.000000000001, 1.4, 250)
            # arr_mM = np.logspace(np.log10(i_mM)/2.,np.log10(i_mM)*2.,50)
            arr_eM = np.logspace( np.log10(i_eM/100.), np.log10(i_eM*4.), 200)
            #
            Arr_gM, Arr_eM = np.meshgrid(arr_gM, arr_eM)
            Zer_gM, Zer_eM = eqs( (Arr_gM, Arr_eM ), gJ, mJ, eJ, gAdiJ, gEJ, mEJ, eEJ, gAdiEJ )
            print("Sol Emin={} Emax={} | Gmin={} Gmax={}".format(Zer_eM.min(),Zer_eM.max(),Zer_gM.min(),Zer_gM.max()))
            fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
            # ax.set_zlim(-1e2,1e2)

            Res = np.copy(Zer_eM)
            Res2 = np.copy(Zer_eM)
            Res[(Res2 > 0)] = np.log10(Res2[(Res2 > 0)])
            Res[(Res2 < 0)] = 0.#-np.log10(-1 * Res2[(Res2 < 0)])
            # if (Res.min()==-1e-2 && 1)
            surf = ax.plot_surface(Arr_gM, np.log10(Arr_eM*cgs.c*cgs.c), Res, cmap='Reds', rstride=4, cstride=4, edgecolor='none',
                                   linewidth=1, antialiased=True, alpha=0.4,
                                   # rstride=100,
                                   # cstride=100
                                   )

            Res = np.copy(Zer_gM)
            Res2 = np.copy(Zer_gM)
            Res[(Res2 > 0)] = np.log10(Res2[(Res2 > 0)])
            Res[(Res2 < 0)] = 0.#-np.log10(-1 * Res2[(Res2 < 0)])
            # if (Res.min()==-1e-2 && 1)
            surf = ax.plot_surface(Arr_gM, np.log10(Arr_eM*cgs.c*cgs.c), Res, cmap='Blues', rstride=4, cstride=1,
                                   edgecolor='none',
                                   linewidth=1, antialiased=True, alpha=0.4
                                   # rstride=100,
                                   # cstride=100
                                   )

            sol = np.array([ i_gM, i_eM ])
            x = eqs( sol,  gJ, mJ, eJ, gAdiJ,  gEJ, mEJ, eEJ, gAdiEJ )
            it = 0
            while ( (abs(x[0]) > 1e-8) | (abs(x[0]) > 1e-8) ):
                it = it + 1
                print("it= {} f(Solution)={}".format(it, eqs( (sol[0], sol[1]), gJ, mJ, eJ, gAdiJ,  gEJ, mEJ, eEJ, gAdiEJ )))
                sol = fsolve( eqs, sol,  # initial guess
                              args=(gJ, mJ, eJ, gAdiJ, gEJ, mEJ, eEJ, gAdiEJ
                                    ), xtol=1e-12, maxfev=111500)
                x = eqs( sol, gJ, mJ, eJ, gAdiJ, gEJ, mEJ, eEJ, gAdiEJ)
            #
            # x = fsolve(eqs,  np.array( [ i_gM, i_eM ] ), # initial guess
            #            args=( gJ, mJ, eJ, gAdiJ,
            #                   gEJ, mEJ, eEJ, gAdiEJ
            #                  ), xtol=1e-08, maxfev=111500
            #            )
            # x = fsolve(eqs,  np.array( [ x[0], x[1] ] ), # initial guess
            #            args=( gJ, mJ, eJ, gAdiJ,
            #                   gEJ, mEJ, eEJ, gAdiEJ
            #                  ), xtol=1e-08, maxfev=111500
            #            )
            # x = fsolve(eqs,  np.array( [ x[0], x[1] ] ), # initial guess
            #            args=( gJ, mJ, eJ, gAdiJ,
            #                   gEJ, mEJ, eEJ, gAdiEJ
            #                  ), xtol=1e-08, maxfev=111500
            #            )
            # gM_, eM_ = x[0], x[1]

            ax.plot([sol[0], sol[0]], np.log10(np.array([sol[1], sol[1]])*cgs.c*cgs.c), zs=[0, 1], lw=9, color="black")

            gM_, eM_ = sol[0], sol[1]
            # eM_ = mJ+mEJ
            print('------------------------------------------------------')
            print("Shells colliding: \n" +
                  "gJ={}, mJ={}, eJ={}, gAdiJ={} \n".format(gJ, mJ, eJ, gAdiJ) +
                  "gEJ={}, mEJ={}, eEJ={}, gAdiEJ={} \n".format(gEJ, mEJ, eEJ, gAdiEJ)+
                  "i_gM={}, i_mM={}, i_eM={} \n".format(i_gM, i_mM, i_eM)+
                  "f_gM={}, f_mM={}, f_eM={} \n".format(gM_, i_mM, eM_),
                  "verification f()={}".format( eqs( (gM_, eM_), gJ, mJ, eJ, gAdiJ,  gEJ, mEJ, eEJ, gAdiEJ ))
                  )
            print('------------------------------------------------------')

            plt.show()

            exit(1);

        out = np.hstack((out_grb,out0,out1,out2)).flatten()
        return out




        # return np.array([dRdt,  # 1 / beta / cgs.c,
        #                  1 / beta / Gamma / cgs.c,
        #                  dGammadR * dRdt,
        #                  dEint2dR * dRdt,
        #                  dthetadR * dRdt,
        #                  dErad2dR * dRdt,
        #                  dEsh2dR * dRdt,
        #                  dEad2dR * dRdt,
        #                  dM2dR * dRdt,
        #                  dttdr * dRdt,
        #                  # -----------------------------
        #                  dRdt_ej_0,  # 1 / beta / cgs.c,
        #                  1 / beta_ej_0 / Gamma_ej_0 / cgs.c,
        #                  dGammadR_ej_0 * dRdt_ej_0,
        #                  dEint2dR_ej_0 * dRdt_ej_0,
        #                  dthetadR_ej_0 * dRdt_ej_0,
        #                  dErad2dR_ej_0 * dRdt_ej_0,
        #                  dEsh2dR_ej_0 * dRdt_ej_0,
        #                  dEad2dR_ej_0 * dRdt_ej_0,
        #                  dM2dR_ej_0 * dRdt_ej_0,
        #                  # -------------------------------
        #                  dRdt_ej_0,  # 1 / beta / cgs.c,
        #                  1 / beta_ej_0 / Gamma_ej_0 / cgs.c,
        #                  dGammadR_ej_0 * dRdt_ej_0,
        #                  dEint2dR_ej_0 * dRdt_ej_0,
        #                  dthetadR_ej_0 * dRdt_ej_0,
        #                  dErad2dR_ej_0 * dRdt_ej_0,
        #                  dEsh2dR_ej_0 * dRdt_ej_0,
        #                  dEad2dR_ej_0 * dRdt_ej_0,
        #                  dM2dR_ej_0 * dRdt_ej_0,
        #                  # ----------------------------------
        #                  dRdt_ej_0,  # 1 / beta / cgs.c,
        #                  1 / beta_ej_0 / Gamma_ej_0 / cgs.c,
        #                  dGammadR_ej_0 * dRdt_ej_0,
        #                  dEint2dR_ej_0 * dRdt_ej_0,
        #                  dthetadR_ej_0 * dRdt_ej_0,
        #                  dErad2dR_ej_0 * dRdt_ej_0,
        #                  dEsh2dR_ej_0 * dRdt_ej_0,
        #                  dEad2dR_ej_0 * dRdt_ej_0,
        #                  dM2dR_ej_0 * dRdt_ej_0
        #                  ])

    # @staticmethod
    # def gammaAdi(Gamma, beta):
    #     """ From Peer 2012 arxiv:1203.5797 """
    #     mom = Gamma * beta
    #     theta = mom / 3. * (mom + 1.07 * mom ** 2.) / (1 + mom + 1.07 * mom ** 2.)
    #     zz = theta / (0.24 + theta)
    #     gamma_adi = (5. - 1.21937 * zz + 0.18203 * zz ** 2. - 0.96583 * zz ** 3. + 2.32513 * zz ** 4. -
    #             2.39332 * zz ** 5. + 1.07136 * zz ** 6.) / 3.
    #     return gamma_adi

    # @staticmethod
    # def gammaAdi(Gamma):
    #     return (4 + 1 / Gamma) / 3.

    # @staticmethod
    # def dthetadr(gammaAdi, Gamma, R, theta, aa=np.nan):
    #     """
    #     Source:
    #     1. / (R * Gamma ** (1. + aa) * theta ** (aa))
    #     is from https://academic.oup.com/mnras/article/421/1/570/990386
    #     vperp / R / Gamma * one_over_beta / cgs.c
    #     is from ???
    #     :param gammaAdi:
    #     :param Gamma:
    #     :param R:
    #     :param one_over_beta:
    #     :param theta:
    #     :param aa:
    #     :return:
    #     """
    #
    #     if theta < np.pi:
    #         if np.isfinite(aa):
    #             return 1. / (R * Gamma ** (1. + aa) * theta ** (aa))
    #         else:
    #             vperp = np.sqrt(gammaAdi * (gammaAdi - 1) * (Gamma - 1) / (1 + gammaAdi * (Gamma - 1))) * cgs.c
    #             one_over_beta = 1. / np.power(1 - np.power(Gamma, -2.), 0.5)
    #             return vperp / R / Gamma * one_over_beta / cgs.c
    #     else:
    #         return 0.

    @staticmethod
    def GammaEff(Gamma, gammaAdi):
        return (gammaAdi * Gamma ** 2. - gammaAdi + 1.) / Gamma

    @staticmethod
    def dGammaEffdGamma(Gamma, gammaAdi):
        return (gammaAdi * Gamma ** 2. + gammaAdi - 1.) / Gamma ** 2
        # return 4. / 3. + 1. / Gamma ** 2. / 3. + 2. / Gamma ** 3. / 3.

    @staticmethod
    def dGammadR_fs(Gamma, gammaAdi, dlnrho1dR, M2, dM2dR, Eint2):
        GammaEff = Nava_fs_rhs.GammaEff(Gamma, gammaAdi)  # , gammaAdi) #(gammaAdi * Gamma ** 2. - gammaAdi + 1.) / Gamma
        dGammaEffdGamma = Nava_fs_rhs.dGammaEffdGamma(Gamma, gammaAdi)  # 4. / 3. + 1. / Gamma ** 2. / 3. + 2. / Gamma ** 3. / 3.
        f_2 = GammaEff * (gammaAdi - 1.) * Eint2 / Gamma
        h_2 = GammaEff * (gammaAdi - 1.) * Eint2 * (dM2dR / M2 - dlnrho1dR)
        dGammadR = -((Gamma - 1.) * (GammaEff + 1.) * dM2dR - h_2) / ((1. + M2) + Eint2 * dGammaEffdGamma + f_2)
        return dGammadR
    @staticmethod
    def dGammadR_fs_withInj(Gamma, gammaAdi, drhodr, M2, dM2dR, Eint2, dEingdR_abs, GammaRho, GammaRel, dGammaRhodR, dGammaRelDGamma):
        xi_inj = 1.
        GammaEff = Nava_fs_rhs.GammaEff(Gamma, gammaAdi); # TODO check if GammaRel should be used for this!!!
        dGammaEffdGamma = Nava_fs_rhs.dGammaEffdGamma(Gamma, gammaAdi);
        num1 = (Gamma - GammaRho + GammaEff * (GammaRel - 1.)) * dM2dR;
        num2 = - GammaEff * (gammaAdi - 1.) * Eint2 * (dM2dR/M2 - drhodr/rho - dGammaRhodR / GammaRho); # - 3.*Eint2/R
        num3 = - (1. - GammaEff / Gamma * xi_inj) * dEingdR_abs * 1;
        denum1 = (1.+M2);
        denum2 = Eint2 * dGammaEffdGamma;
        denom3 = GammaEff * (gammaAdi - 1.) * Eint2 * dGammaRelDGamma / GammaRel;
        dGammadR = -1. * (num1 + num2 + num3) / (denum1+denum2+denom3);
        return dGammadR

class Driver_Nava_FS_t2(Driver_t):

    def __init__(
            self,
            Gamma0=1000.,
            theta0=np.pi / 2.,
            M0=1.e53 * cgs.c ** -2 / 1000.,
            tstart=1.e3, Rstart=1.e18,
            rho0=1.e-2 * cgs.mppme,
            dlnrho1dR_0=0.,

            Gamma0_ej_0=1000.,
            theta0_ej_0=np.pi / 2.,
            M0_ej_0=1.e53 * cgs.c ** -2 / 1000.,
            tstart_ej_0=1.e3, Rstart_ej_0=1.e18,
            rho0_ej_0=1.e-2 * cgs.mppme,
            dlnrho1dR_0_ej_0=0.,

            Gamma0_ej_1=1000.,
            theta0_ej_1=np.pi / 2.,
            M0_ej_1=1.e53 * cgs.c ** -2 / 1000.,
            tstart_ej_1=1.e3, Rstart_ej_1=1.e18,
            rho0_ej_1=1.e-2 * cgs.mppme,
            dlnrho1dR_0_ej_1=0.,

            Gamma0_ej_2=1000.,
            theta0_ej_2=np.pi / 2.,
            M0_ej_2=1.e53 * cgs.c ** -2 / 1000.,
            tstart_ej_2=1.e3, Rstart_ej_2=1.e18,
            rho0_ej_2=1.e-2 * cgs.mppme,
            dlnrho1dR_0_ej_2=0.,

            **kwargs
    ):
        assert (not M0 is None)
        # assert np.isfinite(kwargs["aa"])
        beta0 = get_beta(Gamma0)
        beta0_ej_0 = get_beta(Gamma0_ej_0)
        beta0_ej_1 = get_beta(Gamma0_ej_1)
        beta0_ej_2 = get_beta(Gamma0_ej_2)
        ncells = 1

        # some equations for ODEs
        self.eqs_ode_rhs = {"eq_dthetadr": EqOpts.dthetadr_None,
                            "eq_gammaAdi": EqOpts.gamma_adi_peer}

        Rd = 0#get_Rdec2(E0, rho0 / cgs.mppme, Gamma0)
        Rd_ej_0 = 0#get_Rdec2(E0, rho0 / cgs.mppme, Gamma0)
        Rd_ej_1 = 0#get_Rdec2(E0, rho0 / cgs.mppme, Gamma0)
        Rd_ej_2 = 0#get_Rdec2(E0, rho0 / cgs.mppme, Gamma0)
        Rstart = tstart * get_beta(Gamma0) * cgs.c
        Rstart_ej_0 = tstart_ej_0 * get_beta(Gamma0_ej_0) * cgs.c
        Rstart_ej_1 = tstart_ej_1 * get_beta(Gamma0_ej_1) * cgs.c
        Rstart_ej_2 = tstart_ej_2 * get_beta(Gamma0_ej_2) * cgs.c
        M20 = (2 / 3.) * np.pi * Rstart ** 3. * (1. - np.cos(theta0)) * rho0 / ncells
        M20_ej_0 = (2 / 3.) * np.pi * Rstart_ej_0 ** 3. * (1. - np.cos(theta0_ej_0)) * rho0_ej_0 / ncells
        M20_ej_1 = (2 / 3.) * np.pi * Rstart_ej_1 ** 3. * (1. - np.cos(theta0_ej_1)) * rho0_ej_1 / ncells
        M20_ej_2 = (2 / 3.) * np.pi * Rstart_ej_2 ** 3. * (1. - np.cos(theta0_ej_2)) * rho0_ej_2 / ncells

        self.v_ns_init_vals = [
            "R", "tcomoving", "Gamma", "Eint2", "theta", "Erad2", "Esh2", "Ead2", "M2", "ttin",
            "R_ej_0", "tcomoving_ej_0", "Gamma_ej_0", "Eint2_ej_0", "theta_ej_0", "Erad2_ej_0", "Esh2_ej_0", "Ead2_ej_0", "M2_ej_0", "ttin_ej_0",
            "R_ej_1", "tcomoving_ej_1", "Gamma_ej_1", "Eint2_ej_1", "theta_ej_1", "Erad2_ej_1", "Esh2_ej_1", "Ead2_ej_1", "M2_ej_1", "ttin_ej_1",
            "R_ej_2", "tcomoving_ej_2", "Gamma_ej_2", "Eint2_ej_2", "theta_ej_2", "Erad2_ej_2", "Esh2_ej_2", "Ead2_ej_2", "M2_ej_2", "ttin_ej_2"
        ]
        init_data = {
            "R": Rstart,#Rstart / (beta0 * cgs.c),  # 0
            "tcomoving": Rstart / (beta0 * Gamma0 * cgs.c),  # 1
            "Gamma": Gamma0,  # 2
            "Eint2": (Gamma0 - 1) * M20 / M0,  # 3 # 0.75 *
            "theta": theta0,  # 4
            "Erad2": 0.,  # 5
            "Esh2": 0.,  # 6
            "Ead2": 0.,  # 7
            "M2": M20 / M0,  # 8,
            "ttin": self.init_elapsed_time(Rstart, beta0, Gamma0, **kwargs),

            "R_ej_0": Rstart_ej_0,  # Rstart / (beta0 * cgs.c),  # 0
            "tcomoving_ej_0": Rstart_ej_0 / (beta0_ej_0 * Gamma0_ej_0 * cgs.c),  # 1
            "Gamma_ej_0": Gamma0_ej_0,  # 2
            "Eint2_ej_0": (Gamma0_ej_0 - 1) * M20_ej_0 / M0_ej_0,  # 3 # 0.75 *
            "theta_ej_0": theta0_ej_0,  # 4
            "Erad2_ej_0": 0.,  # 5
            "Esh2_ej_0": 0.,  # 6
            "Ead2_ej_0": 0.,  # 7
            "M2_ej_0": M20_ej_0 / M0_ej_0,  # 8
            "ttin_ej_0": self.init_elapsed_time(Rstart_ej_0, beta0_ej_0, Gamma0_ej_0, **kwargs),

            "R_ej_1": Rstart_ej_1,  # Rstart / (beta0 * cgs.c),  # 0
            "tcomoving_ej_1": Rstart_ej_1 / (beta0_ej_1 * Gamma0_ej_1 * cgs.c),  # 1
            "Gamma_ej_1": Gamma0_ej_1,  # 2
            "Eint2_ej_1": (Gamma0_ej_1 - 1) * M20_ej_1 / M0_ej_1,  # 3 # 0.75 *
            "theta_ej_1": theta0_ej_1,  # 4
            "Erad2_ej_1": 0.,  # 5
            "Esh2_ej_1": 0.,  # 6
            "Ead2_ej_1": 0.,  # 7
            "M2_ej_1": M20_ej_1 / M0_ej_1,  # 8
            "ttin_ej_1": self.init_elapsed_time(Rstart_ej_1, beta0_ej_1, Gamma0_ej_1, **kwargs),

            "R_ej_2": Rstart_ej_2,  # Rstart / (beta0 * cgs.c),  # 0
            "tcomoving_ej_2": Rstart_ej_2 / (beta0_ej_2 * Gamma0_ej_2 * cgs.c),  # 1
            "Gamma_ej_2": Gamma0_ej_2,  # 2
            "Eint2_ej_2": (Gamma0_ej_2 - 1) * M20_ej_2 / M0_ej_2,  # 3 # 0.75 *
            "theta_ej_2": theta0_ej_2,  # 4
            "Erad2_ej_2": 0.,  # 5
            "Esh2_ej_2": 0.,  # 6
            "Ead2_ej_2": 0.,  # 7
            "M2_ej_2": M20_ej_2 / M0_ej_2,  # 8
            "ttin_ej_2": self.init_elapsed_time(Rstart_ej_2, beta0_ej_2, Gamma0_ej_2, **kwargs),

        }
        self.initial_data = self._set_ode_inits(init_data)

        self.pars_ode_rhs = {
            "M0": M0, "theta0": theta0, "Gamma0": Gamma0, "rho": rho0, "dlnrho1dR": dlnrho1dR_0,
            "epsilon_e_rad": kwargs["epsilon_e_rad"], "useSpread": kwargs["useSpread"],
            "adiabLoss": kwargs["adiabLoss"],
            "ncells": kwargs["ncells"], "aa": kwargs["aa"],
            "Rd": Rd, "thetaMax":kwargs["thetaMax"], "thetaE":0., "ctheta_j":kwargs["ctheta_j"],"theta_w_j":kwargs["theta_w_j"],

            "M0_ej_0": M0_ej_0, "theta0_ej_0": theta0_ej_0, "Gamma0_ej_0": Gamma0_ej_0, "rho_ej_0": rho0_ej_0,
            "dlnrho1dR_ej_0": dlnrho1dR_0_ej_0, "Rd_ej_0": Rd_ej_0, "thetaE_ej_0":0., "ctheta_ej_0":kwargs["ctheta_ej_0"],
            "theta_w_ej_0":kwargs["theta_w_ej_0"],
            "M0_ej_1": M0_ej_1, "theta0_ej_1": theta0_ej_1, "Gamma0_ej_1": Gamma0_ej_1, "rho_ej_1": rho0_ej_1,
            "dlnrho1dR_ej_1": dlnrho1dR_0_ej_1, "Rd_ej": Rd_ej_1, "thetaE_ej_1":0., "ctheta_ej_1":kwargs["ctheta_ej_1"],
            "theta_w_ej_1":kwargs["theta_w_ej_1"],
            "M0_ej_2": M0_ej_2, "theta0_ej_2": theta0_ej_2, "Gamma0_ej_2": Gamma0_ej_2, "rho_ej_2": rho0_ej_2,
            "dlnrho1dR_ej_2": dlnrho1dR_0_ej_2, "Rd_ej_2": Rd_ej_2, "thetaE_ej_2":0., "ctheta_ej_2":kwargs["ctheta_ej_2"],
            "theta_w_ej_2":kwargs["theta_w_ej_2"],

            "epsilon_e_rad_ej": kwargs["epsilon_e_rad"], "useSpread_ej": False,
            "adiabLoss_ej": kwargs["adiabLoss"],
            "ncells_ej": kwargs["ncells"], "aa_ej": -1,
            "thetaMax_ej":kwargs["thetaMax"],




            "Rj_arr":[],
            "ctheta_j_arr":[],
            "rho2_j_arr":[],
            "tburst_arr":[],
            "o_sedov": SedovDimensionLess(gamma=5/3),

            "t_ent":0, "is_entered":False, "r_ej_ent":0.,

            "neq" : 10, "use_st_rho" : kwargs["use_st_rho"], "use_ctheta_lim" : kwargs["use_ctheta_lim"],

            "deinjdt" : 0.
        }

        # initialize time elapsed in the comoving frame
        tt0 = self.init_elapsed_time(Rstart, beta0, Gamma0, **kwargs)
        tt0_ej_0 = self.init_elapsed_time(Rstart_ej_0, beta0_ej_0, Gamma0_ej_0, **kwargs)
        tt0_ej_1 = self.init_elapsed_time(Rstart_ej_1, beta0_ej_1, Gamma0_ej_1, **kwargs)
        tt0_ej_2 = self.init_elapsed_time(Rstart_ej_2, beta0_ej_2, Gamma0_ej_2, **kwargs)
        gammaAdi0 = EqOpts.gamma_adi_peer(Gamma0, beta0)
        gammaAdi0_ej_0 = EqOpts.gamma_adi_peer(Gamma0_ej_0, beta0_ej_0)
        gammaAdi0_ej_1 = EqOpts.gamma_adi_peer(Gamma0_ej_1, beta0_ej_1)
        gammaAdi0_ej_2 = EqOpts.gamma_adi_peer(Gamma0_ej_2, beta0_ej_2)
        rho20 = EqOpts.rho2_transrel(Gamma0, beta0, rho0, gammaAdi0)#self.get_rhoprime(rho0, Gamma0)
        rho20_ej_0 = EqOpts.rho2_transrel(Gamma0_ej_0, beta0_ej_0, rho0_ej_0, gammaAdi0_ej_0)#self.get_rhoprime(rho0, Gamma0)
        rho20_ej_1 = EqOpts.rho2_transrel(Gamma0_ej_1, beta0_ej_1, rho0_ej_1, gammaAdi0_ej_1)#self.get_rhoprime(rho0, Gamma0)
        rho20_ej_2 = EqOpts.rho2_transrel(Gamma0_ej_2, beta0_ej_2, rho0_ej_2, gammaAdi0_ej_2)#self.get_rhoprime(rho0, Gamma0)
        self.all_v_ns = self.v_ns_init_vals + ["tburst", "rho", "tt", "thickness", "U_e", "beta", "rho2", "gammaAdi"] + \
                        ["tburst_ej_0", "rho_ej_0", "tt_ej_0", "thickness_ej_0", "U_e_ej_0", "beta_ej_0", "rho2_ej_0", "gammaAdi_ej_0"] + \
                        ["tburst_ej_1", "rho_ej_1", "tt_ej_1", "thickness_ej_1", "U_e_ej_1", "beta_ej_1", "rho2_ej_1", "gammaAdi_ej_1"] + \
                        ["tburst_ej_2", "rho_ej_2", "tt_ej_2", "thickness_ej_2", "U_e_ej_2", "beta_ej_2", "rho2_ej_2", "gammaAdi_ej_2"]
        self.dynamics = np.zeros((1, len(self.all_v_ns)))
        self.dynamics[0, :len(self.initial_data)] = self.initial_data
        self.dynamics[0, len(self.initial_data):] = \
            np.array([tstart, rho0, tt0, 0., 0., beta0, rho20, gammaAdi0] + \
                     [tstart_ej_0, rho0_ej_0, tt0_ej_0, 0., 0., beta0_ej_0, rho20_ej_0, gammaAdi0_ej_0] + \
                     [tstart_ej_1, rho0_ej_1, tt0_ej_1, 0., 0., beta0_ej_1, rho20_ej_1, gammaAdi0_ej_1] + \
                     [tstart_ej_2, rho0_ej_2, tt0_ej_2, 0., 0., beta0_ej_2, rho20_ej_2, gammaAdi0_ej_2])
        self.dynamics[0, self.i_nv("U_e")] = self.get_U_e(idx=0) * cgs.c**2
        self.dynamics[0, self.i_nv("U_e_ej_0")] = self.get_U_e_ej(idx=0,ii=0) * cgs.c**2
        self.dynamics[0, self.i_nv("U_e_ej_1")] = self.get_U_e_ej(idx=0,ii=1) * cgs.c**2
        self.dynamics[0, self.i_nv("U_e_ej_2")] = self.get_U_e_ej(idx=0,ii=2) * cgs.c**2
        self.dynamics[0] = self.apply_units(self.dynamics[0])
        # self.dynamics = np.hstack((self.apply_units(np.copy(self.initial_data)),
        #                            np.array([rho0, tt0, Rstart, 0., 0., beta0])))

        self.rhs = Nava_fs_rhs_t2()

        super(Driver_Nava_FS_t2, self).__init__(tstart, kwargs["ode_rtol"], kwargs["ode_nsteps"])

        self.kwargs = kwargs

    def get_U_e(self, idx=-1):

        rho = self.get("rho")[idx]
        Gamma = self.get("Gamma")[idx]
        M2 = self.get("M2")[idx]
        Eint2 = self.get("Eint2")[idx]

        rhoprim = 4. * rho * Gamma  # comoving density
        V2 = M2 / rhoprim  # comoving volume
        U_e = Eint2 / V2  # comoving energy density (electrons)
        # U_b = eps_b * U_e  # comoving energy density (MF)
        return U_e
    def get_U_e_ej(self, idx=-1, ii=0):

        rho = self.get("rho_ej_{}".format(ii))[idx]
        Gamma = self.get("Gamma_ej_{}".format(ii))[idx]
        M2 = self.get("M2_ej_{}".format(ii))[idx]
        Eint2 = self.get("Eint2_ej_{}".format(ii))[idx]

        rhoprim = 4. * rho * Gamma  # comoving density
        V2 = M2 / rhoprim  # comoving volume
        U_e = Eint2 / V2  # comoving energy density (electrons)
        # U_b = eps_b * U_e  # comoving energy density (MF)
        return U_e
    def apply_units(self, i_res):
        # self.units_dic = {
        #     "M2": M0,
        #     "Eint2": M0 * cgs.c**2,
        #     "Erad2": M0 * cgs.c**2,
        #     "Esh2": M0 * cgs.c**2,
        #     "Ead2": M0 * cgs.c**2,
        # }
        i_res[self.i_nv("M2")] *= self.pars_ode_rhs["M0"]
        i_res[self.i_nv("Eint2")] *= self.pars_ode_rhs["M0"] * cgs.c**2
        i_res[self.i_nv("Erad2")] *= self.pars_ode_rhs["M0"] * cgs.c**2
        i_res[self.i_nv("Esh2")] *= self.pars_ode_rhs["M0"] * cgs.c**2
        i_res[self.i_nv("Ead2")] *= self.pars_ode_rhs["M0"] * cgs.c**2
        for ii in [0,1,2]:

            i_res[self.i_nv(f"M2_ej_{ii}")] *= self.pars_ode_rhs[f"M0_ej_{ii}"]
            i_res[self.i_nv(f"Eint2_ej_{ii}")] *= self.pars_ode_rhs[f"M0_ej_{ii}"] * cgs.c**2
            i_res[self.i_nv(f"Erad2_ej_{ii}")] *= self.pars_ode_rhs[f"M0_ej_{ii}"] * cgs.c**2
            i_res[self.i_nv(f"Esh2_ej_{ii}")] *= self.pars_ode_rhs[f"M0_ej_{ii}"] * cgs.c**2
            i_res[self.i_nv(f"Ead2_ej_{ii}")] *= self.pars_ode_rhs[f"M0_ej_{ii}"] * cgs.c**2

        return i_res


# https://github.com/christophmschaefer/miluphcuda
# miluphcuda/test_cases/sedov/sedovAnalytical.py

class Sedov():
    """
    Analytical solution for the sedov blast wave problem
    """

    def __init__(self, time, r_max, rho0=1.,e0=1.,gamma=5/3.,w=0,n_dim=3):

        # rho0 = 1.0  # 1
        # e0 = 1.0  # 1e5
        # gamma = 5 / 3.  # 1.66 #1.333
        # w = 0  # Power law index
        # n_dim = 3

        self.sol = SedovSolution(e0, rho0, gamma=gamma, w=w, nu=n_dim)
        self.r = np.linspace(0, r_max, 1001)[1:]
        self.t = time

        print("Shock radius: {}".format(self.sol.shock_radius(self.t)))

    def compute(self, y):
        return map(self.determine, ['r', y])

    def determine(self, x):
        if x == 'r':
            return self.r
        elif x == 'velocity':
            return self.sol.velocity(self.r, self.t)
        elif x == 'rho':
            return self.sol.rho(self.r, self.t)
        elif x == 'pressure':
            return self.sol.pressure(self.r, self.t)
        elif x == 'internal_energy':
            return self.sol.internal_energy(self.r, self.t)
        else:
            raise AttributeError("Sedov solution for variable %s not known" % x)
def plot_sedov():
    # sedov = Sedov(0.001, 1)
    t_arr = np.linspace(0.01, 0.5, 5)
    r_arr = np.zeros_like(t_arr)
    v_arr = np.zeros_like(t_arr)
    fig, (ax1, ax2, ax3) = plt.subplots(nrows=3, sharex="all")
    for i, t in enumerate(t_arr):
        sedov = Sedov(t, r_max=1, rho0=1, e0=1)
        r_arr[i] = sedov.sol.shock_radius(t)
        v_arr[i] = sedov.sol.shock_velocity(t)
        ax1.plot(*sedov.compute('rho'), label="rho")
        ax2.plot(*sedov.compute('pressure'), label="pressure")
        ax3.plot(*sedov.compute('internal_energy'), label="internal energy")
        # ax1.semilogx(loc='best')
        # ax2.plot(loc='best')
        # ax3.plot(loc='best')
    plt.xlim(r_arr[0],r_arr[-1])
    plt.show()

    plt.loglog(r_arr, v_arr, '.')
    plt.show()

# https://github.com/pencil-code/pencil-code
# pencil-code/python/pencil/ism_dyn/ism_sedov_taylor.py
def sedov_taylor(*args, **kwargs):
    """
    Compute analytic radial time evolution of SN blast waves for
    comparison with numerical results
    *t_sedov*:
      Time_series object read from the simulation sn_series.dat
    *par*:
      Param object containing the simulation parameters
    *time*:
      list of time in code units
    *nt*:
      Integer size of analytic arrays
    *endt*
      Real end time in code units for the time series
    *dims*:
      Dimension of the simulation default 3D
    *rho0*:
      Ambient ISM density
    """
    st_tmp = SedovTaylor()
    st_tmp.get_st(*args, **kwargs)
    return st_tmp
class SedovTaylor():
    """
    SedovTaylor -- holds blast wave parameters and radial profiles.
    """

    def __init__(self):
        """
        Fill members with default values.
        """

        self.t = []
        self.keys = []

    def get_st(
            self,
            t_sedov=0,
            par=list(),
            time=list(),
            nt=5000,
            startt=0.0,
            endt=0.005,
            dims=3,
            quiet=True,
            rho0=1.6728e-24,
            M0=10,
            lsnowplough=True,
            lcioffi=True,
    ):
        """
        Compute analytic radial time evolution of SN blast waves for
        comparison with numerical results
        *t_sedov*:
          Time_series object read from the simulation sn_series.dat
        *par*:
          Param object containing the simulation parameters
        *time*:
          list of time in code units
        *nt*:
          Integer size of analytic arrays
        *endt*
          Real end time in code units for the time series
        *dims*:
          Dimension of the simulation default 3D
        *rho0*:
          Ambient ISM density
        *lsnowplough:
          Include original snowplough profile
        *lcioffi:
          Include Cioffi et al profile
        """

        if isinstance(par, list):
            unit_length = 3.086e21  # 1pc in cm
            unit_time = 3.1557e16  # Gyr in s
            unit_velocity = unit_length / unit_time  # ~km/s in cm/s
            unit_density = 1.6728e-24  # g/cm3
            unit_energy_density = unit_density * unit_velocity ** 2  # erg/cm3
            unit_energy = unit_energy_density * unit_length ** 3  # erg
            E0 = 1e51 / unit_energy  # erg
            M0 = 10
        else:
            unit_length = par.unit_length
            unit_time = par.unit_time
            unit_velocity = par.unit_velocity
            unit_density = par.unit_density
            unit_energy_density = par.unit_energy_density
            unit_energy = par.unit_energy
            E0 = par.ampl_sn
            if par.lsn_mass:
                M0 = par.mass_sn
            else:
                M0 = 10
        if len(time) > 0:
            time = np.array(time)
        else:
            time = np.linspace(startt, endt, nt) + t_sedov
        xi = 2.026  # McKee/Ostriker 1988
        rho0 /= unit_density
        setattr(self, "unit_density", unit_density)
        setattr(self, "unit_length", unit_length)
        setattr(self, "unit_time", unit_time)
        setattr(self, "unit_velocity", unit_velocity)
        setattr(self, "unit_energy", unit_energy)
        setattr(self, "E0", E0)
        setattr(self, "M0", M0)
        setattr(self, "rho0", rho0)
        setattr(self, "xi", xi)
        setattr(self, "t_sedov", t_sedov)
        setattr(self, "dims", dims)
        setattr(self, "t", time)
        for key in self.__dict__.keys():
            print(key, self.__getattribute__(key))

        # Sedov-Taylor
        rst = (self.xi * self.E0 / self.rho0) ** (1.0 / (self.dims + 2.0)) * self.t ** ( 2.0 / (self.dims + 2.0) )
        rdot = (
                2.0
                / (self.dims + 2.0)
                * (self.xi * self.E0 / self.rho0) ** (1.0 / (self.dims + 2.0))
                * self.t ** (2.0 / (self.dims + 2.0) - 1)
        )
        setattr(self, "stR", rst)
        setattr(self, "stRdot", rdot)
        m_p = 1.67262158e-24
        n0 = self.rho0 * self.unit_density / m_p
        setattr(self, "n0", n0)

        def get_snpl(self, quiet):
            """
            Compute analytic radial time evolution of SN blast waves for
            comparison of snowplough solution with numerical results
            """
            # Woltier transition to momentum conservation
            vtrans = 2.0e2 * self.n0 ** (2.0 / 17.0) * self.E0 ** (1.0 / 17.0)
            ttrans = (
                             2.0
                             / (2.0 + self.dims)
                             * (self.xi * self.E0 / self.rho0) ** (1.0 / (self.dims + 2.0))
                             * vtrans ** (-1)
                     ) ** (1.0 / (1 - 2.0 / (self.dims + 2.0)))
            rtrans = (self.xi * self.E0 / self.rho0) ** (
                    1.0 / (self.dims + 2.0)
            ) * ttrans ** (2.0 / (self.dims + 2.0))
            if not quiet:
                print("Woltier transition to momentum conservation")
                print("ttrans {}, rtrans {}, vtrans {}".format(ttrans, rtrans, vtrans))
            setattr(self, "tWolt", ttrans)
            setattr(self, "rWolt", rtrans)
            setattr(self, "vWolt", vtrans)
            # snowplough
            rsnpl = self.stR.copy()
            vsnpl = self.stRdot.copy()
            isnpl = np.where(self.t > self.tWolt)[0]
            rsnpl[isnpl] = self.rWolt * (
                    8.0 / (self.dims + 2.0) * self.t[isnpl] / self.tWolt
                    - 3.0 / (self.dims + 2.0)
            ) ** (1.0 / 4.0)
            vsnpl[isnpl] = self.vWolt * (
                    8.0 / (self.dims + 2.0) * self.t[isnpl] / self.tWolt
                    - 3.0 / (self.dims + 2.0)
            ) ** (-3.0 / 4.0)
            setattr(self, "snplR", rsnpl)
            setattr(self, "snplRdot", vsnpl)

        def get_cioffi(self, quiet):
            """
            Compute analytic radial time evolution of SN blast waves for
            comparison with numerical results
            """
            vej = (400 * self.E0 / self.M0) ** 0.5
            setattr(self, "vejecta", vej)
            # pressure driven snowplough transition
            tpds = (
                           3.61e-5 * self.E0 ** (3.0 / 14.0) / self.n0 ** (4.0 / 7.0)
                   ) / 2.718281828
            rpds = (self.xi * self.E0 / self.rho0) ** (
                    1.0 / (self.dims + 2.0)
            ) * tpds ** (2.0 / (self.dims + 2.0))
            vpds = (2.0 / (2 + self.dims)) * (self.xi * self.E0 / self.rho0) ** (
                    2.0 / (2 + self.dims) - 1
            )
            if not quiet:
                print("Pressure-driven snowplough transition")
                print("tpds {}, rpds {}, vpds {}".format(tpds, rpds, vpds))
            setattr(self, "tpds", tpds)
            setattr(self, "rpds", rpds)
            setattr(self, "vpds", vpds)
            # momentum conserving snowplough transition
            tmcs = (
                    61
                    * self.vejecta ** 3
                    / self.n0 ** (3.0 / 7.0)
                    / self.E0 ** (3.0 / 14.0)
                    * self.tpds
            )
            rmcs = self.rpds * (4.0 / 3.0 * tmcs / self.tpds - 1.0 / 3.0) ** (
                    3.0 / 10.0
            )
            if not quiet:
                print("Momentum conserving snowplough")
                print("tmcs {}, rmcs {}".format(tmcs, rmcs))
            setattr(self, "tmcs", tmcs)
            setattr(self, "rmcs", rmcs)
            # Cioffi et al
            rcioffi = self.stR.copy()
            vcioffi = self.stRdot.copy()
            icioffi = np.where(self.t > self.tpds)[0]
            rcioffi[icioffi] = self.rpds * (
                    4.0 / 3.0 * self.t[icioffi] / self.tpds - 1.0 / 3.0
            ) ** (3.0 / 10.0)
            vcioffi[icioffi] = (
                    0.3
                    * self.rpds
                    * 4.0
                    / 3.0
                    / self.tpds
                    * (4.0 / 3.0 * self.t[icioffi] / self.tpds - 1.0 / 3.0) ** (-7.0 / 10.0)
            )
            jcioffi = np.where(self.t > self.tmcs)[0]
            rcioffi[jcioffi] = (
                    self.rpds
                    * (
                            4.66
                            * (self.t[jcioffi] - self.tmcs)
                            / self.tpds
                            * (1.0 - 0.779 / (self.tmcs / self.tpds) ** 0.17)
                            + (self.rmcs / self.rpds) ** 4
                    )
                    ** 0.25
            )
            vcioffi[jcioffi] = (
                    0.25
                    * self.rpds
                    * 4.66
                    / self.tpds
                    * (1.0 - 0.779 / (self.tmcs / self.tpds) ** 0.17)
                    * (
                            4.66
                            * (self.t[jcioffi] - self.tmcs)
                            / self.tpds
                            * (1.0 - 0.779 / (self.tmcs / self.tpds) ** 0.17)
                            + (self.rmcs / self.rpds) ** 4
                    )
                    ** (-0.75)
            )
            setattr(self, "CioffiR", rcioffi)
            setattr(self, "CioffiRdot", vcioffi)

        if lsnowplough:
            get_snpl(self, quiet)
        if lcioffi:
            get_cioffi(self, quiet)
        return self

if __name__ == '__main__':

    # run magnetar
    tarrm, ldip, lacc = Magnetar.run()
    print(ldip)
    deinjdt = ldip + lacc


    rho = 1e-2 * cgs.mppme
    RR = np.logspace(14., 23., 1000)
    tarr = tarrm#@np.logspace(-7, 16, 2000)

    # dyn2t = Driver_Nava_FS_t(E0=1e48, Gamma0=get_Gamma(0.7), M0=1e48 / (cgs.c ** 2 * get_Gamma(0.7)), tstart=tarr[0], Rstart = RR[0],
    #                          rho0=rho, useSpread=False, aa=-1., ncells=1,
    #                          adiabLoss=True, epsilon_e_rad=0,
    #                          ode_rtol=1e-9, ode_nsteps=1000,  eq_dthetadr=EqOpts.dthetadr_None, eq_gammaAdi=EqOpts.gamma_adi_peer,
    #                          eq_rhoprime=EqOpts.rho2_transrel,thetaMax=np.pi/2., )
    # # dyn2t_no_a = Driver_Nava_FS_t(E0=1e48, Gamma0=get_Gamma(0.7), M0=1e48 / (cgs.c ** 2 * get_Gamma(0.7)), tstart=tarr[0], Rstart = RR[0],
    # #                               rho0=rho, useSpread=False, aa=-1., ncells=1,
    # #                               adiabLoss=False, epsilon_e_rad=0,
    # #                               ode_rtol=1e-9, ode_nsteps=1000,  eq_dthetadr=EqOpts.dthetadr_None, eq_gammaAdi=EqOpts.gamma_adi_peer,
    # #                               eq_rhoprime=EqOpts.rho2_transrel,thetaMax=np.pi/2., )
    #
    # for i in tqdm(range(1,len(tarr))):
    #     # dyn2.update_ode_pars(rho=rho)
    #     rho, dlnrho1dR = rho_dlnrho1dR(tarr[i], 1e-2, None, None, None, None)
    #     deinjdt_ = 0.
    #     if (tarr[i] > tarrm.min() and tarr[i] < tarrm.max()):
    #         deinjdt_ = np.interp(tarr[i],tarrm,deinjdt)
    #
    #     dyn2t.evolove(tarr[i], rho, dlnrho1dR, deinjdt_)
    #     # dyn2t_no_a.evolove(tarr[i], rho, dlnrho1dR, deinjdt)
    #
    # plt.loglog(dyn2t.get("tburst"),dyn2t.get("Gamma")*get_beta(dyn2t.get("Gamma")))
    # # plt.loglog(dyn2t_no_a.get("tburst"),dyn2t.get("Gamma")*get_beta(dyn2t.get("Gamma")))
    # plt.show()


    # nshells =10
    # nlayers =20
    # ii = 0
    # for ish in range(nshells):
    #     for il in range(nlayers):
    #         x = il + nlayers * ish
    #         print("x={} ii={}".format(x,ii))
    #
    #         ii = ii + 1


    # sedov_taylor()
    # plot_sedov()


    # dyn2 = Driver_Nava_FS(E0=1e51, Gamma0=200, M0=1e50 / (cgs.c ** 2 * 200), Rstart=RR[0], rho0=rho, useSpread=False, aa=-1., ncells=1,
    #                       adiabLoss=True, epsilon_e_rad=0,
    #                       ode_rtol=1e-4, ode_nsteps=1000,  eq_dthetadr=EqOpts.dthetadr_None, eq_gammaAdi=EqOpts.gamma_adi_peer,
    #                       eq_rhoprime=EqOpts.rho2_transrel,thetaMax=np.pi/2., )
    # for i in tqdm(range(1,len(RR))):
    #     # dyn2.update_ode_pars(rho=rho)
    #     rho, dlnrho1dR = rho_dlnrho1dR(RR[i], 1e-2, None, None, None, None)
    #     dyn2.evolove(RR[i], rho, dlnrho1dR)
    #
    #
    # dyn2t = Driver_Nava_FS_t(E0=1e51, Gamma0=200, M0=1e50 / (cgs.c ** 2 * 200), tstart=tarr[0], Rstart = RR[0],
    #                          rho0=rho, useSpread=False, aa=-1., ncells=1,
    #                       adiabLoss=True, epsilon_e_rad=0,
    #                       ode_rtol=1e-4, ode_nsteps=1000,  eq_dthetadr=EqOpts.dthetadr_None, eq_gammaAdi=EqOpts.gamma_adi_peer,
    #                       eq_rhoprime=EqOpts.rho2_transrel,thetaMax=np.pi/2., )
    # for i in tqdm(range(1,len(tarr))):
    #     # dyn2.update_ode_pars(rho=rho)
    #     rho, dlnrho1dR = rho_dlnrho1dR(tarr[i], 1e-2, None, None, None, None)
    #     dyn2t.evolove(tarr[i], rho, dlnrho1dR)

    dyn2t2 = Driver_Nava_FS_t2(E0=1e45, Gamma0=200, M0=1e45 / (cgs.c ** 2 * 200), tstart=tarr[0], Rstart = RR[0], ctheta_j = 0.01, theta_w_j = 0.1,
                               rho0=rho, useSpread=True, aa=-1., ncells=1, theta0=0.16,
                               # ----------
                               Gamma0_ej_0=get_Gamma(0.3), M0_ej_0=1e45 / (cgs.c ** 2 * get_Gamma(0.3)),
                               tstart_ej_0=tarr[0], Rstart_ej_0=RR[0], ctheta_ej_0 = 0.10, theta_w_ej_0 = np.pi/2.,
                               rho0_ej_0=rho,
                               # ----------
                               Gamma0_ej_1=get_Gamma(0.3), M0_ej_1=1e45 / (cgs.c ** 2 * get_Gamma(0.3)),
                               tstart_ej_1=tarr[0], Rstart_ej_1=RR[0], ctheta_ej_1 = 0.10, theta_w_ej_1 = np.pi/2.,
                               rho0_ej_1=rho,
                               # ----------
                               Gamma0_ej_2=get_Gamma(0.3), M0_ej_2=1e45 / (cgs.c ** 2 * get_Gamma(0.3)),
                               tstart_ej_2=tarr[0], Rstart_ej_2=RR[0], ctheta_ej_2 = 0.10, theta_w_ej_2 = np.pi/2.,
                               rho0_ej_2=rho,
                               # ----------
                               useSpread_ej=False, aa_ej=-1., ncells_ej=1,
                               adiabLoss=True, epsilon_e_rad=0,
                               ode_rtol=1e-6, ode_nsteps=1500, thetaMax=np.pi/2.,
                               adiabLoss_ej=True, epsilon_e_rad_ej=0,
                               thetaMax_ej=np.pi / 2.,
                               eq_dthetadr=EqOpts.dthetadr_None, eq_gammaAdi=EqOpts.gamma_adi_peer,
                               eq_rhoprime=EqOpts.rho2_transrel,
                               # -----------------------------
                               use_st_rho=True, use_ctheta_lim=True
                               )
    for i in tqdm(range(1,len(tarr))):
        # dyn2.update_ode_pars(rho=rho)
        deinjdt_ = 0.
        if (tarr[i] > tarrm.min() and tarr[i] < tarrm.max()):
            deinjdt_ = np.interp(tarr[i],tarrm,deinjdt)
        rho, dlnrho1dR = rho_dlnrho1dR(tarr[i], 1e-2, None, None, None, None)
        dyn2t2.evolove(tarr[i], rho, dlnrho1dR, deinjdt_)

    # plt.semilogx(dyn.get("R"), dyn.get("Gamma"), label="P")
    # v_n_x = "R"
    # v_n_y = "beta"
    # plt.loglog(dyn2.get(v_n_x), dyn2.get(v_n_y),  color='blue', label='with raduis')
    # plt.loglog(dyn2t.get(v_n_x), dyn2t.get(v_n_y), color='red', label='with time')
    # plt.loglog(dyn2t2.get(v_n_x), dyn2t2.get(v_n_y), color='green', label='with time')
    # v_n_x = "R_ej"
    # v_n_y = "beta_ej"
    # plt.xlabel(v_n_x)
    # plt.ylabel(v_n_y)
    # plt.loglog(dyn2t2.get(v_n_x), dyn2t2.get(v_n_y), color='lime', label='with time')
    # # plt.loglog(dyn3.get("R"), dyn3.get("Gamma"), label="N2")
    # plt.legend()
    # plt.show()

    fig, axes = plt.subplots(ncols=1,nrows=3,sharex="all")
    for ii, color in zip([0,1,2], ["blue", "green", "red"]):
        gamAdi = EqOpts.gamma_adi_peer(dyn2t2.get(f"Gamma_ej_{ii}"),get_beta(dyn2t2.get(f"Gamma_ej_{ii}")))
        GammaEff = Nava_fs_rhs.GammaEff(dyn2t2.get(f"Gamma_ej_{ii}"),gamAdi)
        Etot = (dyn2t2.get(f"Gamma_ej_{ii}") - 1)* cgs.c ** 2 * \
               (dyn2t2.get(f"M2_ej_{ii}") + 1e48 / (cgs.c ** 2 * dyn2t2.get(f"Gamma_ej_{ii}")[0])) + GammaEff * dyn2t2.get(f"Eint2_ej_{ii}")
        axes[0].plot(dyn2t2.get("tburst"),Etot / 1,label=r"$E_{\rm tot}$",color=color,ls='-')
        axes[0].plot(dyn2t2.get("tburst"),GammaEff * dyn2t2.get(f"Eint2_ej_{ii}") / 1,label=r"$E_{\rm int} \Gamma_{\rm eff}$",color=color,ls='--')
        axes[0].plot(dyn2t2.get("tburst"),( dyn2t2.get(f"Gamma_ej_{ii}") - 1) * cgs.c  ** 2
                     * (dyn2t2.get(f"M2_ej_{ii}")+1e48 / (cgs.c ** 2 * dyn2t2.get(f"Gamma_ej_{ii}")[0])) / 1,label=r"$(\Gamma-1)c M_2$",color=color,ls=':')
        # plt.loglog(dyn2t2.get("R"),cgs.c*np.sqrt(dyn2t2.get("Gamma")**2 - 1)*(dyn2t2.get("M2")+gamAdi*dyn2t2.get("Eint2")/cgs.c/cgs.c) / ,label=r"$P_{\rm tot}$")
        # plt.loglog(dyn2t2.get("R"),dyn2t2.get("Eint2"),label="int")
        # plt.loglog(dyn2t2.get("R"),dyn2t2.get("Esh2"),label="sh")
        # plt.loglog(dyn2t2.get("R"),dyn2t2.get("Eint2")+dyn2t2.get("Esh2"),label="tot")
        axes[1].plot(dyn2t2.get("tburst"),dyn2t2.get(f"Gamma_ej_{ii}")*get_beta(dyn2t2.get(f"Gamma_ej_{ii}")),label=r"Momentum",color=color,ls='-',marker='.')
        axes[2].plot(dyn2t2.get("tburst"),dyn2t2.get(f"R_ej_{ii}"),label=r"Radius",color=color,ls='-',marker='.')
    axes[1].plot(tarrm,np.log10(ldip+lacc),label=r"$log_{10}(L_{\rm dip}+L_{\rm acc})$",color='gray',ls='--')
    for ax in axes:
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.legend(loc="best")
        ax.grid()
    plt.show()




    beta12 = (dyn2t2.get("beta_ej") - dyn2t2.get("beta")) / (1. - dyn2t2.get("beta_ej") * dyn2t2.get("beta"))
    coll_idx = np.argmax( dyn2t2.get("R") < dyn2t2.get("R_ej") )
    cs_jet2 = 5. * dyn2t2.get("beta") ** 2 / 9.

    fig, axes = plt.subplots(figsize=(9,3), ncols=2, nrows=1)
    ax = axes[0]
    ax.plot(dyn2t2.get("tburst"),dyn2t2.get("beta"),ls='-', color='blue',label='Jet BW')
    ax.plot(dyn2t2.get("tburst"), dyn2t2.get("beta_ej"),ls='-', color='red',label='kN BW')
    ax.axvline(x=dyn2t2.get("tburst")[np.argmax( dyn2t2.get("beta") < dyn2t2.get("beta_ej") )], color="gray", linestyle="-",label=r"$R$ at $\beta_{\rm ej}=\beta_{\rm jet}$")
    ax.axvline(x=dyn2t2.get("tburst")[np.argmax( dyn2t2.get("R") < dyn2t2.get("R_ej") )], color="gray", linestyle="--",
               label=r"$R$ at $R_{\rm ej}=R_{\rm jet}$, "
                     r"$\beta_{12}=$"+"{:.2f}".format(beta12[coll_idx])
                     +" $\Gamma_{12}=$"+"{:.2f}".format(get_Gamma(beta12[coll_idx]))
                     +" $c_s=$"+"{:.2f}".format(np.sqrt(cs_jet2[coll_idx]))
                     +" $\mathcal{M}_s=$"+"{:.2f}".format(beta12[coll_idx]/np.sqrt(cs_jet2[coll_idx])))
    ax.set_xlabel(r"time=$\int dr \beta c$ [s]")
    ax.set_ylabel(r"velocity $\beta$ [c] (solid lines)")
    ax.set_xscale("log")
    ax.set_yscale("log")
    # ax.set_xlim(3e17,1e19)
    ax.set_xlim(2e6,2e9)
    ax.set_ylim(8e-2,1.1e0)
    ax.legend(loc="best")
    ax.tick_params(axis='both', which='both', labelleft=True,
                   labelright=False, tick1On=True, tick2On=True,
                   labelsize=12,
                   direction='in',
                   bottom=True, top=True, left=True, right=True)
    ax2 = ax.twinx()
    ax2.plot(dyn2t2.get("tburst"), dyn2t2.get("R_ej"), ls='--', color='lime')
    ax2.plot(dyn2t2.get("tburst"), np.abs(dyn2t2.get("R") - dyn2t2.get("R_ej")), ls='--', color='lime', label=r"$| \Delta R |$")
    ax2.plot(dyn2t2.get("tburst"), dyn2t2.get("R_ej"), ls='--', color='red')
    ax2.legend(loc="best")
    ax2.set_ylabel("Radius [cm] (dashed lines)")
    ax2.set_xscale("log")
    ax2.set_yscale("log")

    ctheta_j = dyn2t2.kwargs["ctheta_j"] + 0.5 * (2.0 * dyn2t2.get("theta") - 2.0 * dyn2t2.kwargs["theta_w_j"])
    ctheta_ej = dyn2t2.kwargs["ctheta_ej"] + 0.5 * (2.0 * dyn2t2.get("theta_ej") - 2.0 * dyn2t2.kwargs["theta_w_ej"])

    axx = axes[1]
    axx.plot(dyn2t2.get("tburst"), ctheta_j, ls='-', color='blue', label=r'Jet')
    axx.plot(dyn2t2.get("tburst"), ctheta_ej, ls='-', color='red', label=r'ejecta')
    axx.legend(loc="best")
    axx.set_xlabel(r"time=$\int dr \beta c$ [s]")
    axx.set_ylabel(r"$\theta_c$")
    axx.set_xscale("log")
    axx.set_yscale("linear")

    axx2 = axx.twinx()
    axx2.plot(dyn2t2.get("tburst"), dyn2t2.get("M2"), ls='--', color='blue')
    axx2.plot(dyn2t2.get("tburst"), dyn2t2.get("M2_ej"), ls='--', color='red')
    axx2.legend(loc="best")
    axx2.set_ylabel("M2 [cm] (dashed lines)")
    axx2.set_xscale("log")
    axx2.set_yscale("log")

    plt.tight_layout()
    plt.savefig("time_intersection2.png",dpi=256)
    plt.show()



