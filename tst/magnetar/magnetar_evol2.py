#!/usr/bin/env python3

import h5py

#
# __author__ = "Michele Ronchi, Vsevolod Nedora"
# __copyright__ = "Copyright 2022"
# __credits__ = ["Michele Ronchi"]
# __license__ = "MIT"
# __maintainer__ = "Michele Ronchi"
# __email__ = "ronchi@ice.csic.es"

'''
This notebook contains the code to perform a parameter study of the spin-period evolution of pulsars 
interacting with supernova fallback disk as in [Ronchi et al. 2022](https://arxiv.org/abs/2201.11704). By using general 
assumptions for the pulsar spin period and magnetic field at birth, initial fallback accretion rates and including 
magnetic field decay, we find that very long spin periods ($100 \, {\rm s}$) can be reached in the presence of strong, 
magnetar-like magnetic fields ($\gtrsim 10^{14} \, {\rm G}$) and moderate initial fallback accretion rates 
($10^{22-27} \, {\rm g \, s^{-1}}$). 
In addition, we study the cases of two recently discovered periodic radio sources, the pulsar 
PSR J0901-4046 [Caleb et al. 2022](https://ui.adsabs.harvard.edu/abs/2022NatAs.tmp..123C/abstract)
($P = 75.9 \, {\rm s}$) and the radio transient GLEAM-X J162759.5-523504.3 
([Hurley-Walker et al. 2022](https://www.nature.com/articles/s41586-021-04272-x)) ($P = 1091 \, {\rm s}$), 
in light of our model. 

'''

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import random

from scipy import interpolate
from scipy.integrate import odeint
from scipy.integrate import ode
from tqdm import tqdm
from scipy.integrate import trapz
from scipy.optimize import curve_fit
from typing import Tuple

# import numpy as np

# Unit conversions.
class const:
    KPC_TO_PC = 1000.0  # Convert from [kpc] to [pc].
    KPC_TO_KM = 3.08567758e16  # Convert from [kpc] to [km].
    KPC_TO_CM = 3.08567758e21  # Convert from [kpc] to [cm].
    YR_TO_S = 3600 * 24 * 365  # Convert from [yr] to [s].
    DEG_TO_MAS = 3600000  # Convert [deg] to [mas].
    DEG_TO_RAD = np.pi / 180  # Convert [deg] to [rad].
    MJY_TO_ERG = 1.0e-26  # Convert [mJy] to [erg cm^-2 s^-1 Hz^-1].
    JY_TO_ERG = 1.0e-23  # Convert [Jy] to [erg cm^-2 s^-1 Hz^-1].

    # Physical constants.

    M_SUN = 2.0e33  # Sun mass [g].
    G = 6.67e-8  # Gravitational constant [cm^3 g^-1 s^-2].
    G_KPC_YR = (
            G / (KPC_TO_CM ** 3) * YR_TO_S ** 2
    )  # Gravitational constant [kpc^3 g^-1 yr^-2].
    C = 2.99792458e10  # Speed of light [cm/s]
    E = 4.80320425e-10  # Electric charge in [statC] = [cm^(3/2)g^(1/2)/s].
    M_E = 9.10938356e-28  # Electron's mass in [g].
    M_P = 1.67262192369e-24  # Proton mass in [g].
    SIGMA_T = (
            8.0 * np.pi / 3.0 * (E ** 2 / (M_E * C ** 2)) ** 2
    )  # Thompson cross section in [cm^2].

    grav = 6.67259e-8

from matplotlib import rcParams
rcParams["mathtext.fontset"] = "stix"
rcParams["font.family"] = "Liberation serif"
# rcParams["font.size"] = "30"
# rcParams['font.weight']='bold'
# rcParams["figure.figsize"] = "8.0, 7.0"
# rcParams["figure.autolayout"] = True
#
# rcParams["axes.linewidth"] = "1.7"
# rcParams["axes.labelpad"] = "15.0"
# rcParams["axes.titlepad"] = "15.0"
#
# rcParams["xtick.direction"] = "in"
# rcParams["xtick.top"] = True
# rcParams["xtick.major.pad"] = "10.0"
# rcParams["xtick.minor.pad"] = "10.0"
# rcParams["xtick.major.size"] = "10.0"
# rcParams["xtick.major.width"] = "1.7"
# rcParams["xtick.minor.size"] = "5.0"
# rcParams["xtick.minor.width"] = "1.7"
# rcParams["xtick.labelsize"] = "30"
#
# rcParams["ytick.direction"] = "in"
# rcParams["ytick.right"] = True
# rcParams["ytick.major.pad"] = "10.0"
# rcParams["ytick.minor.pad"] = "10.0"
# rcParams["ytick.major.size"] = "10.0"
# rcParams["ytick.major.width"] = "1.7"
# rcParams["ytick.minor.size"] = "5.0"
# rcParams["ytick.minor.width"] = "1.7"
# rcParams["ytick.labelsize"] = "30"


''' ---------------------------------------------------------------------------------------------------------------- '''
# TODO Add (i) Bext/Bint for ellipticity; (ii) NS_mass/Radius evolution; (iii) Collapse criterion (iv) does Ldip = f(alpha)?
class MagnetarRonchi22:

    def __init__(self):
        pass

    # -------------------------------| RHS |---------------------
    def __call__(self, t : float, params : np.ndarray, pars : dict):
        """
        Derivative of the magnetic field and spin frequency for the three phases: direct accretion, propeller and
        ejector. See Ronchi et al. (2022), Metzger et al. (2018) for details.

        Args:
            t (float): time variable in [s],
            initial_cond (np.ndarray): current value of the dipolar component of the magnetic field at the
            magnetic pole measured in [G] and value of the spin frequency measured omega in [1/s].
            B0 (float): initial value of the dipolar component of the magnetic field at the
            magnetic pole measured in [G].
            Mdot_d0 (float): accretion rate in [g s^-1].
            alpha(float): index of the power law for the time evolution of the accretion rate.
            t_disrupt(float): time in [s] when the disk is disrupted because r_in > r_out.

        Returns:
            (float): Derivative of the magnetic field and spin frequency for the three phases:
            direct accretion, propeller and ejector [s^-2].
        """
        B = params[0]
        omega = params[1]
        alpha = params[2]
        # mdot = initial_cond[3]

        # magnetic field evolution
        Bdot = MagnetarRonchi22.field_derivative(B=B, B_initial=pars["B0"],
                                                 characteristic_length=pars["characteristic_length"],
                                                 sigma_conductivity=pars["sigma_conductivity"],
                                                 characteristic_electron_dens=pars["characteristic_electron_dens"])
        # evolution of the angle
        alphadot = MagnetarRonchi22.alpha_derivatives(B=B, omega=omega, alpha=alpha,
                                                      NS_radius=pars["NS_radius"],NS_inertia=pars["NS_inertia"],
                                                      k_coefficients=pars["k_coefficients"])

        # viscous timescale
        t_vis = MagnetarRonchi22.t_vis(pars)#(NS_mass=pars["NS_mass"],NS_radius=pars["NS_radius"],
                                            #characteristic_disk_temp=pars["characteristic_disk_temp"],
                                            #circularization_disk_radius=pars["circularization_disk_radius"],
                                            #method_t_vis=pars["method_t_vis"])

        # eddington lum
        Ledd = MagnetarRonchi22.Ledd(NS_mass=pars["NS_mass"], opacity=pars["disk_opacity"])

        # total accretion rate
        Mdot_d = MagnetarRonchi22.disk_accretion_rate_t(Mdot_d0=pars["Mdot_d0"],omega=omega,t=t,t_vis=t_vis,
                                                        disk_acc_rate_pl=pars["disk_acc_rate_pl"],
                                                        method_disk_acc=pars["method_disk_acc"])
        # eddington limited accretion
        Mdot = MagnetarRonchi22.disk_accretion_rate_in_t(
            B=B, omega=omega, Mdot_d=Mdot_d, t=t, t_disrupt=pars["t_disrupt"],
            NS_mass=pars["NS_mass"], NS_radius=pars["NS_radius"], Ledd=Ledd)

        # char. radii
        r_m = MagnetarRonchi22.radius_magnetospheric(B=B, Mdot=Mdot,NS_radius=pars["NS_radius"],NS_mass=pars["NS_mass"])
        r_c = MagnetarRonchi22.radius_corotation(omega=omega,NS_mass=pars["NS_mass"])
        r_lc = MagnetarRonchi22.radius_lc(omega=omega)
        # Compute the inner and outer radius of the disk.
        r_in = np.minimum(r_m, r_lc)
        r_out = MagnetarRonchi22.disk_outer_radius_t(t=t,t_vis=t_vis,
                                                     circularization_disk_radius=pars["circularization_disk_radius"],
                                                     disk_rout_evol_idx=pars["disk_rout_evol_idx"])
        # Check if the field could be buried in the initial stages (before the field decays
        buried = False
        if (r_m < pars["NS_radius"]):
            buried = True
            print(f"buried at t={t}")

        # Check if the disk is disrupted because r_in > r_out and save the time when this happens.
        if (r_in > r_out) and (pars["t_disrupt"] < 0):
            pars["t_disrupt"] = t
            print(f"disk disrupted at t={t}")

        # Following Gibson+ (capping)
        if ((pars["k"]>0)and(r_m > pars["k"] * r_lc)):
            r_m = pars["k"] * r_lc

        if not (np.isfinite(r_m)&np.isfinite(r_c)&np.isfinite(r_lc)):
            raise ValueError()

        # torques
        n_acc = MagnetarRonchi22.Nacc(
            omega=omega,Mdot_din=Mdot,r_m=r_m,r_lc=r_lc,r_corot=r_c,
            NS_radius=pars["NS_radius"],NS_mass=pars["NS_mass"],NS_inertia=pars["NS_inertia"],
            NS_critical_beta=pars["NS_critical_beta"],method_nacc=pars["method_nacc"])
        n_dip = MagnetarRonchi22.Ndip(
            B=B,omega=omega,r_m=r_m,r_lc=r_lc,alpha=alpha,NS_radius=pars["NS_radius"],NS_inertia=pars["NS_inertia"],
            k_coefficients=pars["k_coefficients"],method_ndip=pars["method_ndip"])
        n_grav = MagnetarRonchi22.Ngrav(omega=omega,alpha=alpha,
                                        NS_inertia=pars["NS_inertia"],ns_ellipticity=pars["ns_ellipticity"])

        if not (np.isfinite(n_acc) and np.isfinite(n_dip) and np.isfinite(n_grav) and np.isfinite(omega) and omega>0):
            raise ValueError()

        # luminocities
        ldip = - n_dip * omega
        lprop = - n_acc*omega - const.grav*pars["NS_mass"]*Mdot/r_m # From Gompertz et al. (2014)
        if lprop < 0: lprop = 0.

        # spin-down
        omegadot = (n_dip + n_acc + n_grav) / pars["NS_inertia"]

        # output
        pars["t_vis"] = t_vis
        pars["omegadot"] = omegadot
        pars["Bdot"] = Bdot
        pars["alphadot"] = alphadot
        pars["Ledd"] = Ledd
        pars["Mdot"] = Mdot
        pars["Mdot_d"] = Mdot_d
        pars["r_m"] = r_m
        pars["r_c"] = r_c
        pars["r_lc"] = r_lc
        pars["r_in"] = r_in
        pars["r_out"] = r_out
        pars["n_acc"] = n_acc
        pars["n_dip"] = n_dip
        pars["n_grav"] = n_grav
        pars["ldip"] = ldip
        pars["lacc"] = lprop
        derivatives = np.array([Bdot, omegadot, alphadot])
        return derivatives
    # ---------------------| Utils |----------------------------
    @staticmethod
    def E_rot(omega:float, NS_inertia:float) -> float:
        """
        Rotational energy of the NS
        """
        out=0.5*NS_inertia*omega**2
        return out
    @staticmethod
    def E_bind(NS_mass:float, NS_radius:float) -> float:
        """
        Binding energy of the NS
        Prescription from Lattimer and Prakash (2001)

        """
        num = const.grav*NS_mass
        den = NS_radius*const.C**2-0.5*const.grav*NS_mass
        out = 0.6*NS_mass*const.C**2*num/den
        return out
    @staticmethod
    def Ledd(NS_mass:float,opacity)->float:
        # Eddington luminosity assuming Thompson scattering in [erg / s].
        # eddigton_lum = 4.0 * np.pi * const.G * NS_mass * const.M_P * const.C / const.SIGMA_T
        eddigton_lum = 4.0 * np.pi * const.G * NS_mass * const.M_P * const.C / const.SIGMA_T
        eddigton_lum = 4.0 * np.pi * const.G * NS_mass * const.C / opacity
        return eddigton_lum
    @staticmethod
    def t_vis(pars:dict) -> float:
        NS_mass=pars["NS_mass"]
        NS_radius=pars["NS_radius"]
        method_t_vis=pars["method_t_vis"]
        # typical initial viscousity timescale of the disk in [s] see Menou et al. 2001.
        # method_t_vis = "menou" # "menou" "gompertz" "rrayand"
        if method_t_vis == "menou":
            characteristic_disk_temp=pars["characteristic_disk_temp"]
            circularization_disk_radius=pars["circularization_disk_radius"]
            t_v = 2080.0 * (characteristic_disk_temp / 1.0e6) ** (-1.0) * (circularization_disk_radius / 10 ** 8) ** (1.0 / 2.0) # s
        elif method_t_vis == "gompertz":
            disk_radius = pars["disk_radius"]#5.e8
            disk_alpha = pars["disk_alpha"]#0.1  # disk viscosity parameter
            disk_cs = pars["disk_cs"]#1e7
            t_v = disk_radius / (3. * disk_alpha * disk_cs * disk_radius) # s
        elif method_t_vis == "rrayand": # TODO ADD TEMPERATURE DEPENDENCE HERE!
            disk_radius = pars["disk_radius"]#5.e8
            disk_alpha = pars["disk_alpha"]#0.1  # disk viscosity parameter
            disk_askpect_ratio = pars["disk_askpect_ratio"]#0.3 # Aspect ratio H/R
            disk_H = disk_radius * disk_askpect_ratio
            disk_omega_kep = OmegaKep = (const.grav * NS_mass / NS_radius**3)**0.5
            disc_cs = disk_H * disk_omega_kep * (NS_radius/disk_radius)**1.5
            t_v = disk_radius**2 / (3. * disk_alpha * disc_cs * disk_H) # s
        else:
            raise KeyError()
        return t_v
    @staticmethod
    def omega_keplerian(r: float, NS_mass) -> float:
        """
        Keplerian orbital frequency as a function of the radius.

        Args:
            r (float): distance from the star in [cm].

        Returns:
            (float): keplerian orbital frequency in [1/s].
        """

        omega_k = (const.G * NS_mass / (r ** 3)) ** (1.0 / 2.0)

        return omega_k
    # ----------------------------------------------------------------------
    @staticmethod
    def timescale_Ohm(L: float, sigma: float) -> float:
        """
        Calculating the ohmic diffusion timescale for a given conductivity and characteristic magnetic
        field length scale. Note that for our purposes, we neglect the fact that both quantities can
        vary significantly with depth inside the neutron star, and we simply use effective quantities
        that reflect the ohmic diffusion process.

        Args:
            L (float): characteristic length scale on which the magnetic field varies, measured in [cm].
            sigma (float): conductivity of the dominating dissipative process, measured in [1/s].

        Returns:
            (float): ohmic diffusion timescale in [yr].
        """

        tau_Ohm = 4 * np.pi * sigma * L ** 2 / (const.C ** 2)

        return tau_Ohm
    @staticmethod
    def timescale_Hall(B: float, L: float, n_e: float) -> float:
        """
        Calculating the Hall timescale for a given field strength, characteristic magnetic field length
        scale and electron density. Note that for our purposes, we neglect the fact that all quantities can
        vary significantly with depth inside the neutron star, and we simply use effective quantities
        that reflect the conservative Hall process.
        B will be identified with the initial dipolar magnetic field components at the pulsars' pole.

        Args:
            B (float): (local) magnetic field strength, measured in [G].
            L (float): characteristic length scale on which the magnetic field varies, measured in [cm].
            n_e (float): electron density, measured in [g/cm^3].

        Returns:
            (float): Hall timescale in [yr].
        """

        tau_Hall = 4 * np.pi * const.E * n_e * L ** 2 / (const.C * B)

        return tau_Hall

        # Example of Hall and Ohmic timescale values for B0 = 10^14 G in [yr].
        # print(timescale_Hall(1.e14, L, n_e) / const.YR_TO_S)
        # print(timescale_Ohm(L, sigma) / const.YR_TO_S)
    @staticmethod
    def field_derivative(B: float, B_initial: float,
                         characteristic_length:float, sigma_conductivity:float,characteristic_electron_dens:float) -> float:
        """
        Calculating the change in the magnetic field strength of a pulsar based on a simplified
        differential equation (see eq. (18) of Aguilera et al. (2008)) that captures the
        characteristics of more complicated numerical simulations of pulsar magnetic field
        evolution, i.e., at early timescales the Hall evolution dominates, while at late times
        the exponential magnetic field decay due to Ohmic dissipation kicks in. Note that as
        explained in Aguilera et al. (2008) the Hall timescale corresponds to that of the initial
        field strength.

        Args:
            B (float): pulsar's magnetic field magnitudes evolving with time, measured in [G].
            B_initial (float): pulsar's initial magnetic field magnitudes, measured in [G].

        Returns:
            (float): magnetic field derivatives for a simulated pulsars in [G/yr].
        """

        tau_Ohm = MagnetarRonchi22.timescale_Ohm(characteristic_length, sigma_conductivity)
        tau_Hall = MagnetarRonchi22.timescale_Hall(B_initial, characteristic_length, characteristic_electron_dens)

        B_deriv = -B / tau_Ohm - B ** 2 / (tau_Hall * B_initial)

        return B_deriv
    # -------------------| RADII |----------------------------------------
    @staticmethod
    def radius_magnetospheric(B: float, Mdot: float, NS_radius:float, NS_mass:float) -> float:
        """
        Magnetospheric radius as a function of the magnetic field of the neutron star and the accretion rate.

        Args:
            B (float): value of the dipolar component of the magnetic field at the
            magnetic pole measured in [G].
            Mdot (float): inflow rate in [g s^-1].

        Returns:
            (float): magnetosperic radius in [cm].
        """

        mu = B * NS_radius ** 3 / 2.0
        r_m = 0.5 * (mu ** 4 / (2 * const.G * NS_mass * Mdot ** 2)) ** (1.0 / 7.0)
        # out  = self.mu**(4./7) * (self.gravconst*self.NS_mass)**(-1./7) * Mdot**(-2./7)

        return r_m
    @staticmethod
    def radius_lc(omega: float) -> float:
        """
        Light cylinder radius as a function of the spin frequency.

        Args:
            omega (float): value of the spin frequency measured in [1/s].

        Returns:
            (float): light cylinder radius in [cm].
        """

        r_lc = const.C / omega
        # out = self.lightspeed/Omega

        return r_lc
    @staticmethod
    def radius_corotation(omega: float, NS_mass:float) -> float:
        """
        Corotation radius as a function of the spin frequency.

        Args:
            omega (float): value of the spin frequency measured in [1/s].

        Returns:
            (float): corotation radius in [cm].
        """

        r_cor = (const.G * NS_mass / (omega ** 2)) ** (1.0 / 3.0)
        # out = (self.gravconst * self.NS_mass/ Omega**2)**(1./3)

        return r_cor
    # -------------------| Accretion |------------------------------------
    @staticmethod
    def accretion_rate_Edd_rmag(B: float, NS_radius:float, NS_mass:float, Ledd:float) -> float:
        """
        Eddington limit of the accretion rate at the magnetospheric radius.

        Args:
            B (float): value of the dipolar component of the magnetic field at the magnetic pole measured in [G].

        Returns:
            (float): Eddington limit of the accretion rate at the magnetospheric radius in [g s^-1].
        """

        # Magnetic moment.
        mu = B * NS_radius ** 3 / 2.0

        Mdot_Edd = (Ledd / (const.G * NS_mass)) ** (7.0 / 9.0) \
                   * (mu ** 4.0 / (2.0 * const.G * NS_mass)) ** (1.0 / 9.0)

        return Mdot_Edd
    @staticmethod
    def accretion_rate_Edd_rlc(omega: float,NS_mass:float,Ledd:float) -> float:
        """
        Eddington limit of the accretion rate at the light cylinder radius.

        Args:
            omega (float): value of the spin frequency measured in [1/s].

        Returns:
            (float): Eddington limit of the accretion rate at the light cylinder radius in [g s^-1].
        """

        # Light cylinder radius.
        r_lc = MagnetarRonchi22.radius_lc(omega)

        Mdot_Edd = 2. * Ledd * r_lc / (const.G * NS_mass)

        return Mdot_Edd
    @staticmethod
    def disk_accretion_rate_t(Mdot_d0: float, omega: float, t:float,t_vis:float,
                              disk_acc_rate_pl:float,method_disk_acc:str) -> float:
        """
        Accretion rate inside the disk as a function of time.

        Args:
            Mdot_d0 (float): initial maximum disk accretion rate in [g s^-1].
            alpha (float): index of the power law for the time evolution of the disk accretion rate.
            t (float): time in [s].

        Returns:
            (float): disk accretion rate in [g s^-1] as a function of time.
        """

        # [my] Original
        if (method_disk_acc == "pl"):
            Mdot_d = Mdot_d0 * (1 + t / t_vis) ** (-disk_acc_rate_pl)

        elif (method_disk_acc == "exp"):
            # [my] Accretion rate on the NS
            # [my] Eq. (13) from Zhang and Meszaros 2001
            Mdot_d = Mdot_d0 * np.exp(-t / t_vis)
            mdot_floor=1.e-10
            if hasattr(t, '__len__'):
                Mdot_d[Mdot_d<mdot_floor] = mdot_floor
            else:
                if Mdot_d < mdot_floor: Mdot_d = mdot_floor
        elif (method_disk_acc == "gibson"):
            r_lc = MagnetarRonchi22.radius_lc(omega)
        else:
            raise KeyError()

        return Mdot_d
    @staticmethod
    def disk_outer_radius_t(t:float,t_vis:float,circularization_disk_radius:float,disk_rout_evol_idx:float) -> float:
        """
        Accretion rate inside the disk as a function of time.

        Args:
            t (float): time in [s].

        Returns:
            (float): disk accretion rate in [g s^-1] as a function of time.
        """

        r_out = circularization_disk_radius * (1 + t / t_vis) ** disk_rout_evol_idx

        return r_out
    @staticmethod
    def disk_accretion_rate_in_t( B: float, omega: float, Mdot_d : float, t: float, t_disrupt: float,
                                  NS_mass : float, NS_radius : float, Ledd : float
    ) -> float:
        """
        Accretion rate in the inner boundary of the disk as a function of time.

        Args:
            B (float): value of the dipolar component of the magnetic field at the magnetic pole measured in [G].
            omega (float): value of the spin frequency measured in [1/s].
            Mdot_d0 (float): initial maximum disk accretion rate in [g s^-1].
            alpha (float): index of the power law for the time evolution of the disk accretion rate.
            t (float): time in [s].
            t_disrupt (float): time in [s] when the disk is disrupted because r_in > r_out.

        Returns:
            (float): accretion rate at the inner radius of the disk in [g s^-1] as a function of time.
        """

        # Compute Eddington limits on the accretion rate [g s^-1]
        Mdot_Edd_mag = MagnetarRonchi22.accretion_rate_Edd_rmag(B=B,NS_radius=NS_radius,NS_mass=NS_mass,Ledd=Ledd)
        Mdot_Edd_lc = MagnetarRonchi22.accretion_rate_Edd_rlc(omega=omega,NS_mass=NS_mass,Ledd=Ledd)

        # The disk inner radius is the smaller one between the magnetospheric radius and the light cylinder radius.
        Mdot_Edd = min(Mdot_Edd_mag, Mdot_Edd_lc)

        Mdot_din = 0.0

        # Limit the accretion rate at the inner disk edge to the Eddington limit.
        if Mdot_d < Mdot_Edd:
            Mdot_din = Mdot_d
        elif Mdot_d >= Mdot_Edd:
            Mdot_din = Mdot_Edd
        else:
            raise ValueError()
        # When the disk is disrupted because r_in > r_out set the accretion rate to 0 afterwards.
        if (t_disrupt is not None) and (t_disrupt > 0):
            if t > t_disrupt:
                Mdot_din = 1e-20
        if Mdot_din == 0:
            raise ValueError()


        return Mdot_din

    # ---------------| Inclination angle evolution |------------------------
    @staticmethod
    def alpha_derivatives(B:float, omega:float, alpha:float, NS_radius:float, NS_inertia:float,
                          k_coefficients: list[float, float, float]) -> float:
        '''
        https://iopscience.iop.org/article/10.3847/1538-4357/ab498c/pdf
        TODO find the r_m, r_c, r_lc depedence
        k_coefficients is from Philippov 2014
        Eqs 71 from Pons2019
        '''
        mu = B * NS_radius**3
        # dalphadt = - mu ** 2 * omega ** 2 / const.C ** 3 * k_coefficients[2] * np.sin(alpha) * np.cos(alpha)
        beta = 1.0 / 4.0 * NS_radius ** 6 / (NS_inertia * const.C ** 3)
        beta_1 = beta * k_coefficients[2] * np.sin(alpha) * np.cos(alpha)
        dalphadt = - mu ** 2 * omega ** 2 / const.C ** 3 * beta_1 # TODO doubple check it with Pons
        return dalphadt
    # -----------------------| Torques |-----------------------------------
    @staticmethod
    def Ngrav(omega:float, alpha:float, NS_inertia:float, ns_ellipticity:float) -> float:
        '''
        (2.0 * graviational_constant * moi ** 2 * epsilon_b ** 2 / (5.0 * speed_of_light ** 5)) * omega ** 6 * np.sin(chi)**2 * (
                    1.0 + 15.0 * np.sin(chi)**2)
        :param omega:
        :param alpha:
        :return:
        '''
        # n_grav = - 32./5. \
        #          * 6.67259e-8 \
        #          * NS_inertia**2 \
        #          * ns_ellipticity**2 \
        #          * omega**5 \
        #          / 299792458e2**5
        # TODO check if it is in agreement with above
        beta_1 = np.sin(alpha) ** 2 * (1 + 15. * np.sin(alpha)**2)
        n_grav = - 2/5 * const.G * NS_inertia**2 * ns_ellipticity**2 * omega**5 / const.C**5 * beta_1
        return n_grav
    @staticmethod
    def Ndip(B:float, omega:float, r_m:float, r_lc:float, alpha:float, NS_radius:float,NS_inertia:float,
             k_coefficients: list[float, float, float], method_ndip:str) -> float:
        # Auxiliary quantity beta as defined in eq. (72) of Pons & Vigano (2019).
        beta_ax = 1.0 / 4.0 * NS_radius ** 6 / (NS_inertia * const.C ** 3)
        mu = B * NS_radius**3
        if method_ndip == "original":
            # Eq.70 from Pons2019
            # Incorporate the inclination angle dependence into a constant.
            beta_1 = beta_ax * (k_coefficients[0] + k_coefficients[1] * np.sin(alpha) ** 2)

            if r_m <= NS_radius: # TODO doubple check it with Pons
                n_dip = - (r_lc / NS_radius) ** 2 * beta_1 * B ** 2 * omega ** 3 * NS_inertia
            elif (r_m > NS_radius) & (r_m < r_lc):
                n_dip = - (r_lc / r_m) ** 2 * beta_1 * B ** 2 * omega ** 3 * NS_inertia
            elif r_m >= r_lc:
                n_dip = - beta_1 * B ** 2 * omega ** 3 * NS_inertia
            else:
                raise ValueError()
        elif method_ndip == "gompertz":
            ################################################################
            ## Gompertz uses the disk's alfven radius
            ## in the Bucciantini prescription,
            ## but it should actually be the alfven radius of the NS wind...
            ################################################################
            mu = B * NS_radius**3
            n_dip = - 2./3. * mu**2 * omega**3 / const.C**3 * (r_lc/r_m)**3
        elif method_ndip == "rryand":
            #  Standard dipole spindown, no wind or disk
            mu = B * NS_radius**3
            n_dip = - 1./6. * mu**2 * omega**3 / const.C**3
        else:
            raise KeyError()
        return n_dip
    @staticmethod
    def Nacc(omega:float, Mdot_din:float, r_m:float, r_lc:float, r_corot:float,NS_radius:float,NS_mass:float,NS_inertia:float,
             NS_critical_beta:float, method_nacc:str) -> float:

        omega_k_rns = MagnetarRonchi22.omega_keplerian(r=NS_radius,NS_mass=NS_mass)
        omega_k_rm = MagnetarRonchi22.omega_keplerian(r=r_m,NS_mass=NS_mass)
        omega_k_rlc = MagnetarRonchi22.omega_keplerian(r=r_lc,NS_mass=NS_mass)

        if (method_nacc == "original"):
            if r_m <= NS_radius:
                n_acc = Mdot_din * NS_radius ** 2  * (omega_k_rns - omega)
            elif (r_m > NS_radius) & (r_m < r_lc):
                n_acc = Mdot_din * r_m ** 2  * (omega_k_rm - omega)
            elif r_m >= r_lc:
                n_acc = Mdot_din * r_lc ** 2 * (omega_k_rlc - omega)
            else:
                raise ValueError()

        elif (method_nacc == "rryand"):
            # Accretion torque, taking into account the propeller model
            #         Eq (6-7) of Gompertz et al. 2014
            fastness = (r_m / r_corot)**1.5
            ## Eq. (6)
            n_acc = (1. - fastness) * (const.grav * NS_mass * r_m)**0.5 * Mdot_din
            if r_m<=NS_radius:
                OmegaKep = (const.grav * NS_mass / NS_radius**3)**0.5
                n_acc = ((1. - omega/OmegaKep) * (const.grav*NS_mass*r_m)**0.5 * Mdot_din)

        else:
            raise KeyError()

        # from rayand:
        ###############################################
        ## Check for inhibition by bar-mode instability
        ## with beta = T/|W| parameter (Gompertz 2014)
        ###############################################
        beta_bar_instab = MagnetarRonchi22.E_rot(omega=omega,NS_inertia=NS_inertia) \
                          / abs(MagnetarRonchi22.E_bind(NS_mass=NS_mass,NS_radius=NS_radius))
        if (beta_bar_instab > NS_critical_beta):
            print("beta_bar_instab > NS_critical_beta")
            n_acc *= 0.

        return n_acc

def update_ode_pars(pars:dict, **kwargs):
    if len(kwargs.keys()) > 0:
        for key in kwargs.keys():
            pars[key] = kwargs[key]
    return pars

def plot_ronchi_magnetar(time: np.ndarray, solution : np.ndarray, v_ns : list, pars : dict) -> None:
    # in seconds
    # time *=const.YR_TO_S
    tmin = time.min();#10 / const.YR_TO_S
    t_max = time.max()
    # --------------------------------
    r_m = solution[:, v_ns.index("r_m")]
    r_c = solution[:, v_ns.index("r_c")]
    r_lc = solution[:, v_ns.index("r_lc")]
    # Compute the inner and outer radius of the disk.
    r_in = solution[:, v_ns.index("r_in")]
    r_out = solution[:, v_ns.index("r_out")]
    Mdot_d = solution[:, v_ns.index("Mdot_d")]
    Mdot_din = solution[:, v_ns.index("Mdot")]
    P = 2.0 * np.pi / solution[:, v_ns.index("omega")]
    ldip = solution[:, v_ns.index("ldip")]
    lacc = solution[:, v_ns.index("lacc")]
    # --------------------------------


    # def plot_magnetar():
    # Find the time when the star enters the propeller phase.
    t_prop = time[(r_m > r_c) & (r_m < r_lc)]
    if len(t_prop) == 0 : t_prop = t_max
    else: t_prop = t_prop[0]


    # Plot the spin-period evolution.
    mask1 = r_m >= r_lc
    mask2 = r_m < r_lc

    fig, (ax1, ax2, ax3, ax4) = plt.subplots(
        4, 1, sharex="all", gridspec_kw={"wspace": 0, "hspace": 0.05}, figsize=(6, 10)
    )

    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.set_xlim(tmin, t_max)
    ax1.set_ylabel(r"$\dot{M}$ [g s$^{-1}$]")

    ax2.set_xscale("log")
    ax2.set_yscale("log")
    ax2.set_ylabel(r"$r$ [cm]")

    ax3.set_xscale("log")
    ax3.set_yscale("log")
    ax3.set_xlabel(r"$t$ [yr]")
    ax3.set_ylabel(r"$P$ [s]")

    ax4.set_xscale("log")
    ax4.set_yscale("log")
    ax4.set_ylabel(r"Luminocity [erg s$^{-1}$]")

    ax1.plot(
        time, Mdot_d, linestyle="--", linewidth=1.5, color="black",
    )
    ax1.text(
        tmin*10.4,
        Mdot_d[0] / 20,
        r"$\dot{M}_{\rm disk}$",
        fontsize=14,
        rotation=0,
        color="black",
        )
    ax1.plot(
        time, Mdot_din, linestyle="-", linewidth=2, color="black",
    )
    ax1.text(
        tmin*1.4,
        Mdot_din[0] / 20,
        r"$\dot{M}_{\rm disk, in}$",
        fontsize=14,
        rotation=0,
        color="black",
        )

    # ax1.axvline(
    #     x=t_prop, color="black", linestyle="--", linewidth=1, rasterized=True,
    # )
    # ax1.axvspan(
    #     xmin=10 / const.YR_TO_S,
    #     xmax=t_prop,
    #     facecolor="gold",
    #     alpha=0.3,
    #     rasterized=True,
    # )
    # ax1.axvspan(
    #     xmin=t_prop, xmax=t_max, facecolor="tab:orange", alpha=0.3, rasterized=True,
    # )

    ax2.plot(
        time, r_out, linestyle="-", linewidth=4, color="tab:red",
    )
    ax2.plot(
        time, r_in, linestyle="-", linewidth=3, color="tab:red",
    )
    ax2.plot(
        time, r_lc, linestyle="--", linewidth=1., color="black",
    )
    ax2.plot(
        time, r_c, linestyle="--", linewidth=1., color="black",
    )
    ax2.plot(
        time, r_m, linestyle="-.", linewidth=1.5, color="black",
    )

    # ax2.axvline(
    #     x=t_prop, color="black", linestyle="--", linewidth=1, rasterized=True,
    # )
    # ax2.axvspan(
    #     xmin=10 / const.YR_TO_S,
    #     xmax=t_prop,
    #     facecolor="gold",
    #     alpha=0.3,
    #     rasterized=True,
    # )
    # ax2.axvspan(
    #     xmin=t_prop, xmax=t_max, facecolor="tab:orange", alpha=0.3, rasterized=True,
    # )

    ax2.axhline(
        y=pars["NS_radius"], color="black", linestyle="-", linewidth=1, rasterized=True,
    )
    ax2.axhspan(
        ymin=0.0, ymax=pars["NS_radius"], facecolor="tab:gray", hatch="/", alpha=1, rasterized=True,
    )
    ax2.text(tmin*10, r_m[0] * 1.5, r"$r_{\rm m}$", fontsize=12, rotation=0, color="black")
    ax2.text(tmin*10, r_lc[0] * 2.5, r"$r_{\rm lc}$", fontsize=12, rotation=0, color="black")
    ax2.text(tmin*10, r_c[0] * 2.5, r"$r_{\rm cor}$", fontsize=12, rotation=0, color="black")
    ax2.text(tmin*1, r_in[np.where(time==tmin)[0]]*1.2, r"$r_{\rm in}$", fontsize=12, rotation=30, color="tab:red")
    ax2.text(tmin*1, r_out[np.where(time==tmin)[0]]*1.2, r"$r_{\rm out}$", fontsize=12, rotation=25, color="tab:red")
    ax2.text(
        tmin*100,
        pars["NS_radius"] * 1.5,
        r"Neutron star surface",
        fontsize=14,
        rotation=0,
        color="black",
        )

    ax3.text(tmin*10, 5.0e2, r"Ejector", fontsize=12, rotation=0, color="black")
    ax3.text(t_prop * 10, 5.0, r"Propeller", fontsize=12, rotation=0, color="black")

    ax3.plot(
        np.ma.masked_where(mask1, time),
        np.ma.masked_where(mask1, P),
        color="black",
        linestyle="-",
        linewidth=4,
        rasterized=True,
    )
    ax3.plot(
        np.ma.masked_where(mask2, time),
        np.ma.masked_where(mask2, P),
        color="black",
        linestyle="--",
        linewidth=4,
        rasterized=True,
    )

    # ax3.axvline(
    #     x=t_prop, color="black", linestyle="--", linewidth=1, rasterized=True,
    # )
    # ax3.axvspan(
    #     xmin=10 / const.YR_TO_S,
    #     xmax=t_prop,
    #     facecolor="gold",
    #     alpha=0.3,
    #     rasterized=True,
    # )
    # ax3.axvspan(
    #     xmin=t_prop, xmax=t_max, facecolor="tab:orange", alpha=0.3, rasterized=True,
    # )

    ax4.plot(
        np.ma.masked_where(mask1, time),
        np.ma.masked_where(mask1, ldip),
        color="blue",
        linestyle="-",
        linewidth=4,
        rasterized=True,
    )
    ax4.plot(
        np.ma.masked_where(mask2, time),
        np.ma.masked_where(mask2, ldip),
        color="blue",
        linestyle="--",
        linewidth=4,
        rasterized=True,
    )
    ax4.plot(
        np.ma.masked_where(mask1, time),
        np.ma.masked_where(mask1, lacc),
        color="red",
        linestyle="-",
        linewidth=4,
        rasterized=True,
    )
    ax4.plot(
        np.ma.masked_where(mask2, time),
        np.ma.masked_where(mask2, lacc),
        color="red",
        linestyle="--",
        linewidth=4,
        rasterized=True,
    )
    ax4.plot(
        np.ma.masked_where(mask1, time),
        np.ma.masked_where(mask1, ldip+lacc),
        color="black",
        linestyle="-",
        linewidth=4,
        rasterized=True,
    )
    ax4.plot(
        np.ma.masked_where(mask2, time),
        np.ma.masked_where(mask2, ldip+lacc),
        color="black",
        linestyle="--",
        linewidth=4,
        rasterized=True,
    )

    ax4.text(tmin*10, ldip[0] / 10, r"$L_{\rm dip}$", fontsize=12, rotation=0, color="blue")
    ax4.text(tmin*10, lacc[0] / 700, r"$L_{\rm acc}$", fontsize=12, rotation=0, color="red")
    ax4.text(tmin*10, (ldip+lacc)[0] / 300, r"$L_{\rm tot}$", fontsize=12, rotation=0, color="black")

    for axi in [ax1,ax2,ax3,ax4]:
        axi.axvline(
            x=t_prop, color="black", linestyle="--", linewidth=1, rasterized=True,
        )
        axi.axvspan(
            xmin=10 / const.YR_TO_S,
            xmax=t_prop,
            facecolor="gold",
            alpha=0.3,
            rasterized=True,
        )
        axi.axvspan(
            xmin=t_prop, xmax=t_max, facecolor="tab:orange", alpha=0.3, rasterized=True,
        )
        axi.grid()
        axi.set_xlim(time.min(),time.max())
    ax4.set_xlabel("time [s]")
    fig.savefig("propeller_evolution_ex_Bdecay.pdf", bbox_inches="tight")
    # plt.savefig()
    plt.show()

def run_ronchi_magnetar(time_grid : np.ndarray, pars : dict, v_ns : list) -> np.ndarray:
    """
        Evolving the neutron stars' spin period in case of interaction with a fallback disk.

        Args:
            B0 (float): initial value of the magnetic field magnitude, measured in [G].
            P_initial (float): pulsars' initial rotation periods, measured in [s].
            Mdot_d0 (float): initial accretion rate at the beginning of fall-back in [g s^-1].
            alpha (float): index of the power law for the time evolution of the accretion rate.
            t_max (float): maximum time for the evolution in [yr].

        Returns:
            (np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray,
            np.ndarray, np.ndarray, bool, bool, float): time array in [yr], spin period as a
            function of time in [s], derivative of the spin period in [s/s], magnetic field as
            a function of time in [G],
    """

    def ikey(key:str) -> int:
        return v_ns.index(key)

    # ----------------------------------------------------------
    # Set the initial conditions.
    o = MagnetarRonchi22()
    omega_initial = 2.0 * np.pi / pars["P_in"]
    initial_conditions = np.array([pars["B0"], omega_initial, pars["alpha0"]])
    other_variables = np.zeros(11)
    solution = np.zeros((len(time_grid), len(v_ns)))
    # ---- add initial data and insert it into container
    solution[0, :len(initial_conditions)] = initial_conditions
    _ = o(time_grid[0], initial_conditions, pars)
    for ik in range(len(initial_conditions), len(v_ns)):
        solution[0, ik] = pars[v_ns[ik]]
    # ---- set integrator
    odeinstance = ode(MagnetarRonchi22())
    odeinstance.set_integrator("dop853", rtol=1e-8, nsteps=1000, first_step=time_grid[0])
    odeinstance.set_f_params(pars)
    odeinstance.set_initial_value(initial_conditions,time_grid[0])
    # ---- integrate
    for i in tqdm(range(1,len(time_grid))):
        t = time_grid[i]
        odeinstance.set_f_params(pars)
        # --------------------------------------------
        # check if disk is disrupted and update the pars;
        # solution = np.vstack((solution, np.zeros(len(initial_conditions))))
        i_res = np.copy(odeinstance.integrate(time_grid[i]))
        if (t > time_grid[0]):
            assert np.isfinite(odeinstance.t)
            assert odeinstance.t > time_grid[0]
        solution[i, :len(initial_conditions)] = i_res
        # iB = i_res[0]
        # iOmega = i_res[1]
        # ialpha = i_res[2]
        # compute and store RHS values again for plotting
        _ = o(t, i_res, pars)
        for ik in range(len(initial_conditions), len(v_ns)):
            solution[i, ik] = pars[v_ns[ik]]
        P_t = 2.0 * np.pi / solution[i, ikey("omega")]
        Pdot_t = -(P_t ** 2) / (2.0 * np.pi) * solution[i, ikey("omegadot")]

    # plt.loglog(time_grid, np.abs(solution[:, ikey("omega")]),color='blue')
    # plt.loglog(time_grid, np.abs(solution[:, ikey("n_dip")]),color='blue')
    # plt.loglog(time_grid, np.abs(solution[:, ikey("n_acc")]),color='green')
    # plt.loglog(time_grid, np.abs(solution[:, ikey("n_grav")]),color='red')
    # plt.show()
    return solution


        # # Compute the inner and outer radius of the disk.
        # r_in = np.minimum(r_m, r_lc)
        # r_out = MagnetarRonchi22.disk_outer_radius_t(time_grid)
        #
        # # Check if the field could be buried in the initial stages (before the field decays
        # # on an Ohmic timescale ~10^6 yr).
        # buried = False
        # if np.min(r_m[time_grid / const.YR_TO_S < 10 ^ 6]) < NS_radius:
        #     buried = True
        #     print("buried:", buried)
        #
        # # Check if the disk is disrupted because r_in > r_out and save the time when this happens.
        # disk_disrupted = False
        # t_disrupt = None
        # if any(r_in > r_out):
        #     disk_disrupted = True
        #     t_disrupt = time_grid[r_in > r_out]
        #     t_disrupt = t_disrupt[0]
        #
        #     # Re-run the evolution considering the time when the disk should be disrupted.
        #     evol_output = np.array(
        #         odeint(
        #             Magnetar.derivatives_spin_evolution,
        #             y0=initial_conditions,
        #             t=time_grid,
        #             args=(B0, Mdot_d0, t_disrupt),
        #             tfirst=True,
        #         )
        #     )
        #
        #     print("disk disrupted:", disk_disrupted)
        #
        # # omega_t = np.array(evol_output[:, 1].tolist())
        # # omega_t = np.array(evol_output[:, 1].tolist())
        # # omega_t = np.array(evol_output[:, 1].tolist())
        # B_t = np.array(evol_output[:, 0].tolist())
        # omega_t = np.array(evol_output[:, 1].tolist())
        # alpha_t = np.array(evol_output[:, 2].tolist())
        #
        # P_t = 2 * np.pi / omega_t
        #
        # omegadot_t, n_dip_t, n_acc_t, n_grav_t = \
        #     Magnetar.omegadot_numpy(time_grid, omega_t, B_t, alpha_t, Mdot_d0, t_disrupt)
        # Pdot_t = -(P_t ** 2) / (2.0 * np.pi) * omegadot_t
        #
        # Mdot_din = Magnetar.disk_accretion_rate_in_t_numpy(
        #     B_t, omega_t, Mdot_d0, alpha_t, time_grid, t_disrupt
        # )
        #
        # r_m = Magnetar.radius_magnetospheric(B_t, Mdot_din)
        # r_c = Magnetar.radius_corotation(omega_t)
        # r_lc = Magnetar.radius_lc(omega_t)
        # # n_dip = Magnetar.Nd
        #
        # ldip = Magnetar.Ldip(n_dip_t, omega_t)
        # lacc = Magnetar.Lcc(Mdot_din,n_acc_t,omega_t,r_m)

def main():
    # ----------------- |SET PARAMETERS
    t_max = 1.0e8 # Maximum time reached by the evolution in [yr].
    t_max = t_max * const.YR_TO_S
    time_grid = np.logspace(np.log10(10.0), np.log10(t_max), 2000)
    # --------------------------------------------------------------------------------------------------------------
    pars = dict()
    pars["NS_radius"] = 1.2e6 # Characteristic neutron star radius in [cm].
    pars["NS_mass"] = 1.9 * const.M_SUN # Characteristic neutron star mass in solar masses.
    pars["NS_inertia"] = 2.0 / 5.0 * pars["NS_mass"] * pars["NS_radius"] ** 2 # Canonical neutron star moment of inertia in [g cm^2] assuming a perfect solid sphere.
    pars["NS_critical_beta"] = 0.27 # bar-mode instability criterion (for rryand options only)
    pars["B0"] = 1e15 # NS magnetic field
    pars["ns_ellipticity"] = 0.1 # [my] for GWs only
    pars["alpha0"] = 0.0 # Assume an inclination angle in [rad].
    pars["sigma_conductivity"] = 1e24 # Dominant conductivity based on phonon or impurity scattering, in [1/s]. For details see Cumming et al. (2004) or Gourgouliatos and Cumming (2014).
    pars["characteristic_length"] = 1e5 # Characteristic length scale of the magnetic field in [cm].
    pars["characteristic_electron_dens"] = 1e35 # Characteristic electron density in [g/cm^3].
    # --------------------------------------------------------------------------------------------------------------
    pars["circularization_disk_radius"] = 1.0e9 # Circularization radius of the disk in [cm].
    pars["characteristic_disk_temp"] = 1.0e9 # Initial disk central temperature in [K].
    pars["disk_acc_rate_pl"] = 1.2  # Power-law index for the disk accretion rate decay.
    pars["disk_rout_evol_idx"] = 0.44 # Power-law index for the disk outer radius evolution.
    pars["diskmass"] = 1.e-2 * 1.989e+33
    pars["Mdot0"] = 1.0e29 # 1.0e19,1.0e29; initial accretion rate for the fall-back disk in [g s^-1].
    pars["disk_radius"]=5.e8
    pars["disk_alpha"]=0.1
    pars["disk_askpect_ratio"]=0.3
    pars["disk_opacity"]=1.
    # --------------------------------------------------------------------------------------------------------------
    pars["method_ndip"] = "original" # "original", "gompertz", "rryand"
    pars["method_nacc"] = "original" #
    pars["method_disk_acc"] = "pl" # "exp" "pl"
    pars["method_t_vis"] = "rrayand"# "menou" # "exp" "pl"
    # -------------------------------------------------------------------------------------------------------------
    pars["t_vis"] = 0.
    pars["omegadot"] = 0.
    pars["Bdot"] = 0.
    pars["alphadot"] = 0.
    pars["Ledd"] = 0.
    pars["Mdot_d"] = 0.
    pars["Mdot"] = 0.
    pars["r_m"] = 0.
    pars["r_c"] = 0.
    pars["r_lc"] = 0.
    pars["r_in"] = 0.
    pars["r_out"] = 0.
    pars["n_acc"] = 0.
    pars["n_dip"] = 0.
    pars["n_grav"] = 0.
    pars["ldip"] = 0.
    pars["lacc"] = 0.
    pars["t_disrupt"] = -1.
    # accretion rate0
    method_mdot0 = "disk0" # "given" "disk0"
    if method_mdot0 == "given" :
        pass
    elif method_mdot0 == "disk0" : # [requires disk mass]
        t_v = MagnetarRonchi22.t_vis( pars )
        pars["Mdot0"] = pars["diskmass"] / t_v
        print(f"t_v={t_v} mdot0="+"{:.2e}".format(pars["Mdot0"]))
    else: raise KeyError()
    pars["Mdot_d0"] = pars["Mdot0"]

    # Initial spin period of the neutron star in [s].
    method_P0 = "given" # "given" "crit"
    if method_P0 == "given":
        pars["P_in"] = 0.001#0.01
    elif method_P0 == "crit":
        # NS_critical_beta = 0.27 # bar-mode instability criterion
        Omega0 = np.sqrt(2.*pars["NS_critical_beta"]
                         * MagnetarRonchi22.E_bind(NS_mass=pars["NS_mass"],NS_radius=pars["NS_radius"])
                         /pars["NS_inertia"])
        pars["P_in"] = 2*np.pi/Omega0
    else:
        raise KeyError()

    pars["k"] = -1  #0.9 # # Following Gibson+ (capping) r_m <= k * r_lc

    # ----------------------------------------------------------
    # Dimensionless coefficients k_0, k_1, k_2 for a force-free magnetosphere
    # taken from Spitkovsky (2006) and Philippov et al. (2014).
    # For comparison, in vacuum k_0 = 0 and k_1 = k_2 = 2/3.
    pars["k_coefficients"] = [1.0, 1.0, 1.0]
    # -------------------------------------------------------------------------------------------------------------
    v_ns = ["B","omega","alpha","t_vis","omegadot","Bdot","alphadot","Ledd","Mdot","Mdot_d",
            "r_m","r_c","r_lc","r_in","r_out","n_acc","n_dip","n_grav","ldip","lacc"]

    solution = run_ronchi_magnetar(time_grid, pars, v_ns)
    # prepare evolution file for the C++ afterglow code
    dfile = h5py.File("magnetar_evol.h5", "w")
    dfile.create_dataset(name="time", data=time_grid)
    for i, key in enumerate(v_ns):
        dfile.create_dataset(name=key, data=solution[:, v_ns.index(key)])
    dfile.close()

    # plot the resulted evolution
    plot_ronchi_magnetar(time=time_grid, solution=solution, v_ns=v_ns, pars=pars)

if __name__ == '__main__':
    main()

