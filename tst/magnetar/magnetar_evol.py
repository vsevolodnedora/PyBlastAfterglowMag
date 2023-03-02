#!/usr/bin/env python3
"""propeller_parameter_search_Bdecay.ipynb: a notebook to recreate results of Ronchi et al. (2022)"""
#
# __author__ = "Michele Ronchi"
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




# Characteristic neutron star radius in [cm].
NS_radius = 1.1e6

# Characteristic neutron star mass in solar masses.
NS_mass = 1.4 * const.M_SUN

NS_critical_beta = 0.27 # bar-mode instability criterion (for rryand options only)

# NS magnetic field
B0 = 5.0e15

# Dimensionless coefficients k_0, k_1, k_2 for a force-free magnetosphere
# taken from Spitkovsky (2006) and Philippov et al. (2014).
# For comparison, in vacuum k_0 = 0 and k_1 = k_2 = 2/3.
k_coefficients = [1.0, 1.0, 1.0]

# Canonical neutron star moment of inertia in [g cm^2] assuming a perfect solid sphere.
NS_inertia = 2.0 / 5.0 * NS_mass * NS_radius ** 2

# [my] for GWs only
ns_ellipticity = 0.1

# Auxiliary quantity beta as defined in eq. (72) of Pons & Vigano (2019).
beta_ax = 1.0 / 4.0 * NS_radius ** 6 / (NS_inertia * const.C ** 3)

# Assume an inclination angle in [rad].
alpha0 = 0.2 # TODO why does it circularises wo quickly???

# # Incorporate the inclination angle dependence into a constant.
# beta_1 = beta * (k_coefficients[0] + k_coefficients[1] * np.sin(chi) ** 2)

# Dominant conductivity based on phonon or impurity scattering, in [1/s].
# For details see Cumming et al. (2004) or Gourgouliatos and Cumming (2014).
sigma = 1e24

# Characteristic length scale of the magnetic field in [cm].
L = 1e5

# Characteristic electron density in [g/cm^3].
n_e = 1e35

# Eddington luminosity assuming Thompson scattering in [erg / s].
L_Edd = 4.0 * np.pi * const.G * NS_mass * const.M_P * const.C / const.SIGMA_T

# Circularization radius of the disk in [cm].
r_d = 1.0e8

# Initial disk central temperature in [K].
T_c = 1.0e6


# typical initial viscousity timescale of the disk in [s] see Menou et al. 2001.
method_t_vis = "menou" # "menou" "gompertz" "rrayand"
if method_t_vis == "menou":
    t_v = 2080.0 * (T_c / 1.0e6) ** (-1.0) * (r_d / 10 ** 8) ** (1.0 / 2.0) # s
elif method_t_vis == "gompertz":
    disk_radius = 5.e8
    disk_alpha = 0.1  # disk viscosity parameter
    disk_cs = 1e7
    t_v = disk_radius / (3. * disk_alpha * disk_cs * disk_radius) # s
elif method_t_vis == "rrayand": # TODO ADD TEMPERATURE DEPENDENCE HERE!
    disk_radius = 5.e8
    disk_alpha = 0.1  # disk viscosity parameter
    disk_askpect_ratio = 0.3 # Aspect ratio H/R
    disk_H = disk_radius * disk_askpect_ratio
    disk_omega_kep = OmegaKep = (const.grav * NS_mass / NS_radius**3)**0.5
    disc_cs = disk_H * disk_omega_kep * (NS_radius/disk_radius)**1.5
    t_v = disk_radius**2 / (3. * disk_alpha * disc_cs * disk_H) # s
else:
    raise KeyError()


# accretion rate0
method_mdot0 = "given" # "given" "disk0"
if method_mdot0 == "given" :
    Mdot0 = 1.0e24
elif method_mdot0 == "disk0" : # [requires disk mass]
    diskmass = 1.e-2
    Mdot0 = diskmass / t_v
else: raise KeyError()

t_max = 1.0e8 # Maximum time reached by the evolution in [yr].
t_max = t_max * const.YR_TO_S
time_grid = np.logspace(np.log10(10.0), np.log10(t_max), 2000)


# Power-law index for the disk accretion rate decay.
disk_acc_rate_pl = 1.2

# Power-law index for the disk outer radius evolution.
gamma = 0.44


# Initial spin period of the neutron star in [s].
method_P0 = "given" # "given" "crit"
def E_bind():
    """
    Binding energy of the NS
    Prescription from Lattimer and Prakash (2001)

    """
    num = const.grav*NS_mass
    den = NS_radius*const.C**2-0.5*const.grav*NS_mass
    out = 0.6*NS_mass*const.C**2*num/den
    return out
if method_P0 == "given":
    P_in = 2e-3#0.01
elif method_P0 == "crit":
    # NS_critical_beta = 0.27 # bar-mode instability criterion
    Omega0 = np.sqrt(2.*NS_critical_beta*E_bind()/NS_inertia)
    P_in = 2*np.pi*Omega0
else:
    raise KeyError()

method_ndip = "original" # "original", "gompertz", "rryand"

method_nacc = "original"

def E_rot(omega):
    """
    Rotational energy of the NS
    """
    out=0.5*NS_inertia*omega**2
    return out

# Range of initial accretion rate for the fall-back disk in [g s^-1].
Mdot_d0_min = 1.0e19
Mdot_d0_max = 1.0e29

# Range of initial magnetic fields in [G].
B0_min = 1.0e12
B0_max = 1.0e15

# [my]
method_disk_acc = "pl" # "exp" "pl"


method_omega_crit = "fromTOV" # "given" "fromTOV"
if method_omega_crit == "given":
    OmegaCrit = 1e3
    Pcrit = 2*np.pi*OmegaCrit
elif method_omega_crit == "fromTOV":
    """
    Sun, Zhang & Gao (2017)
    eq. (25)
    
    NS collapse for Omega < Omega_c (P>Pc)
    
    Rem: assume constant NS mass
    
    """
    EOS_Mtov = 2.18 # Msun
    EOS_alpha = 4.678e-10
    EOS_beta = -2.738
    num = NS_mass - EOS_Mtov
    if (num>0):
        den = EOS_alpha * EOS_Mtov
        Pcrit = (num/den)**(1./EOS_beta)
        OmegaCrit = 2*np.pi/Pcrit
    else:
        ## then NS always stable
        Pcrit = np.inf
        OmegaCrit = np.inf
else:
    raise KeyError()

''' ---------------------------------------------------------------------------------------------------------------- '''
# TODO Add (i) Bext/Bint for ellipticity; (ii) NS_mass/Radius evolution; (iii) Collapse criterion (iv) does Ldip = f(alpha)?
class Magnetar:

    def __init__(self):
        pass

    @staticmethod
    def field_derivative(B: float, B_initial: float) -> float:
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

        tau_Ohm = Magnetar.timescale_Ohm(L, sigma)
        tau_Hall = Magnetar.timescale_Hall(B_initial, L, n_e)

        B_deriv = -B / tau_Ohm - B ** 2 / (tau_Hall * B_initial)

        return B_deriv
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
    def radius_magnetospheric(B: float, Mdot: float) -> float:
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
    def radius_corotation(omega: float) -> float:
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

    @staticmethod
    def accretion_rate_Edd_rmag(B: float) -> float:
        """
        Eddington limit of the accretion rate at the magnetospheric radius.

        Args:
            B (float): value of the dipolar component of the magnetic field at the magnetic pole measured in [G].

        Returns:
            (float): Eddington limit of the accretion rate at the magnetospheric radius in [g s^-1].
        """

        # Magnetic moment.
        mu = B * NS_radius ** 3 / 2.0

        Mdot_Edd = (L_Edd / (const.G * NS_mass)) ** (7.0 / 9.0) * (
                mu ** 4.0 / (2.0 * const.G * NS_mass)
        ) ** (1.0 / 9.0)

        return Mdot_Edd

    @staticmethod
    def accretion_rate_Edd_rlc(omega: float) -> float:
        """
        Eddington limit of the accretion rate at the light cylinder radius.

        Args:
            omega (float): value of the spin frequency measured in [1/s].

        Returns:
            (float): Eddington limit of the accretion rate at the light cylinder radius in [g s^-1].
        """

        # Light cylinder radius.
        r_lc = Magnetar.radius_lc(omega)

        Mdot_Edd = 2 * L_Edd * r_lc / (const.G * NS_mass)

        return Mdot_Edd

    @staticmethod
    def disk_accretion_rate_t(Mdot_d0: float, alpha: float, t:float) -> float:
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
            Mdot_d = Mdot_d0 * (1 + t / t_v) ** (-disk_acc_rate_pl)

        elif (method_disk_acc == "exp"):
            # [my] Accretion rate on the NS
            # [my] Eq. (13) from Zhang and Meszaros 2001
            Mdot_d = Mdot_d0 * np.exp(-t / t_v)
            mdot_floor=1.e-10
            if hasattr(t, '__len__'):
                Mdot_d[Mdot_d<mdot_floor] = mdot_floor
            else:
                if Mdot_d < mdot_floor: Mdot_d = mdot_floor
        else:
            raise KeyError()

        return Mdot_d

    @staticmethod
    def disk_outer_radius_t(t:float) -> float:
        """
        Accretion rate inside the disk as a function of time.

        Args:
            t (float): time in [s].

        Returns:
            (float): disk accretion rate in [g s^-1] as a function of time.
        """

        r_out = r_d * (1 + t / t_v) ** gamma

        return r_out

    @staticmethod
    def disk_accretion_rate_in_t(
            B: float,
            omega: float,
            Mdot_d0: float,
            alpha: float,
            t: float,
            t_disrupt: float = None,
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
        Mdot_Edd_mag = Magnetar.accretion_rate_Edd_rmag(B)
        Mdot_Edd_lc = Magnetar.accretion_rate_Edd_rlc(omega)

        # The disk inner radius is the smaller one between the magnetospheric radius and the light cylinder radius.
        Mdot_Edd = min(Mdot_Edd_mag, Mdot_Edd_lc)

        Mdot_d = Magnetar.disk_accretion_rate_t(Mdot_d0, alpha, t)
        Mdot_din = 0.0

        # Limit the accretion rate at the inner disk edge to the Eddington limit.
        if Mdot_d < Mdot_Edd:
            Mdot_din = Mdot_d
        elif Mdot_d >= Mdot_Edd:
            Mdot_din = Mdot_Edd

        # When the disk is disrupted because r_in > r_out set the accretion rate to 0 afterwards.
        if t_disrupt is not None:
            if t > t_disrupt:
                Mdot_din = 0.0

        return Mdot_din
    @staticmethod
    def disk_accretion_rate_in_t_numpy(
            B_arr: np.ndarray,
            omega_arr: np.ndarray,
            Mdot_d0: float,
            alpha_arr: np.ndarray,
            t_arr: np.array,
            t_disrupt: float = None,
    ) -> float:
        """
        Accretion rate in the inner boundary of the disk as a function of time optimized for numpy array.

        Args:
            B (float): value of the dipolar component of the magnetic field at the magnetic pole measured in [G].
            omega (np.ndarray): value of the spin frequency measured in [1/s].
            Mdot_d0 (float): initial maximum fallback rate in [g s^-1].
            alpha (float): index of the power law for the time evolution of the inflow rate.
            t (np.ndarray): time in [s].
            t_disrupt (float): time in [s] when the disk is disrupted because r_in > r_out.

        Returns:
            (np.ndarray): fall back accretion rate in [g s^-1] as a function of time.
        """
        Mdot_din_arr = np.zeros(len(t_arr))
        for i, ti in enumerate(t_arr):
            t = t_arr[i]
            omega = omega_arr[i]
            B = B_arr[i]
            alpha = alpha_arr[i]

            # Compute Eddington limits on the accretion rate [g s^-1].
            Mdot_Edd_mag = Magnetar.accretion_rate_Edd_rmag(B)
            Mdot_Edd_lc = Magnetar.accretion_rate_Edd_rlc(omega)

            # The disk inner radius is the smaller one between the magnetospheric radius and the light cylinder radius.
            Mdot_Edd = np.minimum(Mdot_Edd_mag, Mdot_Edd_lc)

            Mdot_d = Magnetar.disk_accretion_rate_t(Mdot_d0, alpha, t)
            # Mdot_din = np.zeros(len(t))

            # Limit the accretion rate at the inner disk edge to the Eddington limit.
            # Mdot_din[Mdot_d >= Mdot_Edd] = Mdot_Edd[Mdot_d >= Mdot_Edd]
            # Mdot_din[Mdot_d < Mdot_Edd] = Mdot_d[Mdot_d < Mdot_Edd]

            if Mdot_d >= Mdot_Edd: Mdot_din = Mdot_Edd
            elif Mdot_d < Mdot_Edd: Mdot_din = Mdot_d
            else:raise ValueError()

            # When the disk is disrupted because r_in > r_out set the accretion rate to 0 afterwards.
            if t_disrupt is not None:
                if t > t_disrupt: Mdot_din = 0.
                # Mdot_din[t > t_disrupt] = 0.0
            Mdot_din_arr[i] = Mdot_din
        return Mdot_din_arr

    @staticmethod
    def omega_keplerian(r: float) -> float:
        """
        Keplerian orbital frequency as a function of the radius.

        Args:
            r (float): distance from the star in [cm].

        Returns:
            (float): keplerian orbital frequency in [1/s].
        """

        omega_k = (const.G * NS_mass / (r ** 3)) ** (1.0 / 2.0)

        return omega_k
    @staticmethod
    def alpha_derivatives(B, omega, alpha):
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
    @staticmethod
    def Ngrav(omega, alpha):
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
    def Ndip(B, omega, r_m, r_lc, alpha):
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
    def Nacc(omega, Mdot_din, r_m, r_lc, r_corot):

        omega_k_rns = Magnetar.omega_keplerian(NS_radius)
        omega_k_rm = Magnetar.omega_keplerian(r_m)
        omega_k_rlc = Magnetar.omega_keplerian(r_lc)

        if method_nacc == "original":
            if r_m <= NS_radius:
                n_acc = Mdot_din * NS_radius ** 2  * (omega_k_rns - omega)
            elif (r_m > NS_radius) & (r_m < r_lc):
                n_acc = Mdot_din * r_m ** 2  * (omega_k_rm - omega)
            elif r_m >= r_lc:
                n_acc = Mdot_din * r_lc ** 2 * (omega_k_rlc - omega)
            else:
                raise ValueError()

        elif method_nacc == "rryand":
            # Accretion torque, taking into account the propeller model
            #         Eq (6-7) of Gompertz et al. 2014
            fastness = (r_m / r_corot)**1.5
            ## Eq. (6)
            n_acc = (1. - fastness) * (const.grav * NS_mass * r_m)**0.5 * Mdot_din
            if r_m<=NS_radius:
                OmegaKep = (const.grav * NS_mass / NS_radius**3)**0.5
                n_acc = ((1. - omega/OmegaKep) * (const.grav*NS_mass*r_m)**0.5 * Mdot_din)
            ###############################################
            ## Check for inhibition by bar-mode instability
            ## with beta = T/|W| parameter (Gompertz 2014)
            ###############################################
            beta_bar_instab = E_rot(omega)/abs(E_bind())
            if beta_bar_instab > NS_critical_beta:
                n_acc = 0.
        else:
            raise KeyError()
        return n_acc
    @staticmethod
    def derivatives_spin_evolution(
            t: float,
            initial_cond: np.ndarray,
            B0: float,
            Mdot_d0: float,
            t_disrupt: float = None,
    ) -> float:
        """
        Derivative of the magnetic field and spin frequency for the three phases: direct accretion, propeller and
        ejector. See Metzger et al. (2018) for details.

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
        B = initial_cond[0]
        omega = initial_cond[1]
        alpha = initial_cond[2]

        Bdot = Magnetar.field_derivative(B, B0)
        alphadot = Magnetar.alpha_derivatives(B, omega, alpha)

        Mdot_d = Magnetar.disk_accretion_rate_t(Mdot_d0, alpha, t)
        Mdot_din = Magnetar.disk_accretion_rate_in_t(B, omega, Mdot_d0, alpha, t, t_disrupt)

        r_m = Magnetar.radius_magnetospheric(B, Mdot_din)
        r_c = Magnetar.radius_corotation(omega)
        r_lc = Magnetar.radius_lc(omega)

        if not (np.isfinite(r_m)&np.isfinite(r_c)&np.isfinite(r_lc)):
            raise ValueError()

        # omega_k_rns = omega_keplerian(NS_radius)
        # omega_k_rm = omega_keplerian(r_m)
        # omega_k_rlc = omega_keplerian(r_lc)

        # omegadot = 0.0

        n_acc = Magnetar.Nacc(omega,Mdot_din,r_m,r_lc,r_c)
        n_dip = Magnetar.Ndip(B,omega,r_m,r_lc,alpha)
        n_grav = Magnetar.Ngrav(omega,alpha)

        if not (np.isfinite(n_acc) and np.isfinite(n_dip) and np.isfinite(n_grav) and np.isfinite(omega) and omega>0):
            raise ValueError()

        # n_grav = - 32./5. \
        #          * 6.67259e-8 \
        #          * NS_inertia**2 \
        #          * ns_ellipticity**2 \
        #          * omega**5 \
        #          / 299792458e2**5
        # if r_m <= NS_radius:
        #     # Direct accretion phase with likely buried magnetic field.
        #     n_acc = Mdot_din * NS_radius ** 2  * (omega_k_rns - omega)
        #     n_dip = - (r_lc / NS_radius) ** 2 * beta_1 * B ** 2 * omega ** 3 * NS_inertia
        #
        #     # omegadot = (
        #     #         Mdot_din * NS_radius ** 2 / NS_inertia * (omega_k_rns - omega)
        #     #         - (r_lc / NS_radius) ** 2 * beta_1 * B ** 2 * omega ** 3
        #     # )
        #
        # elif (r_m > NS_radius) & (r_m < r_lc):
        #     # Propeller (r_m > r_c) or direct accretion phase (r_m < r_c).
        #     n_acc = Mdot_din * r_m ** 2  * (omega_k_rm - omega)
        #     n_dip = - (r_lc / r_m) ** 2 * beta_1 * B ** 2 * omega ** 3 * NS_inertia
        #
        #     # omegadot = (
        #     #         Mdot_din * r_m ** 2  * (omega_k_rm - omega)
        #     #         - (r_lc / r_m) ** 2 * beta_1 * B ** 2 * omega ** 3 * NS_inertia
        #     # )
        #
        # elif r_m >= r_lc:
        #     # Ejector phase.
        #     n_acc = Mdot_din * r_lc ** 2 * (omega_k_rlc - omega)
        #     n_dip = - beta_1 * B ** 2 * omega ** 3 * NS_inertia
        #
        #     # Ejector phase.
        #     # omegadot = (
        #     #         Mdot_din * r_lc ** 2 / NS_inertia * (omega_k_rlc - omega)
        #     #         - beta_1 * B ** 2 * omega ** 3
        #     # )

        omegadot = (n_dip + n_acc + n_grav) / NS_inertia

        derivatives = np.array([Bdot, omegadot, alphadot])

        return derivatives


    # def omegadot_numpy_old(
    #         t: np.ndarray,
    #         omega: np.ndarray,
    #         B: np.ndarray,
    #         Mdot_d0: float,
    #         alpha: float,
    #         t_disrupt: float = None,
    # ) -> float:
    #     """
    #     Derivative of the spin frequency for the three phases: direct accretion, propeller and ejector
    #     optimized for numpy arrays. See Metzger et al. (2018) for details.
    #
    #     Args:
    #         t (np.ndarray): time variable in [s],
    #         omega (np.ndarray): value of the spin frequency measured in [1/s].
    #         B (np.ndarray): value of the dipolar component of the magnetic field at the
    #         magnetic pole measured in [G] as a function of time.
    #         Mdot (float): accretion rate in [g s^-1].
    #         alpha (float): index of the power law for the time evolution of the accretion rate.
    #         t_disrupt (float): time in [s] when the disk is disrupted because r_in > r_out.
    #
    #     Returns:
    #         (np.ndarray): Derivative of the spin frequency for the three phases:
    #         direct accretion, propeller and ejector [s^-2].
    #     """
    #
    #     Mdot_d = disk_accretion_rate_t(Mdot_d0, alpha, t)
    #     Mdot_din = disk_accretion_rate_in_t_numpy(B, omega, Mdot_d0, alpha, t, t_disrupt)
    #
    #     r_m = radius_magnetospheric(B, Mdot_din)
    #     r_c = radius_corotation(omega)
    #     r_lc = radius_lc(omega)
    #
    #     omegadot = np.zeros(len(t))
    #
    #     omega_k_rns = omega_keplerian(NS_radius)
    #     omega_k_rm = omega_keplerian(r_m)
    #     omega_k_rlc = omega_keplerian(r_lc)
    #
    #     cond_acc = r_m <= NS_radius
    #     cond_prop_acc = (r_m > NS_radius) & (r_m < r_lc)
    #     cond_ej = r_m >= r_lc
    #
    #     n_acc = np.zeros(len(t))
    #     n_dip = np.zeros(len(t))
    #     n_grav = - 32./5. \
    #              * 6.67259e-8 \
    #              * NS_inertia**2 \
    #              * ns_ellipticity**2 \
    #              * omega**5 \
    #              / 299792458e2**5
    #
    #     n_acc[cond_acc] = Mdot_din[cond_acc] * NS_radius ** 2 * (omega_k_rns - omega[cond_acc])
    #     n_dip[cond_acc] = - (r_lc[cond_acc] / NS_radius) ** 2 * beta_1 * B[cond_acc] ** 2 * omega[cond_acc] ** 3 * NS_inertia
    #
    #     # omegadot[cond_acc] = (
    #     #         Mdot_din[cond_acc]
    #     #         * NS_radius ** 2
    #     #         / NS_inertia
    #     #         * (omega_k_rns - omega[cond_acc])
    #     #         - (r_lc[cond_acc] / NS_radius) ** 2
    #     #         * beta_1
    #     #         * B[cond_acc] ** 2
    #     #         * omega[cond_acc] ** 3
    #     # )
    #
    #     n_acc[cond_prop_acc] = Mdot_din[cond_prop_acc] * r_m[cond_prop_acc] ** 2 * (omega_k_rm[cond_prop_acc] - omega[cond_prop_acc])
    #     n_dip[cond_prop_acc] = - (r_lc[cond_prop_acc] / r_m[cond_prop_acc]) ** 2 * beta_1 * B[cond_prop_acc] ** 2 * omega[cond_prop_acc] ** 3 * NS_inertia
    #
    #     # omegadot[cond_prop_acc] = (
    #     #         Mdot_din[cond_prop_acc]
    #     #         * r_m[cond_prop_acc] ** 2
    #     #         / NS_inertia
    #     #         * (omega_k_rm[cond_prop_acc] - omega[cond_prop_acc])
    #     #         - (r_lc[cond_prop_acc] / r_m[cond_prop_acc]) ** 2
    #     #         * beta_1
    #     #         * B[cond_prop_acc] ** 2
    #     #         * omega[cond_prop_acc] ** 3
    #     # )
    #
    #     n_acc[cond_ej] = Mdot_din[cond_ej] * r_lc[cond_ej] ** 2 * (omega_k_rlc[cond_ej] - omega[cond_ej])
    #     n_dip[cond_ej] = - beta_1 * B[cond_ej] ** 2 * omega[cond_ej] ** 3 * NS_inertia
    #
    #     # omegadot[cond_ej] = (
    #     #         Mdot_din[cond_ej]
    #     #         * r_lc[cond_ej] ** 2
    #     #         / NS_inertia
    #     #         * (omega_k_rlc[cond_ej] - omega[cond_ej])
    #     #         - beta_1 * B[cond_ej] ** 2 * omega[cond_ej] ** 3
    #     # )
    #
    #     omegadot = (n_dip + n_acc + n_grav) / NS_inertia
    #
    #     return omegadot
    @staticmethod
    def omegadot_numpy(
            t: np.ndarray,
            omega_arr: np.ndarray,
            B_arr: np.ndarray,
            alpha_arr :np.ndarray,
            Mdot_d0: float,
            t_disrupt: float = None,
    ):
        """
        Derivative of the spin frequency for the three phases: direct accretion, propeller and ejector
        optimized for numpy arrays. See Metzger et al. (2018) for details.

        Args:
            t (np.ndarray): time variable in [s],
            omega (np.ndarray): value of the spin frequency measured in [1/s].
            B (np.ndarray): value of the dipolar component of the magnetic field at the
            magnetic pole measured in [G] as a function of time.
            Mdot (float): accretion rate in [g s^-1].
            alpha (float): index of the power law for the time evolution of the accretion rate.
            t_disrupt (float): time in [s] when the disk is disrupted because r_in > r_out.

        Returns:
            (np.ndarray): Derivative of the spin frequency for the three phases:
            direct accretion, propeller and ejector [s^-2].
        """

        omegadot = np.zeros(len(t))
        n_acc = np.zeros(len(t))
        n_dip = np.zeros(len(t))
        n_grav = np.zeros(len(t))
        for i, ti in enumerate(t):

            omega = omega_arr[i]
            B = B_arr[i]
            alpha = alpha_arr[i]

            Mdot_d = Magnetar.disk_accretion_rate_t(Mdot_d0, alpha, ti)
            Mdot_din = Magnetar.disk_accretion_rate_in_t(B, omega, Mdot_d0, alpha, ti, t_disrupt)

            r_m = Magnetar.radius_magnetospheric(B, Mdot_din)
            r_c = Magnetar.radius_corotation(omega)
            r_lc = Magnetar.radius_lc(omega)

            # omegadot = np.zeros(len(t))

            # omega_k_rns = omega_keplerian(NS_radius)
            # omega_k_rm = omega_keplerian(r_m)
            # omega_k_rlc = omega_keplerian(r_lc)

            # cond_acc = r_m <= NS_radius
            # cond_prop_acc = (r_m > NS_radius) & (r_m < r_lc)
            # cond_ej = r_m >= r_lc
            #
            # n_grav = Ngrav(omega)
            # if cond_acc:
            #     n_acc = Mdot_din * NS_radius ** 2 * (omega_k_rns - omega)
            #     n_dip = (r_lc / NS_radius) ** 2 * beta_1 * B ** 2 * omega ** 3 * NS_inertia
            # elif cond_prop_acc:
            #     n_acc = Mdot_din * r_m ** 2 * (omega_k_rm - omega)
            #     n_dip = - (r_lc / r_m) ** 2 * beta_1 * B ** 2 * omega ** 3 * NS_inertia
            # elif cond_ej:
            #     n_acc = Mdot_din * r_lc ** 2 * (omega_k_rlc - omega)
            #     n_dip = - beta_1 * B ** 2 * omega ** 3 * NS_inertia
            # else:
            #     raise ValueError()

            n_acc[i] = Magnetar.Nacc(omega,Mdot_din,r_m,r_lc,r_c)
            n_dip[i] = Magnetar.Ndip(B,omega,r_m,r_lc,alpha)
            n_grav[i] = Magnetar.Ngrav(omega,alpha)
            omegadot[i] = (n_dip[i] + n_acc[i] + n_grav[i]) / NS_inertia

        # n_dip = np.zeros(len(t))
        # n_grav = - 32./5. \
        #          * 6.67259e-8 \
        #          * NS_inertia**2 \
        #          * ns_ellipticity**2 \
        #          * omega**5 \
        #          / 299792458e2**5
        #
        # n_acc[cond_acc] = Mdot_din[cond_acc] * NS_radius ** 2 * (omega_k_rns - omega[cond_acc])
        # n_dip[cond_acc] = - (r_lc[cond_acc] / NS_radius) ** 2 * beta_1 * B[cond_acc] ** 2 * omega[cond_acc] ** 3 * NS_inertia
        #
        # # omegadot[cond_acc] = (
        # #         Mdot_din[cond_acc]
        # #         * NS_radius ** 2
        # #         / NS_inertia
        # #         * (omega_k_rns - omega[cond_acc])
        # #         - (r_lc[cond_acc] / NS_radius) ** 2
        # #         * beta_1
        # #         * B[cond_acc] ** 2
        # #         * omega[cond_acc] ** 3
        # # )
        #
        # n_acc[cond_prop_acc] = Mdot_din[cond_prop_acc] * r_m[cond_prop_acc] ** 2 * (omega_k_rm[cond_prop_acc] - omega[cond_prop_acc])
        # n_dip[cond_prop_acc] = - (r_lc[cond_prop_acc] / r_m[cond_prop_acc]) ** 2 * beta_1 * B[cond_prop_acc] ** 2 * omega[cond_prop_acc] ** 3 * NS_inertia
        #
        # # omegadot[cond_prop_acc] = (
        # #         Mdot_din[cond_prop_acc]
        # #         * r_m[cond_prop_acc] ** 2
        # #         / NS_inertia
        # #         * (omega_k_rm[cond_prop_acc] - omega[cond_prop_acc])
        # #         - (r_lc[cond_prop_acc] / r_m[cond_prop_acc]) ** 2
        # #         * beta_1
        # #         * B[cond_prop_acc] ** 2
        # #         * omega[cond_prop_acc] ** 3
        # # )
        #
        # n_acc[cond_ej] = Mdot_din[cond_ej] * r_lc[cond_ej] ** 2 * (omega_k_rlc[cond_ej] - omega[cond_ej])
        # n_dip[cond_ej] = - beta_1 * B[cond_ej] ** 2 * omega[cond_ej] ** 3 * NS_inertia
        #
        # # omegadot[cond_ej] = (
        # #         Mdot_din[cond_ej]
        # #         * r_lc[cond_ej] ** 2
        # #         / NS_inertia
        # #         * (omega_k_rlc[cond_ej] - omega[cond_ej])
        # #         - beta_1 * B[cond_ej] ** 2 * omega[cond_ej] ** 3
        # # )
        #
        # omegadot = (n_dip + n_acc + n_grav) / NS_inertia

        return (omegadot, n_dip, n_acc, n_grav)

    @staticmethod
    def Ldip(n_dip, omega):
        """
        Dipole spindown luminosity, for a general
        time evolution of the NS angular velocity
        """
        ldip = - n_dip * omega
        return ldip
    @staticmethod
    def Lcc(mdot, n_acc, omega, r_mag):
        """
        Propeller luminosity, taking into account
        positive and negative torques due to the
        interaction with the accretion disk
        From Gompertz et al. (2014)
        """

        ### intermediate variables
        # Mdot = self.Accretion_rate(T)
        # Nacc = self.Torque_accretion(T,Omega)
        # rmag = self.Magnetospheric_radius(T,Omega)

        ### output
        lprop = - n_acc*omega - const.grav*NS_mass*mdot/r_mag
        if hasattr(omega,'__len__'):
            lprop[lprop<0.] = 0.
        else:
            if lprop < 0: lprop = 0.
        return lprop

    @staticmethod
    def rotational_evolution(
            B0: float, P_initial: float, Mdot_d0: float, alpha0: float
    ) -> Tuple[
        np.ndarray,
        np.ndarray,
        np.ndarray,
        np.ndarray,
        np.ndarray,
        np.ndarray,
        np.ndarray,
        np.ndarray,
        np.ndarray,
        bool,
        bool,
        float,
    ]:
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

        # Define the time grid in s.
        # t_max = t_max * const.YR_TO_S
        # time_grid = np.logspace(np.log10(10.0), np.log10(t_max), 2000)

        # Set the initial conditions.
        omega_initial = 2.0 * np.pi / P_initial

        initial_conditions = np.array([B0, omega_initial, alpha0])

        # Evolve in time.
        evol_output = np.array(
            odeint(
                Magnetar.derivatives_spin_evolution,
                y0=initial_conditions,
                t=time_grid,
                args=(B0, Mdot_d0,),
                tfirst=True,
            )
        )

        # Extract the evolution output for B and omega.
        B_t = np.array(evol_output[:, 0].tolist())
        omega_t = np.array(evol_output[:, 1].tolist())
        alpha_t = np.array(evol_output[:, 2].tolist())

        Mdot_d = Magnetar.disk_accretion_rate_t(Mdot_d0, alpha_t, time_grid)
        Mdot_din = Magnetar.disk_accretion_rate_in_t_numpy(B_t, omega_t, Mdot_d0, alpha_t, time_grid)

        r_m = Magnetar.radius_magnetospheric(B_t, Mdot_din)
        r_lc = Magnetar.radius_lc(omega_t)

        # Compute the inner and outer radius of the disk.
        r_in = np.minimum(r_m, r_lc)
        r_out = Magnetar.disk_outer_radius_t(time_grid)

        # Check if the field could be buried in the initial stages (before the field decays
        # on an Ohmic timescale ~10^6 yr).
        buried = False
        if np.min(r_m[time_grid / const.YR_TO_S < 10 ^ 6]) < NS_radius:
            buried = True
            print("buried:", buried)

        # Check if the disk is disrupted because r_in > r_out and save the time when this happens.
        disk_disrupted = False
        t_disrupt = None
        if any(r_in > r_out):
            disk_disrupted = True
            t_disrupt = time_grid[r_in > r_out]
            t_disrupt = t_disrupt[0]

            # Re-run the evolution considering the time when the disk should be disrupted.
            evol_output = np.array(
                odeint(
                    Magnetar.derivatives_spin_evolution,
                    y0=initial_conditions,
                    t=time_grid,
                    args=(B0, Mdot_d0, t_disrupt),
                    tfirst=True,
                )
            )

            print("disk disrupted:", disk_disrupted)

        # omega_t = np.array(evol_output[:, 1].tolist())
        # omega_t = np.array(evol_output[:, 1].tolist())
        # omega_t = np.array(evol_output[:, 1].tolist())
        B_t = np.array(evol_output[:, 0].tolist())
        omega_t = np.array(evol_output[:, 1].tolist())
        alpha_t = np.array(evol_output[:, 2].tolist())

        P_t = 2 * np.pi / omega_t

        omegadot_t, n_dip_t, n_acc_t, n_grav_t = \
            Magnetar.omegadot_numpy(time_grid, omega_t, B_t, alpha_t, Mdot_d0, t_disrupt)
        Pdot_t = -(P_t ** 2) / (2.0 * np.pi) * omegadot_t

        Mdot_din = Magnetar.disk_accretion_rate_in_t_numpy(
            B_t, omega_t, Mdot_d0, alpha_t, time_grid, t_disrupt
        )

        r_m = Magnetar.radius_magnetospheric(B_t, Mdot_din)
        r_c = Magnetar.radius_corotation(omega_t)
        r_lc = Magnetar.radius_lc(omega_t)
        # n_dip = Magnetar.Nd

        ldip = Magnetar.Ldip(n_dip_t, omega_t)
        lacc = Magnetar.Lcc(Mdot_din,n_acc_t,omega_t,r_m)

        return (
            time_grid / const.YR_TO_S,
            P_t,
            Pdot_t,
            B_t,
            alpha_t,
            r_m,
            r_c,
            r_lc,
            Mdot_d,
            Mdot_din,
            ldip,
            lacc,
            buried,
            disk_disrupted,
            t_disrupt,
        )

    @staticmethod
    def run():
        time, P, Pdot, B, alpha, r_m, r_c, r_lc, Mdot_d, Mdot_din, ldip, lacc, _, _, _ \
            = Magnetar.rotational_evolution( B0, P_in, Mdot0, alpha0 )
        return (time, ldip, lacc)

def run_magnetar():
    '''
    ## Example of spin-period evolution
    Example of neutron star spin-down in presence of a fallback disk on a timescale of $10^7 \, {\rm yr}$.
    The top panel shows the total disk accretion $\dot{M}_{\rm d}$ and the
    Eddington-limited accretion rate at the inner radius $\dot{M}_{\rm d, in}$.
    The middle panel illustrates the evolution of the three critical radii $r_{\rm cor}$,
    $r_{\rm m}$, $r_{\rm lc}$. The evolution of the disk's inner and outer radii is highlighted
    by the gray lines. The bottom panel shows the resulting time evolution of the spin period.
    We assume an initial spin period $P_0 = 10 \, {\rm ms}$, initial magnetic field
    $B_0 = 4 \times 10^{14} \, {\rm G}$ and an initial disk accretion rate $\dot{M}_{\rm d,0} = 10^{24} \,
    {\rm g \, s^{-1}}$. We also highlight the duration of the ejector and propeller phases by
    shading the background in blue and red, respectively. See fig. 2 in
    [Ronchi et al. 2022](https://arxiv.org/abs/2201.11704).

    '''

    time, P, Pdot, B, alpha, r_m, r_c, r_lc, Mdot_d, Mdot_din, ldip, lacc, _, _, _ \
        = Magnetar.rotational_evolution( B0, P_in, Mdot0, alpha0 )

    print(alpha)
    # Compute the inner and outer radius of the disk.
    r_in = np.minimum(r_m, r_lc)
    r_out = Magnetar.disk_outer_radius_t(time * const.YR_TO_S)




# def plot_magnetar():
    # Find the time when the star enters the propeller phase.
    t_prop = time[(r_m > r_c) & (r_m < r_lc)]
    if len(t_prop) == 0 : t_prop = t_max
    else: t_prop = t_prop[0]

    # Compute the inner and outer radius of the disk.
    r_in = np.minimum(r_m, r_lc)
    r_out = Magnetar.disk_outer_radius_t(time * const.YR_TO_S)

    # Plot the spin-period evolution.
    mask1 = r_m >= r_lc
    mask2 = r_m < r_lc

    fig, (ax1, ax2, ax3, ax4) = plt.subplots(
        4, 1, sharex="all", gridspec_kw={"wspace": 0, "hspace": 0.05}, figsize=(6, 10)
    )

    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.set_xlim(10 / const.YR_TO_S, t_max)
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
        1.0e-6,
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
        1.0e-6,
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
        y=NS_radius, color="black", linestyle="-", linewidth=1, rasterized=True,
    )
    ax2.axhspan(
        ymin=0.0, ymax=NS_radius, facecolor="tab:gray", hatch="/", alpha=1, rasterized=True,
    )
    ax2.text(1.0e-6, r_m[0] * 1.5, r"$r_{\rm m}$", fontsize=12, rotation=0, color="black")
    ax2.text(1.0e-6, r_lc[0] / 2.5, r"$r_{\rm lc}$", fontsize=12, rotation=0, color="black")
    ax2.text(1.0e-6, r_c[0] / 2.5, r"$r_{\rm cor}$", fontsize=12, rotation=0, color="black")
    ax2.text(2, r_in[900], r"$r_{\rm in}$", fontsize=12, rotation=30, color="tab:red")
    ax2.text(0.3, r_out[1000], r"$r_{\rm out}$", fontsize=12, rotation=25, color="tab:red")
    ax2.text(
        1.0e-2,
        NS_radius * 1.5,
        r"Neutron star surface",
        fontsize=14,
        rotation=0,
        color="black",
        )

    ax3.text(1.0e-6, 5.0e2, r"Ejector", fontsize=12, rotation=0, color="black")
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

    ax4.text(1.0e-6, ldip[0] / 10, r"$L_{\rm dip}$", fontsize=12, rotation=0, color="blue")
    ax4.text(1.0e-6, lacc[0] / 700, r"$L_{\rm acc}$", fontsize=12, rotation=0, color="red")
    ax4.text(1.0e-6, (ldip+lacc)[0] / 300, r"$L_{\rm tot}$", fontsize=12, rotation=0, color="black")

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

    ax1.grid()
    ax2.grid()
    ax3.grid()
    ax4.grid()
    ax1.set_xlim(time.min(),time.max())
    ax2.set_xlim(time.min(),time.max())
    ax3.set_xlim(time.min(),time.max())
    ax4.set_xlim(time.min(),time.max())
    fig.savefig("propeller_evolution_ex_Bdecay.pdf", bbox_inches="tight")
    # plt.savefig()
    plt.show()

if __name__ == '__main__':
    run_magnetar()

