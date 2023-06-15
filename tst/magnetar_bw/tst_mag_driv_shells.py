from scipy.interpolate import RegularGridInterpolator, interp1d
from scipy.integrate import quad, cumtrapz
import numpy as np
import astropy.constants as cc
import astropy.units as uu
from collections import namedtuple

ev_to_erg = 1.60218e-12
speed_of_light = cc.c.cgs.value
planck = cc.h.cgs.value
proton_mass = cc.m_p.cgs.value
electron_mass = cc.m_e.cgs.value
mev_cgs = uu.MeV.cgs.scale
qe = 4.80320425e-10
solar_mass = cc.M_sun.cgs.value
sigma_sb = cc.sigma_sb.cgs.value
sigma_T = cc.sigma_T.cgs.value
radiation_constant = sigma_sb*4 / speed_of_light
boltzmann_constant = cc.k_B.cgs.value
km_cgs = uu.km.cgs.scale
day_to_s = 86400
au_cgs = uu.au.cgs.scale
solar_radius = cc.R_sun.cgs.value
graviational_constant = cc.G.cgs.value
angstrom_cgs = uu.Angstrom.cgs.scale
ang_to_hz = 2.9979245799999995e+18
speed_of_light_si = cc.c.si.value
solar_radius = cc.R_sun.cgs.value
jansky = 1e-26 # W m^-2 Hz^-1

def interpolated_barnes_and_kasen_thermalisation_efficiency(mej, vej):
    """
    Uses Barnes+2016 and interpolation to calculate the r-process thermalisation efficiency
    depending on the input mass and velocity
    :param mej: ejecta mass in solar masses
    :param vej: initial ejecta velocity as a fraction of speed of light
    :return: av, bv, dv constants in the thermalisation efficiency equation Eq 25 in Metzger 2017
    """
    v_array = np.array([0.1, 0.2, 0.3])
    mass_array = np.array([1.0e-3, 5.0e-3, 1.0e-2, 5.0e-2])
    a_array = np.asarray([[2.01, 4.52, 8.16], [0.81, 1.9, 3.2],
                          [0.56, 1.31, 2.19], [.27, .55, .95]])
    b_array = np.asarray([[0.28, 0.62, 1.19], [0.19, 0.28, 0.45],
                          [0.17, 0.21, 0.31], [0.10, 0.13, 0.15]])
    d_array = np.asarray([[1.12, 1.39, 1.52], [0.86, 1.21, 1.39],
                          [0.74, 1.13, 1.32], [0.6, 0.9, 1.13]])
    a_func = RegularGridInterpolator((mass_array, v_array), a_array, bounds_error=False, fill_value=None)
    b_func = RegularGridInterpolator((mass_array, v_array), b_array, bounds_error=False, fill_value=None)
    d_func = RegularGridInterpolator((mass_array, v_array), d_array, bounds_error=False, fill_value=None)

    av = a_func([mej, vej])[0]
    bv = b_func([mej, vej])[0]
    dv = d_func([mej, vej])[0]
    return av, bv, dv

def electron_fraction_from_kappa(kappa):
    """
    Uses interpolation from Tanaka+19 to calculate
    the electron fraction based on the temperature independent gray opacity
    :param kappa: temperature independent gray opacity
    :return: electron_fraction
    """

    kappa_array = np.array([1, 3, 5, 20, 30])
    ye_array = np.array([0.4,0.35,0.25,0.2, 0.1])
    kappa_func = interp1d(kappa_array, y=ye_array)
    electron_fraction = kappa_func(kappa)
    return electron_fraction

def _general_metzger_magnetar_driven_kilonova_model(
        time, mej_tot, vej_tot, beta_vel_dist_const, kappa, magnetar_luminosity,
        use_gamma_ray_opacity,
        **kwargs):
    """
    :param time: time array to evaluate model on in source frame in seconds
    :param redshift: redshift
    :param mej_tot: ejecta mass in solar masses
    :param vej_tot: minimum initial velocity
    :param beta_vel_dist_const: velocity power law slope (M=v^-beta)
    :param kappa: opacity
    :param magnetar_luminosity: evaluated magnetar luminosity in source frame
    :param pair_cascade_switch: whether to account for pair cascade losses
    :param use_gamma_ray_opacity: whether to use gamma ray opacity to calculate thermalisation efficiency
    :param kwargs: Additional parameters
    :param ejecta albedo: ejecta albedo; default is 0.5
    :param pair_cascade_fraction: fraction of magnetar luminosity lost to pair cascades; default is 0.05
    :param kappa_gamma: Gamma-ray opacity for leakage efficiency, only used if use_gamma_ray_opacity = True
    :param thermalisation_efficiency: magnetar thermalisation efficiency only used if use_gamma_ray_opacity = False
    :param neutron_precursor_switch: whether to have neutron precursor emission, default True
    :param pair_cascade_switch: whether to account for pair cascade losses, default is True
    :param magnetar_heating: whether magnetar heats all layers or just the bottom layer.
    :param vmax: maximum initial velocity of mass layers, default is 0.7c
    :return: named tuple with 'lorentz_factor', 'bolometric_luminosity', 'temperature',
                'r_photosphere', 'kinetic_energy_updated','erad_total', 'thermalisation_efficiency_total_ej'
    """
    pair_cascade_switch = kwargs.get('pair_cascade_switch', True)
    ejecta_albedo = kwargs.get('ejecta_albedo', 0.5)
    pair_cascade_fraction = kwargs.get('pair_cascade_fraction', 0.01)
    neutron_precursor_switch = kwargs.get('neutron_precursor_switch', True)
    magnetar_heating = kwargs.get('magnetar_heating', 'first_layer')
    vmax = kwargs.get('vmax', 0.7)

    tdays = time/day_to_s
    time_len = len(time)
    mass_len = 200

    # set up kilonova physics
    av, bv, dv = interpolated_barnes_and_kasen_thermalisation_efficiency(mej_tot, vej_tot)
    # thermalisation from Barnes+16
    e_th = 0.36 * (np.exp(-av * tdays) + np.log1p(2.0 * bv * tdays ** dv) / (2.0 * bv * tdays ** dv))
    electron_fraction = electron_fraction_from_kappa(kappa)
    t0 = 1.3 #seconds
    sig = 0.11  #seconds
    tau_neutron = 900  #seconds

    # convert to astrophysical units [my: TOTAL EJECTA MASS AND VELOCITY]
    m_ej_tot_0 = mej_tot * solar_mass
    vej_tot_0 = vej_tot * speed_of_light
    ek_total_0 = 0.5 * m_ej_tot_0 * vej_tot_0 ** 2 #[my: Total Ek of the entire ejecta]

    # set up mass and velocity layers [MY: CREATE VELOCITY SHELLS FROM TOTAL USING SLOPE]
    vmin = vej_tot
    vel = np.linspace(vmin, vmax, mass_len) # [My: arrange velocity [v_0=min... ... .. v_n=max]]
    # [My: arrange masses in these shells [m_0... ... ... m_n]] m ~ v^-beta
    mass_shells = mej_tot * (vel / vmin) ** (-beta_vel_dist_const)
    #
    vel_mass_shells = vel * speed_of_light

    # set up arrays
    time_array = np.tile(time, (mass_len, 1))
    e_th_array = np.tile(e_th, (mass_len, 1))
    edotr = np.zeros((mass_len, time_len)) # nuclear heating in ejecta

    time_mask = time > t0
    time_1 = time_array[:, time_mask]
    time_2 = time_array[:, ~time_mask]
    edotr[:, time_mask] = 2.1e10 * e_th_array[:, time_mask] * ((time_1/ (3600. * 24.)) ** (-1.3))
    edotr[:, ~time_mask] = 4.0e18 * (0.5 - (1. / np.pi) * np.arctan((time_2 - t0) / sig)) ** (1.3) * e_th_array[:,~time_mask]

    # [my: array of magnetar spin-down luminocity]
    lsd = magnetar_luminosity

    # set up empty arrays
    ek_mass_shells = np.zeros((mass_len, time_len))
    lum_rad = np.zeros((mass_len, time_len))
    qdot_rp = np.zeros((mass_len, time_len))
    td_v = np.zeros((mass_len, time_len))
    tau = np.zeros((mass_len, time_len))
    v_photosphere = np.zeros(time_len)
    v0_array = np.zeros(time_len)
    qdot_magnetar = np.zeros(time_len)
    r_photosphere = np.zeros(time_len)

    if neutron_precursor_switch == True:
        neutron_mass = 1e-8 * solar_mass # [my: Fixed fraction of the neutrons]
        # [my: is it from korobkin???]
        neutron_mass_fraction = 1 - 2*electron_fraction * 2 * np.arctan(neutron_mass / mass_shells / solar_mass) / np.pi
        rprocess_mass_fraction = 1.0 - neutron_mass_fraction
        initial_neutron_mass_fraction_array = np.tile(neutron_mass_fraction, (time_len, 1)).T
        rprocess_mass_fraction_array = np.tile(rprocess_mass_fraction, (time_len, 1)).T
        neutron_mass_fraction_array = initial_neutron_mass_fraction_array*np.exp(-time_array / tau_neutron)
        edotn = 3.2e14 * neutron_mass_fraction_array
        edotn = edotn * neutron_mass_fraction_array
        edotr = edotn + edotr
        kappa_n = 0.4 * (1.0 - neutron_mass_fraction_array - rprocess_mass_fraction_array)
        kappa = kappa * rprocess_mass_fraction_array
        kappa = kappa_n + kappa

    dt = np.diff(time)
    dm = np.abs(np.diff(mass_shells))

    #initial conditions
    ek_mass_shells[:, 0] = 0.5 * mass_shells*vel_mass_shells**2
    lum_rad[:, 0] = 0
    qdot_rp[:, 0] = 0
    kinetic_energy_updated = ek_total_0

    # solve ODE using euler method for all mass shells v
    for ii in range(time_len - 1):
        # # evolve the velocity due to pdv work of central shell of mass M and thermal energy Ev0
        kinetic_energy_updated = kinetic_energy_updated + (ek_mass_shells[0, ii] / time[ii]) * dt[ii]
        # kinetic_energy_updated = kinetic_energy_updated + (np.sum(ek_mass_shells[:, ii]) / time[ii]) * dt[ii]
        vej_tot_0 = (2 * kinetic_energy_updated / m_ej_tot_0) ** 0.5
        v0_array[ii] = vej_tot_0
        vel_mass_shells = vej_tot_0 * (mass_shells / (mej_tot)) ** (-1 / beta_vel_dist_const)
        vel_mass_shells[vel_mass_shells > 3e10] = speed_of_light

        if use_gamma_ray_opacity:
            kappa_gamma = kwargs["kappa_gamma"]
            prefactor = 3 * kappa_gamma * mej_tot / (4 * np.pi * vej_tot ** 2)
            thermalisation_efficiency_total_ej = 1 - np.exp(-prefactor * time[ii] ** -2)
        else:
            thermalisation_efficiency_total_ej = kwargs["thermalisation_efficiency"]
        qdot_magnetar[ii] = thermalisation_efficiency_total_ej * lsd[ii]

        if magnetar_heating == 'all_layers':
            if neutron_precursor_switch:
                td_v[:-1, ii] = (kappa[:-1, ii] * mass_shells[:-1] * solar_mass * 3) \
                                / (4 * np.pi * vel_mass_shells[:-1] * speed_of_light * time[ii] * beta_vel_dist_const)
            else:
                td_v[:-1, ii] = (kappa * mass_shells[:-1] * solar_mass * 3) \
                                / (4 * np.pi * vel_mass_shells[:-1] * speed_of_light * time[ii] * beta_vel_dist_const)
            lum_rad[:-1, ii] = ek_mass_shells[:-1, ii] \
                               / (td_v[:-1, ii] + time[ii] * (vel_mass_shells[:-1] / speed_of_light))
            ek_mass_shells[:-1, ii + 1] = ek_mass_shells[:-1, ii] \
                                          + (qdot_magnetar[ii]
                                            + edotr[:-1, ii] * dm * solar_mass
                                            - (ek_mass_shells[:-1, ii] / time[ii])
                                            - lum_rad[:-1, ii]
                                            ) * dt[ii]

        # first mass layer
        # only bottom layer i.e., 0'th mass layer gets magnetar contribution
        if magnetar_heating == 'first_layer':
            # [My: is it td_v[] diffusion timescale???]
            if neutron_precursor_switch:
                td_v[0, ii] = (kappa[0, ii] * mass_shells[0] * solar_mass * 3) \
                            / (4 * np.pi * vel_mass_shells[0] * speed_of_light * time[ii] * beta_vel_dist_const)
                td_v[1:-1, ii] = (kappa[1:-1, ii] * mass_shells[1:-1] * solar_mass * 3) \
                               / (4 * np.pi * vel_mass_shells[1:-1] * speed_of_light * time[ii] * beta_vel_dist_const)
            else:
                td_v[0, ii] = (kappa * mass_shells[0] * solar_mass * 3) \
                            / (4 * np.pi * vel_mass_shells[0] * speed_of_light * time[ii] * beta_vel_dist_const)
                td_v[1:-1, ii] = (kappa * mass_shells[1:-1] * solar_mass * 3) \
                               / (4 * np.pi * vel_mass_shells[1:-1] * speed_of_light * time[ii] * beta_vel_dist_const)

            # [my: radiative luminocity from each shell]
            lum_rad[0, ii] = ek_mass_shells[0, ii] \
                             / (td_v[0, ii] + time[ii] * (vel_mass_shells[0] / speed_of_light))
            ek_mass_shells[0, ii + 1] = ek_mass_shells[0, ii] + \
                                        (qdot_magnetar[ii] # [my: thermalized magnetar emission]
                                            + edotr[0, ii] * dm[0] * solar_mass # [my: thermalized r-process heating]
                                            - (ek_mass_shells[0, ii] / time[ii]) # [my: ??ADIABATIC??]
                                            - lum_rad[0, ii] # [my: energy lost to thermal radiation]
                                         ) * dt[ii] # [my: dE/dt * dt]

            # other layers [my: radiative luminocity from each shell]
            lum_rad[1:-1, ii] = ek_mass_shells[1:-1, ii] \
                                / (td_v[1:-1, ii] + time[ii] * (vel_mass_shells[1:-1] / speed_of_light))
            ek_mass_shells[1:-1, ii + 1] = ek_mass_shells[1:-1, ii] + \
                                           (edotr[1:-1, ii] * dm[1:] * solar_mass # [my: thermalized r-process heating]
                                            - (ek_mass_shells[1:-1, ii] / time[ii]) # [my: ??ADIABATIC??]
                                            - lum_rad[1:-1, ii] # [my: energy lost to thermal radiation]
                                          ) * dt[ii] # [my: dE/dt * dt]


        # [My: tau[] is the optical depth?]
        if neutron_precursor_switch:
            tau[:-1, ii] = (mass_shells[:-1] * solar_mass * kappa[:-1, ii] # [My: Only -1 element!!! ]
                         / (4 * np.pi * (time[ii] * vel_mass_shells[:-1]) ** 2))
        else:
            tau[:-1, ii] = (mass_shells[:-1] * solar_mass * kappa  #
                         / (4 * np.pi * (time[ii] * vel_mass_shells[:-1]) ** 2))

        # [My: photosphere location based on tau = 1 shell]
        tau[mass_len - 1, ii] = tau[mass_len - 2, ii] # [my: ?????????]
        photosphere_index = np.argmin(np.abs(tau[:, ii] - 1))
        v_photosphere[ii] = vel_mass_shells[photosphere_index]
        r_photosphere[ii] = v_photosphere[ii] * time[ii] # r = v * t

    bolometric_luminosity = np.sum(lum_rad, axis=0)

    if pair_cascade_switch == True:
        tlife_t = (0.6/(1 - ejecta_albedo))\
                  *(pair_cascade_fraction/0.1)**0.5 \
                  * (lsd/1.0e45)**0.5 \
                  * (vej_tot_0/(0.3*speed_of_light))**(0.5) \
                  * (time/day_to_s)**(-0.5)
        bolometric_luminosity = bolometric_luminosity / (1.0 + tlife_t)

    # black body spectrum
    temperature = (bolometric_luminosity / (4.0 * np.pi * (r_photosphere) ** (2.0) * sigma_sb)) ** (0.25)

    dynamics_output = namedtuple('dynamics_output', ['lorentz_factor', 'bolometric_luminosity', 'temperature',
                                                     'r_photosphere', 'kinetic_energy_updated','erad_total',
                                                     'thermalisation_efficiency_total_ej'])
    gamma_beta = v0_array/speed_of_light
    lorentz_factor = 1/(np.sqrt(1 - gamma_beta**2))
    dynamics_output.lorentz_factor = lorentz_factor
    dynamics_output.bolometric_luminosity = bolometric_luminosity
    dynamics_output.temperature = temperature
    dynamics_output.r_photosphere = r_photosphere
    dynamics_output.kinetic_energy = (lorentz_factor - 1)*m_ej_tot_0*speed_of_light**2 # total Ek
    dynamics_output.erad_total = np.trapz(bolometric_luminosity, x=time)
    dynamics_output.thermalisation_efficiency = qdot_magnetar/lsd
    return dynamics_output

def _evolving_gw_and_em_magnetar(time, bint, bext, p0, chi0, radius, moi, **kwargs):
    """
    Assumes a combination of GW and EM spin down with a constant spin-magnetic field inclination angle.
    Only EM contributes to observed emission.

    :param time: time in source frame in seconds (must be a large array as this function is semianalytic)
    :param bint: internal magnetic field in G
    :param bext: external magnetic field in G
    :param p0: spin period in s
    :param chi0: initial inclination angle
    :param radius: radius of NS in cm
    :param moi: moment of inertia of NS
    :param kwargs: None
    :return: luminosity
    """
    epsilon_b = -3e-4 * (bint / bext) ** 2 * (bext / 1e16) ** 2
    omega_0 = 2.0 * np.pi / p0
    erot = 0.5 * moi * omega_0**2

    dt = time[1:] - time[:-1]
    omega = np.zeros_like(time)
    chi = chi0

    omega[0] = omega_0

    for i in range(len(time) - 1):
        omega[i + 1] = omega[i] + dt[i] * (
                -(bext**2*radius**6/(moi*speed_of_light**3))*omega[i]**3*(1. + np.sin(chi)) - (
                2.0*graviational_constant*moi*epsilon_b**2/(5.0*speed_of_light**5)) * omega[i]**5*np.sin(chi)**2*(1.0+15.0*np.sin(chi)**2))

    Edot_d = (bext ** 2 * radius ** 6 / (4*speed_of_light ** 3)) * omega ** 4 * (1.0 + np.sin(chi)**2)
    Edot_gw = (2.0 * graviational_constant * moi ** 2 * epsilon_b ** 2 / (5.0 * speed_of_light ** 5)) * omega ** 6 * np.sin(chi)**2 * (
            1.0 + 15.0 * np.sin(chi)**2)

    Ed = cumtrapz(Edot_d, x=time)
    Egw = cumtrapz(Edot_gw, x=time)

    En_t =  3.5e50*(bint/1e17)**2*(radius/1.5e6)**3
    En_p = 5.5e47 * (bext / 1e14) ** 2 * (radius / 1.5e6) ** 3

    output = namedtuple('output', ['e_gw', 'e_em', 'tsd', 'epsilon_b', 'e_magnetic', 'Edot_d', 'Edot_gw', 'erot'])
    output.e_gw = Egw[-1]
    output.e_em = Ed[-1]
    output.erot = erot
    period = p0
    output.tsd = 2.4 * (period/1e-3)**2 *((bext/1e14)**(2) + 7.2*(bint/1e16)**4*(period/1e-3)**(-2))**(-1) * (60*60)
    output.epsilon_b = epsilon_b
    output.e_magnetic = En_t + En_p
    output.Edot_d = Edot_d
    output.Edot_gw = Edot_gw
    return output


@citation_wrapper('https://ui.adsabs.harvard.edu/abs/2022arXiv220514159S/abstract')
def general_metzger_magnetar_driven_evolution(time, redshift, mej, vej, beta, kappa_r, logbint,
                                              logbext, p0, chi0, radius, logmoi, kappa_gamma, **kwargs):
    """
    :param time: observer frame time in days
    :param redshift: redshift
    :param mej: ejecta mass in solar masses
    :param vej: minimum initial velocity
    :param beta: velocity power law slope (M=v^-beta)
    :param kappa_r: opacity
    :param logbint: log10 internal magnetic field in G
    :param logbext: log10 external magnetic field in G
    :param p0: spin period in s
    :param chi0: initial inclination angle
    :param radius: radius of NS in KM
    :param logmoi: log10 moment of inertia of NS
    :param kappa_gamma: gamma-ray opacity used to calculate magnetar thermalisation efficiency
    :param kwargs: Additional parameters
    :param ejecta albedo: ejecta albedo; default is 0.5
    :param pair_cascade_fraction: fraction of magnetar luminosity lost to pair cascades; default is 0.05
    :param neutron_precursor_switch: whether to have neutron precursor emission, default true
    :param pair_cascade_switch: whether to account for pair cascade losses, default is True
    :param magnetar_heating: whether magnetar heats all layers or just the bottom layer. default first layer only
    :param vmax: maximum initial velocity of mass layers, default is 0.7c
    :param frequency: Required if output_format is 'flux_density'.
        frequency to calculate - Must be same length as time array or a single number).
    :param bands: Required if output_format is 'magnitude' or 'flux'.
    :param output_format: 'flux_density', 'magnitude', 'spectra', 'flux', 'sncosmo_source'
    :param lambda_array: Optional argument to set your desired wavelength array (in Angstroms) to evaluate the SED on.
    :return: set by output format - 'flux_density', 'magnitude', 'spectra', 'flux', 'sncosmo_source'
    """
    use_gamma_ray_opacity = True
    kwargs['use_relativistic_blackbody'] = False

    time_temp = np.geomspace(1e-4, 1e7, 500) #in source frame
    bint = 10 ** logbint
    bext = 10 ** logbext
    radius = radius * km_cgs
    moi = 10 ** logmoi
    dl = cosmo.luminosity_distance(redshift).cgs.value
    output = _evolving_gw_and_em_magnetar(time=time_temp, bint=bint, bext=bext, p0=p0, chi0=chi0, radius=radius, moi=moi)
    magnetar_luminosity = output.Edot_d
    output = _general_metzger_magnetar_driven_kilonova_model(time=time_temp, mej=mej, vej=vej, beta=beta, kappa=kappa_r,
                                                             magnetar_luminosity=magnetar_luminosity,
                                                             use_gamma_ray_opacity=use_gamma_ray_opacity,
                                                             kappa_gamma=kappa_gamma, **kwargs)
    time_obs = time
    if kwargs['output_format'] == 'flux_density':
        return _process_flux_density(dl=dl, output=output, redshift=redshift,
                                     time=time, time_temp=time_temp, **kwargs)
    else:
        return _processing_other_formats(dl=dl, output=output, redshift=redshift,
                                         time_obs=time_obs, time_temp=time_temp, **kwargs)
