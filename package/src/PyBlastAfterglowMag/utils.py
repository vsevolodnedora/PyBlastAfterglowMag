import copy
import os.path
import numpy as np
import h5py
import hashlib

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

get_beta = lambda Gamma: np.sqrt(1. - np.power(Gamma, -2))
get_Gamma = lambda beta: np.float64(np.sqrt(1. / (1. - np.float64(beta) ** 2.)))

def get_Beta(Gamma):
    return (1. / Gamma) * np.sqrt((Gamma - 1.) * (Gamma + 1.))

def MomFromGam(gam):
    return np.sqrt(gam*gam - 1.)

def GamFromMom(mom):
    return np.sqrt(1.0+mom*mom)

def BetFromMom(mom):
    return mom / GamFromMom(mom)




def latex_float(f, format="{0:.2g}"):
    float_str = format.format(f)
    if "e" in float_str:
        base, exponent = float_str.split("e")
        if (base == "1"):
            return r"10^{{{0}}}".format(int(exponent))
        return r"{0} \times 10^{{{1}}}".format(base, int(exponent))
    else:
        return float_str

def make_prefix_for_pars(pars, keys=("n_ism", "theta_obs", "eps_e", "eps_b", "p", "eps_t", "nlayers", "d_l")):
    fname = ""
    for key in keys:
        if ((key == "theta_obs")or(key=="theta_c")or(key=="theta_w")):
            val = "{:.1f}".format( float(pars[key]/np.pi*180.) )
        elif (key == "d_l"):
            val = pars[key] / 1.e9 / cgs.pc
        else:
            val = pars[key]
        if len(str(val)) > 7:
            val = "{:.5f}".format(val)
        val = str(val).replace(".", "")
        fname+=key.replace("_","")
        fname+=val
        fname+="_"
    return fname

def find_nearest_index(array, value):
    ''' Finds index of the value in the array that is the closest to the provided one '''
    idx = (np.abs(array - value)).argmin()
    return idx

def make_hash(workdir : str, grb_skymap_fpath : str, kn1_skymap_fpath : str, kn2_skymap_fpath : str, time : float, freq  : float):
    assert os.path.isdir(workdir)
    if grb_skymap_fpath is None : grb_skymap_fpath = "none"
    if kn1_skymap_fpath is None : kn1_skymap_fpath = "none"
    if kn2_skymap_fpath is None : kn2_skymap_fpath = "none"
    time = ("{:.1f}".format(np.float(time))).replace(".", "") if (not type(time) is str) else time
    freq = ("{:.1f}".format(np.float(freq))).replace(".", "") if (not type(freq) is str) else freq
    s = "KEY--workdir=" + workdir.replace('/', "_") \
        + "--grb_skymap_fname=" + grb_skymap_fpath.replace('/', "_") \
        + "--kn1_skymap_fname=" + kn1_skymap_fpath.replace('/', "_") \
        + "--kn2_skymap_fname=" + kn2_skymap_fpath.replace('/', "_") \
        + "--time={}".format(time) \
        + "--freq={}".format(freq) \
        + "--END"
    hash = hashlib.sha1(str.encode(s)).hexdigest()
    return (hash, s)

def d2d(default: dict, new: dict):
    default_ = copy.deepcopy(default)
    for key, new_dict_ in new.items():
        if not key in default_.keys():
            default_[key] = {}
        for key_, val_ in new_dict_.items():
            default_[key][key_] = val_
    return default_