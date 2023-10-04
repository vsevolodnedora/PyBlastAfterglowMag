import numpy as np


class SettingsGaussian:
    name_1 = "Tophat"
    figname = "gaussian"
    structure = {
        "struct":"gaussian",
        "Eiso_c":1.e54, "Gamma0c": 1000., "M0c": -1.,
        "theta_c": 0.1, "theta_w": 0.2, "nlayers_pw": 150, "nlayers_a": 52
    }
    pars = {
        "skymap_times" : np.array([1.,10.,40.,100.,200.,300.,400.,600.,800.,1000.,1300.,1600.,2000.,2500.,3000.,4000.,5000.,7000,10000]),
        "obs_freq":3e9,
        "n_ism": 1e-4,      "eps_e": 0.1,
        "d_l": 3.09e26,     "eps_b": 0.01,
        "z": 0.028,         "p": 2.2,
        "theta_obs": 0.9,#1.57,#0.785,

        "ntb": 5000
    }


    pars_fsrs = {
        "skymap_times" : np.array([1.,10.,40.,100.,200.,300.,400.,600.,800.,1000.,1300.,1600.,2000.,2500.,3000.,4000.,5000.,7000,10000]),
        "obs_freq":3e9,
        "n_ism": 1e-1,      "eps_e": 0.1,   "eps_e_rs": 0.2,
        "d_l": 3.09e26,     "eps_b": 0.01,  "eps_b_rs": 0.02,
        "z": 0.028,         "p": 2.2,       "p_rs": 2.4,
        "theta_obs": 0.9,

        "ntb":3000
    }

    opts_a = {"method_eats":"adaptive", "method_spread":"AFGPY"}
    opts_pw = {"method_eats":"piece-wise", "method_spread":"AA"}

    # opts_1 = {"rtol": 1e-6, "ntb":10000, "iout": 10}
    # opts_a_1 = {"method_eats":"adaptive", "method_spread":"AFGPY", "rhs_type":"grb_fs", "do_rs": "no"}

    resolutions_compare_plot = \
        ((50,100,150,250), # pw
         # ((22,20,40),(22,40,80),(22,80,160),(22,160,320)), # a
         # ((22,20,40),(42,20,40),(62,20,40),(82,20,40)), # a
         ((42,20,40),(42,40,60),(42,60,80),(42,80,100)),
         ('red','orange','yellow', 'cyan'))


    resolutions_pw=((50,100,150,200),
                    ('red','orange','yellow', 'cyan'),
                    ("PW",))

    resolutions_a=(#((22,20,40),(42,20,40),(62,20,40),(82,20,40)),
                    # ((42,20,40),(42,40,60),(42,60,80),(42,80,100)),
                    # ((12,10,10),(12,10,10),(12,10,10),(12,10,10)),
                    # ((22,20,40),(32,20,40),(42,20,40),(52,20,40)), # a
                    ((22,50,0),(32,50,0),(42,50,0),(52,50,0)),#(32,50,60),(42,50,80),(52,50,100)), # a
                   ('red','orange','yellow', 'cyan'),
                   ("22-40-50","32-30-40","A161","A201"))



class SettingsGRB170917A:
    figname = "170817A"
    strucure = {
        "struct":"gaussian",
        "Eiso_c":1.e52, "Gamma0c": 300., "M0c": -1.,
        "theta_c": 0.085, "theta_w": 0.2618, "nlayers_pw":150,"nlayers_a": 50
    }
    pars = {
        "skymap_times" : np.array([1.,10.,40.,100.,200.,300.,400.,600.,800.,1000.,1300.,1600.,2000.,2500.,3000.,4000.,5000.,7000,10000]),
        "obs_freq":3e9,
        "n_ism":0.00031,    "eps_e":0.0708,
        "d_l":1.27e+26,     "eps_b": 0.0052,
        "z": 0.0099,        "p":2.16,
        "theta_obs": 0.3752,
        "ntb": 5000, "iout": 5
    }
    opts_a = {"method_eats":"adaptive", "method_spread":"AFGPY"}
    opts_pw = {"method_eats":"piece-wise", "method_spread":"AA"}