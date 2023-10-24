import numpy as np

class SettingsTopHat:
    name_1 = "Tophat"
    figname = "tophat"
    structure = {
        "struct":"tophat",
        "Eiso_c":1.e53, "Gamma0c": 1000., "M0c": -1.,
        "theta_c": 0.1, "theta_w": 0.1, "nlayers_pw": 150, "nlayers_a": 1
    }
    pars = {
        "skymap_times" : np.array([1.,10.,40.,100.,200.,300.,400.,600.,800.,1000.,1300.,1600.,2000.,2500.,3000.,4000.,5000.,7000,10000]),
        "obs_freq":3e9,
        "n_ism": 1e1,      "eps_e": 0.1,
        "d_l": 3.09e26,     "eps_b": 0.01,
        "z": 0.028,         "p": 2.2,
        "theta_obs": 1.57,#0.785,#0.785,
        "nsublayers":1
    }
    opts_a = {"method_eats": "adaptive", "method_spread": "AFGPY"}
    opts_pw = {"method_eats": "piece-wise", "method_spread": "AA"}

    resolutions_compare_plot = \
        ((80,120,160,200), # pw
         ((1,10,20),(1,20,40),(1,40,80),(1,80,120)), # a
         ('red','orange','yellow', 'cyan'))

    resolutions_pw=((50,100,150,200),
                    ('red','orange','yellow', 'cyan'),
                    ("PW",))

    resolutions_a=(((1,20,0),(1,40,0),(1,60,0),(1,80,0)), # a
                   ('red','orange','yellow', 'cyan'),
                   ("A81","A121","A161","A201"))