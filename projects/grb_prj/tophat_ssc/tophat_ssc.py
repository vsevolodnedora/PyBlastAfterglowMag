import os
import numpy as np

try:
    import package.src.PyBlastAfterglowMag as PBA
except ImportError:
    try:
        import PyBlastAfterglowMag as PBA
    except:
        raise ImportError("Cannot import PyBlastAfterglowMag")

try:
    import afterglowpy as grb
except:
    afterglowpy = False
    print("Error! could not import afteglowpy")



def main():
    structure = {
        "struct":"tophat",
        "Eiso_c":1.e53, "Gamma0c": 1000., "M0c": -1.,
        "theta_c": 0.1, "theta_w": 0.1, "nlayers_pw": 150, "nlayers_a": 1
    }
    pars = {
        "skymap_times" : np.array([1.,10.,40.,100.,200.,300.,400.,600.,800.,1000.,1300.,1600.,2000.,2500.,3000.,4000.,5000.,7000,10000]),
        "obs_freq":3e9,
        "n_ism": 1e-2,      "eps_e": 0.1,
        "d_l": 3.09e26,     "eps_b": 0.01,
        "z": 0.028,         "p": 2.2,
        "theta_obs": 1.57,#0.785,#0.785,
        "nsublayers":1
    }
    opts_a = {"method_eats": "adaptive", "method_spread": "AFGPY"}
    opts_pw = {"method_eats": "piece-wise", "method_spread": "AA"}

    run = PBA.wrappers.CasesFS(default_parfile_fpath=os.getcwd()+'/'+"parfile_def.par",workingdir=os.getcwd()+'/')
    pba = run.run_a(pars=pars, struct=structure,opts=opts_a,opts_grb=opts_a)
    run.plot_generic()




if __name__ == '__main__':
    main()