# import pytest
import os
# import path


try:
    import afterglowpy as grb
except:
    afterglowpy = False
    print("Error! could not import afteglowpy")

try:
    import package.src.PyBlastAfterglowMag as PBA
except ImportError:
    try:
        import PyBlastAfterglowMag as PBA
    except:
        raise ImportError("Cannot import PyBlastAfterglowMag")

curdir = os.getcwd() + "/structured" + '/' #"/home/vsevolod/Work/GIT/GitHub/PyBlastAfterglow_dev/PyBlastAfterglow/src/PyBlastAfterglow/tests/dyn/"

from .settings import SettingsGaussian

def main_structured_skymap():
    grb = PBA.wrappers.CasesFS(default_parfile_fpath=curdir+"parfile_def.par",
                      workingdir=curdir+"output/")
    tsk = SettingsGaussian()
    PBA.wrappers.runset_for_skymap(tsk=tsk,
                      default_parfile_fpath=curdir+"parfile_def.par",
                      workingdir=curdir+"output/",
                      figdir=curdir+"figs/")

def main():
    main_structured_skymap()
    # grb.plot_generic(struct=tsk.structure, pars=tsk.pars, opts_a=tsk.opts_a, title=tsk.figname)
    # grb.paper_plot_compare_spreading(struct=tsk.structure, pars=tsk.pars, opts_a=tsk.opts_a, title=tsk.figname)
    # grb.compare_skymaps_3d_theta_im_max(struct=tsk.structure, pars=tsk.pars, opts_a=tsk.opts_a, title=tsk.figname,
    #                                     theta_maxs=(
    #                                         (1.57,), #0.4, .9, 1.2,  1.57),
    #                                         # (21,121,121,121,121),
    #                                         ((20,20),),#121,121,121,121),
    #                                         ('red','orange','yellow', 'cyan', 'lime')))

if __name__ == '__main__':
    main()
    exit(0)
