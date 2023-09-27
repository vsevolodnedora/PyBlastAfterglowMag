# import PyBlastAfterglowMag

# import pytest
import os
# import path

# from PyBlastAfterglowMag import BPA_METHODS as PBA
# from package.src.PyBlastAfterglowMag.interface import BPA_METHODS as PBA
# from package.src.PyBlastAfterglowMag

# try:
#     from PyBlastAfterglowMag.interface import *
#     from PyBlastAfterglowMag.utils import *
#     from PyBlastAfterglowMag.id_maker_analytic import *
#
# except ImportError:
#     try:
#         from package.src.PyBlastAfterglowMag.interface import modify_parfile_par_opt
#         from package.src.PyBlastAfterglowMag.interface import PyBlastAfterglow
#         from package.src.PyBlastAfterglowMag.interface import (distribute_and_parallel_run, get_str_val, set_parlists_for_pars)
#         from package.src.PyBlastAfterglowMag.utils import (latex_float, cgs, get_beta, get_Gamma)
#         from package.src.PyBlastAfterglowMag.id_maker_analytic import prepare_grb_ej_id_1d, prepare_grb_ej_id_2d
#     except ImportError:
#         raise ImportError("Cannot import PyBlastAfterglowMag")


try:
    import afterglowpy as grb
except:
    afterglowpy = False
    print("Error! could not import afteglowpy")

curdir = os.getcwd() + "/structured" + '/' #"/home/vsevolod/Work/GIT/GitHub/PyBlastAfterglow_dev/PyBlastAfterglow/src/PyBlastAfterglow/tests/dyn/"
# parfiledir = os.getcwd().replace("structured","")
# paperfigdir = "/home/vsevolod/Work/GIT/overleaf/grb_afg/figs/"

# directory = Path(__file__)
# sys.path.append(directory.parent.parent)
# directory = os.path.dirname(os.path.abspath("__file__"))
# sys.path.append(os.path.dirname(os.path.dirname(directory)))

# from projects.grb_prj.grb_tst_tools import *
from grb_tst_tools import *
from .settings import SettingsGaussian

def main_structured():
    grb = TestCasesFS(default_parfile_fpath=curdir+"parfile_def.par",
                      workingdir=curdir+"output/")
    tsk = SettingsGaussian()
    runset_for_skymap(tsk=tsk,
                      default_parfile_fpath=curdir+"parfile_def.par",
                      workingdir=curdir+"output/",
                      figdir=curdir+"figs/")

def main():
    main_structured()
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
