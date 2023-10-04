import os
from matplotlib import cm
import os
curdir = os.getcwd() + '/'
from settings import SettingsGaussian,SettingsGRB170917A

try:
    import package.src.PyBlastAfterglowMag as PBA
except ImportError:
    try:
        import PyBlastAfterglowMag as PBA
    except:
        raise ImportError("Cannot import PyBlastAfterglowMag")

def main():
    tsk = SettingsGRB170917A()
    grb = PBA.wrappers.CasesFS(default_parfile_fpath=curdir+"parfile_def.par", workingdir=curdir+"output/")
    grb.plot_3d_stacked_skymaps(struct=tsk.strucure, pars=tsk.pars, opts_a=tsk.opts_a, opts_pw=tsk.opts_pw, title=tsk.figname,
                                figpath=curdir+"figs/"+"skymaps_3D_170817A", show_fig=True, save_pdf=True)

if __name__ == '__main__':
    main()