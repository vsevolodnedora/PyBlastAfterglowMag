import numpy as np
import os
import h5py
from scipy import ndimage, interpolate
import copy
import subprocess
from shutil import copyfile
from multiprocessing import Pool

from .utils import cgs, get_beta, find_nearest_index

''' 
pip uninstall --no-cache-dir PyBlastAfterglowMag & pip install .
'''

def read_parfile(workingdir, fname="parfile.par", comment="#", sep1="", sep2=""):
    par_sep = "* Parameters"
    opt_sep = "* Settings"
    par_var_sep = " = "
    lines = []
    pars = {}
    opts = {}
    reading_pars = False
    reading_opts = False
    with open(workingdir+fname) as file:
        lines = [line.rstrip() for line in file]
    skip = True
    for iline, line in enumerate(lines):
        if (line == sep1):
            skip = False
        if (line == sep2):
            skip = True
        if (not skip):
            if (len(line) == 1 or line == '' or line[0] == comment):
                continue
            if (line == par_sep):
                reading_pars = True
                reading_opts = False
                continue
            if (line == opt_sep):
                reading_pars = False
                reading_opts = True
                continue
            if (reading_pars):
                if (line.__contains__("#")):
                    line = line.split("#")[0]
                if (not line.__contains__(par_var_sep)):
                    raise KeyError("line does not contain par-var separator {} : {}".format(par_var_sep, line))
                name=line.split(par_var_sep)[0]
                val=line.split(par_var_sep)[-1]
                val = np.double(val)
                pars[name] = val
            if (reading_opts):
                if (line.__contains__("#")):
                    line = line.split("#")[0]
                if (not line.__contains__(par_var_sep)):
                    raise KeyError("line does not contain par-var separator {} : {}".format(par_var_sep, line))
                name=line.split(par_var_sep)[0]
                val=line.split(par_var_sep)[-1]
                opts[name] = val
            if (iline==len(lines)-1):
                return ({}, {})
    # if (len(pars.keys())==0):
    #     raise ValueError("Empty pars.dict something went wrong.")
    # if (len(opts.keys())==0):
    #     raise ValueError("empty opts.dict: something went wrong.")
    return (pars, opts)

def modify_parfile(newpars : dict, newopts : dict, workingdir, comment="#",
                   fname="parfile.par", newfname="parfile2.par", sep1="", sep2=""):
    par_sep = "* Parameters"
    opt_sep = "* Settings"
    par_var_sep = " = "
    new_lines = []
    reading_pars = False
    reading_opts = False
    with open(workingdir+fname) as file:
        lines = [line.rstrip() for line in file]
    skip = True
    for iline, line in enumerate(lines):
        new_lines.append(line)
        if (line == sep1):
            skip = False
        if (line == sep2):
            skip = True
        if (not skip):
            if (len(line) == 1 or line == '' or line[0] == comment):
                continue
            if (line == par_sep):
                reading_pars = True
                reading_opts = False
                continue
            if (line == opt_sep):
                reading_pars = False
                reading_opts = True
                continue
            if (reading_pars):
                for key in newpars.keys():
                    if ((line.split(par_var_sep)[0]) == key):
                        _comment = ""
                        if (line.__contains__("#")):
                            _comment = " # "+ line.split(" # ")[-1]
                        new_lines[iline] = line.replace(line.split(par_var_sep)[-1],str(newpars[key])) + _comment
                        i = 1
            if (reading_opts):
                for key in newopts.keys():
                    if ((line.split(par_var_sep)[0]) == key):
                        _comment = ""
                        if (line.__contains__("#")):
                            _comment = " # "+ line.split(" # ")[-1]
                        new_lines[iline] = line.replace(line.split(par_var_sep)[-1],str(newopts[key])) + _comment
    # for line in new_lines :
    #     line += "\n"
    with open(newfname, 'w') as f:
        for line in new_lines:
            f.write(f"{line}\n")
    print("saved {}".format(newfname))


def modify_parfile_par_opt(workingdir : str, part : str, newpars : dict, newopts : dict,
                           parfile="parfile.par",newparfile="parfile2.par",keep_old=True):
    if (keep_old and parfile == newparfile):
        raise NameError("Cannot keep old parfile if the new name is the same as old")
    if not (os.path.isfile(workingdir+parfile)):
        raise FileNotFoundError("parfile {} not found".format(workingdir+parfile))
    copyfile(workingdir+parfile,workingdir+"tmp_{}".format(parfile))
    # -----------------------------------------------------------------------
    if (part=="main"):
        sep1="# -------------------------- main ---------------------------"
        sep2="# --------------------------- END ---------------------------"
    elif (part=="grb"):
        sep1="# ---------------------- GRB afterglow ----------------------"
        sep2="# --------------------------- END ---------------------------"
    elif (part=="kn"):
        sep1="# ----------------------- kN afterglow ----------------------"
        sep2="# --------------------------- END ---------------------------"
    elif (part=="magnetar"):
        sep1="# ------------------------ Magnetar -------------------------"
        sep2="# --------------------------- END ---------------------------"
    else:
        raise NameError("no part: {} in parfile".format(part))
    # -------------------------------------------------------------------------
    modify_parfile(newpars=newpars,newopts=newopts,workingdir=workingdir,comment="#",
                   fname="tmp_{}".format(parfile),
                   newfname="tmp_mod_{}".format(newparfile),
                   sep1=sep1,
                   sep2=sep2)
    if not keep_old:
        os.remove(workingdir+parfile)
    copyfile(workingdir+"tmp_mod_{}".format(newparfile),workingdir+newparfile)
    os.remove(workingdir+"tmp_{}".format(parfile))
    os.remove(workingdir+"tmp_mod_{}".format(newparfile))


class PBA_BASE:
    '''
        Pars initial data and parameters into C++ code interface
        Load result files
    '''
    def __init__(self, workingdir, readparfileforpaths=True,parfile="pafile.par"):

        if not os.path.isdir(workingdir):
            raise IOError("Working directory not found {}".format(workingdir))
        self.parfile = parfile
        self.workingdir = workingdir
        self.res_dir_mag = workingdir
        self.res_dir_kn = workingdir
        self.res_dir_grb = workingdir
        self.grb_prefix = "grb_"
        self.kn_prefix = "kn_"
        # ------------------ MAIN
        if readparfileforpaths:
            self.main_pars, self.main_opts = self.read_main_part_parfile( self.parfile)
        # ------------------ MAGNETAR
        self.fpath_mag = None
        if readparfileforpaths:
            self.mag_pars, self.mag_opts = self.read_magnetar_part_parfile( self.parfile)
        else:
            self.fpath_mag = self.res_dir_mag + "magnetar.h5"
        # ------------------ GRB
        self.fpath_grb_dyn = None
        self.fpath_grb_spec = None
        self.fpath_grb_light_curve = None
        self.fpath_grb_sky_map = None
        if readparfileforpaths:
            self.grb_pars, self.grb_opts = self.read_grb_part_parfile( self.parfile)
        else:
            self.fpath_grb_dyn = self.res_dir_kn + self.grb_prefix + "dynamics_layers.h5"
            self.fpath_grb_spec = self.res_dir_kn + self.grb_prefix + "spectra.h5"
            self.fpath_grb_light_curve = self.res_dir_kn + self.grb_prefix + "lightcurves_layers.h5"
            self.fpath_grb_sky_map = self.res_dir_kn + self.grb_prefix + "skymap.h5"
        # ------------------ EJECTA
        self.fpath_kn_dyn = None
        self.fpath_kn_spec = None
        self.fpath_kn_light_curve = None
        self.fpath_kn_sky_map = None
        if readparfileforpaths:
            self.kn_pars, self.kn_opts = self.read_kn_part_parfile( self.parfile)
        else:
            self.fpath_kn_dyn = self.res_dir_kn + self.kn_prefix + "dynamics_layers.h5"
            self.fpath_kn_spec = self.res_dir_kn + self.kn_prefix + "spectra_layers.h5"
            self.fpath_kn_light_curve = self.res_dir_kn + self.kn_prefix + "lightcurves_layers.h5"
            self.fpath_kn_sky_map = self.res_dir_kn + self.kn_prefix + "skymap.h5"
        # ------------------- MAGNETAR
        # TODO
        # -----------------
        self.mag_dfile = None
        self.grb_dyn_dfile = None
        self.grb_spec_dfile = None
        self.grb_lc_dfile = None
        self.grb_skymap_dfile = None
        self.kn_dyn_dfile = None
        self.kn_spec_dfile = None
        self.kn_lc_dfile = None
        self.kn_skymap_dfile = None
        # ----------------

    def clear(self):
        # self.overwrite = True
        if (not self.mag_dfile is None):
            self.mag_dfile.close()
            self.mag_dfile = None
        if (not self.grb_dyn_dfile is None):
            self.grb_dyn_dfile.close()
            self.grb_dyn_dfile = None
        if (not self.grb_spec_dfile is None):
            self.grb_spec_dfile.close()
            self.grb_spec_dfile = None
        if (not self.grb_lc_dfile is None):
            self.grb_lc_dfile.close()
            self.grb_lc_dfile = None
        if (not self.grb_skymap_dfile is None):
            self.grb_skymap_dfile.close()
            self.grb_skymap_dfile = None
        if (not self.kn_dyn_dfile is None):
            self.kn_dyn_dfile.close()
            self.kn_dyn_dfile = None
        if (not self.kn_spec_dfile is None):
            self.kn_spec_dfile.close()
            self.kn_spec_dfile = None
        if (not self.kn_lc_dfile is None):
            self.kn_lc_dfile.close()
            self.kn_lc_dfile = None
        if (not self.kn_skymap_dfile is None):
            self.kn_skymap_dfile.close()
            self.kn_skymap_dfile = None
    def read_main_part_parfile(self, parfile="parfile.par"):
        main_pars, main_opts = read_parfile(workingdir=self.workingdir,fname=parfile,comment="#",
                                          sep1="# -------------------------- main ---------------------------",
                                          sep2="# --------------------------- END ---------------------------")
        return (main_pars,main_opts)
    def read_magnetar_part_parfile(self, parfile="parfile.par"):
        mag_pars, mag_opts = read_parfile(workingdir=self.workingdir,fname=parfile,comment="#",
                                            sep1="# ------------------------ Magnetar -------------------------",
                                            sep2="# --------------------------- END ---------------------------")
        if "fname_mag" in mag_opts.keys(): self.fpath_mag = self.res_dir_mag + mag_opts["fname_mag"]
        return (mag_pars, mag_opts)
    def read_grb_part_parfile(self,parfile="parfile.par"):
        grb_pars, grb_opts = read_parfile(workingdir=self.workingdir, fname=parfile,comment="#",
                                          sep1="# ---------------------- GRB afterglow ----------------------",
                                          sep2="# --------------------------- END ---------------------------")
        if "fname_dyn" in grb_opts.keys(): self.fpath_grb_dyn = self.res_dir_grb + grb_opts["fname_dyn"]
        if "fname_spec" in grb_opts.keys(): self.fpath_grb_spec = self.res_dir_grb + grb_opts["fname_spec"]
        if "fname_light_curve" in grb_opts.keys(): self.fpath_grb_light_curve = self.res_dir_grb + grb_opts["fname_light_curve"]
        if "fname_sky_map" in grb_opts.keys(): self.fpath_grb_sky_map = self.res_dir_grb + grb_opts["fname_sky_map"]
        return (grb_pars,grb_opts)
    def read_kn_part_parfile(self,parfile="parfile.par"):
        kn_pars, kn_opts = read_parfile(workingdir=self.workingdir, fname=parfile,comment="#",
                                          sep1="# ----------------------- kN afterglow ----------------------",
                                          sep2="# --------------------------- END ---------------------------")
        if "fname_dyn" in kn_opts.keys(): self.fpath_kn_dyn = self.res_dir_kn + kn_opts["fname_dyn"]
        if "fname_spec" in kn_opts.keys(): self.fpath_kn_spec = self.res_dir_kn + kn_opts["fname_spec"]
        if "fname_light_curve" in kn_opts.keys(): self.fpath_kn_light_curve = self.res_dir_kn + kn_opts["fname_light_curve"]
        if "fname_sky_map" in kn_opts.keys(): self.fpath_kn_sky_map = self.res_dir_kn + kn_opts["fname_sky_map"]
        return (kn_pars,kn_opts)
    # def modify_parfile(self, base_parfile,
    #                    new_main_pars:dict, new_main_opts:dict,
    #                    new_grb_pars:dict, new_grb_opts:dict,
    #                    new_kn_pars:dict, new_kn_opts:dict):
    #     if(len(new_main_pars.keys())+len(new_main_opts.keys())>0):
    #         modify_parfile(newpars=new_main_pars,newopts=new_main_opts,workingdir=self.workingdir,
    #                        comment="#",fname=de,
    #                        newfname="parfile2.par",
    #                        sep1="# -------------------------- main ---------------------------",
    #                        sep2="# --------------------------- END ---------------------------")
    #         copyfile(self.workingdir+"parfile.par",self.workingdir+"old_parfile.par")
    #         copyfile(self.workingdir+"parfile2.par",self.workingdir+"parfile.par")
    #         os.remove(self.workingdir+"parfile2.par")
    #
    #
    #
    #
    #     if (reload_parfile): self.main_pars, self.main_opts = self.read_main_part_parfile()


    # def modify_main_part_parfile(self, newpars : dict, newopts : dict, parfile="parfile.par",newparfile="parfile2.par",keep_old=True):
    #     # if save_old:
    #     #     copyfile(self.workingdir+parfile,self.workingdir+"old_{}".format(parfile))
    #     if (keep_old and parfile == newparfile):
    #         raise NameError("Cannot keep old parfile if the new name is the same as old")
    #     copyfile(self.workingdir+parfile,self.workingdir+"tmp_{}".format(parfile))
    #     modify_parfile(newpars=newpars,newopts=newopts,workingdir=self.workingdir,comment="#",
    #                    fname="tmp_{}".format(parfile),
    #                    newfname="tmp_mod_{}".format(newparfile),
    #                    sep1="# -------------------------- main ---------------------------",
    #                    sep2="# --------------------------- END ---------------------------")
    #     if not keep_old:
    #         os.remove(self.workingdir+parfile)
    #     copyfile(self.workingdir+"tmp_mod_{}".format(newparfile),self.workingdir+newparfile)
    #     os.remove(self.workingdir+"tmp_{}".format(parfile))
    #     os.remove(self.workingdir+"tmp_mod_{}".format(newparfile))
    #     # if not save_old:
    #     #     os.remove(self.workingdir+"parfile2.par")
    #     # copyfile(self.workingdir+"parfile.par",self.workingdir+"old_parfile.par")
    #     # copyfile(self.workingdir+"parfile2.par",self.workingdir+"parfile.par")
    #     # os.remove(self.workingdir+"parfile2.par")
    #     # if (reload_parfile): self.main_pars, self.main_opts = self.read_main_part_parfile()
    # def modify_grb_part_parfile(self, newpars : dict, newopts : dict,parfile="parfile.par",newparfile="parfile2.par",keep_old=True):
    #     if (keep_old and parfile == newparfile):
    #         raise NameError("Cannot keep old parfile if the new name is the same as old")
    #     copyfile(self.workingdir+parfile,self.workingdir+"tmp_{}".format(parfile))
    #     modify_parfile(newpars=newpars,newopts=newopts,workingdir=self.workingdir,comment="#",
    #                    fname="tmp_{}".format(parfile),
    #                    newfname="tmp_mod_{}".format(newparfile),
    #                    sep1="# ---------------------- GRB afterglow ----------------------",
    #                    sep2="# --------------------------- END ---------------------------")
    #     if not keep_old:
    #         os.remove(self.workingdir+parfile)
    #     copyfile(self.workingdir+"tmp_mod_{}".format(newparfile),self.workingdir+newparfile)
    #     os.remove(self.workingdir+"tmp_{}".format(parfile))
    #     os.remove(self.workingdir+"tmp_mod_{}".format(newparfile))
    #     # copyfile(self.workingdir+"parfile.par",self.workingdir+"old_parfile.par")
    #     # copyfile(self.workingdir+"parfile2.par",self.workingdir+"parfile.par")
    #     # os.remove(self.workingdir+"parfile2.par")
    #     # if (reload_parfile): self.grb_pars, self.grb_opts = self.read_grb_part_parfile()
    # def modify_mag_part_parfile(self, newpars : dict, newopts : dict,parfile="parfile.par",newparfile="parfile2.par",keep_old=True):
    #     modify_parfile(newpars=newpars,newopts=newopts,workingdir=self.workingdir,comment="#",fname=parfile,
    #                    newfname="parfile2.par",
    #                    sep1="# ------------------------ Magnetar -------------------------",
    #                    sep2="# --------------------------- END ---------------------------")
    #     copyfile(self.workingdir+"parfile.par",self.workingdir+"old_parfile.par")
    #     copyfile(self.workingdir+"parfile2.par",self.workingdir+"parfile.par")
    #     os.remove(self.workingdir+"parfile2.par")
    #     # if (reload_parfile): self.grb_pars, self.grb_opts = self.read_grb_part_parfile()
    # def modify_kn_part_parfile(self, newpars : dict, newopts : dict,parfile="parfile.par",newparfile="parfile2.par",keep_old=True):
    #     modify_parfile(newpars=newpars,newopts=newopts,workingdir=self.workingdir,comment="#",fname=parfile,
    #                    newfname="parfile2.par",
    #                    sep1="# ----------------------- kN afterglow ----------------------",
    #                    sep2="# --------------------------- END ---------------------------")
    #     copyfile(self.workingdir+"parfile.par",self.workingdir+"old_parfile.par")
    #     copyfile(self.workingdir+"parfile2.par",self.workingdir+"parfile.par")
    #     os.remove(self.workingdir+"parfile2.par")
    #     # if (reload_parfile): self.kn_pars, self.kn_opts = self.read_kn_part_parfile()
    def reload_parfile(self):
        self.main_pars, self.main_opts = self.read_main_part_parfile()
        self.grb_pars, self.grb_opts = self.read_grb_part_parfile()
        self.kn_pars, self.kn_opts = self.read_kn_part_parfile()
    def run(self):
        # this mess is because I did not figure out how $PATH thing works...
        curdir = os.getcwd()
        pbadir = curdir.split("PyBlastAfterglowMag")[0]
        path_to_cpp_executable = pbadir+"PyBlastAfterglowMag"+"/src/pba.out"
        # print(os.getcwd())
        # os.chdir("../../../src/")
        # path_to_executable = "pba.out"
        if not os.path.isfile(path_to_cpp_executable):
            raise IOError("pba.out executable is not found: {}".format(path_to_cpp_executable))
        # subprocess.call(path_to_executable, input="")
        print("{} {} {}".format(path_to_cpp_executable, self.workingdir, self.parfile))
        # subprocess.run(path_to_cpp_executable, input=self.workingdir)
        subprocess.check_call([path_to_cpp_executable, self.workingdir, self.parfile])
    ''' --------- magnetar ---------- '''
    def _check_if_loaded_mag_obj(self):
        if (self.fpath_mag is None):
            raise IOError("self.fpath_mag is not set")
        if (self.mag_dfile is None):
            self.mag_dfile = h5py.File(self.fpath_mag)
    ''' --------- jet (only layers) -------------- '''
    # jet dynamics
    def _ckeck_if_loaded_j_dyn_obj(self):
        if (self.fpath_grb_dyn is None):
            raise IOError("self.fpath_grb_dyn is not set")
        if (self.grb_dyn_dfile is None):
            self.grb_dyn_dfile = h5py.File(self.fpath_grb_dyn)
    # jet spectrum
    def _check_if_loaded_j_spec(self):
        if (self.fpath_grb_spec is None):
            raise IOError("self.fpath_grb_spec is not set")
        if (self.grb_spec_dfile is None):
            self.grb_spec_dfile = h5py.File(self.fpath_grb_spec)
    # jet lightcurves
    def _check_if_loaded_jet_lc(self):
        if (self.fpath_grb_light_curve is None):
            raise IOError("self.fpath_grb_light_curve is not set")
        if (self.grb_lc_dfile is None):
            self.grb_lc_dfile = h5py.File(self.fpath_grb_light_curve)
    # jet skymaps
    def _check_if_loaded_jet_skymap(self):
        if (self.fpath_grb_sky_map is None):
            raise IOError("self.fpath_grb_sky_map is not set")
        if (self.grb_skymap_dfile is None):
            self.grb_skymap_dfile = h5py.File(self.fpath_grb_sky_map)

    ''' -------- ejecta (shells and layers) -------- '''

    # ejecta dynamics
    def _ckeck_if_loaded_ej_dyn_obj(self):
        if (self.fpath_kn_dyn is None):
            raise IOError("self.fpath_kn_dyn is not set")
        if (self.kn_dyn_dfile is None):
            self.kn_dyn_dfile = h5py.File(self.fpath_kn_dyn)
    # ejecta spectrum
    def _check_if_loaded_ej_spec(self):
        if (self.fpath_kn_spec is None):
            raise IOError("self.fpath_kn_spec is not set")
        if (self.kn_spec_dfile is None):
            self.kn_spec_dfile = h5py.File(self.fpath_kn_spec)
    # ejecta lightcurves
    def _check_if_loaded_ej_lc(self):
        if (self.fpath_kn_light_curve is None):
            raise IOError("self.fpath_kn_light_curve is not set")
        if (self.kn_lc_dfile is None):
            self.kn_lc_dfile = h5py.File(self.fpath_kn_light_curve)
    # ejecta skymaps
    def _check_if_loaded_ej_skymap(self):
        if (self.fpath_kn_sky_map is None):
            raise IOError("self.fpath_kn_sky_map is not set")
        if (self.kn_skymap_dfile is None):
            try:
                self.kn_skymap_dfile = h5py.File(self.fpath_kn_sky_map)
            except OSError:
                raise OSError("failed to load the file: \n {}".format(self.fpath_kn_sky_map))

def interp(xxs, yys, fluxes, x_grid, y_grid, method='linear'):
    return interpolate.griddata(np.vstack((xxs, yys)).T, fluxes, np.array((x_grid, y_grid)).T, method=method, fill_value=0.)

def compute_position_of_the_flux_centroid(all_x_arrs, all_y_arrs, all_z_arrs, d_l):
    rad2mas = 2.062648062470964e+08
    xcs = np.average(np.concatenate(all_x_arrs), weights=np.concatenate(all_z_arrs))
    xcs_m = xcs# * rad2mas / d_l
    ycs = np.average(np.concatenate(all_y_arrs), weights=np.concatenate(all_z_arrs))
    ycs_m = ycs# * rad2mas / d_l
    return(xcs_m, ycs_m) # in milli-arc-seconds

def lateral_distributions(grid_x, grid_y, image, collapse_axis='y'):
    # if collapse_axis == "y":
    #     return (grid_x[:, 0], latAvDist, latMaxDist)
    # else:
    #     return (grid_y[0, :], latAvDist, latMaxDist)
    if collapse_axis == 'x':
        points = image.shape[1]
        latAvDist = np.array([image[:, ii].mean() for ii in range(points)])
        # latAvDist2 = np.array([np.sum(image[:-1, ii]*np.diff(grid_y[:,ii]))/np.sum(np.diff(grid_y[:,ii])) for ii in range(points)])
        # latAvDist2 = np.array([np.sum(image[:-1, ii]*np.diff(np.abs(grid_y[:,ii])))/np.sum(np.diff(np.abs(grid_y[:,ii]))) for ii in range(points)])

        # latAvDist2 = np.zeros(points)
        # for ii in range(points):
        #     diff = grid_x[:, ii]
        #     dx = np.diff(diff)[0]
        #     arr = image[:, ii]
        #     latAvDist2[ii] = np.trapz(arr, diff, dx=dx)
        latAvDist2 = np.array( [ np.trapz(image[:, ii], grid_x[:, ii], dx=np.diff(grid_x[:, ii])[0] ) for ii in range(points) ] )

        # latAvDist2 = np.array([np.trapz(image[:, ii], grid_y[:,ii]) for ii in range(points)]) # , dx=np.abs(np.diff(grid_y[:,ii])[0])
        latMaxDist = np.array([image[:, ii].max() for ii in range(points)])

    elif collapse_axis == 'y':
        points = image.shape[0]
        latAvDist = np.array([image[ii, :].mean() for ii in range(points)])
        # latAvDist2 = np.array([np.sum(image[ii, :-1]*np.diff(grid_x[ii,:]))/np.sum(np.diff(grid_x[ii,:])) for ii in range(points)])
        # latAvDist2 = np.array([np.sum(image[ii, :-1]*np.diff(np.abs(grid_x[ii,:])))/np.sum(np.diff(np.abs(grid_x[ii,:]))) for ii in range(points)])
        latAvDist2 = np.array( [ np.trapz(image[ii, :], grid_y[ii, :], dx=np.diff(grid_y[ii, :])[0] ) for ii in range(points) ] )

        latMaxDist = np.array([image[ii, :].max() for ii in range(points)])

    else:
        raise KeyError()

    return (latAvDist, latAvDist2, latMaxDist)

def image_slice(fluxes, xxs, yys, position, nn=50, axis='y', inter='linear', fac=1):
    nn = complex(0, nn)
    if axis == 'y':
        grid_x, grid_y = np.mgrid[position:position:1j, yys.min():yys.max():nn]
    else:
        grid_x, grid_y = np.mgrid[xxs.min():xxs.max():nn, position:position:1j]

    slice = interpolate.gdd(np.array([xxs, yys]).T * fac, fluxes,
                            (grid_x * fac, grid_y * fac), method=inter, fill_value=0)

    return grid_x, grid_y, slice

class BPA_METHODS(PBA_BASE):
    '''
        Process output_uniform_grb files: load, extract for a specific way
    '''
    def __init__(self, workingdir, readparfileforpaths=True, parfile="parfile.par"):
        super().__init__(workingdir=workingdir,readparfileforpaths=readparfileforpaths,parfile=parfile)

    # magnetar
    def get_mag_obj(self):
        self._check_if_loaded_mag_obj()
        return self.mag_dfile


    # jet dynamics

    def get_jet_dyn_obj(self, ilayer=None):
        self._ckeck_if_loaded_j_dyn_obj()
        if (ilayer is None):
            return self.grb_dyn_dfile
        else:
            return self.grb_dyn_dfile["layer={}".format(ilayer)]

    def get_jet_dyn_arr(self, v_n, ilayer=0):
        self._ckeck_if_loaded_j_dyn_obj()
        dfile = self.grb_dyn_dfile
        layer = "layer={}".format(ilayer)
        if (not layer in list(dfile.keys())):
            raise NameError("Layer {} (aka '{}') is not in the jet dyn. file.\n Available: {}"
                            .format(ilayer, layer, dfile.keys()))
        if (not v_n in dfile[layer].keys()):
            raise NameError("v_n {} is not in the jet dyn. dfile[{}].keys() \n Avaialble: {}"
                            .format(v_n, layer, dfile[layer].keys()))
        return np.array(dfile[layer][v_n])

    # jet spectrum

    def get_jet_spec_obj(self):
        self._check_if_loaded_j_spec()
        return self.grb_spec_dfile

    def get_jet_spec_2d_arr(self, v_n="em", ilayer=0):
        dfile = self.get_jet_spec_obj()
        layer = "layer={}".format(ilayer)
        if (not layer in list(dfile.keys())):
            raise NameError("Layer {} (aka '{}') is not in the jet dyn. file.\n Available: {}"
                            .format(ilayer, layer, dfile.keys()))
        if (not v_n in dfile[layer].keys()):
            raise NameError("v_n {} is not in the jet dyn. dfile[{}].keys() \n Avaialble: {}"
                            .format(v_n, layer, dfile[layer].keys()))
        return np.array(dfile[layer][v_n])

    def get_jet_spec_times(self):
        dfile = self.get_jet_spec_obj()
        return np.array(dfile["times"])

    def get_jet_spec_freqs(self):
        dfile = self.get_jet_spec_obj()
        return np.array(dfile["freqs"])

    def get_jet_spec_1d_arr(self, freq=None, time=None, v_n="em", ilayer=0):
        arr = self.get_jet_spec_2d_arr(v_n=v_n, ilayer=ilayer)
        freqs = self.get_jet_spec_freqs()
        times = self.get_jet_spec_times()
        if (not freq is None):
            if (freq < freqs.min() or freq > freqs.max()):
                raise ValueError("Freq={} (jet comov cpec) is not in the limit of freqs[{}, {}]"
                                 .format(freq, freqs.min(), freqs.max()))
            if (not freq in freqs):
                idx = find_nearest_index(freqs, freq)
            else:
                idx = int(np.where(freqs == freq)[0])
            arr_ = arr[idx, :]
            assert len(arr_) == len(times)
            return arr_
        if (not times is None):
            if (time < times.min() or time > times.max()):
                raise ValueError("Time={} (jet comov cpec) is not in the limit of times[{}, {}]"
                                 .format(time, times.min(), times.max()))
            if (not time in times):
                idx = find_nearest_index(time, times)
            else:
                idx = int(np.where(time == times)[0][0])
            arr_ = arr[:, idx]
            assert len(arr_) == len(freqs)
            return arr_

    # jet lightcurves

    def get_jet_lc_obj(self):
        self._check_if_loaded_jet_lc()
        return self.grb_lc_dfile

    def get_jet_lc_times(self):
        dfile = self.get_jet_lc_obj()
        # print(dfile.keys())
        # print(dfile.attrs())

        return np.array(dfile["times"])

    def get_jet_lc_freqs(self):
        dfile = self.get_jet_lc_obj()
        return np.array(dfile["freqs"])

    def get_jet_lc_totalflux(self, freq):
        dfile = self.get_jet_lc_obj()
        nlayers = int(dfile.attrs["nlayers"])
        times = self.get_jet_lc_times()
        freqs = self.get_jet_lc_freqs()
        try:
            key = str("totalflux at freq={:.4e}".format(freq)).replace('.', ',')
            if not key in dfile.keys():
                raise NameError("Not found: {} among keys:{}".format(key, dfile.keys()))
        except NameError:
            key = str("totalflux at freq={:.4e}".format(freq))
            if not key in dfile.keys():
                raise NameError("Not found in lc_fname={}; key={}, among keys:{}"
                                .format(self.fpath_grb_light_curve,
                                        key, [key for key in dfile.keys() if key.__contains__("totalflux")]))
        except:
            raise NameError()
        return np.array(dfile[key])

    def get_jet_lc(self, freq=None, ilayer=None, v_n=None):
        dfile = self.get_jet_lc_obj()
        nlayers = int(dfile.attrs["nlayers"])
        times = self.get_jet_lc_times()
        freqs = self.get_jet_lc_freqs()
        if (v_n is None): v_n = list(dfile["layer=0"].keys())[0]  # Default key

        if (freq is None):
            # spectum
            if (ilayer is None):
                # fluxes2d = np.zeros((len(freqs), len(times)))
                # for ish in range(nshells):
                #     for il in range(nlayers):
                #         fluxes2d += np.array(dfile["shell={} layer={}".format(ish, il)][v_n])# [freq,time]
                # return fluxes2d
                fluxes2d = []
                for ifreq in freqs:
                    fluxes2d.append(self.get_jet_lc_totalflux(freq=ifreq))  # [freq,time]
                fluxes2d = np.reshape(np.array(fluxes2d), (len(freqs), len(times)))
                return fluxes2d
            elif (not ilayer is None):
                fluxes2d = np.array(dfile["layer={}".format(ilayer)][v_n])  # [freq,time]
                return fluxes2d
            else:
                raise NameError()
        else:
            # light curves
            if (not freq in self.get_jet_lc_freqs()):
                raise ValueError("freq:{} is not in ej_lc Given:{}".format(freq, self.get_jet_lc_freqs()))
            ifreq = find_nearest_index(self.get_jet_lc_freqs(), freq)
            if (ilayer is None):
                return self.get_jet_lc_totalflux(freq=freq)
            elif (not ilayer is None):
                fluxes1d = np.array(dfile["layer={}".format(ilayer)][v_n][ifreq])  # [freq,time]
                return fluxes1d
            else:
                raise NameError()

    def _alpha_jet(self, freq, freqm1, ilayer, v_n):
        values_i = self.get_jet_lc(freq=freq, ilayer=ilayer, v_n=v_n)
        values_im1 = self.get_jet_lc(freq=freqm1, ilayer=ilayer, v_n=v_n)
        ffreq_i = np.full(values_i.shape, freq)
        ffreq_im1 = np.full(values_im1.shape, freqm1)
        num = np.log10(values_i) - np.log10(values_im1)
        denum = np.log10(ffreq_i) - np.log10(ffreq_im1)
        values = 1 * num / denum
        return values

    def get_jet_lc_spec_idx(self, freq, ishell=None, ilayer=None):
        freqs = self.get_jet_lc_freqs()
        times = self.get_jet_lc_times()
        if ((freq is None) and (ilayer is None)):
            arr = []
            for ifreq in range(1, len(freqs)):
                arr.append(self._alpha_jet(freq=freqs[ifreq], freqm1=freqs[ifreq - 1], ilayer=ilayer, v_n=None))
            return np.reshape(np.array(arr), newshape=(len(freqs) - 1, len(times)))
        else:
            idx = find_nearest_index(freqs, freq)
            return self._alpha(freq=freqs[idx], freqm1=freqs[idx - 1], ishell=ishell, ilayer=ilayer, v_n=None)

    def get_jet_lc_temp_idx(self, freq, ishell, ilayer):
        freqs = self.get_jet_lc_freqs()
        t = self.get_jet_lc_times()
        Lnu = self.get_jet_lc(freq=freq, ilayer=ilayer, v_n=None)
        val = np.zeros_like(Lnu)
        val[1:-1] = np.log10(Lnu[2:] / Lnu[:-2]) / np.log10(t[2:] / t[:-2])
        val[0] = np.log10(Lnu[1] / Lnu[0]) / np.log10(t[1] / t[0])
        val[-1] = np.log10(Lnu[-1] / Lnu[-2]) / np.log10(t[-1] / t[-2])
        return val

    # jet skymaps

    def get_jet_skymap_obj(self):
        self._check_if_loaded_jet_skymap()
        return self.grb_skymap_dfile

    def get_jet_skymap_times(self):
        self._check_if_loaded_jet_skymap()
        return np.array(self.grb_skymap_dfile["times"])

    def get_jet_skymap_freqs(self):
        self._check_if_loaded_jet_skymap()
        return np.array(self.grb_skymap_dfile["freqs"])

    def get_jet_skymap_cm(self, all_xrs, all_yrs, all_zz):
        dfile = self.get_jet_skymap_obj()
        xc_m, yc_m = compute_position_of_the_flux_centroid(all_xrs, all_yrs, all_zz, float(dfile.attrs["d_l"]))
        return (xc_m, yc_m)

    # def get_jet_skymap(self, time=None, freq=None, verbose=False, remove_mu=False, int_or_hist="int"):
    #     ish = 0
    #     nx = 2000
    #     ny = 1000
    #     min_mu = 1e-5  # limit for mu that if too small blow up the intensity for some reason...
    #     min_mu_frac = 0.85
    #     # limit_data = False
    #     times = self.get_jet_skymap_times()
    #     freqs = self.get_jet_skymap_freqs()
    #     dfile = self.get_jet_skymap_obj()
    #     # nshells = int(dfile.attrs["nshells"])
    #     if ((not time is None) and (not time in times)):
    #         raise ValueError("time={} is not in the list for skypams={}".format(time, times))
    #     if ((not freq is None) and (not freq in freqs)):
    #         raise ValueError("freq={} is not in the list for skypams={}".format(freq, freqs))
    #     if ((not time is None) and (not freq is None)):
    #         print(dfile.keys())
    #         key = "time={:.4e} freq={:.4e}".format(time, freq)
    #         try:
    #             ddfile = dfile[key]
    #         except KeyError:
    #             try:
    #                 ddfile = dfile[key.replace(".", ",")]
    #             except KeyError:
    #                     raise KeyError("Failed to find ky for jet skymap : {} Avaialble:\n{}".format(key, dfile.keys()))
    #         r_i = np.array(ddfile["r"][ish])
    #         mu_i = np.array(ddfile["mu"][ish])
    #         xrs_i = np.array(ddfile["xrs"][ish])
    #         yrs_i = np.array(ddfile["yrs"][ish])
    #         int_i = np.array(ddfile["intensity"][ish])  # * dfile.attrs["d_L"] ** 2
    #         gam_i = np.array(ddfile["gamma"][ish])
    #         B_i = np.array(ddfile["B"][ish])
    #         tb_i = np.array(ddfile["tburst"][ish])
    #         theta_i = np.array(ddfile["theta"][ish])
    #         theta_ij = np.array(ddfile["theta_j"][ish])
    #         theta0_i = np.array(ddfile["theta0"][ish])
    #         phi_i = np.array(ddfile["phi"][ish])
    #
    #         # remove nans
    #         # int_i = _int_i[(np.isfinite(_xrs_i))&(_xrs_i > 1.)&(_yrs_i > 1.)]
    #         # yrs_i = _yrs_i[(np.isfinite(_xrs_i))&(_xrs_i > 1.)&(_yrs_i > 1.)]
    #         # xrs_i = _xrs_i[(np.isfinite(_xrs_i))&(_xrs_i > 1.)&(_yrs_i > 1.)]
    #
    #
    #         # print(np.arccos(mu_i[np.abs(mu_i) < 1e-4]))
    #         print("Min={} max={}".format(np.min(np.abs(mu_i)), np.max(np.abs(mu_i))))
    #
    #         # plt.plot(mu_i, theta0_i, marker='.', ls='none')
    #         # plt.semilogy(mu_i, int_i, marker='.', ls='none')
    #         # plt.show()
    #
    #         if verbose:
    #
    #
    #             _mu = np.arccos(mu_i)
    #             x = np.exp(1 + abs(mu_i).max() / -abs(mu_i))
    #             idx1 = abs(mu_i) < min_mu
    #             int_i2 = np.copy(int_i)  * abs(mu_i) #/ (1. - 2.*np.sin(_mu/2)) #/ mu_i# (1.-mu_i*get_beta(gam_i))/(mu_i-get_beta(gam_i))#np.abs(mu_i) / mu_i**2
    #             #fac = np.exp(1 / (1 / abs(mu_i))-1)
    #             # int_i2[abs(mu_i) < 1e-3] *= 1e-3
    #             #int_i2[abs(mu_i) < 1e-1] *= abs(mu_i)[abs(mu_i) < 1e-1]
    #             # int_i2[abs(mu_i) < 1e-5] *= 1e-5
    #             # int_i2[abs(mu_i) < 1e-6] *= 1e-6
    #             #idx1 = phi_i == 0.
    #             idx1 = abs(mu_i) < min_mu  # TODO this is overritten!
    #             # int_i[idx1] =
    #             # print("len(gam[idx1])={}".format(gam_i[idx1]))
    #             # idx = int_i > int_i.max() * min_mu_frac
    #             idx = idx1
    #             if verbose:
    #                 if (len(mu_i[idx1]) != len(mu_i[idx])):
    #                     print('\t', mu_i[idx1])
    #                     print('\t', mu_i[idx])
    #                     # exit(1)
    #             if verbose:
    #                 # if len(idx) < 2:
    #                 #     print("No excess = {}".format(ii))
    #                 print("len(gam[idx])={}".format(gam_i[idx]))
    #                 print("Gamma={}".format(gam_i[idx]))
    #                 print("beta={}".format(get_beta(gam_i[idx])))
    #                 print("mu_i={}".format(mu_i[idx]))
    #                 print("r_i={}".format(r_i[idx]))
    #                 print("B={}".format(B_i[idx]))
    #                 print("tb_i={}".format(tb_i[idx]))
    #             print(np.abs(mu_i).min())
    #             print(theta_i[np.abs(mu_i) < 1e-5])
    #             print(phi_i[np.abs(mu_i) < 1e-5])
    #
    #             fig, axes = plt.subplots(figsize=(4.6+2, 3.2), ncols=4, nrows=1)
    #             ax = axes[0]
    #             ax.plot(mu_i, int_i / int_i.max(), '.')
    #             ax.plot(mu_i, int_i2 / int_i.max(), 'o', zorder=-1)
    #             if len(mu_i[idx]) > 0: ax.plot(mu_i[idx], int_i[idx] / int_i[idx].max(), 'x', color='red')
    #             # plt.semilogy(mu_i,int_i2/int_i2.max(), 'o', color='red')
    #             ax.set_xlabel(r"$\mu=\sin(\theta_{\rm obs}) \sin(\theta) \sin(\phi) + \cos(\theta_{\rm obs}) \cos(\theta)$")
    #             ax.set_ylabel(r"$I/I_{\rm max}$")
    #             ax.set_xscale("linear")
    #             ax.set_yscale("log")
    #
    #             ax = axes[1]
    #             ax.plot(np.sin(phi_i), int_i / int_i.max(), '.')
    #             ax.plot(np.sin(phi_i), int_i2 / int_i.max(), 'o', zorder=-1)
    #             if len(mu_i[idx]) > 0: ax.plot(np.sin(phi_i)[idx], int_i[idx] / int_i[idx].max(), 'x', color='red')
    #             # plt.semilogy(mu_i,int_i2/int_i2.max(), 'o', color='red')
    #             ax.set_xlabel(r"$\phi$")
    #             ax.set_ylabel(r"$I/I_{\rm max}$")
    #             ax.set_xscale("linear")
    #             ax.set_yscale("log")
    #
    #             ax = axes[2]
    #             ax.plot(np.sin(theta_i), int_i / int_i.max(), '.')
    #             ax.plot(np.sin(theta_i), int_i2 / int_i2.max(), 'o')
    #             if len(mu_i[idx]) > 0: ax.plot(np.sin(theta_i)[idx], int_i[idx] / int_i[idx].max(), 'x', color='red')
    #             # plt.semilogy(mu_i,int_i2/int_i2.max(), 'o', color='red')
    #             ax.set_xlabel(r"$\sin\theta$")
    #             ax.set_ylabel(r"$I/I_{\rm max}$")
    #             ax.set_xscale("linear")
    #             ax.set_yscale("log")
    #
    #             ax = axes[3]
    #             ax.plot(np.cos(theta_i), int_i / int_i.max(), '.')
    #             ax.plot(np.cos(theta_i), int_i2 / int_i.max(), 'o')
    #             if len(mu_i[idx]) > 0: ax.plot(np.cos(theta_i)[idx], int_i[idx] / int_i[idx].max(), 'x', color='red')
    #             # plt.semilogy(mu_i,int_i2/int_i2.max(), 'o', color='red')
    #             ax.set_xlabel(r"$\cos\theta$")
    #             ax.set_ylabel(r"$I/I_{\rm max}$")
    #             ax.set_xscale("linear")
    #             ax.set_yscale("log")
    #
    #             plt.tight_layout()
    #             plt.show()
    #
    #         # plt.plot(theta0_i, mu_i, marker='.', ls='none')
    #         # plt.show()
    #
    #
    #         if remove_mu:
    #             int_i *= np.abs(mu_i)
    #
    #         # int_i[~np.isfinite(int_i)] = 0.0
    #         # int_i[idx] = 0.
    #         # int_i = int_i2
    #
    #         if int_or_hist == "int":
    #             nx = np.complex(0, nx)
    #             ny = np.complex(0, ny)
    #             grid_x, grid_y = np.mgrid[xrs_i.min()*1.2:xrs_i.max()*1.2:nx, yrs_i.min()*1.2:yrs_i.max()*1.2:ny]
    #             i_zz = interp(xrs_i, yrs_i, int_i, grid_x, grid_y, 'linear')
    #         elif int_or_hist == "both":
    #
    #             nx = 800
    #             ny = 400
    #
    #             nx = np.complex(0, nx)
    #             ny = np.complex(0, ny)
    #             grid_x, grid_y = np.mgrid[xrs_i.min():xrs_i.max():nx,
    #                              yrs_i.min():yrs_i.max():ny]
    #             i_zz = interp(xrs_i, yrs_i, int_i, grid_x, grid_y, 'linear')
    #
    #
    #             nx = np.complex(0, nx + 1)
    #             ny = np.complex(0, ny + 1)
    #             edges_x = np.mgrid[xrs_i.min() * 1.2:xrs_i.max() * 1.2:nx]
    #             edges_y = np.mgrid[yrs_i.min() * 1.2:yrs_i.max() * 1.2:ny]
    #             # grid_x, grid_y = np.mgrid[xrs_i.min() * 1.2:xrs_i.max() * 1.2:nx,
    #             #                  yrs_i.min() * 1.2:yrs_i.max() * 1.2:ny]
    #             i_zz, _ = np.histogramdd(tuple([grid_x.flatten(), grid_y.flatten()]), bins=tuple([edges_x, edges_y]), weights=i_zz.T.flatten())
    #             grid_x = 0.5 * (edges_x[1:] + edges_x[:-1])
    #             grid_y = 0.5 * (edges_y[1:] + edges_y[:-1])
    #             print(i_zz.shape)
    #         else:
    #
    #
    #             nx = 200
    #             ny = 100
    #             nx = np.complex(0, nx+1)
    #             ny = np.complex(0, ny+1)
    #             edges_x = np.mgrid[xrs_i.min() * 1.2:xrs_i.max() * 1.2:nx]
    #             edges_y = np.mgrid[yrs_i.min() * 1.2:yrs_i.max() * 1.2:ny]
    #             # grid_x, grid_y = np.mgrid[xrs_i.min() * 1.2:xrs_i.max() * 1.2:nx,
    #             #                  yrs_i.min() * 1.2:yrs_i.max() * 1.2:ny]
    #             i_zz, _ = np.histogramdd(tuple([xrs_i, yrs_i]), bins=tuple([edges_x, edges_y]), weights=int_i)
    #             grid_x = 0.5 * (edges_x[1:] + edges_x[:-1])
    #             grid_y = 0.5 * (edges_y[1:] + edges_y[:-1])
    #             print(i_zz.shape)
    #
    #
    #         # i_zz = ndimage.uniform_filter(i_zz, )
    #         # i_zz = ndimage.filters.gaussian_filter(i_zz, [10,10], mode='reflect')
    #
    #         return (grid_x, grid_y, i_zz, xrs_i, yrs_i, int_i)
    #
    #         #
    #         # else:
    #         #     # assess what is the combined grid extend (for all images)
    #         #     xmin_neg, xmin_pos, xmax = [], [], []
    #         #     ymin_neg, ymin_pos, ymax = [], [], []
    #         #     i_min, i_max = [], []
    #         #     i_shells = []
    #         #     for ish in range(nshells):
    #         #         xrs_i = np.array(ddfile["xrs"][ish])
    #         #         yrs_i = np.array(ddfile["yrs"][ish])
    #         #         int_i = np.array(ddfile["intensity"][ish])
    #         #         if (np.sum(int_i) == 0):
    #         #             continue
    #         #         i_shells.append(ish)
    #         #         xmin_neg.append(xrs_i.min())
    #         #         xmin_pos.append(xrs_i[xrs_i > 0].min())
    #         #         xmax.append(xrs_i.max())
    #         #         ymin_neg.append(yrs_i.min())
    #         #         ymin_pos.append(yrs_i[yrs_i > 0].min())
    #         #         ymax.append(yrs_i.max())
    #         #         i_min.append(int_i.min())
    #         #         i_max.append(int_i.max())
    #         #     xmin_neg = np.array(xmin_neg)
    #         #     xmin_pos = np.array(xmin_pos)
    #         #     xmax = np.array(xmax)
    #         #     ymin_neg = np.array(ymin_neg)
    #         #     ymin_pos = np.array(ymin_pos)
    #         #     ymax = np.array(ymax)
    #         #     i_min = np.array(i_min)
    #         #     i_max = np.array(i_max)
    #         #     if verbose:
    #         #         for i in range(len(i_shells)):
    #         #             print(
    #         #                 "\tshell={} xmin_neg={:.2e} xmin_pos={:.2e} xmax={:.2e}".format(i, xmin_neg[i], xmin_pos[i],
    #         #                                                                                 xmax[i]))
    #         #             print(
    #         #                 "\t         ymin_neg={:.2e} ymin_pos={:.2e} ymax={:.2e}".format(i, ymin_neg[i], ymin_pos[i],
    #         #                                                                                 ymax[i]))
    #         #         print("\toverall X min(xmin_pos) = {:.2e} max(xmax) = {:.2e}".format(xmin_pos.min(), xmax.max()))
    #         #         print("\toverall Y min(ymin_pos) = {:.2e} max(ymax) = {:.2e}".format(ymin_pos.min(), ymax.max()))
    #         #         print("\toverall I min(i_min)    = {:.2e} max(ymax) = {:.2e}".format(i_min.min(), i_max.max()))
    #         #     x_grid = np.concatenate((-1. * np.logspace(np.log10(xmin_pos.min()), np.log10(xmax.max()), nx)[::-1],
    #         #                              np.logspace(np.log10(xmin_pos.min()), np.log10(xmax.max()), nx)))
    #         #     y_grid = np.concatenate((-1. * np.logspace(np.log10(ymin_pos.min()), np.log10(ymax.max()), ny)[::-1],
    #         #                              np.logspace(np.log10(ymin_pos.min()), np.log10(ymax.max()), ny)))
    #         #     if verbose:
    #         #         print("\tGrid xmin_neg={:.2e} xmin_pos={:.2e} xmax={:.2e}".format(x_grid.min(),
    #         #                                                                           x_grid[x_grid > 0].min(),
    #         #                                                                           x_grid.max()))
    #         #         print("\tGrid ymin_neg={:.2e} ymin_pos={:.2e} ymax={:.2e}".format(y_grid.min(),
    #         #                                                                           y_grid[y_grid > 0].min(),
    #         #                                                                           y_grid.max()))
    #         #     xx_grid, yy_grid = np.meshgrid(x_grid, y_grid)
    #         #     zz = np.zeros_like((xx_grid))
    #         #     # interpolate onto the grid that covers all images
    #         #     all_xrs, all_yrs, all_zz = [], [], []
    #         #     for ii, ish in enumerate(i_shells):  # range(len(i_shells))
    #         #         if verbose: print("Pocessing: shell={} [{}/{}]".format(ish, ii, len(i_shells)))
    #         #         r_i = np.array(ddfile["r"][ish])
    #         #         mu_i = np.array(ddfile["mu"][ish])
    #         #         xrs_i = np.array(ddfile["xrs"][ish])
    #         #         yrs_i = np.array(ddfile["yrs"][ish])
    #         #         int_i = np.array(ddfile["intensity"][ish])  # * dfile.attrs["d_L"] ** 2
    #         #         gam_i = np.array(ddfile["gamma"][ish])
    #         #         B_i = np.array(ddfile["B"][ish])
    #         #         tb_i = np.array(ddfile["tburst"][ish])
    #         #         theta_i = np.array(ddfile["theta"][ish])
    #         #         if limit_data:
    #         #             idx1 = abs(mu_i) < min_mu  # TODO this is overritten!
    #         #             # int_i[idx1] =
    #         #             # print("len(gam[idx1])={}".format(gam_i[idx1]))
    #         #             idx = int_i > int_i.max() * min_mu_frac
    #         #             idx = idx1
    #         #             if verbose:
    #         #                 if (len(mu_i[idx1]) != len(mu_i[idx])):
    #         #                     print('\t', mu_i[idx1])
    #         #                     print('\t', mu_i[idx])
    #         #                     # exit(1)
    #         #             if verbose:
    #         #                 if len(idx) < 2:
    #         #                     print("No excess = {}".format(ii))
    #         #                 print("len(gam[idx])={}".format(gam_i[idx]))
    #         #                 print("Gamma={}".format(gam_i[idx]))
    #         #                 print("beta={}".format(get_beta(gam_i[idx])))
    #         #                 print("mu_i={}".format(mu_i[idx]))
    #         #                 print("r_i={}".format(r_i[idx]))
    #         #                 print("B={}".format(B_i[idx]))
    #         #                 print("tb_i={}".format(tb_i[idx]))
    #         #             # print(gam_i[idx], Beta(gam_i[idx]))
    #         #
    #         #             # int_i[mu_i < 1e-3] = 0.
    #         #             # int_i2 = int_i * abs(mu_i)# * np.sqrt(1 - mu_i**2)#/ ((1-mu_i)*(1+mu_i))
    #         #             # plt.figure(figsize=(4.6,3.2))
    #         #             # plt.semilogy(mu_i,int_i/int_i.max(), '.')
    #         #             # plt.semilogy(mu_i[idx],int_i[idx]/int_i[idx].max(), 'x', color='red')
    #         #             # plt.semilogy(mu_i,int_i2/int_i2.max(), 'o', color='red')
    #         #             # plt.xlabel(r"$\mu=\sin(\theta_{\rm obs}) \sin(\theta) \sin(\phi) + \cos(\theta_{\rm obs}) \cos(\theta)$")
    #         #             # plt.ylabel(r"$I/I_{\rm max}$")
    #         #             # plt.tight_layout()
    #         #             # plt.show()
    #         #             int_i[idx] = 0.
    #         #         # int_i /= abs(mu_i)[abs(mu_i)>0.01]
    #         #         i_zz = interp(xrs_i, yrs_i, int_i, xx_grid, yy_grid, 'linear')
    #         #         zz += i_zz
    #         #         all_xrs.append(xrs_i)
    #         #         all_yrs.append(yrs_i)
    #         #         all_zz.append(int_i)
    #         #     if verbose:
    #         #         print("\tAfter interpolation: I [{:.2e}, {:.2e}] Sum = {:.2e} [Total expected {:.2e}]"
    #         #               .format(zz.min(), zz.max(), np.sum(zz), float(
    #         #             np.array(dfile["totalflux at freq={:.4e}".format(freq)])[find_nearest_index(times, time)])))
    #         #     return (xx_grid, yy_grid, zz, all_xrs, all_yrs, all_zz)
    #     else:
    #         raise NotImplementedError("Not implemented")
    def get_jet_skymap(self, time=None, freq=None, verbose=False, remove_mu=False, plots=False, renormalize=True):
        ish = 0
        # nx = 1000
        # ny = 500
        min_mu = 1e-5  # limit for mu that if too small blow up the intensity for some reason...
        # min_mu_frac = 0.85
        # limit_data = False
        d_l = float(self.get_jet_skymap_obj().attrs["d_l"])
        # theta_obs = float(self.get_jet_skymap_obj().attrs["thetaObs"])
        times = self.get_jet_skymap_times()
        freqs = self.get_jet_skymap_freqs()
        dfile = self.get_jet_skymap_obj()
        # nshells = int(dfile.attrs["nshells"])
        if ((not time is None) and (not time in times)):
            raise ValueError("time={} ({} days) is not in the list for skypams={} ({} days)"
                             .format(time, time/cgs.day, times, times/cgs.day))
        if ((not freq is None) and (not freq in freqs)):
            raise ValueError("freq={} is not in the list for skypams={}".format(freq, freqs))
        if ((not time is None) and (not freq is None)):
            # print(dfile.keys())
            key = "time={:.4e} freq={:.4e}".format(time, freq)
            try:
                ddfile = dfile[key]
            except KeyError:
                try:
                    ddfile = dfile[key.replace(".", ",")]
                except KeyError:
                    raise KeyError("Failed to find ky for jet skymap : {} Avaialble:\n{}".format(key, dfile.keys()))
            # r_i = np.array(ddfile["r"][ish])
            mu_i = np.array(ddfile["mu"][ish])
            xrs_i = np.array(ddfile["xrs"][ish]) * cgs.rad2mas / d_l  # m -> mas
            yrs_i = np.array(ddfile["yrs"][ish]) * cgs.rad2mas / d_l  # m -> mas
            int_i = np.array(ddfile["intensity"][ish]) * (d_l ** 2 / cgs.rad2mas ** 2)  # -> mJy / mas^2
            # gam_i = np.array(ddfile["gamma"][ish])
            # B_i = np.array(ddfile["B"][ish])
            # tb_i = np.array(ddfile["tburst"][ish])
            # theta_i = np.array(ddfile["theta"][ish])
            # theta_ij = np.array(ddfile["theta_j"][ish])
            # theta0_i = np.array(ddfile["theta0"][ish])
            # phi_i = np.array(ddfile["phi"][ish])
            # tt_i = np.array(ddfile["tt"][ish])

            # xrs =  -1. * np.cos(theta_obs) * np.sin(theta_i) * np.sin(phi_i) + np.sin(theta_obs) * np.cos(theta_i)
            # yrs = np.sin(theta_i) * np.cos(phi_i)

            # xrs_i = xrs * r_i * cgs.rad2mas / d_l
            # yrs_i = yrs * r_i * cgs.rad2mas / d_l

            # xrs_i[np.isnan(int_i)] = 0.
            # yrs_i[np.isnan(int_i)] = 0.
            # int_i[np.isnan(int_i)] = 0.
            # int_i[phi_i < np.pi ] = 0.
            # yrs_i[phi_i < np.pi] = 0.
            # xrs_i[phi_i < np.pi] = 0.

            if plots:
                layers = int(self.get_jet_dyn_obj().attrs["nlayers"])
                plt.semilogx(self.get_jet_dyn_arr("tt", ilayer=0) / cgs.day, self.get_jet_dyn_arr("ctheta", ilayer=0),
                             label="nl={}".format(0))
                plt.semilogx(self.get_jet_dyn_arr("tt", ilayer=layers - 1) / cgs.day,
                             self.get_jet_dyn_arr("ctheta", ilayer=layers - 1), label="nl={}".format(layers - 1))
                plt.semilogx(self.get_jet_dyn_arr("tt", ilayer=layers - int(layers / 2)) / cgs.day,
                             self.get_jet_dyn_arr("ctheta", ilayer=layers - int(layers / 2)),
                             label="nl={}".format(layers - int(layers / 2)))
                plt.xlabel("Time [days]")
                plt.grid()
                plt.legend()
                plt.show()

                fig = plt.figure()
                cmap = cm.get_cmap('viridis')
                my_norm = LogNorm(int_i.max() * 1e-2, int_i.max())
                ax = fig.add_subplot(projection='3d')
                # ax.scatter(xrs_i.flatten(), yrs_i.flatten(), ((theta_i-theta0_i)*180/np.pi).flatten(),  c=cmap(my_norm(int_i.flatten())))
                ax.scatter(xrs_i.flatten(), yrs_i.flatten(), np.log10(r_i).flatten(), c=cmap(my_norm(int_i.flatten())))
                # ax.scatter(xrs_i.flatten(), yrs_i.flatten(),mu_i.flatten(),  c=cmap(my_norm(int_i.flatten())))
                ax.set_xlabel('X Label')
                ax.set_ylabel('Y Label')
                ax.set_zlabel('I Label')
                plt.show()

            # remove nans
            # int_i = _int_i[(np.isfinite(_xrs_i))&(_xrs_i > 1.)&(_yrs_i > 1.)]
            # yrs_i = _yrs_i[(np.isfinite(_xrs_i))&(_xrs_i > 1.)&(_yrs_i > 1.)]
            # xrs_i = _xrs_i[(np.isfinite(_xrs_i))&(_xrs_i > 1.)&(_yrs_i > 1.)]

            # print(np.arccos(mu_i[np.abs(mu_i) < 1e-4]))
            # print("Jet skymap min={} max={}".format(np.min(np.abs(mu_i)), np.max(np.abs(mu_i))))

            # plt.plot(mu_i, theta0_i, marker='.', ls='none')
            # plt.semilogy(mu_i, int_i, marker='.', ls='none')
            # plt.show()

            if verbose:

                _mu = np.arccos(mu_i)
                x = np.exp(1 + abs(mu_i).max() / -abs(mu_i))
                idx1 = abs(mu_i) < min_mu
                int_i2 = np.copy(int_i) * abs(
                    mu_i)  # / (1. - 2.*np.sin(_mu/2)) #/ mu_i# (1.-mu_i*get_beta(gam_i))/(mu_i-get_beta(gam_i))#np.abs(mu_i) / mu_i**2
                # fac = np.exp(1 / (1 / abs(mu_i))-1)
                # int_i2[abs(mu_i) < 1e-3] *= 1e-3
                # int_i2[abs(mu_i) < 1e-1] *= abs(mu_i)[abs(mu_i) < 1e-1]
                # int_i2[abs(mu_i) < 1e-5] *= 1e-5
                # int_i2[abs(mu_i) < 1e-6] *= 1e-6
                # idx1 = phi_i == 0.
                idx1 = abs(mu_i) < min_mu  # TODO this is overritten!
                # int_i[idx1] =
                # print("len(gam[idx1])={}".format(gam_i[idx1]))
                # idx = int_i > int_i.max() * min_mu_frac
                idx = idx1
                if verbose:
                    if (len(mu_i[idx1]) != len(mu_i[idx])):
                        print('\t', mu_i[idx1])
                        print('\t', mu_i[idx])
                        # exit(1)
                if verbose:
                    # if len(idx) < 2:
                    #     print("No excess = {}".format(ii))
                    print("len(gam[idx])={}".format(gam_i[idx]))
                    print("Gamma={}".format(gam_i[idx]))
                    print("beta={}".format(get_beta(gam_i[idx])))
                    print("mu_i={}".format(mu_i[idx]))
                    print("r_i={}".format(r_i[idx]))
                    print("B={}".format(B_i[idx]))
                    print("tb_i={}".format(tb_i[idx]))
                print(np.abs(mu_i).min())
                print(theta_i[np.abs(mu_i) < 1e-5])
                print(phi_i[np.abs(mu_i) < 1e-5])

                fig, axes = plt.subplots(figsize=(4.6 + 2, 3.2), ncols=4, nrows=1)
                ax = axes[0]
                ax.plot(mu_i, int_i / int_i.max(), '.')
                ax.plot(mu_i, int_i2 / int_i.max(), 'o', zorder=-1)
                if len(mu_i[idx]) > 0: ax.plot(mu_i[idx], int_i[idx] / int_i[idx].max(), 'x', color='red')
                # plt.semilogy(mu_i,int_i2/int_i2.max(), 'o', color='red')
                ax.set_xlabel(
                    r"$\mu=\sin(\theta_{\rm obs}) \sin(\theta) \sin(\phi) + \cos(\theta_{\rm obs}) \cos(\theta)$")
                ax.set_ylabel(r"$I/I_{\rm max}$")
                ax.set_xscale("linear")
                ax.set_yscale("log")

                ax = axes[1]
                ax.plot(np.sin(phi_i), int_i / int_i.max(), '.')
                ax.plot(np.sin(phi_i), int_i2 / int_i.max(), 'o', zorder=-1)
                if len(mu_i[idx]) > 0: ax.plot(np.sin(phi_i)[idx], int_i[idx] / int_i[idx].max(), 'x', color='red')
                # plt.semilogy(mu_i,int_i2/int_i2.max(), 'o', color='red')
                ax.set_xlabel(r"$\phi$")
                ax.set_ylabel(r"$I/I_{\rm max}$")
                ax.set_xscale("linear")
                ax.set_yscale("log")

                ax = axes[2]
                ax.plot(np.sin(theta_i), int_i / int_i.max(), '.')
                ax.plot(np.sin(theta_i), int_i2 / int_i2.max(), 'o')
                if len(mu_i[idx]) > 0: ax.plot(np.sin(theta_i)[idx], int_i[idx] / int_i[idx].max(), 'x', color='red')
                # plt.semilogy(mu_i,int_i2/int_i2.max(), 'o', color='red')
                ax.set_xlabel(r"$\sin\theta$")
                ax.set_ylabel(r"$I/I_{\rm max}$")
                ax.set_xscale("linear")
                ax.set_yscale("log")

                ax = axes[3]
                ax.plot(np.cos(theta_i), int_i / int_i.max(), '.')
                ax.plot(np.cos(theta_i), int_i2 / int_i.max(), 'o')
                if len(mu_i[idx]) > 0: ax.plot(np.cos(theta_i)[idx], int_i[idx] / int_i[idx].max(), 'x', color='red')
                # plt.semilogy(mu_i,int_i2/int_i2.max(), 'o', color='red')
                ax.set_xlabel(r"$\cos\theta$")
                ax.set_ylabel(r"$I/I_{\rm max}$")
                ax.set_xscale("linear")
                ax.set_yscale("log")

                plt.tight_layout()
                plt.show()

            # plt.plot(theta0_i, mu_i, marker='.', ls='none')
            # plt.show()

            if remove_mu:
                print("Removing 'mu' from jet skymap")
                int_i *= np.abs(mu_i)

            # int_i[~np.isfinite(int_i)] = 0.0
            # int_i[idx] = 0.
            # int_i = int_i2

            if renormalize:
                print("Renormalizing jet skymap to total flux")
                # Problem: 'all_fluxes_jet' decreased for increasing nlayers; Solution: Below
                all_fluxes_jet = np.array(int_i)
                fnus = self.get_jet_skymap_totfluxes(freq=freq)
                fnu = fnus[find_nearest_index(self.get_jet_skymap_times(), time)]
                delta_x = xrs_i.max() - xrs_i.min()
                delta_y = yrs_i.max() - yrs_i.min()
                dfnu = fnu / (delta_x * delta_y)
                if verbose:
                    print("Jet")
                    print("\tfnu = {:.2e} ".format(fnu))
                    print("\tall_x = [{:.2e}, {:.2e}]".format(np.array(xrs_i).min(), np.array(xrs_i).max()))
                    print("\tall_y = [{:.2e}, {:.2e}]".format(np.array(yrs_i).min(), np.array(yrs_i).max()))
                    print("\tDelta_x = {:.2f}, Delta_y = {:.2f}]".format(delta_x, delta_y))
                    print("\tFnu/mas^2 = {:.2e} mJy/mas^2".format(dfnu))
                int_i = (all_fluxes_jet / all_fluxes_jet.max()) * dfnu

            return (xrs_i, yrs_i, int_i)

            #
            # else:
            #     # assess what is the combined grid extend (for all images)
            #     xmin_neg, xmin_pos, xmax = [], [], []
            #     ymin_neg, ymin_pos, ymax = [], [], []
            #     i_min, i_max = [], []
            #     i_shells = []
            #     for ish in range(nshells):
            #         xrs_i = np.array(ddfile["xrs"][ish])
            #         yrs_i = np.array(ddfile["yrs"][ish])
            #         int_i = np.array(ddfile["intensity"][ish])
            #         if (np.sum(int_i) == 0):
            #             continue
            #         i_shells.append(ish)
            #         xmin_neg.append(xrs_i.min())
            #         xmin_pos.append(xrs_i[xrs_i > 0].min())
            #         xmax.append(xrs_i.max())
            #         ymin_neg.append(yrs_i.min())
            #         ymin_pos.append(yrs_i[yrs_i > 0].min())
            #         ymax.append(yrs_i.max())
            #         i_min.append(int_i.min())
            #         i_max.append(int_i.max())
            #     xmin_neg = np.array(xmin_neg)
            #     xmin_pos = np.array(xmin_pos)
            #     xmax = np.array(xmax)
            #     ymin_neg = np.array(ymin_neg)
            #     ymin_pos = np.array(ymin_pos)
            #     ymax = np.array(ymax)
            #     i_min = np.array(i_min)
            #     i_max = np.array(i_max)
            #     if verbose:
            #         for i in range(len(i_shells)):
            #             print(
            #                 "\tshell={} xmin_neg={:.2e} xmin_pos={:.2e} xmax={:.2e}".format(i, xmin_neg[i], xmin_pos[i],
            #                                                                                 xmax[i]))
            #             print(
            #                 "\t         ymin_neg={:.2e} ymin_pos={:.2e} ymax={:.2e}".format(i, ymin_neg[i], ymin_pos[i],
            #                                                                                 ymax[i]))
            #         print("\toverall X min(xmin_pos) = {:.2e} max(xmax) = {:.2e}".format(xmin_pos.min(), xmax.max()))
            #         print("\toverall Y min(ymin_pos) = {:.2e} max(ymax) = {:.2e}".format(ymin_pos.min(), ymax.max()))
            #         print("\toverall I min(i_min)    = {:.2e} max(ymax) = {:.2e}".format(i_min.min(), i_max.max()))
            #     x_grid = np.concatenate((-1. * np.logspace(np.log10(xmin_pos.min()), np.log10(xmax.max()), nx)[::-1],
            #                              np.logspace(np.log10(xmin_pos.min()), np.log10(xmax.max()), nx)))
            #     y_grid = np.concatenate((-1. * np.logspace(np.log10(ymin_pos.min()), np.log10(ymax.max()), ny)[::-1],
            #                              np.logspace(np.log10(ymin_pos.min()), np.log10(ymax.max()), ny)))
            #     if verbose:
            #         print("\tGrid xmin_neg={:.2e} xmin_pos={:.2e} xmax={:.2e}".format(x_grid.min(),
            #                                                                           x_grid[x_grid > 0].min(),
            #                                                                           x_grid.max()))
            #         print("\tGrid ymin_neg={:.2e} ymin_pos={:.2e} ymax={:.2e}".format(y_grid.min(),
            #                                                                           y_grid[y_grid > 0].min(),
            #                                                                           y_grid.max()))
            #     xx_grid, yy_grid = np.meshgrid(x_grid, y_grid)
            #     zz = np.zeros_like((xx_grid))
            #     # interpolate onto the grid that covers all images
            #     all_xrs, all_yrs, all_zz = [], [], []
            #     for ii, ish in enumerate(i_shells):  # range(len(i_shells))
            #         if verbose: print("Pocessing: shell={} [{}/{}]".format(ish, ii, len(i_shells)))
            #         r_i = np.array(ddfile["r"][ish])
            #         mu_i = np.array(ddfile["mu"][ish])
            #         xrs_i = np.array(ddfile["xrs"][ish])
            #         yrs_i = np.array(ddfile["yrs"][ish])
            #         int_i = np.array(ddfile["intensity"][ish])  # * dfile.attrs["d_L"] ** 2
            #         gam_i = np.array(ddfile["gamma"][ish])
            #         B_i = np.array(ddfile["B"][ish])
            #         tb_i = np.array(ddfile["tburst"][ish])
            #         theta_i = np.array(ddfile["theta"][ish])
            #         if limit_data:
            #             idx1 = abs(mu_i) < min_mu  # TODO this is overritten!
            #             # int_i[idx1] =
            #             # print("len(gam[idx1])={}".format(gam_i[idx1]))
            #             idx = int_i > int_i.max() * min_mu_frac
            #             idx = idx1
            #             if verbose:
            #                 if (len(mu_i[idx1]) != len(mu_i[idx])):
            #                     print('\t', mu_i[idx1])
            #                     print('\t', mu_i[idx])
            #                     # exit(1)
            #             if verbose:
            #                 if len(idx) < 2:
            #                     print("No excess = {}".format(ii))
            #                 print("len(gam[idx])={}".format(gam_i[idx]))
            #                 print("Gamma={}".format(gam_i[idx]))
            #                 print("beta={}".format(get_beta(gam_i[idx])))
            #                 print("mu_i={}".format(mu_i[idx]))
            #                 print("r_i={}".format(r_i[idx]))
            #                 print("B={}".format(B_i[idx]))
            #                 print("tb_i={}".format(tb_i[idx]))
            #             # print(gam_i[idx], Beta(gam_i[idx]))
            #
            #             # int_i[mu_i < 1e-3] = 0.
            #             # int_i2 = int_i * abs(mu_i)# * np.sqrt(1 - mu_i**2)#/ ((1-mu_i)*(1+mu_i))
            #             # plt.figure(figsize=(4.6,3.2))
            #             # plt.semilogy(mu_i,int_i/int_i.max(), '.')
            #             # plt.semilogy(mu_i[idx],int_i[idx]/int_i[idx].max(), 'x', color='red')
            #             # plt.semilogy(mu_i,int_i2/int_i2.max(), 'o', color='red')
            #             # plt.xlabel(r"$\mu=\sin(\theta_{\rm obs}) \sin(\theta) \sin(\phi) + \cos(\theta_{\rm obs}) \cos(\theta)$")
            #             # plt.ylabel(r"$I/I_{\rm max}$")
            #             # plt.tight_layout()
            #             # plt.show()
            #             int_i[idx] = 0.
            #         # int_i /= abs(mu_i)[abs(mu_i)>0.01]
            #         i_zz = interp(xrs_i, yrs_i, int_i, xx_grid, yy_grid, 'linear')
            #         zz += i_zz
            #         all_xrs.append(xrs_i)
            #         all_yrs.append(yrs_i)
            #         all_zz.append(int_i)
            #     if verbose:
            #         print("\tAfter interpolation: I [{:.2e}, {:.2e}] Sum = {:.2e} [Total expected {:.2e}]"
            #               .format(zz.min(), zz.max(), np.sum(zz), float(
            #             np.array(dfile["totalflux at freq={:.4e}".format(freq)])[find_nearest_index(times, time)])))
            #     return (xx_grid, yy_grid, zz, all_xrs, all_yrs, all_zz)
        else:
            raise NotImplementedError("Not implemented")

    def get_jet_skymap_totfluxes(self, freq, time=None):

        self._check_if_loaded_jet_skymap()
        if time is None:
            return np.array(self.grb_skymap_dfile["totalflux at freq={:.4e}".format(freq)])
        else:
            arr = np.array(self.grb_skymap_dfile["totalflux at freq={:.4e}".format(freq)])
            times = self.get_jet_skymap_times()
            if len(times) > 1:
                assert self.get_jet_skymap_times().min() < time < self.get_jet_skymap_times().max()
                val = arr[find_nearest_index(self.get_jet_skymap_totfluxes(freq=freq, time=None), time)]
                return val
            else:
                return arr[0]

    ''' -------- ejecta (shells and layers) -------- '''

    # ejecta dynamics

    def get_ej_dyn_obj(self):
        self._ckeck_if_loaded_ej_dyn_obj()
        return self.kn_dyn_dfile

    def get_ej_dyn_arr(self, v_n, ishell=None, ilayer=None):
        self._ckeck_if_loaded_ej_dyn_obj()
        obj_dyn = self.get_ej_dyn_obj()
        nlayers = int(obj_dyn.attrs["nlayers"])
        nshells = int(obj_dyn.attrs["nshells"])
        dfile = self.kn_dyn_dfile
        if ((ishell is None) and (ilayer is None)):
            raise ValueError("both ishell and ilayer cannot be None for dynamical data")
        elif ((not ishell is None) and (ilayer is None)):
            arr = []
            for il in range(nlayers):
                arr.append(self.get_ej_dyn_arr(v_n, ishell=ishell, ilayer=il))
            arr = np.reshape(np.array(arr), newshape=(nlayers, len(arr[0])))
            return arr
        elif ((ishell is None) and (not ilayer is None)):
            arr = []
            for ish in range(nshells):
                arr.append(self.get_ej_dyn_arr(v_n, ishell=ish, ilayer=ilayer))
            arr = np.reshape(np.array(arr), newshape=(nshells, len(arr[0])))
            return arr
        elif ((not ishell is None) and (not ilayer is None)):
            layer = "shell={} layer={}".format(ishell, ilayer)
            if (not layer in list(dfile.keys())):
                raise NameError("Layer {} (aka '{}') is not in the ejecta dyn. file.\n Available: {}"
                                .format(ilayer, layer, dfile.keys()))
            if (not v_n in dfile[layer].keys()):
                raise NameError("v_n {} is not in the ejecta dyn. dfile[{}].keys() \n Avaialble: {}"
                                .format(v_n, layer, dfile[layer].keys()))
            return np.array(dfile[layer][v_n])
        else:
            raise NameError()

    def get_ej_dyn_1d_arr_layers(self, v_n="em", idx=0, ishell=None, ilayer=None):
        obj_dyn = self.get_ej_dyn_obj()
        nlayers = int(obj_dyn.attrs["nlayers"])
        nshells = int(obj_dyn.attrs["nshells"])
        x_arr = []
        if (not ishell is None and ilayer is None):
            for il in range(nlayers):
                x_arr.append(self.get_ej_dyn_arr(v_n, ishell=ishell, ilayer=il)[idx])
        else:
            for ish in range(nshells):
                x_arr.append(self.get_ej_dyn_arr(v_n, ishell=ish, ilayer=ilayer)[idx])
        return np.array(x_arr)

    # ejecta spectrum

    def get_ej_spec_obj(self):
        self._check_if_loaded_ej_spec()
        return self.kn_spec_dfile

    def get_ej_spec_2d_arr(self, v_n="em", ishell=0, ilayer=0):
        dfile = self.get_ej_spec_obj()
        layer = "shell={} layer={}".format(ishell, ilayer)
        if (not layer in list(dfile.keys())):
            raise NameError("Layer {} (aka '{}') is not in the ejecta comov.spec. file.\n Available: {}"
                            .format(ilayer, layer, dfile.keys()))
        if (not v_n in dfile[layer].keys()):
            raise NameError("v_n {} is not in the ejecta comov.spec. dfile[{}].keys() \n Avaialble: {}"
                            .format(v_n, layer, dfile[layer].keys()))
        return np.array(dfile[layer][v_n])

    def get_ej_spec_times(self):
        dfile = self.get_ej_spec_obj()
        return np.array(dfile["times"])

    def get_ej_spec_freqs(self):
        dfile = self.get_ej_spec_obj()
        return np.array(dfile["freqs"])

    def get_ej_spec_1d_arr(self, freq=None, time=None, v_n="em", ishell=0, ilayer=0):
        arr = self.get_ej_spec_2d_arr(v_n=v_n, ishell=ishell, ilayer=ilayer)
        freqs = self.get_ej_spec_freqs()
        times = self.get_ej_spec_times()
        if (not freq is None):
            if (freq < freqs.min() or freq > freqs.max()):
                raise ValueError("Freq={} (jet comov cpec) is not in the limit of freqs[{}, {}]"
                                 .format(freq, freqs.min(), freqs.max()))
            if (not freq in freqs):
                idx = find_nearest_index(freqs, freq)
            else:
                idx = int(np.where(freqs == freq)[0])
            arr_ = arr[idx, :]
            assert len(arr_) == len(times)
            return arr_
        if (not times is None):
            if (time < times.min() or time > times.max()):
                raise ValueError("Time={} (jet comov cpec) is not in the limit of times[{}, {}]"
                                 .format(time, times.min(), times.max()))
            if (not time in times):
                idx = find_nearest_index(time, times)
            else:
                idx = int(np.where(time == times)[0][0])
            arr_ = arr[:, idx]
            assert len(arr_) == len(freqs)
            return arr_

    def get_ej_spec_2d_arr_layers(self, freq=None, time=None, v_n="em", ishell=None, ilayer=None):
        obj_dyn = self.get_ej_dyn_obj()
        arr = []
        x_arr = []
        # x_arr = self.get_ej_dyn_arr("ctheta",ishell=ishell,ilayer=0)
        if (not ishell is None and ilayer is None):
            y_arr = self.get_ej_dyn_arr("R", ishell=ishell, ilayer=0)
            # freqs = self.get_ej_spec_freqs()
            # times = self.get_ej_spec_times()
            nlayers = int(obj_dyn.attrs["nlayers"])
            for il in range(nlayers):
                x_arr.append(self.get_ej_dyn_arr("ctheta", ishell=ishell, ilayer=il)[0])
                arr.append(self.get_ej_spec_1d_arr(freq=freq, time=time, v_n=v_n, ishell=ishell, ilayer=il))
            x_arr = np.array(x_arr)
            arr = np.reshape(np.array(arr), (len(x_arr), len(y_arr)))
        else:
            y_arr = self.get_ej_dyn_arr("R", ishell=0, ilayer=ilayer)
            nshells = int(obj_dyn.attrs["nshells"])
            for ish in range(nshells):
                x_arr.append(self.get_ej_dyn_arr("Gamma", ishell=ish, ilayer=ilayer)[0])
                arr.append(self.get_ej_spec_1d_arr(freq=freq, time=time, v_n=v_n, ishell=ish, ilayer=ilayer))
            x_arr = np.array(x_arr)
            arr = np.reshape(np.array(arr), (len(x_arr), len(y_arr)))
        return arr

    # ejecta lightcurves

    def get_ej_lc_obj(self):
        self._check_if_loaded_ej_lc()
        return self.kn_lc_dfile

    def get_ej_lc_times(self):
        dfile = self.get_ej_lc_obj()
        # print(dfile.keys())
        # print(dfile.attrs())

        return np.array(dfile["times"])

    def get_ej_lc_freqs(self):
        dfile = self.get_ej_lc_obj()
        return np.array(dfile["freqs"])

    def get_ej_lc_totalflux(self, freq):
        dfile = self.get_ej_lc_obj()
        nlayers = int(dfile.attrs["nlayers"])
        nshells = int(dfile.attrs["nshells"])
        times = self.get_ej_lc_times()
        freqs = self.get_ej_lc_freqs()
        try:
            key = str("totalflux at freq={:.4e}".format(freq)).replace('.', ',')
            if not key in dfile.keys():
                raise NameError("Not found: {} among keys:{}".format(key, [key for key in dfile.keys() if key.__contains__("totalflux")]))
        except NameError:
            key = str("totalflux at freq={:.4e}".format(freq))
            if not key in dfile.keys():
                raise NameError("Not found for ej. lightcurve: {} among keys:{} \n ejecta_prefix:{}"
                                .format(key, [key for key in dfile.keys() if key.__contains__("totalflux")], self.ejecta_prefix))
        except:
            raise NameError()
        return np.array(dfile[key])

    def get_ej_lc(self, freq=None, ishell=None, ilayer=None, v_n=None):
        dfile = self.get_ej_lc_obj()
        nlayers = int(dfile.attrs["nlayers"])
        nshells = int(dfile.attrs["nshells"])
        times = self.get_ej_lc_times()
        freqs = self.get_ej_lc_freqs()
        if (v_n is None): v_n = list(dfile["shell=0 layer=0"].keys())[0]  # Default key

        if (freq is None):
            # spectum
            if ((ishell is None) and (ilayer is None)):
                # fluxes2d = np.zeros((len(freqs), len(times)))
                # for ish in range(nshells):
                #     for il in range(nlayers):
                #         fluxes2d += np.array(dfile["shell={} layer={}".format(ish, il)][v_n])# [freq,time]
                # return fluxes2d
                fluxes2d = []
                for ifreq in freqs:
                    fluxes2d.append(self.get_ej_lc_totalflux(freq=ifreq))  # [freq,time]
                fluxes2d = np.reshape(np.array(fluxes2d), (len(freqs), len(times)))
                return fluxes2d
            elif ((ishell is None) and (not ilayer is None)):
                fluxes2d = np.zeros((len(freqs), len(times)))
                for ish in range(nshells):
                    fluxes2d += np.array(dfile["shell={} layer={}".format(ish, ilayer)][v_n])  # [freq,time]
                return fluxes2d
            elif ((not ishell is None) and (ilayer is None)):
                fluxes2d = np.zeros((len(freqs), len(times)))
                for il in range(nlayers):
                    fluxes2d += np.array(dfile["shell={} layer={}".format(ishell, il)][v_n])  # [freq,time]
                return fluxes2d
            elif ((not ishell is None) and (not ilayer is None)):
                fluxes2d = np.array(dfile["shell={} layer={}".format(ishell, ilayer)][v_n])  # [freq,time]
                return fluxes2d
            else:
                raise NameError()
        else:
            # light curves
            if (not freq in self.get_ej_lc_freqs()):
                raise ValueError("freq:{} is not in ej_lc Given:{}".format(freq, self.get_ej_lc_freqs()))
            ifreq = find_nearest_index(self.get_ej_lc_freqs(), freq)
            if ((ishell is None) and (ilayer is None)):
                return self.get_ej_lc_totalflux(freq=freq)
            elif ((ishell is None) and (not ilayer is None)):
                if (ilayer > nlayers - 1):
                    raise ValueError("Layer={} is > nlayers={}".format(ilayer, nlayers - 1))
                # fluxes1d = np.zeros_like(times)
                # for ish in range(nshells):
                #     fluxes1d += np.array(dfile[v_n]["shell={} layer={}".format(ish, ilayer)][ifreq])  # [freq,time]
                # return fluxes1d
                fluxes2d = []
                for ish in range(nshells):
                    fluxes2d.append(np.array(dfile["shell={} layer={}".format(ish, ilayer)][v_n][ifreq]))  # [freq,time]
                return np.reshape(fluxes2d, newshape=(nshells, len(times)))
            elif ((not ishell is None) and (ilayer is None)):
                if (ishell > nshells - 1):
                    raise ValueError("Shell={} is > nshell={}".format(ishell, nshells - 1))
                # fluxes1d = np.zeros_like(times)
                # for il in range(nlayers):
                #     fluxes1d += np.array(dfile[v_n]["shell={} layer={}".format(ishell, il)][ifreq])  # [freq,time]
                # return fluxes1d
                fluxes2d = []
                for il in range(nlayers):
                    fluxes2d.append(np.array(dfile["shell={} layer={}".format(ishell, il)][v_n][ifreq]))  # [freq,time]
                return np.reshape(fluxes2d, newshape=(nlayers, len(times)))
            elif ((not ishell is None) and (not ilayer is None)):
                fluxes1d = np.array(dfile["shell={} layer={}".format(ishell, ilayer)][v_n][ifreq])  # [freq,time]
                return fluxes1d
            else:
                raise NameError()

    def _alpha(self, freq, freqm1, ishell, ilayer, v_n):
        values_i = self.get_ej_lc(freq=freq, ishell=ishell, ilayer=ilayer, v_n=v_n)
        values_im1 = self.get_ej_lc(freq=freqm1, ishell=ishell, ilayer=ilayer, v_n=v_n)
        ffreq_i = np.full(values_i.shape, freq)
        ffreq_im1 = np.full(values_im1.shape, freqm1)
        num = np.log10(values_i) - np.log10(values_im1)
        denum = np.log10(ffreq_i) - np.log10(ffreq_im1)
        values = 1 * num / denum
        return values

    def get_ej_lc_spec_idx(self, freq, ishell=None, ilayer=None, freq1=None, freq2=None):
        freqs = self.get_ej_lc_freqs()
        times = self.get_ej_lc_times()
        if ((freq is None) and (ishell is None) and (ilayer is None)):
            if (not freq1 is None) or (not freq2 is None):
                raise KeyError("For 2d spectral index freq1 and freq2 should be None (for now)")
            arr = []
            for ifreq in range(1, len(freqs)):
                arr.append( self._alpha(freq=freqs[ifreq], freqm1=freqs[ifreq - 1], ishell=ishell, ilayer=ilayer, v_n=None) )
            return np.reshape(np.array(arr), newshape=(len(freqs) - 1, len(times)))
        else:
            idx = find_nearest_index(freqs, freq)
            _freq1 = freqs[idx] if freq1 is None else freq1
            _freq2 = freqs[idx - 1] if freq2 is None else freq2
            return self._alpha(freq=_freq1, freqm1=_freq2, ishell=ishell, ilayer=ilayer, v_n=None)

    def get_ej_lc_temp_idx(self, freq, ishell, ilayer):
        freqs = self.get_ej_lc_freqs()
        t = self.get_ej_lc_times()
        Lnu = self.get_ej_lc(freq=freq, ishell=ishell, ilayer=ilayer, v_n=None)
        val = np.zeros_like(Lnu)
        val[1:-1] = np.log10(Lnu[2:] / Lnu[:-2]) / np.log10(t[2:] / t[:-2])
        val[0] = np.log10(Lnu[1] / Lnu[0]) / np.log10(t[1] / t[0])
        val[-1] = np.log10(Lnu[-1] / Lnu[-2]) / np.log10(t[-1] / t[-2])
        return val

    # ejecta skymaps

    def get_ej_skymap_obj(self):
        self._check_if_loaded_ej_skymap()
        return self.kn_skymap_dfile

    def get_ej_skymap_times(self):
        self._check_if_loaded_ej_skymap()
        # print(self.ej_skymap.keys())
        return np.array(self.kn_skymap_dfile["times"])

    def get_ej_skymap_freqs(self):
        self._check_if_loaded_ej_skymap()
        return np.array(self.kn_skymap_dfile["freqs"])

    def get_ej_skymap_totfluxes(self, freq, shell=None, time=None):
        self._check_if_loaded_ej_skymap()
        if time is None:
            if (shell is None):
                return np.array(self.kn_skymap_dfile["totalflux at freq={:.4e}".format(freq)])
            else:
                return np.array(self.kn_skymap_dfile["totalflux at freq={:.4e} shell={}".format(freq, shell)])
        else:
            if (shell is None):
                arr = np.array(self.kn_skymap_dfile["totalflux at freq={:.4e}".format(freq)])
                if (len(arr)==1):
                    return arr[0]
                assert self.get_ej_skymap_times().min() < time < self.get_ej_skymap_times().max()
                val = arr[find_nearest_index(self.get_ej_skymap_totfluxes(freq=freq, shell=None, time=None), time)]
                return val
            else:
                raise KeyError("Not finished...")

    def get_ej_skymap_cm(self, all_xrs, all_yrs, all_zz):
        dfile = self.get_ej_skymap_obj()
        # _x = np.concatenate(all_xrs)
        # _y = np.concatenate(all_yrs)
        # _z = np.concatenate(all_zz)
        xc_m, yc_m = compute_position_of_the_flux_centroid(all_xrs, all_yrs, all_zz, float(dfile.attrs["d_L"]))
        return (xc_m, yc_m)

    @staticmethod
    def combine_images(xs, ys, datas, verbose=False, hist_or_int="int", shells=False, nx=200, ny=100, extend=2,
                       retrun_edges=False, edges_x=None, edges_y=None):

        if shells:
            assert len(xs) == len(ys)
            assert len(xs) == len(datas)
            nshells = len(xs)

            xmin_neg, xmin_pos, xmax, xmin = [], [], [], []
            ymin_neg, ymin_pos, ymax, ymin = [], [], [], []
            i_min, i_max = [], []
            i_shells = []
            for ish in range(nshells):
                xrs_i = np.array(xs[ish])
                yrs_i = np.array(ys[ish])
                int_i = np.array(datas[ish])
                if (np.sum(int_i) == 0):
                    continue
                i_shells.append(ish)
                xmin_neg.append(xrs_i[xrs_i < 0].min())
                xmin_pos.append(xrs_i[xrs_i > 0].min())
                xmax.append(xrs_i.max())
                xmin.append(xrs_i.min())
                ymin_neg.append(yrs_i[yrs_i < 0].min())
                ymin_pos.append(yrs_i[yrs_i > 0].min())
                ymax.append(yrs_i.max())
                ymin.append(yrs_i.min())
                i_min.append(int_i.min())
                i_max.append(int_i.max())
            xmin_neg = np.array(xmin_neg)
            xmin_pos = np.array(xmin_pos)
            xmax = np.array(xmax)
            xmin = np.array(xmin)
            ymin_neg = np.array(ymin_neg)
            ymin_pos = np.array(ymin_pos)
            ymax = np.array(ymax)
            ymin = np.array(ymin)
            i_min = np.array(i_min)
            i_max = np.array(i_max)
            if verbose:
                for i in range(len(i_shells)):
                    print("\tshell={} xmin_neg={:.2e} xmin_pos={:.2e} xmax={:.2e}".format(i, xmin_neg[i], xmin_pos[i],
                                                                                          xmax[i]))
                    print("\t         ymin_neg={:.2e} ymin_pos={:.2e} ymax={:.2e}".format(i, ymin_neg[i], ymin_pos[i],
                                                                                          ymax[i]))
                print("\toverall X min(xmin_pos) = {:.2e} max(xmax) = {:.2e}".format(xmin_pos.min(), xmax.max()))
                print("\toverall Y min(ymin_pos) = {:.2e} max(ymax) = {:.2e}".format(ymin_pos.min(), ymax.max()))
                print("\toverall I min(i_min)    = {:.2e} max(ymax) = {:.2e}".format(i_min.min(), i_max.max()))
            x_grid = np.concatenate((-1. * np.logspace(np.log10(xmin_pos.min()), np.log10(xmax.max()), nx)[::-1],
                                     np.logspace(np.log10(xmin_pos.min()), np.log10(xmax.max()), nx)))
            y_grid = np.concatenate((-1. * np.logspace(np.log10(ymin_pos.min()), np.log10(ymax.max()), ny)[::-1],
                                     np.logspace(np.log10(ymin_pos.min()), np.log10(ymax.max()), ny)))

            if edges_x is None: edges_x = np.linspace(xmin.min()*extend, xmax.max()*extend, num=nx)
            if edges_y is None: edges_y = np.linspace(ymin.min()*extend, ymax.max()*extend, num=ny)

            if verbose:
                print("\tGrid xmin_neg={:.2e} xmin_pos={:.2e} xmax={:.2e}".format(x_grid.min(),
                                                                                  x_grid[x_grid > 0].min(),
                                                                                  x_grid.max()))
                print("\tGrid ymin_neg={:.2e} ymin_pos={:.2e} ymax={:.2e}".format(y_grid.min(),
                                                                                  y_grid[y_grid > 0].min(),
                                                                                  y_grid.max()))
            xx_grid, yy_grid = np.meshgrid(x_grid, y_grid)
            if hist_or_int == "int":
                zz = np.zeros_like((xx_grid))
            else:
                zz = np.zeros((len(edges_x) - 1, len(edges_y) - 1))
            # interpolate onto the grid that covers all images
            all_xrs, all_yrs, all_zz = [], [], []
            for ii, ish in enumerate(i_shells):  # range(len(i_shells))
                if verbose: print("Pocessing: shell={} [{}/{}]".format(ish, ii, len(i_shells)))
                xrs_i = np.array(xs[ish])
                yrs_i = np.array(ys[ish])
                int_i = np.array(datas[ish])  # * dfile.attrs["d_L"] ** 2

                if hist_or_int == "int":
                    # nx = np.complex(0, nx)
                    # ny = np.complex(0, ny)
                    # grid_x, grid_y = np.mgrid[xrs_i.min() * 1.2:xrs_i.max() * 1.2:nx,
                    #                  yrs_i.min() * 1.2:yrs_i.max() * 1.2:ny]
                    i_zz = interp(xrs_i, yrs_i, int_i, xx_grid, yy_grid, 'linear')
                    # return (grid_x, grid_y, i_zz, xrs_i, yrs_i, int_i)
                else:
                    # nx = 100
                    # ny = 100
                    # nx = np.complex(0, nx + 1)
                    # ny = np.complex(0, ny + 1)
                    # edges_x = np.mgrid[xrs_i.min() * 1.2:xrs_i.max() * 1.2:nx]
                    # edges_y = np.mgrid[yrs_i.min() * 1.2:yrs_i.max() * 1.2:ny]
                    # grid_x, grid_y = np.mgrid[xrs_i.min() * 1.2:xrs_i.max() * 1.2:nx,
                    #                  yrs_i.min() * 1.2:yrs_i.max() * 1.2:ny]
                    i_zz, _ = np.histogramdd(tuple([xrs_i, yrs_i]), bins=tuple([edges_x, edges_y]), weights=int_i)
                    # grid_x = 0.5 * (edges_x[1:] + edges_x[:-1])
                    # grid_y = 0.5 * (edges_y[1:] + edges_y[:-1])
                    # return (grid_x, grid_y, i_zz, xrs_i, yrs_i, int_i)

                zz += i_zz
                # all_xrs.append(xrs_i)
                # all_yrs.append(yrs_i)
                # all_zz.append(int_i)
                #
                # if hist_or_int == "int":
                #     i_zz = interp(xrs_i, yrs_i, int_i, xx_grid, yy_grid, 'linear')
                #     zz += i_zz
                #     all_xrs.append(xrs_i)
                #     all_yrs.append(yrs_i)
                #     all_zz.append(int_i)
                # else:
                #     nx = 2000
                #     ny = 1000
                #     nx = np.complex(0, nx + 1)
                #     ny = np.complex(0, ny + 1)
                #     edges_x = np.mgrid[x_grid.min():x_grid.max():nx]
                #     edges_y = np.mgrid[y_grid.min():yrs_i.max():ny]
                #     # grid_x, grid_y = np.mgrid[xrs_i.min() * 1.2:xrs_i.max() * 1.2:nx,
                #     #                  yrs_i.min() * 1.2:yrs_i.max() * 1.2:ny]
                #     i_zz, _ = np.histogramdd(tuple([xrs_i, yrs_i]), bins=tuple([edges_x, edges_y]), weights=int_i)
            xx_grid = 0.5 * (edges_x[1:] + edges_x[:-1])
            yy_grid = 0.5 * (edges_y[1:] + edges_y[:-1])
            #     return (grid_x, grid_y, i_zz, xrs_i, yrs_i, int_i)

            if retrun_edges:
                return (xx_grid, yy_grid, zz, edges_x, edges_y)
            else:
                return (xx_grid, yy_grid, zz)
        else:

            xs = np.concatenate(xs)
            ys = np.concatenate(ys)
            datas = np.concatenate(datas)

            if hist_or_int == "int":
                nx = np.complex(0, nx)
                ny = np.complex(0, ny)
                grid_x, grid_y = np.mgrid[xs.min()*extend:xs.max()*extend:nx, ys.min()*extend:ys.max()*extend:ny]
                i_zz = interp(xs, ys, datas, grid_x, grid_y, 'linear')
            elif hist_or_int == "both":

                nx = np.complex(0, nx)
                ny = np.complex(0, ny)
                grid_x, grid_y = np.mgrid[xs.min():xs.max():nx, ys.min():ys.max():ny]
                i_zz = interp(xs, ys, datas, grid_x, grid_y, 'linear')

                nx = np.complex(0, nx + 1)
                ny = np.complex(0, ny + 1)
                edges_x = np.mgrid[xs.min():xs.max():nx]
                edges_y = np.mgrid[ys.min():ys.max():ny]
                # grid_x, grid_y = np.mgrid[xrs_i.min() * 1.2:xrs_i.max() * 1.2:nx,
                #                  yrs_i.min() * 1.2:yrs_i.max() * 1.2:ny]
                i_zz, _ = np.histogramdd(tuple([grid_x.flatten(), grid_y.flatten()]), bins=tuple([edges_x, edges_y]),
                                         weights=i_zz.T.flatten())
                grid_x = 0.5 * (edges_x[1:] + edges_x[:-1])
                grid_y = 0.5 * (edges_y[1:] + edges_y[:-1])
                # print(i_zz.shape)
            else:

                # nx = 200
                # ny = 100
                nx = np.complex(0, nx + 1)
                ny = np.complex(0, ny + 1)
                edges_x = np.mgrid[xs.min()*extend:xs.max()*extend:nx]
                edges_y = np.mgrid[ys.min()*extend:ys.max()*extend:ny]
                # grid_x, grid_y = np.mgrid[xrs_i.min() * 1.2:xrs_i.max() * 1.2:nx,
                #                  yrs_i.min() * 1.2:yrs_i.max() * 1.2:ny]
                i_zz, _ = np.histogramdd(tuple([xs, ys]), bins=tuple([edges_x, edges_y]), weights=datas)
                grid_x = 0.5 * (edges_x[1:] + edges_x[:-1])
                grid_y = 0.5 * (edges_y[1:] + edges_y[:-1])
                print(i_zz.shape)

            # i_zz = ndimage.uniform_filter(i_zz, )
            # i_zz = ndimage.filters.gaussian_filter(i_zz, [10,10], mode='reflect')

            return (grid_x, grid_y, i_zz)
        #
        # assert len(xs) == len(ys)
        # assert len(xs) == len(datas)
        # nshells = len(xs)
        # nx = 1000
        # ny = 1000
        # xmin_neg, xmin_pos, xmax, xmin = [], [], [], []
        # ymin_neg, ymin_pos, ymax, ymin = [], [], [], []
        # i_min, i_max = [], []
        # i_shells = []
        # for ish in range(nshells):
        #     xrs_i = np.array(xs[ish])
        #     yrs_i = np.array(ys[ish])
        #     int_i = np.array(datas[ish])
        #     if (np.sum(int_i) == 0):
        #         continue
        #     i_shells.append(ish)
        #     xmin_neg.append(xrs_i[xrs_i < 0].min())
        #     xmin_pos.append(xrs_i[xrs_i > 0].min())
        #     xmax.append(xrs_i.max())
        #     xmin.append(xrs_i.min())
        #     ymin_neg.append(yrs_i[yrs_i < 0].min())
        #     ymin_pos.append(yrs_i[yrs_i > 0].min())
        #     ymax.append(yrs_i.max())
        #     ymin.append(yrs_i.min())
        #     i_min.append(int_i.min())
        #     i_max.append(int_i.max())
        # xmin_neg = np.array(xmin_neg)
        # xmin_pos = np.array(xmin_pos)
        # xmax = np.array(xmax)
        # ymin_neg = np.array(ymin_neg)
        # ymin_pos = np.array(ymin_pos)
        # ymax = np.array(ymax)
        # i_min = np.array(i_min)
        # i_max = np.array(i_max)
        # if verbose:
        #     for i in range(len(i_shells)):
        #         print("\tshell={} xmin_neg={:.2e} xmin_pos={:.2e} xmax={:.2e}".format(i, xmin_neg[i], xmin_pos[i],
        #                                                                               xmax[i]))
        #         print("\t         ymin_neg={:.2e} ymin_pos={:.2e} ymax={:.2e}".format(i, ymin_neg[i], ymin_pos[i],
        #                                                                               ymax[i]))
        #     print("\toverall X min(xmin_pos) = {:.2e} max(xmax) = {:.2e}".format(xmin_pos.min(), xmax.max()))
        #     print("\toverall Y min(ymin_pos) = {:.2e} max(ymax) = {:.2e}".format(ymin_pos.min(), ymax.max()))
        #     print("\toverall I min(i_min)    = {:.2e} max(ymax) = {:.2e}".format(i_min.min(), i_max.max()))
        # x_grid = np.concatenate((-1. * np.logspace(np.log10(xmin_pos.min()), np.log10(xmax.max()), nx)[::-1],
        #                          np.logspace(np.log10(xmin_pos.min()), np.log10(xmax.max()), nx)))
        # y_grid = np.concatenate((-1. * np.logspace(np.log10(ymin_pos.min()), np.log10(ymax.max()), ny)[::-1],
        #                          np.logspace(np.log10(ymin_pos.min()), np.log10(ymax.max()), ny)))
        #
        # edges_x = np.linspace(xmin.min(), xmax.max(), num=nx)
        # edges_y = np.linspace(ymin.min(), ymax.max(), num=ny)
        #
        # if verbose:
        #     print("\tGrid xmin_neg={:.2e} xmin_pos={:.2e} xmax={:.2e}".format(x_grid.min(),
        #                                                                       x_grid[x_grid > 0].min(),
        #                                                                       x_grid.max()))
        #     print("\tGrid ymin_neg={:.2e} ymin_pos={:.2e} ymax={:.2e}".format(y_grid.min(),
        #                                                                       y_grid[y_grid > 0].min(),
        #                                                                       y_grid.max()))
        # xx_grid, yy_grid = np.meshgrid(x_grid, y_grid)
        # if hist_or_int == "int":
        #     zz = np.zeros_like((xx_grid))
        # else:
        #     zz = np.zeros((len(edges_x) - 1, len(edges_y) - 1))
        # # interpolate onto the grid that covers all images
        # all_xrs, all_yrs, all_zz = [], [], []
        # for ii, ish in enumerate(i_shells):  # range(len(i_shells))
        #     if verbose: print("Pocessing: shell={} [{}/{}]".format(ish, ii, len(i_shells)))
        #     xrs_i = np.array(xs[ish])
        #     yrs_i = np.array(ys[ish])
        #     int_i = np.array(datas[ish])  # * dfile.attrs["d_L"] ** 2
        #
        #     if hist_or_int == "int":
        #         # nx = np.complex(0, nx)
        #         # ny = np.complex(0, ny)
        #         # grid_x, grid_y = np.mgrid[xrs_i.min() * 1.2:xrs_i.max() * 1.2:nx,
        #         #                  yrs_i.min() * 1.2:yrs_i.max() * 1.2:ny]
        #         i_zz = interp(xrs_i, yrs_i, int_i, xx_grid, yy_grid, 'linear')
        #         # return (grid_x, grid_y, i_zz, xrs_i, yrs_i, int_i)
        #     else:
        #         # nx = 100
        #         # ny = 100
        #         # nx = np.complex(0, nx + 1)
        #         # ny = np.complex(0, ny + 1)
        #         # edges_x = np.mgrid[xrs_i.min() * 1.2:xrs_i.max() * 1.2:nx]
        #         # edges_y = np.mgrid[yrs_i.min() * 1.2:yrs_i.max() * 1.2:ny]
        #         # grid_x, grid_y = np.mgrid[xrs_i.min() * 1.2:xrs_i.max() * 1.2:nx,
        #         #                  yrs_i.min() * 1.2:yrs_i.max() * 1.2:ny]
        #         i_zz, _ = np.histogramdd(tuple([xrs_i, yrs_i]), bins=tuple([edges_x, edges_y]), weights=int_i)
        #         # grid_x = 0.5 * (edges_x[1:] + edges_x[:-1])
        #         # grid_y = 0.5 * (edges_y[1:] + edges_y[:-1])
        #         # return (grid_x, grid_y, i_zz, xrs_i, yrs_i, int_i)
        #
        #     zz += i_zz
        #     all_xrs.append(xrs_i)
        #     all_yrs.append(yrs_i)
        #     all_zz.append(int_i)
        #     #
        #     # if hist_or_int == "int":
        #     #     i_zz = interp(xrs_i, yrs_i, int_i, xx_grid, yy_grid, 'linear')
        #     #     zz += i_zz
        #     #     all_xrs.append(xrs_i)
        #     #     all_yrs.append(yrs_i)
        #     #     all_zz.append(int_i)
        #     # else:
        #     #     nx = 2000
        #     #     ny = 1000
        #     #     nx = np.complex(0, nx + 1)
        #     #     ny = np.complex(0, ny + 1)
        #     #     edges_x = np.mgrid[x_grid.min():x_grid.max():nx]
        #     #     edges_y = np.mgrid[y_grid.min():yrs_i.max():ny]
        #     #     # grid_x, grid_y = np.mgrid[xrs_i.min() * 1.2:xrs_i.max() * 1.2:nx,
        #     #     #                  yrs_i.min() * 1.2:yrs_i.max() * 1.2:ny]
        #     #     i_zz, _ = np.histogramdd(tuple([xrs_i, yrs_i]), bins=tuple([edges_x, edges_y]), weights=int_i)
        #     #     grid_x = 0.5 * (edges_x[1:] + edges_x[:-1])
        #     #     grid_y = 0.5 * (edges_y[1:] + edges_y[:-1])
        #     #     return (grid_x, grid_y, i_zz, xrs_i, yrs_i, int_i)
        # # if verbose:
        #     # print("\tAfter interpolation: I [{:.2e}, {:.2e}] Sum = {:.2e} [Total expected {:.2e}]"
        #     #       .format(zz.min(), zz.max(), np.sum(zz), float(
        #     #     np.array(dfile["totalflux at freq={:.4e}".format(freq)])[find_nearest_index(times, time)])))

        return (xx_grid, yy_grid, zz, all_xrs, all_yrs, all_zz)

    def get_combained_ej_skymaps_adjusted_to_other(self, time, freq, other_pb_instance, nx=100, ny=50):

        all_x, all_y, all_fluxes \
            = self.get_ej_skymap(time=time * cgs.day, freq=freq, verbose=False, remove_mu=True)
        # xcs_m, ycs_m = pb.get_ej_skymap_cm(all_x, all_y, all_fluxes)
        int_x, int_y, int_zz, edges_x, edges_y = self.combine_images(all_x, all_y, all_fluxes,
                                                                     hist_or_int="hist", shells=True, nx=nx,
                                                                     ny=ny, retrun_edges=True)
        all_x_m1, all_y_m1, all_fluxes_m1 \
            = other_pb_instance.get_ej_skymap(time=time * cgs.day, freq=freq, verbose=False, remove_mu=True)
        _, _, int_zz_m1 = other_pb_instance.combine_images(all_x_m1, all_y_m1, all_fluxes_m1,
                                                           hist_or_int="hist", shells=True, nx=nx,
                                                           ny=ny, edges_x=edges_x, edges_y=edges_y)

        return (int_x, int_y, int_zz, int_zz_m1)

    def get_combined_ej_spectral_map(self, time, freq, nx=100, ny=50, extend=2):
        freqs = self.get_ej_skymap_freqs()
        assert len(freqs) > 2
        idx = find_nearest_index(freqs, freq)
        assert idx != len(freqs)-1
        all_x, all_y, all_fluxes \
            = self.get_ej_skymap(time=time * cgs.day, freq=freqs[idx], verbose=False, remove_mu=True)
        # xcs_m, ycs_m = pb.get_ej_skymap_cm(all_x, all_y, all_fluxes)
        int_x, int_y, int_zz, edges_x, edges_y = self.combine_images(all_x, all_y, all_fluxes,
                                                                     hist_or_int="hist", shells=True, nx=nx,
                                                                     ny=ny, extend=extend, retrun_edges=True)

        all_x_m1, all_y_m1, all_fluxes_m1 \
            = self.get_ej_skymap(time=time * cgs.day, freq=freqs[idx - 1], verbose=False, remove_mu=True)
        _, _, int_zz_m1 = self.combine_images(all_x_m1, all_y_m1, all_fluxes_m1,
                                              hist_or_int="hist", shells=True, nx=nx,
                                              ny=ny, extend=extend, edges_x=edges_x, edges_y=edges_y)

        all_x_p1, all_y_p1, all_fluxes_p1 \
            = self.get_ej_skymap(time=time * cgs.day, freq=freqs[idx + 1], verbose=False, remove_mu=True)
        _, _, int_zz_p1 = self.combine_images(all_x_p1, all_y_p1, all_fluxes_p1,
                                              hist_or_int="hist", shells=True, nx=nx,
                                              ny=ny, extend=extend, edges_x=edges_x, edges_y=edges_y)

        # ffreq_i = np.full(int_zz.shape, freqs[idx])
        ffreq_ip1 = np.full(int_zz_p1.shape, freqs[idx + 1])
        ffreq_im1 = np.full(int_zz_m1.shape, freqs[idx - 1])
        ffreq = np.full(int_zz.shape, freqs[idx])
        # num = np.log10(int_zz) - np.log10(int_zz_m1)
        num = np.log10(int_zz_p1) - np.log10(int_zz_m1)
        # denum = np.log10(ffreq) - np.log10(ffreq_im1)
        denum = np.log10(ffreq_ip1) - np.log10(ffreq_im1)
        # int_zz = num / denum
        int_zz = 1.0 * num / denum / 2.0

        return (int_x, int_y, int_zz)

    def get_combined_jet_spectral_map(self, time, freq, nx=100, ny=50, extend=2):

        # TODO NOT TESTED

        freqs = self.get_jet_skymap_freqs()
        assert len(freqs) > 2
        idx = find_nearest_index(freqs, freq)
        assert idx != len(freqs)-1
        all_x, all_y, all_fluxes \
            = self.get_jet_skymap(time=time * cgs.day, freq=freqs[idx], verbose=False, remove_mu=True)
        # xcs_m, ycs_m = pb.get_ej_skymap_cm(all_x, all_y, all_fluxes)
        int_x, int_y, int_zz, edges_x, edges_y = self.combine_images(all_x, all_y, all_fluxes,
                                                                     hist_or_int="hist", shells=False, nx=nx,
                                                                     ny=ny, extend=extend, retrun_edges=True)

        all_x_m1, all_y_m1, all_fluxes_m1 \
            = self.get_jet_skymap(time=time * cgs.day, freq=freqs[idx - 1], verbose=False, remove_mu=True)
        _, _, int_zz_m1 = self.combine_images(all_x_m1, all_y_m1, all_fluxes_m1,
                                              hist_or_int="hist", shells=False, nx=nx,
                                              ny=ny, extend=extend, edges_x=edges_x, edges_y=edges_y)

        all_x_p1, all_y_p1, all_fluxes_p1 \
            = self.get_jet_skymap(time=time * cgs.day, freq=freqs[idx + 1], verbose=False, remove_mu=True)
        _, _, int_zz_p1 = self.combine_images(all_x_p1, all_y_p1, all_fluxes_p1,
                                              hist_or_int="hist", shells=False, nx=nx,
                                              ny=ny, extend=extend, edges_x=edges_x, edges_y=edges_y)

        # ffreq_i = np.full(int_zz.shape, freqs[idx])
        ffreq_ip1 = np.full(int_zz_p1.shape, freqs[idx + 1])
        ffreq_im1 = np.full(int_zz_m1.shape, freqs[idx - 1])
        num = np.log10(int_zz_p1) - np.log10(int_zz_m1)
        denum = np.log10(ffreq_ip1) - np.log10(ffreq_im1)
        int_zz = 1. * num / denum / 2.

        return (int_x, int_y, int_zz)

    def get_combined_kn_grb_spectral_map(self, time, freq, nx=100, ny=50, extend=2):

        # TODO NOT TESTED

        raise ValueError(" I AM NOT DONE!!! ")

        # freqs = self.get_ej_skymap_freqs()
        # assert len(freqs) > 2
        # idx = find_nearest_index(freqs, freq)
        # assert idx != len(freqs)-1
        # all_x, all_y, all_fluxes \
        #     = self.get_jet_skymap(time=time * cgs.day, freq=freqs[idx], verbose=False, remove_mu=True)
        # # xcs_m, ycs_m = pb.get_ej_skymap_cm(all_x, all_y, all_fluxes)
        # int_x, int_y, int_zz, edges_x, edges_y = self.combine_images(all_x, all_y, all_fluxes,
        #                                                            hist_or_int="hist", shells=False, nx=nx,
        #                                                            ny=ny, extend=extend, retrun_edges=True)
        #
        # all_x_m1, all_y_m1, all_fluxes_m1 \
        #     = self.get_jet_skymap(time=time * cgs.day, freq=freqs[idx - 1], verbose=False, remove_mu=True)
        # _, _, int_zz_m1 = self.combine_images(all_x_m1, all_y_m1, all_fluxes_m1,
        #                                     hist_or_int="hist", shells=False, nx=nx,
        #                                     ny=ny, extend=extend, edges_x=edges_x, edges_y=edges_y)
        #
        # all_x_p1, all_y_p1, all_fluxes_p1 \
        #     = self.get_jet_skymap(time=time * cgs.day, freq=freqs[idx + 1], verbose=False, remove_mu=True)
        # _, _, int_zz_p1 = self.combine_images(all_x_p1, all_y_p1, all_fluxes_p1,
        #                                     hist_or_int="hist", shells=False, nx=nx,
        #                                     ny=ny, extend=extend, edges_x=edges_x, edges_y=edges_y)
        #
        # # ffreq_i = np.full(int_zz.shape, freqs[idx])
        # ffreq_ip1 = np.full(int_zz_p1.shape, freqs[idx + 1])
        # ffreq_im1 = np.full(int_zz_m1.shape, freqs[idx - 1])
        # num = np.log10(int_zz_p1) - np.log10(int_zz_m1)
        # denum = np.log10(ffreq_ip1) - np.log10(ffreq_im1)
        # int_zz = 1. * num / denum / 2.
        #
        # return (int_x, int_y, int_zz)

    def get_ej_skymap(self, time=None, freq=None, ishell=None, verbose=False, remove_mu=False, renormalize=True):

        # nx = 200
        # ny = 100
        # min_mu = 1e-4 # limit for mu that if too small blow up the intensity for some reason...
        # min_mu_frac = 0.85
        # # limit_data = True
        times = self.get_ej_skymap_times()
        freqs = self.get_ej_skymap_freqs()
        dfile = self.get_ej_skymap_obj()
        nshells = int(dfile.attrs["nshells"])
        d_l = float(self.get_ej_skymap_obj().attrs["d_L"])
        if ((not time is None) and (not time in times)):
            raise ValueError(
                "time={} day is not in the list for skypams={} days".format(time / cgs.day, times / cgs.day))
        if ((not freq is None) and (not freq in freqs)):
            raise ValueError("freq={} is not in the list for skypams={}".format(freq, freqs))
        if ((not ishell is None) and (ishell > nshells - 1)):
            raise ValueError("shell={} is beying skymap nshells={}".format(ishell, nshells))
        if ((not time is None) and (not freq is None)):
            # print(dfile.keys())
            ddfile = dfile["time={:.4e} freq={:.4e}".format(time, freq)]
            if (not ishell is None):
                r_i = np.array(ddfile["r"][ishell])
                mu_i = np.array(ddfile["mu"][ishell])
                xrs_i = np.array(ddfile["xrs"][ishell]) * cgs.rad2mas / d_l  # m -> mas
                yrs_i = np.array(ddfile["yrs"][ishell]) * cgs.rad2mas / d_l  # m -> mas
                int_i = np.array(ddfile["intensity"][ishell]) * (d_l ** 2 / cgs.rad2mas ** 2)  # -> mJy / mas^2
                gam_i = np.array(ddfile["gamma"][ishell])
                B_i = np.array(ddfile["B"][ishell])
                tb_i = np.array(ddfile["tburst"][ishell])
                theta_i = np.array(ddfile["theta"][ishell])
                phi_i = np.array(ddfile["phi"][ishell])

                if remove_mu:
                    print("Removing 'mu' from ejecta skymap")
                    int_i *= np.abs( mu_i)  # TODO I was produced as F / (R^2 abs(mu)), where abs(mu)->0 and I->inf. Problem!!!

                if renormalize:
                    print("Renormalizing ejecta skymap (shell my shell separately)")
                    fnus = self.get_ej_skymap_totfluxes(freq=freq, shell=ishell)
                    fnu = fnus[find_nearest_index(self.get_ej_skymap_times(), time)]
                    all_fluxes_arr = np.array(int_i)
                    delta_x = np.array(xrs_i).max() - np.array(xrs_i).min()
                    delta_y = np.array(yrs_i).max() - np.array(yrs_i).min()
                    dfnu = fnu / (delta_x * delta_y)
                    if verbose:
                        print("Ejecta shell {}".format(ishell))
                        print("\tFnu/mas^2 = {:.2e} mJy/mas^2".format(dfnu))
                        print("\tall_x = [{:.2e}, {:.2e}]".format(np.array(xrs_i).min(), np.array(xrs_i).max()))
                        print("\tall_y = [{:.2e}, {:.2e}]".format(np.array(yrs_i).min(), np.array(yrs_i).max()))
                        print("\tDelta_x = {:.2f}, Delta_y = {:.2f}]".format(delta_x, delta_y))
                        print("\tFnu/mas^2 = {:.2e} mJy/mas^2".format(dfnu))
                    int_i = (all_fluxes_arr / all_fluxes_arr.max()) * dfnu
                    fnus_tot = int_i + fnus

                return (xrs_i, yrs_i, int_i)

                # plt.figure(figsize=(4.6, 3.2))
                # plt.semilogy(mu_i, int_i / int_i.max(), '.')
                # # plt.semilogy(mu_i, int_i / int_i.max(), 'x', color='red')
                # # plt.semilogy(mu_i,int_i2/int_i2.max(), 'o', color='red')
                # plt.xlabel(
                #     r"$\mu=\sin(\theta_{\rm obs}) \sin(\theta) \sin(\phi) + \cos(\theta_{\rm obs}) \cos(\theta)$")
                # plt.ylabel(r"$I/I_{\rm max}$")
                # plt.tight_layout()
                # plt.show()

                # nx = np.complex(0, nx)
                # ny = np.complex(0, ny)
                # grid_x, grid_y = np.mgrid[xrs_i.min():xrs_i.max():nx, yrs_i.min():yrs_i.max():ny]
                # i_zz = interp(xrs_i, yrs_i, int_i, grid_x, grid_y, 'linear')
                # return (grid_x, grid_y, i_zz, xrs_i, yrs_i, int_i)

                # if int_or_hist == "int":
                #     nx = np.complex(0, nx)
                #     ny = np.complex(0, ny)
                #     grid_x, grid_y = np.mgrid[xrs_i.min() * 1.2:xrs_i.max() * 1.2:nx,
                #                      yrs_i.min() * 1.2:yrs_i.max() * 1.2:ny]
                #     i_zz = interp(xrs_i, yrs_i, int_i, grid_x, grid_y, 'linear')
                #     return (grid_x, grid_y, i_zz, xrs_i, yrs_i, int_i)
                # else:
                #     nx = 100
                #     ny = 100
                #     nx = np.complex(0, nx + 1)
                #     ny = np.complex(0, ny + 1)
                #     edges_x = np.mgrid[xrs_i.min() * 1.2:xrs_i.max() * 1.2:nx]
                #     edges_y = np.mgrid[yrs_i.min() * 1.2:yrs_i.max() * 1.2:ny]
                #     # grid_x, grid_y = np.mgrid[xrs_i.min() * 1.2:xrs_i.max() * 1.2:nx,
                #     #                  yrs_i.min() * 1.2:yrs_i.max() * 1.2:ny]
                #     i_zz, _ = np.histogramdd(tuple([xrs_i, yrs_i]), bins=tuple([edges_x, edges_y]), weights=int_i)
                #     grid_x = 0.5 * (edges_x[1:] + edges_x[:-1])
                #     grid_y = 0.5 * (edges_y[1:] + edges_y[:-1])
                #     return (grid_x, grid_y, i_zz, xrs_i, yrs_i, int_i)
                #     # print(i_zz.shape)
            else:
                i_shells = []
                # total_fluxes = []
                for ish in range(nshells):
                    xrs_i = np.array(ddfile["xrs"][ish]) * cgs.rad2mas / d_l  # m -> mas
                    yrs_i = np.array(ddfile["yrs"][ish]) * cgs.rad2mas / d_l  # m -> mas
                    int_i = np.array(ddfile["intensity"][ish]) * (d_l ** 2 / cgs.rad2mas ** 2)  # -> mJy / mas^2
                    if (np.sum(int_i) == 0):
                        continue
                    i_shells.append(ish)
                if remove_mu:
                    print("Removing 'mu' from ejecta skymap")
                all_xrs, all_yrs, all_zz = [], [], []
                for ii, ish in enumerate(i_shells):
                    mu_i = np.array(ddfile["mu"][ish])
                    xrs_i = np.array(ddfile["xrs"][ish]) * cgs.rad2mas / d_l  # m -> mas
                    yrs_i = np.array(ddfile["yrs"][ish]) * cgs.rad2mas / d_l  # m -> mas
                    int_i = np.array(ddfile["intensity"][ish]) * ( d_l ** 2 / cgs.rad2mas ** 2 )  # * dfile.attrs["d_L"] ** 2
                    if remove_mu:
                        # idx1 = abs(mu_i) < min_mu # TODO this is overritten!
                        # int_i[idx1] =
                        # print("len(gam[idx1])={}".format(gam_i[idx1]))
                        # idx = int_i > int_i.max() * min_mu_frac
                        # idx = idx1
                        # if verbose:
                        #     if (len(mu_i[idx1]) != len(mu_i[idx])):
                        #         print('\t', mu_i[idx1])
                        #         print('\t', mu_i[idx])
                        #         # exit(1)
                        # if verbose:
                        #     if len(idx) < 2:
                        #         print("No excess = {}".format(ii))
                        #     print("len(gam[idx])={}".format(gam_i[idx]))
                        #     print("Gamma={}".format(gam_i[idx]))
                        #     print("beta={}".format(get_beta(gam_i[idx])))
                        #     print("mu_i={}".format(mu_i[idx]))
                        #     print("r_i={}".format(r_i[idx]))
                        #     print("B={}".format(B_i[idx]))
                        #     print("tb_i={}".format(tb_i[idx]))
                        # print(np.abs(mu_i).min())
                        # print(theta_i[np.abs(mu_i) < 1e-5])
                        # print(phi_i[np.abs(mu_i) < 1e-5])
                        # print(np.abs(mu_i).min())

                        # int_i[mu_i < 1e-3] = 0.
                        int_i *= abs( mu_i)  # TODO I was produced as F / (R^2 abs(mu)), where abs(mu)->0 and I->inf. Problem!!!
                        # * np.sqrt(1 - mu_i**2)#/ ((1-mu_i)*(1+mu_i))
                        # plt.figure(figsize=(4.6,3.2))
                        # plt.semilogy(mu_i,int_i/int_i.max(), '.')
                        # if len(mu_i[idx]) > 0: plt.semilogy(mu_i[idx],int_i[idx]/int_i[idx].max(), 'x', color='red')
                        # # plt.semilogy(mu_i,int_i2/int_i2.max(), 'o', color='red')
                        # plt.xlabel(r"$\mu=\sin(\theta_{\rm obs}) \sin(\theta) \sin(\phi) + \cos(\theta_{\rm obs}) \cos(\theta)$")
                        # plt.ylabel(r"$I/I_{\rm max}$")
                        # plt.tight_layout()
                        # plt.show()
                        # int_i[idx] = 0.
                    all_xrs.append(xrs_i)
                    all_yrs.append(yrs_i)
                    all_zz.append(int_i)
                # Problem: changing nlayers changes the image/image, Fnu per pixel; Solution:
                if renormalize:
                    print("Renormalizing ejecta skymap (shell by shell separately)")
                    fnus_tot = np.zeros_like(self.get_ej_skymap_times())
                    for i, i_ish in enumerate(i_shells):
                        fnus = self.get_ej_skymap_totfluxes(freq=freq, shell=i_ish)
                        fnu = fnus[find_nearest_index(self.get_ej_skymap_times(), time)]
                        all_fluxes_arr = np.array(all_zz[i])
                        delta_x = np.array(all_xrs[i]).max() - np.array(all_xrs[i]).min()
                        delta_y = np.array(all_yrs[i]).max() - np.array(all_yrs[i]).min()
                        dfnu = fnu / (delta_x * delta_y)
                        if verbose:
                            print("SHELL {}".format(i_ish))
                            print("\tfnu = {:.2e} ".format(fnu))
                            print("\tall_x = [{:.2e}, {:.2e}]".format(np.array(all_xrs).min(), np.array(all_xrs).max()))
                            print("\tall_y = [{:.2e}, {:.2e}]".format(np.array(all_yrs).min(), np.array(all_yrs).max()))
                            print("\tDelta_x = {:.2f}, Delta_y = {:.2f}]".format(delta_x, delta_y))
                            print("\tFnu/mas^2 = {:.2e} mJy/mas^2".format(dfnu))
                        all_zz[i] = (all_fluxes_arr / all_fluxes_arr.max()) * dfnu
                        fnus_tot = fnus_tot + fnus

                return (all_xrs, all_yrs, all_zz)

                # assess what is the combined grid extend (for all images)
                # xmin_neg, xmin_pos, xmax, xmin = [], [], [], []
                # ymin_neg, ymin_pos, ymax, ymin = [], [], [], []
                # i_min, i_max = [], []
                # i_shells = []
                # for ish in range(nshells):
                #     xrs_i = np.array(ddfile["xrs"][ish]) * cgs.rad2mas / d_l  # m -> mas
                #     yrs_i = np.array(ddfile["yrs"][ish]) * cgs.rad2mas / d_l  # m -> mas
                #     int_i = np.array(ddfile["intensity"][ish]) * (d_l ** 2 / cgs.rad2mas ** 2)  # -> mJy / mas^2
                #     if (np.sum(int_i) == 0):
                #         continue
                #     i_shells.append(ish)
                #     xmin_neg.append(xrs_i[xrs_i < 0].min())
                #     xmin_pos.append(xrs_i[xrs_i > 0].min())
                #     xmax.append(xrs_i.max())
                #     xmin.append(xrs_i.min())
                #     ymin_neg.append(yrs_i[yrs_i < 0].min())
                #     ymin_pos.append(yrs_i[yrs_i > 0].min())
                #     ymax.append(yrs_i.max())
                #     ymin.append(yrs_i.min())
                #     i_min.append(int_i.min())
                #     i_max.append(int_i.max())
                # xmin_neg = np.array(xmin_neg)
                # xmin_pos = np.array(xmin_pos)
                # xmax = np.array(xmax)
                # xmin = np.array(xmin)
                # ymin_neg = np.array(ymin_neg)
                # ymin_pos = np.array(ymin_pos)
                # ymax = np.array(ymax)
                # ymin = np.array(ymin)
                # i_min = np.array(i_min)
                # i_max = np.array(i_max)
                #
                # if verbose:
                #     for i in range(len(i_shells)):
                #         print("\tshell={} xmin_neg={:.2e} xmin_pos={:.2e} xmax={:.2e}".format(i, xmin_neg[i], xmin_pos[i], xmax[i]))
                #         print("\t         ymin_neg={:.2e} ymin_pos={:.2e} ymax={:.2e}".format(i, ymin_neg[i], ymin_pos[i],  ymax[i]))
                #     print("\toverall X min(xmin_pos) = {:.2e} max(xmax) = {:.2e}".format(xmin_pos.min(), xmax.max()))
                #     print("\toverall Y min(ymin_pos) = {:.2e} max(ymax) = {:.2e}".format(ymin_pos.min(), ymax.max()))
                #     print("\toverall I min(i_min)    = {:.2e} max(ymax) = {:.2e}".format(i_min.min(), i_max.max()))
                #
                # edges_x = np.linspace(xmin.min(), xmax.max(), num=nx)
                # edges_y = np.linspace(ymin.min(), ymax.max(), num=ny)
                #
                #
                # x_grid = np.concatenate((-1. * np.logspace(np.log10(xmin_pos.min()), np.log10(xmax.max()), nx)[::-1],
                #                          np.logspace(np.log10(xmin_pos.min()), np.log10(xmax.max()), nx)))
                # y_grid = np.concatenate((-1. * np.logspace(np.log10(ymin_pos.min()), np.log10(ymax.max()), ny)[::-1],
                #                          np.logspace(np.log10(ymin_pos.min()), np.log10(ymax.max()), ny)))
                #
                # # edges_x = np.concatenate((-1. * np.logspace(np.log10(xmin_pos.min()), np.log10(xmax.max()), 101)[::-1],
                # #                          np.logspace(np.log10(xmin_pos.min()), np.log10(xmax.max()), nx)))
                # # edges_y = np.concatenate((-1. * np.logspace(np.log10(ymin_pos.min()), np.log10(ymax.max()), 101)[::-1],
                # #                          np.logspace(np.log10(ymin_pos.min()), np.log10(ymax.max()), ny)))
                # # plt.plot(edges_x, edges_y, marker='.', ls='none')
                # # plt.show()
                #
                # if verbose:
                #     print("\tGrid xmin_neg={:.2e} xmin_pos={:.2e} xmax={:.2e}".format(x_grid.min(),
                #                                                                     x_grid[x_grid > 0].min(),
                #                                                                     x_grid.max()))
                #     print("\tGrid ymin_neg={:.2e} ymin_pos={:.2e} ymax={:.2e}".format(y_grid.min(),
                #                                                                     y_grid[y_grid > 0].min(),
                #                                                                     y_grid.max()))
                # xx_grid, yy_grid = np.meshgrid(x_grid, y_grid)
                # if int_or_hist == "int":
                #     zz = np.zeros_like((xx_grid))
                # else:
                #     zz = np.zeros((len(edges_x)-1,len(edges_y)-1))
                # # interpolate onto the grid that covers all images
                # all_xrs, all_yrs, all_zz = [], [], []
                # for ii, ish in enumerate(i_shells):  # range(len(i_shells))
                #     if verbose: print("Pocessing: shell={} [{}/{}]".format(ish, ii, len(i_shells)))
                #     r_i = np.array(ddfile["r"][ish])
                #     mu_i = np.array(ddfile["mu"][ish])
                #     xrs_i = np.array(ddfile["xrs"][ish]) * cgs.rad2mas / d_l  # m -> mas
                #     yrs_i = np.array(ddfile["yrs"][ish]) * cgs.rad2mas / d_l  # m -> mas
                #     int_i = np.array(ddfile["intensity"][ish]) * (d_l ** 2 / cgs.rad2mas ** 2) # * dfile.attrs["d_L"] ** 2
                #     gam_i = np.array(ddfile["gamma"][ish])
                #     B_i = np.array(ddfile["B"][ish])
                #     tb_i = np.array(ddfile["tburst"][ish])
                #     theta_i = np.array(ddfile["theta"][ish])
                #     theta_ij = np.array(ddfile["theta_j"][ish])
                #     theta0_i = np.array(ddfile["theta0"][ish])
                #     phi_i = np.array(ddfile["phi"][ish])
                #
                #
                #     # plt.plot(mu_i, theta0_i, marker='.', ls='none')
                #     # plt.show()
                #     # int_i[mu_i < 5e-2] = 0.
                #     # int_i *= abs(mu_i)
                #
                #     if remove_mu:
                #         # idx1 = abs(mu_i) < min_mu # TODO this is overritten!
                #         # int_i[idx1] =
                #         # print("len(gam[idx1])={}".format(gam_i[idx1]))
                #         # idx = int_i > int_i.max() * min_mu_frac
                #         # idx = idx1
                #         # if verbose:
                #         #     if (len(mu_i[idx1]) != len(mu_i[idx])):
                #         #         print('\t', mu_i[idx1])
                #         #         print('\t', mu_i[idx])
                #         #         # exit(1)
                #         # if verbose:
                #         #     if len(idx) < 2:
                #         #         print("No excess = {}".format(ii))
                #         #     print("len(gam[idx])={}".format(gam_i[idx]))
                #         #     print("Gamma={}".format(gam_i[idx]))
                #         #     print("beta={}".format(get_beta(gam_i[idx])))
                #         #     print("mu_i={}".format(mu_i[idx]))
                #         #     print("r_i={}".format(r_i[idx]))
                #         #     print("B={}".format(B_i[idx]))
                #         #     print("tb_i={}".format(tb_i[idx]))
                #         # print(np.abs(mu_i).min())
                #         # print(theta_i[np.abs(mu_i) < 1e-5])
                #         # print(phi_i[np.abs(mu_i) < 1e-5])
                #         # print(np.abs(mu_i).min())
                #
                #         # int_i[mu_i < 1e-3] = 0.
                #         int_i *= abs(mu_i) # TODO I was produced as F / (R^2 abs(mu)), where abs(mu)->0 and I->inf. Problem!!!
                #         # * np.sqrt(1 - mu_i**2)#/ ((1-mu_i)*(1+mu_i))
                #         # plt.figure(figsize=(4.6,3.2))
                #         # plt.semilogy(mu_i,int_i/int_i.max(), '.')
                #         # if len(mu_i[idx]) > 0: plt.semilogy(mu_i[idx],int_i[idx]/int_i[idx].max(), 'x', color='red')
                #         # # plt.semilogy(mu_i,int_i2/int_i2.max(), 'o', color='red')
                #         # plt.xlabel(r"$\mu=\sin(\theta_{\rm obs}) \sin(\theta) \sin(\phi) + \cos(\theta_{\rm obs}) \cos(\theta)$")
                #         # plt.ylabel(r"$I/I_{\rm max}$")
                #         # plt.tight_layout()
                #         # plt.show()
                #         # int_i[idx] = 0.
                #     # int_i /= abs(mu_i)[abs(mu_i)>0.01]
                #     # else:
                #     #     int_i2 = int_i
                #
                #     if int_or_hist == "int":
                #         # nx = np.complex(0, nx)
                #         # ny = np.complex(0, ny)
                #         # grid_x, grid_y = np.mgrid[xrs_i.min() * 1.2:xrs_i.max() * 1.2:nx,
                #         #                  yrs_i.min() * 1.2:yrs_i.max() * 1.2:ny]
                #         i_zz = interp(xrs_i, yrs_i, int_i, xx_grid, yy_grid, 'linear')
                #         #return (grid_x, grid_y, i_zz, xrs_i, yrs_i, int_i)
                #     else:
                #         # nx = 100
                #         # ny = 100
                #         # nx = np.complex(0, nx + 1)
                #         # ny = np.complex(0, ny + 1)
                #         # edges_x = np.mgrid[xrs_i.min() * 1.2:xrs_i.max() * 1.2:nx]
                #         # edges_y = np.mgrid[yrs_i.min() * 1.2:yrs_i.max() * 1.2:ny]
                #         # grid_x, grid_y = np.mgrid[xrs_i.min() * 1.2:xrs_i.max() * 1.2:nx,
                #         #                  yrs_i.min() * 1.2:yrs_i.max() * 1.2:ny]
                #         i_zz, _ = np.histogramdd(tuple([xrs_i, yrs_i]), bins=tuple([edges_x, edges_y]), weights=int_i)
                #         # grid_x = 0.5 * (edges_x[1:] + edges_x[:-1])
                #         # grid_y = 0.5 * (edges_y[1:] + edges_y[:-1])
                #         #return (grid_x, grid_y, i_zz, xrs_i, yrs_i, int_i)
                #
                #     # i_zz = interp(xrs_i, yrs_i, int_i, xx_grid, yy_grid, 'linear')
                #     zz += i_zz
                #     all_xrs.append(xrs_i)
                #     all_yrs.append(yrs_i)
                #     all_zz.append(int_i)
                # if verbose:
                #     print("\tAfter interpolation: I [{:.2e}, {:.2e}] Sum = {:.2e} [Total expected {:.2e}]"
                #           .format(zz.min(), zz.max(), np.sum(zz), float(np.array(dfile["totalflux at freq={:.4e}".format(freq)])[find_nearest_index(times, time)])))
                #
                # all_xrs = np.concatenate(all_xrs)
                # all_yrs = np.concatenate(all_yrs)
                # all_zz = np.concatenate(all_zz)
                # if int_or_hist == "int":
                #
                #     return (xx_grid, yy_grid, zz, all_xrs, all_yrs, all_zz)
                # else:
                #     grid_x = 0.5 * (edges_x[1:] + edges_x[:-1])
                #     grid_y = 0.5 * (edges_y[1:] + edges_y[:-1])
                #     return (grid_x, grid_y, zz, all_xrs, all_yrs, all_zz)
        else:
            raise NotImplementedError("Not implemented")



    def get_ej_skymap_spec_idx(self, time=None, freq=None, ishell=None, verbose=False, remove_mu=False):

        freqs = self.get_ej_skymap_freqs()
        assert len(freqs) > 1
        idx = find_nearest_index(freqs, freq)
        dfile = self.get_ej_skymap_obj()
        nshells = int(dfile.attrs["nshells"])

        all_xrs, all_yrs, all_zz = \
            self.get_ej_skymap(time=time, freq=freqs[idx], ishell=ishell, verbose=verbose, remove_mu=remove_mu)

        all_xrs_m1, all_yrs_m1, all_zz_m1 = \
            self.get_ej_skymap(time=time, freq=freqs[idx], ishell=ishell, verbose=verbose, remove_mu=remove_mu)

        assert len(all_zz) == len(all_zz_m1)

        all_values = []
        for i in range(len(all_zz_m1)):
            ffreq_i = np.full(all_zz[i].shape, freqs[idx])
            ffreq_im1 = np.full(all_zz_m1[i].shape, freqs[idx - 1])
            num = np.log10(all_zz[i]) - np.log10(all_zz_m1[i])
            denum = np.log10(ffreq_i) - np.log10(ffreq_im1)
            values = 1 * num / denum
            all_values.append(values)

        return (all_xrs, all_yrs, all_values)

    # def get_ej_skymap(self, time=None, freq=None, ishell=None, verbose=False, remove_mu=False, int_or_hist="int"):
    #
    #
    #     nx = 200
    #     ny = 100
    #     min_mu = 1e-4 # limit for mu that if too small blow up the intensity for some reason...
    #     min_mu_frac = 0.85
    #     # limit_data = True
    #     times = self.get_ej_skymap_times()
    #     freqs = self.get_ej_skymap_freqs()
    #     dfile = self.get_ej_skymap_obj()
    #     nshells = int(dfile.attrs["nshells"])
    #     d_l = float(self.get_ej_skymap_obj().attrs["d_L"])
    #     if ((not time is None)and(not time in times)):
    #         raise ValueError("time={} is not in the list for skypams={}".format(time, times))
    #     if((not freq is None)and(not freq in freqs)):
    #         raise ValueError("freq={} is not in the list for skypams={}".format(freq, freqs))
    #     if ((not ishell is None)and(ishell > nshells-1)):
    #         raise ValueError("shell={} is beying skymap nshells={}".format(ishell, nshells))
    #     if ((not time is None) and (not freq is None)):
    #         print(dfile.keys())
    #         ddfile = dfile["time={:.4e} freq={:.4e}".format(time, freq)]
    #         if (not ishell is None):
    #             r_i = np.array(ddfile["r"][ishell])
    #             mu_i = np.array(ddfile["mu"][ishell ])
    #             xrs_i = np.array(ddfile["xrs"][ishell]) * cgs.rad2mas / d_l # m -> mas
    #             yrs_i = np.array(ddfile["yrs"][ishell]) * cgs.rad2mas / d_l # m -> mas
    #             int_i = np.array(ddfile["intensity"][ishell]) * (d_l ** 2 / cgs.rad2mas ** 2) #  -> mJy / mas^2
    #             gam_i = np.array(ddfile["gamma"][ishell])
    #             B_i = np.array(ddfile["B"][ishell])
    #             tb_i = np.array(ddfile["tburst"][ishell])
    #             theta_i = np.array(ddfile["theta"][ishell])
    #             phi_i = np.array(ddfile["phi"][ishell])
    #
    #             if remove_mu:
    #                 int_i *= np.abs(mu_i) # TODO I was produced as F / (R^2 abs(mu)), where abs(mu)->0 and I->inf. Problem!!!
    #
    #             # plt.figure(figsize=(4.6, 3.2))
    #             # plt.semilogy(mu_i, int_i / int_i.max(), '.')
    #             # # plt.semilogy(mu_i, int_i / int_i.max(), 'x', color='red')
    #             # # plt.semilogy(mu_i,int_i2/int_i2.max(), 'o', color='red')
    #             # plt.xlabel(
    #             #     r"$\mu=\sin(\theta_{\rm obs}) \sin(\theta) \sin(\phi) + \cos(\theta_{\rm obs}) \cos(\theta)$")
    #             # plt.ylabel(r"$I/I_{\rm max}$")
    #             # plt.tight_layout()
    #             # plt.show()
    #
    #             # nx = np.complex(0, nx)
    #             # ny = np.complex(0, ny)
    #             # grid_x, grid_y = np.mgrid[xrs_i.min():xrs_i.max():nx, yrs_i.min():yrs_i.max():ny]
    #             # i_zz = interp(xrs_i, yrs_i, int_i, grid_x, grid_y, 'linear')
    #             # return (grid_x, grid_y, i_zz, xrs_i, yrs_i, int_i)
    #
    #             if int_or_hist == "int":
    #                 nx = np.complex(0, nx)
    #                 ny = np.complex(0, ny)
    #                 grid_x, grid_y = np.mgrid[xrs_i.min() * 1.2:xrs_i.max() * 1.2:nx,
    #                                  yrs_i.min() * 1.2:yrs_i.max() * 1.2:ny]
    #                 i_zz = interp(xrs_i, yrs_i, int_i, grid_x, grid_y, 'linear')
    #                 return (grid_x, grid_y, i_zz, xrs_i, yrs_i, int_i)
    #             else:
    #                 nx = 100
    #                 ny = 100
    #                 nx = np.complex(0, nx + 1)
    #                 ny = np.complex(0, ny + 1)
    #                 edges_x = np.mgrid[xrs_i.min() * 1.2:xrs_i.max() * 1.2:nx]
    #                 edges_y = np.mgrid[yrs_i.min() * 1.2:yrs_i.max() * 1.2:ny]
    #                 # grid_x, grid_y = np.mgrid[xrs_i.min() * 1.2:xrs_i.max() * 1.2:nx,
    #                 #                  yrs_i.min() * 1.2:yrs_i.max() * 1.2:ny]
    #                 i_zz, _ = np.histogramdd(tuple([xrs_i, yrs_i]), bins=tuple([edges_x, edges_y]), weights=int_i)
    #                 grid_x = 0.5 * (edges_x[1:] + edges_x[:-1])
    #                 grid_y = 0.5 * (edges_y[1:] + edges_y[:-1])
    #                 return (grid_x, grid_y, i_zz, xrs_i, yrs_i, int_i)
    #                 # print(i_zz.shape)
    #
    #         else:
    #             # assess what is the combined grid extend (for all images)
    #             xmin_neg, xmin_pos, xmax, xmin = [], [], [], []
    #             ymin_neg, ymin_pos, ymax, ymin = [], [], [], []
    #             i_min, i_max = [], []
    #             i_shells = []
    #             for ish in range(nshells):
    #                 xrs_i = np.array(ddfile["xrs"][ish]) * cgs.rad2mas / d_l  # m -> mas
    #                 yrs_i = np.array(ddfile["yrs"][ish]) * cgs.rad2mas / d_l  # m -> mas
    #                 int_i = np.array(ddfile["intensity"][ish]) * (d_l ** 2 / cgs.rad2mas ** 2)  # -> mJy / mas^2
    #                 if (np.sum(int_i) == 0):
    #                     continue
    #                 i_shells.append(ish)
    #                 xmin_neg.append(xrs_i[xrs_i < 0].min())
    #                 xmin_pos.append(xrs_i[xrs_i > 0].min())
    #                 xmax.append(xrs_i.max())
    #                 xmin.append(xrs_i.min())
    #                 ymin_neg.append(yrs_i[yrs_i < 0].min())
    #                 ymin_pos.append(yrs_i[yrs_i > 0].min())
    #                 ymax.append(yrs_i.max())
    #                 ymin.append(yrs_i.min())
    #                 i_min.append(int_i.min())
    #                 i_max.append(int_i.max())
    #             xmin_neg = np.array(xmin_neg)
    #             xmin_pos = np.array(xmin_pos)
    #             xmax = np.array(xmax)
    #             xmin = np.array(xmin)
    #             ymin_neg = np.array(ymin_neg)
    #             ymin_pos = np.array(ymin_pos)
    #             ymax = np.array(ymax)
    #             ymin = np.array(ymin)
    #             i_min = np.array(i_min)
    #             i_max = np.array(i_max)
    #
    #             if verbose:
    #                 for i in range(len(i_shells)):
    #                     print("\tshell={} xmin_neg={:.2e} xmin_pos={:.2e} xmax={:.2e}".format(i, xmin_neg[i], xmin_pos[i], xmax[i]))
    #                     print("\t         ymin_neg={:.2e} ymin_pos={:.2e} ymax={:.2e}".format(i, ymin_neg[i], ymin_pos[i],  ymax[i]))
    #                 print("\toverall X min(xmin_pos) = {:.2e} max(xmax) = {:.2e}".format(xmin_pos.min(), xmax.max()))
    #                 print("\toverall Y min(ymin_pos) = {:.2e} max(ymax) = {:.2e}".format(ymin_pos.min(), ymax.max()))
    #                 print("\toverall I min(i_min)    = {:.2e} max(ymax) = {:.2e}".format(i_min.min(), i_max.max()))
    #
    #             edges_x = np.linspace(xmin.min(), xmax.max(), num=nx)
    #             edges_y = np.linspace(ymin.min(), ymax.max(), num=ny)
    #
    #
    #             x_grid = np.concatenate((-1. * np.logspace(np.log10(xmin_pos.min()), np.log10(xmax.max()), nx)[::-1],
    #                                      np.logspace(np.log10(xmin_pos.min()), np.log10(xmax.max()), nx)))
    #             y_grid = np.concatenate((-1. * np.logspace(np.log10(ymin_pos.min()), np.log10(ymax.max()), ny)[::-1],
    #                                      np.logspace(np.log10(ymin_pos.min()), np.log10(ymax.max()), ny)))
    #
    #             # edges_x = np.concatenate((-1. * np.logspace(np.log10(xmin_pos.min()), np.log10(xmax.max()), 101)[::-1],
    #             #                          np.logspace(np.log10(xmin_pos.min()), np.log10(xmax.max()), nx)))
    #             # edges_y = np.concatenate((-1. * np.logspace(np.log10(ymin_pos.min()), np.log10(ymax.max()), 101)[::-1],
    #             #                          np.logspace(np.log10(ymin_pos.min()), np.log10(ymax.max()), ny)))
    #             # plt.plot(edges_x, edges_y, marker='.', ls='none')
    #             # plt.show()
    #
    #             if verbose:
    #                 print("\tGrid xmin_neg={:.2e} xmin_pos={:.2e} xmax={:.2e}".format(x_grid.min(),
    #                                                                                 x_grid[x_grid > 0].min(),
    #                                                                                 x_grid.max()))
    #                 print("\tGrid ymin_neg={:.2e} ymin_pos={:.2e} ymax={:.2e}".format(y_grid.min(),
    #                                                                                 y_grid[y_grid > 0].min(),
    #                                                                                 y_grid.max()))
    #             xx_grid, yy_grid = np.meshgrid(x_grid, y_grid)
    #             if int_or_hist == "int":
    #                 zz = np.zeros_like((xx_grid))
    #             else:
    #                 zz = np.zeros((len(edges_x)-1,len(edges_y)-1))
    #             # interpolate onto the grid that covers all images
    #             all_xrs, all_yrs, all_zz = [], [], []
    #             for ii, ish in enumerate(i_shells):  # range(len(i_shells))
    #                 if verbose: print("Pocessing: shell={} [{}/{}]".format(ish, ii, len(i_shells)))
    #                 r_i = np.array(ddfile["r"][ish])
    #                 mu_i = np.array(ddfile["mu"][ish])
    #                 xrs_i = np.array(ddfile["xrs"][ish]) * cgs.rad2mas / d_l  # m -> mas
    #                 yrs_i = np.array(ddfile["yrs"][ish]) * cgs.rad2mas / d_l  # m -> mas
    #                 int_i = np.array(ddfile["intensity"][ish]) * (d_l ** 2 / cgs.rad2mas ** 2) # * dfile.attrs["d_L"] ** 2
    #                 gam_i = np.array(ddfile["gamma"][ish])
    #                 B_i = np.array(ddfile["B"][ish])
    #                 tb_i = np.array(ddfile["tburst"][ish])
    #                 theta_i = np.array(ddfile["theta"][ish])
    #                 theta_ij = np.array(ddfile["theta_j"][ish])
    #                 theta0_i = np.array(ddfile["theta0"][ish])
    #                 phi_i = np.array(ddfile["phi"][ish])
    #
    #
    #                 # plt.plot(mu_i, theta0_i, marker='.', ls='none')
    #                 # plt.show()
    #                 # int_i[mu_i < 5e-2] = 0.
    #                 # int_i *= abs(mu_i)
    #
    #                 if remove_mu:
    #                     # idx1 = abs(mu_i) < min_mu # TODO this is overritten!
    #                     # int_i[idx1] =
    #                     # print("len(gam[idx1])={}".format(gam_i[idx1]))
    #                     # idx = int_i > int_i.max() * min_mu_frac
    #                     # idx = idx1
    #                     # if verbose:
    #                     #     if (len(mu_i[idx1]) != len(mu_i[idx])):
    #                     #         print('\t', mu_i[idx1])
    #                     #         print('\t', mu_i[idx])
    #                     #         # exit(1)
    #                     # if verbose:
    #                     #     if len(idx) < 2:
    #                     #         print("No excess = {}".format(ii))
    #                     #     print("len(gam[idx])={}".format(gam_i[idx]))
    #                     #     print("Gamma={}".format(gam_i[idx]))
    #                     #     print("beta={}".format(get_beta(gam_i[idx])))
    #                     #     print("mu_i={}".format(mu_i[idx]))
    #                     #     print("r_i={}".format(r_i[idx]))
    #                     #     print("B={}".format(B_i[idx]))
    #                     #     print("tb_i={}".format(tb_i[idx]))
    #                     # print(np.abs(mu_i).min())
    #                     # print(theta_i[np.abs(mu_i) < 1e-5])
    #                     # print(phi_i[np.abs(mu_i) < 1e-5])
    #                     # print(np.abs(mu_i).min())
    #
    #                     # int_i[mu_i < 1e-3] = 0.
    #                     int_i *= abs(mu_i) # TODO I was produced as F / (R^2 abs(mu)), where abs(mu)->0 and I->inf. Problem!!!
    #                     # * np.sqrt(1 - mu_i**2)#/ ((1-mu_i)*(1+mu_i))
    #                     # plt.figure(figsize=(4.6,3.2))
    #                     # plt.semilogy(mu_i,int_i/int_i.max(), '.')
    #                     # if len(mu_i[idx]) > 0: plt.semilogy(mu_i[idx],int_i[idx]/int_i[idx].max(), 'x', color='red')
    #                     # # plt.semilogy(mu_i,int_i2/int_i2.max(), 'o', color='red')
    #                     # plt.xlabel(r"$\mu=\sin(\theta_{\rm obs}) \sin(\theta) \sin(\phi) + \cos(\theta_{\rm obs}) \cos(\theta)$")
    #                     # plt.ylabel(r"$I/I_{\rm max}$")
    #                     # plt.tight_layout()
    #                     # plt.show()
    #                     # int_i[idx] = 0.
    #                 # int_i /= abs(mu_i)[abs(mu_i)>0.01]
    #                 # else:
    #                 #     int_i2 = int_i
    #
    #                 if int_or_hist == "int":
    #                     # nx = np.complex(0, nx)
    #                     # ny = np.complex(0, ny)
    #                     # grid_x, grid_y = np.mgrid[xrs_i.min() * 1.2:xrs_i.max() * 1.2:nx,
    #                     #                  yrs_i.min() * 1.2:yrs_i.max() * 1.2:ny]
    #                     i_zz = interp(xrs_i, yrs_i, int_i, xx_grid, yy_grid, 'linear')
    #                     #return (grid_x, grid_y, i_zz, xrs_i, yrs_i, int_i)
    #                 else:
    #                     # nx = 100
    #                     # ny = 100
    #                     # nx = np.complex(0, nx + 1)
    #                     # ny = np.complex(0, ny + 1)
    #                     # edges_x = np.mgrid[xrs_i.min() * 1.2:xrs_i.max() * 1.2:nx]
    #                     # edges_y = np.mgrid[yrs_i.min() * 1.2:yrs_i.max() * 1.2:ny]
    #                     # grid_x, grid_y = np.mgrid[xrs_i.min() * 1.2:xrs_i.max() * 1.2:nx,
    #                     #                  yrs_i.min() * 1.2:yrs_i.max() * 1.2:ny]
    #                     i_zz, _ = np.histogramdd(tuple([xrs_i, yrs_i]), bins=tuple([edges_x, edges_y]), weights=int_i)
    #                     # grid_x = 0.5 * (edges_x[1:] + edges_x[:-1])
    #                     # grid_y = 0.5 * (edges_y[1:] + edges_y[:-1])
    #                     #return (grid_x, grid_y, i_zz, xrs_i, yrs_i, int_i)
    #
    #                 # i_zz = interp(xrs_i, yrs_i, int_i, xx_grid, yy_grid, 'linear')
    #                 zz += i_zz
    #                 all_xrs.append(xrs_i)
    #                 all_yrs.append(yrs_i)
    #                 all_zz.append(int_i)
    #             if verbose:
    #                 print("\tAfter interpolation: I [{:.2e}, {:.2e}] Sum = {:.2e} [Total expected {:.2e}]"
    #                       .format(zz.min(), zz.max(), np.sum(zz), float(np.array(dfile["totalflux at freq={:.4e}".format(freq)])[find_nearest_index(times, time)])))
    #
    #             all_xrs = np.concatenate(all_xrs)
    #             all_yrs = np.concatenate(all_yrs)
    #             all_zz = np.concatenate(all_zz)
    #             if int_or_hist == "int":
    #
    #                 return (xx_grid, yy_grid, zz, all_xrs, all_yrs, all_zz)
    #             else:
    #                 grid_x = 0.5 * (edges_x[1:] + edges_x[:-1])
    #                 grid_y = 0.5 * (edges_y[1:] + edges_y[:-1])
    #                 return (grid_x, grid_y, zz, all_xrs, all_yrs, all_zz)
    #     else:
    #         raise NotImplementedError("Not implemented")
    @staticmethod
    def get_skymap_lat_dist(all_x, all_y, all_fluxes, collapse_axis="y", fac=1.0, nx=1000, ny=1000, extend=2):

        all_x = np.concatenate(all_x)
        all_y = np.concatenate(all_y)
        all_fluxes = np.concatenate(all_fluxes)

        # nx, ny = len(all_x), len(all_y)
        nx = np.complex(0, nx + 1)
        ny = np.complex(0, ny + 1)
        # edges_x = np.mgrid[all_x.min() * 1.0:all_x.max() * 1.0:nx]
        # edges_y = np.mgrid[all_y.min() * 1.0:all_y.max() * 1.0:ny]
        #
        # grid_x = 0.5 * (edges_x[1:] + edges_x[:-1])
        # grid_y = 0.5 * (edges_y[1:] + edges_y[:-1])

        # grid_X, grid_Y = np.meshgrid(grid_x,grid_y,indexing="ij")

        # image = interp(all_x, all_y, all_fluxes, grid_X, grid_Y, method='linear')
        fac = 1.0
        grid_x, grid_y = np.mgrid[all_x.min()*extend:all_x.max()*extend:nx,
                         all_y.min()*extend:all_y.max()*extend:ny]

        # nx = np.complex(0, nx + 1)
        # ny = np.complex(0, ny + 1)
        # edges_x = np.mgrid[all_x.min() * extend:all_y.max() * extend:nx]
        # edges_y = np.mgrid[all_x.min() * extend:all_y.max() * extend:ny]

        image = interpolate.griddata(np.array([all_x, all_y]).T * fac, all_fluxes,
                                     (grid_x * fac, grid_y * fac), method='linear', fill_value=0)
        latAvDist, latAvDist2, latMaxDist = lateral_distributions(grid_x, grid_y, image, collapse_axis=collapse_axis)
        if collapse_axis == "y":
            return (grid_x[:, 0], latAvDist, latAvDist2, latMaxDist)
        else:
            return (grid_y[0, :], latAvDist, latAvDist2, latMaxDist)

        # int_x, int_y, int_zz, all_x, all_y, all_fluxes \
        #     = pb.get_ej_skymap(time=time * cgs.day, freq=freq, verbose=False, remove_mu=True, int_or_hist="hist")
        if axis == "x":
            # nx = 200
            # ny = 100
            nx = np.complex(0, nx + 1)
            ny = np.complex(0, ny + 1)
            edges_x = np.mgrid[all_x.min() * 1.2:all_x.max() * 1.2:nx]
            edges_y = np.mgrid[all_y.min() * 1.2:all_y.max() * 1.2:ny]
            # grid_x, grid_y = np.mgrid[xrs_i.min() * 1.2:xrs_i.max() * 1.2:nx,
            #                  yrs_i.min() * 1.2:yrs_i.max() * 1.2:ny]
            i_zz, _ = np.histogramdd(tuple([all_x, all_y]), bins=tuple([edges_x, edges_y]), weights=all_fluxes)
            grid_x = 0.5 * (edges_x[1:] + edges_x[:-1])
            grid_y = 0.5 * (edges_y[1:] + edges_y[:-1])
            # tmp = []
            # for i in range(len(grid_x)):
            #     mask = i_zz[i, :] >  1e-3*i_zz[i, :].max()
            #     tmp.append(np.sum( (i_zz[i, :] * np.abs(np.diff(edges_y)))[mask] )  / np.sum(np.diff(edges_y)[mask]) )
            #     # tmp.append(np.average())
            # y_sum_ii = np.array(tmp)

            y_sum_ii = np.sum(i_zz * np.diff(edges_y), axis=1) / np.sum(np.diff(edges_y))
            max_ii = np.max(y_sum_ii)
            x1 = grid_x[np.argmin(y_sum_ii < max_ii * fac)]
            x2 = grid_x[::-1][np.argmin(y_sum_ii[::-1] < max_ii * fac)]
            # assert x2 >= x1
            return (x1, x2, grid_x, y_sum_ii)
        elif axis == "z":
            # nx = 200
            # ny = 100
            nx = np.complex(0, nx + 1)
            ny = np.complex(0, ny + 1)
            edges_x = np.mgrid[all_x.min() * 1.2:all_x.max() * 1.2:nx]
            edges_y = np.mgrid[all_y.min() * 1.2:all_y.max() * 1.2:ny]
            # grid_x, grid_y = np.mgrid[xrs_i.min() * 1.2:xrs_i.max() * 1.2:nx,
            #                  yrs_i.min() * 1.2:yrs_i.max() * 1.2:ny]
            i_zz, _ = np.histogramdd(tuple([all_x, all_y]), bins=tuple([edges_x, edges_y]), weights=all_fluxes)
            grid_x = 0.5 * (edges_x[1:] + edges_x[:-1])
            grid_y = 0.5 * (edges_y[1:] + edges_y[:-1])
            tmp = []
            for i in range(len(grid_y)):
                tmp.append(np.sum(i_zz[:, i] * np.diff(edges_x)) / np.sum(np.diff(edges_x)))
            tmp = np.array(tmp)

            x_sum_ii = np.sum(i_zz, axis=0)
            max_ii = np.max(x_sum_ii)
            y1 = grid_y[np.argmin(x_sum_ii < max_ii * fac)]
            y2 = grid_y[::-1][np.argmin(x_sum_ii[::-1] < max_ii * fac)]
            assert y2 > y1
            return (y1, y2, grid_y, x_sum_ii)
        else:
            raise KeyError()

    #
    #
    # def _get_lc_for_freq(self, freq, times, freqs):
    #     _times, _fluxes = [], []
    #     for i in range(len(times)):
    #         if (times[i] == freq):
    #             _times.append(times[i])
    #             _fluxes.append(freqs[i])
    #     _times = np.array(_times)
    #     _fluxes = np.array(_fluxes)
    #     if (len(_times) < 1):
    #         raise ValueError("Freq={} not found in data. Avaialable={}".format(freq, set(freqs)))
    #     return (times, freqs)#(_times, _fluxes / 1e3 / 1e23)  # [s, mJy -> CGS]
    # def _check_if_loaded_ej_lc(self):
    #     if (self.ej_lc is None):
    #         self.ej_lc = h5py.File(self.res_dir + self.ejecta_prefix + self.def_name_lc)
    #
    # def get_ej_lc_2d_arr(self, freq=None, ishell=None, ilayer=None):
    #     dfile = self.get_ej_lc_obj()
    #     v_n = list( dfile.keys() )[0]
    #     times_ejecta = np.array(dfile["times"])
    #     freqs_ejecta = np.array(dfile["freqs"])
    #     fluxes_ejecta = np.array(dfile[v_n]["totalflux"])
    #     if (freq is None and ishell is None and ilayer is None):
    #         return fluxes_ejecta / 1e3 / 1e23
    #     elif (ishell is None and ilayer is None):
    #         return self._get_lc_for_freq(freq, times_ejecta, freqs_ejecta)
    #
    #     layer = "shell={} layer={}".format(ishell, ilayer)
    #     if (not layer in list(dfile.keys())):
    #         raise NameError("Layer {} (aka '{}') is not in the ejecta comov.spec. file.\n Available: {}"
    #                         .format(ilayer, layer, dfile.keys()))
    #     if (not v_n in dfile[layer].keys()):
    #         raise NameError("v_n {} is not in the ejecta comov.spec. dfile[{}].keys() \n Avaialble: {}"
    #                         .format(v_n, layer, dfile[layer].keys()))
    #     if (freq is None):
    #         return times_ejecta, fluxes_ejecta / 1e3 / 1e23  # [s, mJy -> CGS]
    #     else:
    #         _times, _fluxes = [], []
    #         for i in range(len(times_ejecta)):
    #             if (freqs_ejecta[i] == freq):
    #                 _times.append(times_ejecta[i])
    #                 _fluxes.append(fluxes_ejecta[i])
    #         _times = np.array(_times)
    #         _fluxes = np.array(_fluxes)
    #         if (len(_times) < 1):
    #             raise ValueError("Freq={} not found in data. Avaialable={}".format(freq, set(freqs_ejecta)))
    #         return (_times, _fluxes / 1e3 / 1e23)  # [s, mJy -> CGS]
    #
    #     return np.array(dfile[layer][v_n])
    @staticmethod
    def get_skymap_fwhm(grid, latAvDist, cc, fac=0.5):
        assert len(grid) == len(latAvDist)
        val = latAvDist[find_nearest_index(grid, cc)]
        val =  np.interp(cc, grid, latAvDist)# latAvDist[find_nearest_index(grid, cc)]
        max_val = np.max(latAvDist)
        x1 = grid[np.argmin(latAvDist < val * fac)]
        # x1 = np.interp(val * fac, latAvDist, grid) #  grid[np.argmin(latAvDist < val * fac)]
        # plt.plot(grid, latAvDist)
        # plt.axhline(y=val * fac)
        # plt.show()
        x2 = grid[::-1][np.argmin(latAvDist[::-1] < val * fac)]
        return (x1, x2)

    @staticmethod
    def smooth_interpolated_skymap_with_gaussian_kernel(i_zz, type="uniform", sigma=3):
        if type=="uniform":
            i_zz = ndimage.uniform_filter(i_zz, size=sigma )
        elif type == "gaussian":
            if i_zz.ndim == 1:
                i_zz = ndimage.filters.gaussian_filter(i_zz, sigma=[sigma], mode='reflect')
            elif i_zz.ndim == 2:
                i_zz = ndimage.filters.gaussian_filter(i_zz, sigma=[sigma,sigma], mode='reflect')
        else:
            raise KeyError("smooth type is not recognize: {}".format(type))
        return i_zz

def tmp_for_file2(time, freq, grp, settings, plot_dict, **kwargs):

    '''
    RETURN :: (xc, yc, x1, x2, y1, y2, int_x, int_y, int_zz)
    '''

    if(len(settings["kn_grb_skymap"].keys()) > 0):
        if (len(settings["kn_skymap"].keys()) == 0):
            raise KeyError("For 'kn_grb_skymap' kn_skymap has to be set")
        if (len(settings["grb_skymap"].keys()) == 0):
            raise KeyError("For 'kn_grb_skymap' grb_skymap has to be set")

    xc, yc, x1, x2, y1, y2, int_x, int_y, int_zz = 0., 0., 0., 0., 0., 0., None, None, None

    if not "precompute" in settings.keys():
        print("setting 'precompute' to False")
        settings["precompute"] = False

    ''' ---- size/cm of the ejecta ---- '''
    pb = BPA_METHODS()
    pb.res_dir_kn = plot_dict["res_dir_kn"]
    pb.ejecta_prefix = plot_dict["ejecta_prefix1"]
    pb.res_dir_grb = plot_dict["res_dir_grb"]
    pb.jet_prefix = plot_dict["jet_prefix"]

    pb_w = BPA_METHODS()
    pb_w.res_dir_kn = plot_dict["res_dir_kn"]
    pb_w.ejecta_prefix = plot_dict["ejecta_prefix2"]
    pb_w.res_dir_grb = plot_dict["res_dir_grb"]
    pb_w.jet_prefix = plot_dict["jet_prefix"]

    if type(time) is str:
        if time == "tp1":
            t_lc_max = pb.get_ej_lc_times()[find_nearest_index(pb.get_ej_lc_totalflux(freq), pb.get_ej_lc_totalflux(freq).max())]
            t_sm_where_lc_max = pb.get_ej_skymap_times()[find_nearest_index(pb.get_ej_skymap_times(), t_lc_max)]
            # idx =find_nearest_index(pb.get_ej_skymap_totfluxes(freq), pb.get_ej_lc_totalflux(freq).max())
            # time = pb.get_ej_skymap_times()[find_nearest_index(pb.get_ej_skymap_totfluxes(freq), pb.get_ej_lc_totalflux(freq).max() )]
            # time_lc = pb.get_ej_lc_times()[find_nearest_index(pb.get_ej_lc_totalflux(freq),pb.get_ej_lc_totalflux(freq).max())]
            print("Selected nearest time {} [d] while the light curve peak time is {} [d]".format(t_sm_where_lc_max/cgs.day,t_lc_max/cgs.day))
        elif time == "tp2":
            t_lc_max = pb_w.get_ej_lc_times()[ find_nearest_index(pb_w.get_ej_lc_totalflux(freq), pb_w.get_ej_lc_totalflux(freq).max())]
            t_sm_where_lc_max = pb_w.get_ej_skymap_times()[find_nearest_index(pb_w.get_ej_skymap_times(), t_lc_max)]
            # time = pb.get_ej_skymap_times()[ find_nearest_index(pb.get_ej_skymap_totfluxes(freq), pb.get_ej_lc_totalflux(freq).max())]
            # time_lc = pb.get_ej_lc_times()[ find_nearest_index(pb.get_ej_lc_totalflux(freq), pb.get_ej_lc_totalflux(freq).max())]
            print("Selected nearest time {} [d] while the light curve peak time is {} [d]".format(t_sm_where_lc_max/cgs.day,t_lc_max/cgs.day))
        elif time == "tpgrb":
            t_lc_max = pb.get_jet_lc_times()[find_nearest_index(pb.get_jet_lc_totalflux(freq), pb.get_jet_lc_totalflux(freq).max())]
            t_sm_where_lc_max = pb.get_jet_skymap_times()[find_nearest_index(pb.get_jet_skymap_times(), t_lc_max)]
            # time = pb.get_ej_skymap_times()[ find_nearest_index(pb.get_ej_skymap_totfluxes(freq), pb.get_jet_lc_totalflux(freq).max())]
            # time_lc = pb.get_ej_lc_times()[ find_nearest_index(pb.get_ej_lc_totalflux(freq), pb.get_jet_lc_totalflux(freq).max())]
            print("Selected nearest time {} [d] while the light curve peak time is {} [d]".format(t_sm_where_lc_max/cgs.day,t_lc_max/cgs.day))
        else:
            raise KeyError("time can be a number, 'tp1', 'tp2' or 'tpj'. Given:{}".format(time))

        time = t_sm_where_lc_max / cgs.day

    if (len(settings["kn_skymap"].keys()) > 0) or (len(settings["kn_grb_skymap"].keys()) > 0):
        tmp = copy.deepcopy(settings["kn_skymap"])
        key = "kn nx={} ny={}".format(tmp["hist_nx"], tmp["hist_ny"])
        if settings["precompute"]:
            all_x, all_y, all_fluxes \
                = pb.get_ej_skymap(time=time * cgs.day, freq=freq, verbose=False, remove_mu=True)

            if tmp["spec"]:
                int_x, int_y, int_zz = pb.get_combined_ej_spectral_map(time=time, freq=freq,
                                                                       nx=tmp["hist_nx"], ny=tmp["hist_ny"], extend=2)
            else:
                int_x, int_y, int_zz = pb.combine_images(all_x, all_y, all_fluxes, hist_or_int="hist",
                                                         shells=True, nx=tmp["hist_nx"], ny=tmp["hist_ny"], extend=2)

            grid_y, _, i_zz_y, _ = pb.get_skymap_lat_dist(all_x, all_y, all_fluxes, collapse_axis="x",
                                                          nx=tmp["hist_nx"], ny=tmp["hist_ny"])
            grid_x, _, i_zz_x, _ = pb.get_skymap_lat_dist(all_x, all_y, all_fluxes, collapse_axis="y",
                                                          nx=tmp["hist_nx"], ny=tmp["hist_ny"])
            xc, yc = pb.get_ej_skymap_cm(all_x, all_y, all_fluxes)
            #
            if key in grp.keys():
                raise KeyError("key: {} already exists: {}".format(key, grp.keys()))
            grp_kn = grp.create_group(key)
            grp_kn.create_dataset("int_x", data=int_x)
            grp_kn.create_dataset("int_y", data=int_y)
            grp_kn.create_dataset("int_zz", data=int_zz)
            grp_kn.create_dataset("grid_y", data=grid_y)
            grp_kn.create_dataset("i_zz_y", data=i_zz_y)
            grp_kn.create_dataset("grid_x", data=grid_x)
            grp_kn.create_dataset("i_zz_x", data=i_zz_x)
            grp_kn.attrs.create("xc", data=xc)
            grp_kn.attrs.create("yc", data=yc)
        else:
            if not key in grp.keys():
                raise KeyError("key={} is not in the list for a group: {}".format(key, grp.keys()))
            grp_kn = grp[key]
            int_x = np.array(np.array(grp_kn["int_x"]))
            int_y = np.array(np.array(grp_kn["int_y"]))
            int_zz = np.array(np.array(grp_kn["int_zz"]))  # plotting
            grid_y = np.array(np.array(grp_kn["grid_y"]))
            i_zz_y = np.array(np.array(grp_kn["i_zz_y"]))
            grid_x = np.array(np.array(grp_kn["grid_x"]))
            i_zz_x = np.array(np.array(grp_kn["i_zz_x"]))
            xc = float(grp_kn.attrs["xc"])
            yc = float(grp_kn.attrs["yc"])
        i_zz_x /= i_zz_x.max()
        i_zz_y /= i_zz_y.max()
        if (int_zz.max() == int_zz.min()):
            raise ValueError("i_zz_x.max() == i_zz_x.min() = {} \n plot_dic=\n{}".format(int_zz.max(), plot_dict))
        y1, y2 = pb.get_skymap_fwhm(grid_y, i_zz_y, yc)
        x1, x2 = pb.get_skymap_fwhm(grid_x, i_zz_x, xc)
        if (not (x2 > x1)) and (x2!=0) and (x1!=0) :
            raise ValueError("kn_skymap failed computing width x2={} <= x1={} Should be other way around"
                             .format(x2, x1))
        # if len(tmp["cm"].keys()) > 0: ax.plot(xc, yc, **tmp["cm"])  # color='red', marker="o")
        # if len(tmp["ysize"].keys()) > 0: ax.errorbar([xc, xc], [y1, y2], xerr=[int_x.max() / 10, int_x.max() / 10],
        #                                              **tmp["ysize"])
        # if len(tmp["xsize"].keys()) > 0: ax.errorbar([x1, x2], [yc, yc], yerr=[int_y.max() / 10, int_y.max() / 10],
        #                                              **tmp["xsize"])
        # # --------------------
        # if len(tmp["pcolormesh"].keys()) > 0:
        #     # levels = np.geomspace(int_ration[~np.isnan(int_ration)].min(), int_ration[~np.isnan(int_ration)].max(), 50)
        #     int_zz = int_zz.T
        #     int_zz[~np.isfinite(int_zz)] = 0.
        #     if (len(tmp["smooth"].keys()) > 0):
        #         int_zz = pb.smooth_interpolated_skymap_with_gaussian_kernel(i_zz=int_zz, type=tmp["smooth"]["type"])
        #     im = plot_pcolomesh(ax=ax, task_dict=tmp["pcolormesh"], int_x=int_x, int_y=int_y, int_zz=int_zz)

        # return (xc, yc, x1, x2, y1, y2, int_x, int_y, int_zz)
        if (len(kwargs.keys())>0) and ("kn_skymap" in kwargs.keys()):
            return (xc, yc, x1, x2, y1, y2, int_x, int_y, int_zz)

    ''' ---- size/cm of the grb ---- '''
    if  ("grb_skymap" in settings.keys()) and (len(settings["grb_skymap"].keys()) > 0) or (len(settings["kn_grb_skymap"].keys()) > 0):
        tmp = copy.deepcopy(settings["grb_skymap"])
        key = "grb nx={} ny={}".format(tmp["hist_nx"], tmp["hist_ny"])
        if settings["precompute"]:
            all_x_jet, all_y_jet, all_fluxes_jet \
                = pb.get_jet_skymap(time=time * cgs.day, freq=freq, verbose=False, remove_mu=True)
            if tmp["spec"]:
                int_x_j, int_y_j, int_zz_j = pb.get_combined_jet_spectral_map(time=time, freq=freq,
                                                                              nx=tmp["hist_nx"], ny=tmp["hist_ny"], extend=2)
            else:
                int_x_j, int_y_j, int_zz_j = pb.combine_images([all_x_jet], [all_y_jet], [all_fluxes_jet],
                                                               hist_or_int="hist", shells=False,
                                                               nx=tmp["hist_nx"], ny=tmp["hist_ny"], extend=2)
            grid_y_j, _, i_zz_y_j, _ = pb.get_skymap_lat_dist([all_x_jet], [all_y_jet], [all_fluxes_jet],
                                                              collapse_axis="x", nx=tmp["hist_nx"], ny=tmp["hist_ny"])
            grid_x_j, _, i_zz_x_j, _ = pb.get_skymap_lat_dist([all_x_jet], [all_y_jet], [all_fluxes_jet],
                                                              collapse_axis="y", nx=tmp["hist_nx"], ny=tmp["hist_ny"])
            xc_m_j, yc_m_j = pb.get_ej_skymap_cm([all_x_jet], [all_y_jet], [all_fluxes_jet])
            #
            if key in grp.keys():
                raise KeyError("key: {} already exists: {}".format(key, grp.keys()))
            grp_j = grp.create_group(key)
            grp_j.create_dataset("int_x", data=int_x_j)
            grp_j.create_dataset("int_y", data=int_y_j)
            grp_j.create_dataset("int_zz", data=int_zz_j)
            grp_j.create_dataset("grid_y", data=grid_y_j)
            grp_j.create_dataset("i_zz_y", data=i_zz_y_j)
            grp_j.create_dataset("grid_x", data=grid_x_j)
            grp_j.create_dataset("i_zz_x", data=i_zz_x_j)
            grp_j.attrs.create("xc", data=xc_m_j)
            grp_j.attrs.create("yc", data=yc_m_j)
        else:
            grp_j = grp[key]
            int_x_j = np.array(np.array(grp_j["int_x"]))
            int_y_j = np.array(np.array(grp_j["int_y"]))
            int_zz_j = np.array(np.array(grp_j["int_zz"]))  # plotting
            grid_y_j = np.array(np.array(grp_j["grid_y"]))
            i_zz_y_j = np.array(np.array(grp_j["i_zz_y"]))
            grid_x_j = np.array(np.array(grp_j["grid_x"]))
            i_zz_x_j = np.array(np.array(grp_j["i_zz_x"]))
            xc_m_j = float(grp_j.attrs["xc"])
            yc_m_j = float(grp_j.attrs["yc"])
        if (int_zz_j.max() == int_zz_j.min()):
            raise ValueError("i_zz_x.max() == i_zz_x.min() = {} \n plot_dic=\n{}".format(int_zz_j.max(), plot_dict))
        x1_j, x2_j = pb.get_skymap_fwhm(grid_x_j, i_zz_x_j, xc_m_j)
        y1_j, y2_j = pb.get_skymap_fwhm(grid_y_j, i_zz_y_j, yc_m_j)
        if (not (x2_j > x1_j)) and (x2_j != 0) and (x1_j != 0):
            raise ValueError("kn_skymap failed computing width x2={} <= x1={} Should be other way around"
                             .format(x2_j, x1_j))
        # if len(tmp["cm"].keys()) > 0: ax.plot(xc_m_j, yc_m_j, **tmp["cm"])  # color='red', marker="o")
        # if len(tmp["ysize"].keys()) > 0: ax.errorbar([xc_m_j, xc_m_j], [y1_j, y2_j],
        #                                              xerr=[int_x_j.max() / 10, int_x_j.max() / 10], **tmp["ysize"])
        # if len(tmp["xsize"].keys()) > 0: ax.errorbar([x1_j, x2_j], [yc_m_j, yc_m_j],
        #                                              yerr=[int_y_j.max() / 10, int_y_j.max() / 10], **tmp["xsize"])
        # # --------------------
        # if len(tmp["pcolormesh"].keys()) > 0:
        #     int_zz_j = int_zz_j.T
        #     int_zz_j[~np.isfinite(int_zz_j)] = 0.
        #     im = plot_pcolomesh(ax=ax, task_dict=tmp["pcolormesh"], int_x=int_x_j, int_y=int_y_j, int_zz=int_zz_j)
        xc, yc, x1, x2, y1, y2, int_x, int_y, int_zz = \
            xc_m_j, yc_m_j, x1_j, x2_j, y1_j, y2_j, int_x_j, int_y_j, int_zz_j
        # return (xc_m_j, yc_m_j, x1_j, x2_j, y1_j, y2_j, int_x_j, int_y_j, int_zz_j)
        if (len(kwargs.keys())>0) and ("grb_skymap" in kwargs.keys()):
            return (xc, yc, x1, x2, y1, y2, int_x, int_y, int_zz)

    ''' --- size/cm of the kn with grb (using previous two) ---'''
    if ("kn_grb_skymap" in settings.keys()) and len(settings["kn_grb_skymap"].keys()) > 0:
        assert (len(settings["grb_skymap"].keys()) > 0)
        assert (len(settings["kn_skymap"].keys()) > 0)
        tmp = copy.deepcopy(settings["kn_grb_skymap"])
        key = "kn_and_grb nx={} ny={}".format(tmp["hist_nx"], tmp["hist_ny"])
        if settings["precompute"]:
            all_x_pj = [all_x_jet] + all_x
            all_y_pj = [all_y_jet] + all_y
            all_fluxes_pj = [all_fluxes_jet] + all_fluxes
            if tmp["spec"]:
                int_x_pj, int_y_pj, int_zz_pj = pb.get_combined_kn_grb_spectral_map(time=time, freq=freq,
                                                                                    nx=tmp["hist_nx"], ny=tmp["hist_ny"], extend=2)
            else:
                int_x_pj, int_y_pj, int_zz_pj = \
                    pb.combine_images(all_x_pj, all_y_pj, all_fluxes_pj, hist_or_int="hist", shells=False,
                                      nx=tmp["hist_nx"], ny=tmp["hist_ny"], extend=2)
            # _, _, int_zz_wjet_nojet = pb.combine_images(all_x_pj, all_y_pj, all_fluxes_pj, hist_or_int="hist",
            # shells=False, nx=tmp["hist_nx"], ny=tmp["hist_ny"])
            xc_m_pj, yc_m_pj = pb.get_jet_skymap_cm(all_x_pj, all_y_pj, all_fluxes_pj)
            grid_y_pj, _, i_zz_y_pj, _ = pb.get_skymap_lat_dist(all_x_pj, all_y_pj, all_fluxes_pj,
                                                                collapse_axis="x", nx=tmp["hist_nx"], ny=tmp["hist_ny"])
            grid_x_pj, _, i_zz_x_pj, _ = pb.get_skymap_lat_dist(all_x_pj, all_y_pj, all_fluxes_pj,
                                                                collapse_axis="y", nx=tmp["hist_nx"], ny=tmp["hist_ny"])
            if key in grp.keys():
                raise KeyError("key: {} already exists: {}".format(key, grp.keys()))
            grp_kn_plus_grb = grp.create_group(key)
            # grp_kn_plus_grb = grp.create_group("kn_and_grb nx={} ny={}".format(tmp["hist_nx"], tmp["hist_ny"]))
            grp_kn_plus_grb.create_dataset("int_x", data=int_x_pj)
            grp_kn_plus_grb.create_dataset("int_y", data=int_y_pj)
            grp_kn_plus_grb.create_dataset("int_zz", data=int_zz_pj)
            grp_kn_plus_grb.create_dataset("grid_y", data=grid_y_pj)
            grp_kn_plus_grb.create_dataset("i_zz_y", data=i_zz_y_pj)
            grp_kn_plus_grb.create_dataset("grid_x", data=grid_x_pj)
            grp_kn_plus_grb.create_dataset("i_zz_x", data=i_zz_x_pj)
            grp_kn_plus_grb.attrs.create("xc", data=xc_m_pj)
            grp_kn_plus_grb.attrs.create("yc", data=yc_m_pj)
        else:
            grp_kn_plus_grb = grp[key]
            int_x_pj = np.array(np.array(grp_kn_plus_grb["int_x"]))
            int_y_pj = np.array(np.array(grp_kn_plus_grb["int_y"]))
            int_zz_pj = np.array(np.array(grp_kn_plus_grb["int_zz"]))  # plotting
            grid_y_pj = np.array(np.array(grp_kn_plus_grb["grid_y"]))
            i_zz_y_pj = np.array(np.array(grp_kn_plus_grb["i_zz_y"]))
            grid_x_pj = np.array(np.array(grp_kn_plus_grb["grid_x"]))
            i_zz_x_pj = np.array(np.array(grp_kn_plus_grb["i_zz_x"]))
            xc_m_pj = float(grp_kn_plus_grb.attrs["xc"])
            yc_m_pj = float(grp_kn_plus_grb.attrs["yc"])
        if (int_zz_pj.max() == int_zz_pj.min()):
            raise ValueError("i_zz_x.max() == i_zz_x.min() = {} \n plot_dic=\n{}".format(int_zz_pj.max(), plot_dict))
        x1_pj, x2_pj = pb.get_skymap_fwhm(grid_x_pj, i_zz_x_pj, xc_m_pj)
        y1_pj, y2_pj = pb.get_skymap_fwhm(grid_y_pj, i_zz_y_pj, yc_m_pj)
        if (not (x2_pj > x1_pj)) and (x2_pj != 0) and (x1_pj != 0):
            raise ValueError("kn_grb_skymap failed computing width x2={} <= x1={} Should be other way around"
                             .format(x2_pj, x1_pj))
        # if len(tmp["cm"].keys()) > 0: ax.plot(xc_m_pj, yc_m_pj, **tmp["cm"])
        # if len(tmp["ysize"].keys()) > 0: ax.errorbar([xc_m_pj, xc_m_pj], [y1_pj, y2_pj],
        #                                              xerr=[int_x_pj.max() / 10, int_x_pj.max() / 10], **tmp["ysize"])
        # if len(tmp["xsize"].keys()) > 0: ax.errorbar([x1_pj, x2_pj], [yc_m_pj, yc_m_pj],
        #                                              yerr=[int_y_pj.max() / 10, int_y_pj.max() / 10], **tmp["xsize"])
        # # --------------------
        # if len(tmp["pcolormesh"].keys()) > 0:
        #     # levels = np.geomspace(int_ration[~np.isnan(int_ration)].min(), int_ration[~np.isnan(int_ration)].max(), 50)
        #     int_zz_pj = int_zz_pj.T
        #     int_zz_pj[~np.isfinite(int_zz_pj)] = 0.
        #     im = plot_pcolomesh(ax=ax, task_dict=tmp["pcolormesh"], int_x=int_x_pj, int_y=int_y_pj, int_zz=int_zz_pj)

        xc, yc, x1, x2, y1, y2, int_x, int_y, int_zz = \
            xc_m_pj, yc_m_pj, x1_pj, x2_pj, y1_pj, y2_pj, int_x_pj, int_y_pj, int_zz_pj

        if (len(kwargs.keys())>0) and ("kn_grb_skymap" in kwargs.keys()):
            return (xc, yc, x1, x2, y1, y2, int_x, int_y, int_zz)
        # return (xc_m_pj, yc_m_pj, x1_pj, x2_pj, y1_pj, y2_pj, int_x_pj, int_y_pj, int_zz_pj)

    ''' ---- size/cm of the ejecta with interaction ---- '''
    if ("kn_w_skymap" in settings.keys()) and len(settings["kn_w_skymap"].keys()) > 0:
        tmp = copy.deepcopy(settings["kn_w_skymap"])
        key = "kn_w nx={} ny={}".format(tmp["hist_nx"], tmp["hist_ny"])
        if settings["precompute"]:
            all_x_w, all_y_w, all_fluxes_w \
                = pb_w.get_ej_skymap(time=time * cgs.day, freq=freq, verbose=False, remove_mu=True)
            int_x_w, int_y_w, int_zz_w = pb_w.combine_images(all_x_w, all_y_w, all_fluxes_w, hist_or_int="hist",
                                                             shells=True, nx=tmp["hist_nx"], ny=tmp["hist_ny"],
                                                             extend=2)
            grid_y_w, _, i_zz_y_w, _ = pb_w.get_skymap_lat_dist(all_x_w, all_y_w, all_fluxes_w, collapse_axis="x",
                                                                nx=tmp["hist_nx"], ny=tmp["hist_ny"])
            grid_x_w, _, i_zz_x_w, _ = pb_w.get_skymap_lat_dist(all_x_w, all_y_w, all_fluxes_w, collapse_axis="y",
                                                                nx=tmp["hist_nx"], ny=tmp["hist_ny"])
            xc_w, yc_w = pb_w.get_ej_skymap_cm(all_x_w, all_y_w, all_fluxes_w)
            #
            if key in grp.keys():
                raise KeyError("key: {} already exists: {}".format(key, grp.keys()))
            grp_kn_w = grp.create_group(key)
            # grp_kn_w = grp.create_group("kn_w nx={} ny={}".format(tmp["hist_nx"], tmp["hist_ny"]))
            grp_kn_w.create_dataset("int_x", data=int_x_w)
            grp_kn_w.create_dataset("int_y", data=int_y_w)
            grp_kn_w.create_dataset("int_zz", data=int_zz_w)
            grp_kn_w.create_dataset("grid_y", data=grid_y_w)
            grp_kn_w.create_dataset("i_zz_y", data=i_zz_y_w)
            grp_kn_w.create_dataset("grid_x", data=grid_x_w)
            grp_kn_w.create_dataset("i_zz_x", data=i_zz_x_w)
            grp_kn_w.attrs.create("xc", data=xc_w)
            grp_kn_w.attrs.create("yc", data=yc_w)
        else:
            grp_kn_w = grp[key]
            int_x_w = np.array(np.array(grp_kn_w["int_x"]))
            int_y_w = np.array(np.array(grp_kn_w["int_y"]))
            int_zz_w = np.array(np.array(grp_kn_w["int_zz"]))  # plotting
            grid_y_w = np.array(np.array(grp_kn_w["grid_y"]))
            i_zz_y_w = np.array(np.array(grp_kn_w["i_zz_y"]))
            grid_x_w = np.array(np.array(grp_kn_w["grid_x"]))
            i_zz_x_w = np.array(np.array(grp_kn_w["i_zz_x"]))
            xc_w = float(grp_kn_w.attrs["xc"])
            yc_w = float(grp_kn_w.attrs["yc"])
        if (int_zz_w.max() == int_zz_w.min()):
            raise ValueError("i_zz_x.max() == i_zz_x.min() = {} \n plot_dic=\n{}".format(int_zz_w.max(), plot_dict))
        i_zz_x_w /= i_zz_x_w.max()
        i_zz_y_w /= i_zz_y_w.max()
        y1_w, y2_w = pb_w.get_skymap_fwhm(grid_y_w, i_zz_y_w, yc_w)
        x1_w, x2_w = pb_w.get_skymap_fwhm(grid_x_w, i_zz_x_w, xc_w)
        if (not (x2_w > x1_w)) and (x2_w != 0) and (x1_w != 0):
            raise ValueError("kn_w_skymap failed computing width x2={} <= x1={} Should be other way around"
                             .format(x2_w, x1_w))
        # if len(tmp["cm"].keys()) > 0: ax.plot(xc_w, yc_w, **tmp["cm"])
        # if len(tmp["ysize"].keys()) > 0: ax.errorbar([xc_w, xc_w], [y1_w, y2_w],
        #                                              xerr=[int_x_w.max() / 10, int_x_w.max() / 10], **tmp["ysize"])
        # if len(tmp["xsize"].keys()) > 0: ax.errorbar([x1_w, x2_w], [yc_w, yc_w],
        #                                              yerr=[int_y_w.max() / 10, int_y_w.max() / 10], **tmp["xsize"])
        # # --------------------
        # if len(tmp["pcolormesh"].keys()) > 0:
        #     int_zz_w = int_zz_w.T
        #     int_zz_w[~np.isfinite(int_zz_w)] = 0.
        #     im = plot_pcolomesh(ax=ax, task_dict=tmp["pcolormesh"], int_x=int_x_w, int_y=int_y_w, int_zz=int_zz_w)

        xc, yc, x1, x2, y1, y2, int_x, int_y, int_zz = \
            xc_w, yc_w, x1_w, x2_w, y1_w, y2_w, int_x_w, int_y_w, int_zz_w
        # return (xc_w, yc_w, x1_w, x2_w, y1_w, y2_w, int_x_w, int_y_w, int_zz_w)
        if (len(kwargs.keys())>0) and ("kn_w_skymap" in kwargs.keys()):
            return (xc, yc, x1, x2, y1, y2, int_x, int_y, int_zz)

    ''' --- relation --- '''
    if ("kn_skymap_ratio" in settings.keys()) and len(settings["kn_skymap_ratio"].keys()) > 0:
        assert (len(settings["kn_skymap"].keys()) > 0)
        assert (len(settings["kn_w_skymap"].keys()) > 0)
        tmp = copy.deepcopy(settings["kn_skymap_ratio"])
        key = "kn_ratio nx={} ny={}".format(tmp["hist_nx"], tmp["hist_ny"])
        if settings["precompute"]:
            int_x_w, int_y_w, int_zz, int_zz_w = \
                pb.get_combained_ej_skymaps_adjusted_to_other(time=time, freq=freq, other_pb_instance=pb_w,
                                                              nx=tmp["hist_nx"], ny=tmp["hist_ny"])
            if key in grp.keys():
                raise KeyError("key: {} already exists: {}".format(key, grp.keys()))
            kn_adjusted = grp.create_group(key)
            # kn_adjusted = grp.create_group("kn_adjusted nx={} ny={}".format(tmp["hist_nx"], tmp["hist_ny"]))
            kn_adjusted.create_dataset("int_x", data=int_x_w)
            kn_adjusted.create_dataset("int_y", data=int_y_w)
            kn_adjusted.create_dataset("int_zz", data=int_zz)
            kn_adjusted.create_dataset("int_zz_w", data=int_zz_w)
        else:
            kn_adjusted = grp[key]
            int_x_w = np.array(np.array(kn_adjusted["int_x"]))
            int_y_w = np.array(np.array(kn_adjusted["int_y"]))
            int_zz = np.array(np.array(kn_adjusted["int_zz"]))
            int_zz_w = np.array(np.array(kn_adjusted["int_zz_w"]))
        if (int_zz_w.max() == int_zz_w.min()):
            raise ValueError("i_zz_x.max() == i_zz_x.min() = {} \n plot_dic=\n{}".format(int_zz_w.max(), plot_dict))
        int_ration = int_zz_w / int_zz
        int_ration = int_ration.T
        int_ration[~np.isfinite(int_ration)] = 1.

        # norm = Normalize(int_ration[np.isfinite(int_ration)].min(), int_ration[np.isfinite(int_ration)].max())
        # # levels = np.geomspace(int_ration[~np.isnan(int_ration)].min(), int_ration[~np.isnan(int_ration)].max(), 50)
        # cmap = cm.get_cmap('viridis')
        # int_ration = int_ration.T
        # int_ration[~np.isfinite(int_ration)] = 0.
        # im = ax.pcolormesh(int_x_w, int_y_w, int_ration, norm=norm, cmap=cmap, alpha=1.0)
        # im.set_rasterized(True)
        # ax.set_facecolor(cmap(0.0))

        # vmin = tmp["pcolormesh"]["vmin"]
        # vmax = tmp["pcolormesh"]["vmax"]
        # vcenter = tmp["pcolormesh"]["vcenter"]
        # levels = MaxNLocator(nbins=40).tick_values(int_ration.min(), int_ration.max())
        # # levels = MaxNLocator(nbins=40).tick_values(-5, 1)
        # cmap = plt.get_cmap(tmp["pcolormesh"]["cmap"])
        # norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
        # norm = colors.TwoSlopeNorm(vmin=vmin, vcenter=vcenter, vmax=vmax)
        # # im = aux_ax2.contourf(theta, r, values, levels=np.logspace(np.log10(1e-4),np.log10(1e0),20), cmap=cmap, norm=norm)
        # # im = aux_ax2.contourf(theta, r, values, levels=np.linspace(val.min(),val.max(),20), cmap=cmap, norm=norm)
        # im = ax.pcolormesh(int_x_w, int_y_w, int_ration, cmap=cmap, norm=norm, alpha=1.0)
        # im.set_rasterized(True)

        # im = plot_pcolomesh(ax=ax, task_dict=tmp["pcolormesh"], int_x=int_x_w, int_y=int_y_w, int_zz=int_ration)
        int_x, int_y, int_zz = \
            int_x_w, int_y_w, int_ration
        # return (None, None, None, None, None, None, int_x_w, int_y_w, int_ration)
        if (len(kwargs.keys())>0) and ("kn_skymap_ratio" in kwargs.keys()):
            return (xc, yc, x1, x2, y1, y2, int_x, int_y, int_zz)

    return (xc, yc, x1, x2, y1, y2, int_x, int_y, int_zz)

''' parallel runs '''

class PBA_PARALLEL:
    def __init__(self, workind_dir, list_parfile_fnames):
        assert os.path.isdir(workind_dir)
        assert len(list_parfile_fnames) > 0
        for fname in list_parfile_fnames:
            assert os.path.isfile(workind_dir+fname)
        self.working_dir = workind_dir
        self.list_fnames = list_parfile_fnames
    def __call__(self, idx):
        if (idx > (len(self.list_fnames) - 1)):
            raise ValueError("parfiles give {} while index requiested {}".format(len(self.list_fnames), idx))
        if not (os.path.isfile(self.working_dir+self.list_fnames[idx])):
            raise FileNotFoundError("parfile not found {}".format(self.working_dir+self.list_fnames[idx]))
        pba = BPA_METHODS(workingdir=self.working_dir, readparfileforpaths=True, parfile=self.list_fnames[idx])
        pba.run()
        pba.clear()

def distribute_and_run(working_dir:str,list_parfiles:list, n_cpu:int):
    assert len(list_parfiles) > 0
    pba_parallel = PBA_PARALLEL(workind_dir=working_dir, list_parfile_fnames=list_parfiles)
    if (n_cpu == 1):
        for i, pars in enumerate(list_parfiles):
            pba_parallel(i)
    else:
        if (n_cpu is None):
            ncpus = os.cpu_count()
        else:
            ncpus = int(n_cpu)
        try:
            pool = Pool(ncpus)  # on 8 processors
            pool.map(pba_parallel, range(len(list_parfiles)))  # distribute tasks for each clkass
            # pool.map(pb, *[(pars, kn_pars, grb_pars) for pars, kn_pars, grb_pars in zip(pars_list, grb_pars_list, kn_pars_list)])  # distribute tasks for each clkass
        finally:  # To make sure processes are closed in the end, even if errors happen
            pool.close()
            pool.join()
    print("Finished 'distribute_and_run()'")

