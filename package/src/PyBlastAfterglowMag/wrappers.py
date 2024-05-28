import os,shutil,copy,subprocess

from .interface import PyBlastAfterglow
from .parfile_tools import read_parfile,create_parfile
from .skymap_process import ProcessRawSkymap

from .id_analytic import JetStruct
from .id_david import prepare_kn_ej_id_2d
from .id_kenta import EjectaData

def run_grb(working_dir: str, P: dict, run: bool = True, process_skymaps: bool = True,
            loglevel:str="info", path_to_cpp:str or None=None) -> PyBlastAfterglow:
    """
            conf = {"nx": 64, "ny": 32, "extend_grid": 2, "fwhm_fac": 0.5, "lat_dist_method": "integ",
                "intp_filter": {"type": None, "sigma": 2, "mode": 'reflect'},  # "gaussian"
                "hist_filter": {"type": None, "sigma": 2, "mode": 'reflect'}}
    :param working_dir:
    :param struct:
    :param P:
    :param run:
    :return:
    """
    # clean he temporary direcotry
    if run and os.path.isdir(working_dir):
        shutil.rmtree(working_dir)
    if not os.path.isdir(working_dir):
        os.mkdir(working_dir)

    # generate initial data for blast waves
    P = copy.deepcopy(P)
    struct = copy.deepcopy(P["grb"]["structure"])
    del P["grb"]["structure"]
    if ("corr_fpath_david" in struct.keys()):
        prepare_kn_ej_id_2d(
            nlayers=None if not "n_layers_pw" in struct.keys() else struct["n_layers_pw"],
            corr_fpath=struct["corr_fpath_david"],outfpath=working_dir+"id_pw.h5", dist="pw"
        )
    else:
        pba_id = JetStruct(
            n_layers_pw=81 if not "n_layers_pw" in struct.keys() else struct["n_layers_pw"],
            n_layers_a=(1 if struct["struct"] == "tophat" else
                        (21 if not "n_layers_a" in struct.keys() else struct["n_layers_a"])))

        # save piece-wise EATS ID
        id_dict, id_pars = pba_id.get_1D_id(pars=struct, type="piece-wise")
        pba_id.save_1d_id(id_dict=id_dict, id_pars=id_pars, outfpath=working_dir + "id_pw.h5")

        # save adaptive EATS ID
        id_dict, id_pars = pba_id.get_1D_id(pars=struct, type="adaptive")
        pba_id.save_1d_id(id_dict=id_dict, id_pars=id_pars, outfpath=working_dir + "id_a.h5")

    # get eats type for the grb (it defines initial data and the model itself)
    if not "eats_type" in P["grb"].keys():
        type = "a"
    else:
        type = P["grb"]["eats_type"]
        del P["grb"]["eats_type"]
    P["grb"]["fname_ejecta_id"] = "id_a.h5" if type == "a" else "id_pw.h5"
    P["grb"]["method_eats"] = "piece-wise" if type == "pw" else "adaptive"
    if (struct["struct"]=="tophat"): P["grb"]["nsublayers"] = 35 # for skymap resolution
    grb_skymap_config = None
    if "skymap_conf" in P["grb"].keys():
        grb_skymap_config = copy.deepcopy(P["grb"]["skymap_conf"])
        del P["grb"]["skymap_conf"]

    # create new parfile for the simulation
    create_parfile(working_dir=working_dir, P=P)

    # instantiate PyBlastAfterglow
    pba = PyBlastAfterglow(workingdir=working_dir)

    # run the code with given parfile
    if run:
        # this mess is because I did not figure out how $PATH thing works...
        # curdir = os.getcwd()
        if path_to_cpp is None:
            path_to_cpp = str(__file__).split("package")[0]+"src/pba.out" # todo make it proper through bashrc and setting a path
        # pbadir = curdir.split("PyBlastAfterglowMag")[0]
        # path_to_cpp_executable = pbadir+"PyBlastAfterglowMag"+"/src/pba.out"
        # print(os.getcwd())
        # os.chdir("../../../src/")
        # path_to_executable = "pba.out"
        if not os.path.isfile(path_to_cpp):
            raise IOError("executable is not found: {}".format(path_to_cpp))
        # subprocess.call(path_to_executable, input="")
        # print("{} {} {} {}".format(path_to_cpp_executable, self.workingdir, self.parfile, self.loglevel))
        # subprocess.run(path_to_cpp_executable, input=self.workingdir)
        subprocess.check_call([path_to_cpp, working_dir, pba.parfile, loglevel])

    # process skymap
    if (process_skymaps and pba.GRB.opts["do_skymap"] == "yes"):
        prep = ProcessRawSkymap(conf=grb_skymap_config, verbose=False)
        prep.process_singles(infpaths=working_dir + "raw_skymap_*.h5",
                             outfpath=pba.GRB.fpath_sky_map,
                             remove_input=False)

    return pba

def run_kn(working_dir: str, struct: dict, P: dict, run: bool = True, process_skymaps: bool = True,
            loglevel:str="info", path_to_cpp:str or None=None) -> PyBlastAfterglow:
    """
            conf = {"nx": 64, "ny": 32, "extend_grid": 2, "fwhm_fac": 0.5, "lat_dist_method": "integ",
                "intp_filter": {"type": None, "sigma": 2, "mode": 'reflect'},  # "gaussian"
                "hist_filter": {"type": None, "sigma": 2, "mode": 'reflect'}}
    :param working_dir:
    :param struct:
    :param P:
    :param type:
    :param run:
    :return:
    """
    # clean he temporary direcotry
    if run and os.path.isdir(working_dir):
        shutil.rmtree(working_dir)
    if not os.path.isdir(working_dir):
        os.mkdir(working_dir)

    # generate initial data for blast waves
    struct = copy.deepcopy(struct)
    if ("corr_fpath_david" in struct.keys()):
        prepare_kn_ej_id_2d(
            nlayers=None if not "n_layers_pw" in struct.keys() else struct["n_layers_pw"],
            corr_fpath=struct["corr_fpath_david"],outfpath=working_dir+"id_pw.h5", dist="pw"
        )
    else:
        raise KeyError()

    # create new parfile
    P = copy.deepcopy(P)
    P["kn"]["fname_ejecta_id"] = "id_a.h5" if type == "a" else "id_pw.h5"
    P["kn"]["method_eats"] = "piece-wise" if type == "pw" else "adaptive"
    if (struct["struct"]=="tophat"): P["kn"]["nsublayers"] = 35 # for skymap resolution
    grb_skymap_config = None
    if "skymap_conf" in P["kn"].keys():
        grb_skymap_config = copy.deepcopy(P["kn"]["skymap_conf"])
        del P["kn"]["skymap_conf"]

    # create new parfile for the simulation
    create_parfile(working_dir=working_dir, P=P)

    # instantiate PyBlastAfterglow
    pba = PyBlastAfterglow(workingdir=working_dir)

    # run the code with given parfile
    if run:
        # this mess is because I did not figure out how $PATH thing works...
        # curdir = os.getcwd()
        if path_to_cpp is None:
            path_to_cpp = str(__file__).split("package")[0]+"src/pba.out" # todo make it proper through bashrc and setting a path
        # pbadir = curdir.split("PyBlastAfterglowMag")[0]
        # path_to_cpp_executable = pbadir+"PyBlastAfterglowMag"+"/src/pba.out"
        # print(os.getcwd())
        # os.chdir("../../../src/")
        # path_to_executable = "pba.out"
        if not os.path.isfile(path_to_cpp):
            raise IOError("executable is not found: {}".format(path_to_cpp))
        # subprocess.call(path_to_executable, input="")
        # print("{} {} {} {}".format(path_to_cpp_executable, self.workingdir, self.parfile, self.loglevel))
        # subprocess.run(path_to_cpp_executable, input=self.workingdir)
        subprocess.check_call([path_to_cpp, working_dir, pba.parfile, loglevel])

    # process skymap
    if (process_skymaps and pba.KN.opts["do_skymap"] == "yes"):
        prep = ProcessRawSkymap(conf=grb_skymap_config, verbose=False)
        prep.process_singles(infpaths=working_dir + "raw_skymap_*.h5",
                             outfpath=pba.KN.fpath_sky_map,
                             remove_input=False)

    return pba
