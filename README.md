# PyBlastAfterglow

> The code is under development. If there is a bug or missing feature, feel free to _create an issue_ in the repository

PyBlastAfterglow is a C++ code with python interface designed to model electromagnetic signatures (light curves, spectra, sky maps) of gamma-ray bursts and kilonova afterglows. The model approximates the hydrodynamics of the fluid flow by treating each fluid element as an independent element, a system comprised of forward shock, contact discontinuity, and reverse shock (a blastwave).  

At these hydrodynamic shocks, collisionless shocks form where electrons that cross shocks' front are accelerated and seed magnetic fields are amplified. The code _numerically_ follows the evolution of electron distribution in shocks' downstream, accounting for radiation and adiabatic cooling. Afterwards, the comoving synchrotron and synchrotron self-Compton (SSC) emission are computed numerically. Synchrotron-self absorption and pair-production are included in an approximate form. The observed emission is obtained via equal-time arrival surface integration. The code is modular and designed so that new physics (e.g., new radiation methods, energy injection) can be easily incorporated.  

The code also includes several formulations of the same physics for testing, comparison, and analysis. For example analytical model of electron distribution and synchrotron radiation, only forward shock dynamics, two methods to discretize a jet into lateral layers (see paper).  

_In summation:_ PyBlastAfterglow is designed to stand in between full numeric radiation relativistic hydrodynamics codes (e.g., GAMMA, JET) and semi-analytic codes (e.g., afterglowpy, jetsimpy). It has more physics than the latter and is faster than the former. 

## Repository structure

* `/analysis/` contains the simulation analysis
* `/data/` contains tabulated data used by the code (e.g., EBL table)
* `/notebooks/` contains .ipynb files with examples on how to run/use the code
* `/package/` contains python library that allows to easily launch runs and process simulation output
* `/sandbox/` contains experiments with the code
* `/src/` contains C++ source code (header files only)
* `/tst/` contains test to compare code output with other codes (afterglowpy and jetsimpy)

## Installation

1. Install the C++ source code.
   To install C++ code base, `g++` with lhdf5 is required.  
   C++ requirements are noted in `Makefile`. Then, execute   
   `cd src/`  
   `make main`  
   or  
   `cd src/`  
   `make alterantive`  
   depending on how h5lib was installed, it might be tricky to compile.  
   Check if `/src/pba.out` was successfully generated.

2. Install the Python interface for the code input/output.
   Python requirements are noted in `package/requirements.py`.
   A separate conda environment is advised. Installation of the package is done via `pip` as:  
   `cd package/`  
   `pip uninstall --no-cache-dir PyBlastAfterglowMag & pip install .`  
   Check if code runs by executing one of the notebooks in the  `/notebooks/`

## Usage

PyBlastAfterglow consists of two main components: _core numerical model_ written in C++ 
and located in `/src/` and a _Python package_ located in `/package/`. 

The _core numerical model_ is a self-contained simulation model that is independent from 
Python package and can be used by (i) compiling C++ code and (ii) using text parfile `parfile.par`
that contains model settings. Once the code is compled and parfile is prepared. The code can be launched as: 
```bash
pba.out ./path/to/working/dir/ parfile.par loglevel_name
```
where `loglevel_name` can be one of the following: [debug, info, warn, err]. 

The code will then generate a set of output HDF5 files in the `working directory`. 
What files are generated depends on the settings specified in the parfile.  

_Python package_ provides a simplified way to interact with the _compiled_ code by (i) generating parfile automatically
(ii) providing interface to the code output and connecting them to simple plotting methods. 


Example:

```python
import numpy as np, os, matplotlib.pyplot as plt
from PyBlastAfterglowMag.wrappers import run_grb
from PyBlastAfterglowMag.utils import cgs # constants

# set MAIN model parameters (independent of GRB ejecta)
main_pars=dict(
    d_l= 1.27e+26,              # luminocity distance to the source [cm]
    z = 0.0099,                 # redshift of the source (used in Doppler shifring and EBL table)
    n_ism=np.power(10, -1.60),  # ISM density [cm^-3] (assuming uniform)
    theta_obs=np.deg2rad(20.8), # observer angle [rad] (from pol to jet axis)  
    rtol=5e-7,                  # relative tolerance for adaptive quadrature that compute observed emission
    lc_freqs='array logspace 1e8 1e29 96', # frequencies for light curve calculation
    lc_times='array logspace 3e3 1e10 128', # times for light curve calculation
    tb0=3e2,tb1=1e14,ntb=1000, # burster frame time grid boundary, resolution, for the simulation
)

# set GRB ejecta lateral structure
struct = dict(
    struct="gaussian",          # type of the structure tophat or gaussian
    Eiso_c=np.power(10, 54.1),  # isotropic equivalent energy of the burst 
    Gamma0c=300.,               # lorentz factor of the core of the jet 
    M0c=-1.,                    # mass of the ejecta (if -1 -- infer from Eiso_c and Gamma0c)
    theta_c=np.deg2rad(3.50),   # half-opening angle of the core of the jet
    theta_w=np.deg2rad(25.),    # half-opening angle of the winds of the jet
    n_layers_a=21               # resolution of the jet (number of individual blastwaves)
)

# Set GRB ejecta properties
grb_pars=dict(
    structure=struct,               # set structure of the ejecta
    eps_e_fs=np.power(10,-3.42),    # [microphysics] - FS - fraction energy in electrons
    eps_b_fs=np.power(10,-4.02),    # [microphysics] - FS - fraction energy in magnetic fields
    p_fs=2.10, #2.16,               # [microphysics ]- FS - slope of the injection electron spectrum
    # Settings:
    do_lc='yes',                    # [task] - whether to compute and save light curves
    save_spec='yes',                # [task] - whether to compute and save comoving spectra
    method_synchrotron_fs = "CSYN", # [numeric] method for the synchrotron radiation
    method_ele_fs = "numeric",      # [numeric] use analytical or numeric electron distribution evolution model
)

# Run GRB afterglow simulation
pba = run_grb(
    working_dir=os.getcwd() + '/tmp_gauss/', # directory to save/load from simulation data
    P={
        'main':main_pars, 
        'grb':grb_pars
    },   
    run=True,                # run code itself (if False, it will try to load results)
    path_to_cpp="/home/vsevolod/Work/GIT/GitHub/PyBlastAfterglowMag/src/pba.out",
    loglevel="err",
    process_skymaps=False
)
# get the class that contains GRB ejecta methods 
ejecta = pba.GRB

# access and plot the result
times = ejecta.get_lc_times()
freqs = (3e9, 3.80e+14, 2.41e+17, 1e28)
fig, axes = plt.subplots(nrows=1, ncols=len(freqs), figsize=(12,3), layout='constrained')
for i, freq in enumerate(freqs):
    lc = ejecta.get_lc(freq=freq)
    axes[i].plot(times/cgs.day, lc,color='black')
    axes[i].set_xscale('log')
    axes[i].set_yscale('log')
    axes[i].set_ylim(1e-4*np.max(lc), 10*np.max(lc))
    axes[i].set_title(f"log(freq)={np.log10(freq)}")
    axes[i].set_xlabel("time [day]")
    axes[i].grid(ls=':')
axes[0].set_ylabel("Flux density [mJy]")
axes[-1].legend()
plt.show()
```

In the example above, first the required parameters are set. Parameters are split between:  
- (i) main parameters that define general properties of the source
- (ii) GRB jet lateral structure (that in this case is set by an analytical formula)  
- (iii) GRB parameters that define the ejecta properties and microphysics parameters 

Then the `run_grb()` function creates `parfile.par` in the `working_dir` and if `run=True` it 
runs the PyBlastAfterglow executable at `path_to_cpp` and if not it just instantiates the class that 
has methods to load and process output of the code. The class has a class variable `pba.GRB` that is of 
the type `Ejecta` and contains methods to load and process PyBlastAfterglow output HDF5 files. 

In the example above the `ejecta=pba.GRB` is used to get light curve times and flux densities as `get_lc()`

### Additional Information

- Fill list of parameters that C++ code expects in a parfile is given in `parfile_parameter_list.py` file in the parckage folder
- The parfile is created by the `create_parfile()` function in the `parfile_tools.py` file in the package
- Examples on how to use the code are given in `/notebooks/`
- Comparisons between the light curves generated with PyBlastAfterglow and other code are presented in `/tst/`

## Acknowledgments

The code is inspired by [Gavin. P.~Lamb](https://doi.org/10.1093/mnras/stab2879) 
and [J.~Ryan](https://iopscience.iop.org/article/10.3847/1538-4357/ab93cf) afterglow models.  
Numerical microphsyics module for electron distribution evolution is inspired by [Miceli & Nava paper](https://www.mdpi.com/2075-4434/10/3/66).  


## Attribution  

If you use the code, please consider citing  
```
@ARTICLE{2024arXiv240916852N,
       author = {{Nedora}, Vsevolod and {Crosato Menegazzi}, Ludovica and {Peretti}, Enrico and {Dietrich}, Tim and {Shibata}, Masaru},
        title = "{Multi-physics framework for fast modeling of gamma-ray burst afterglows}",
      journal = {arXiv e-prints},
     keywords = {Astrophysics - High Energy Astrophysical Phenomena},
         year = 2024,
        month = sep,
          eid = {arXiv:2409.16852},
        pages = {arXiv:2409.16852},
          doi = {10.48550/arXiv.2409.16852},
archivePrefix = {arXiv},
       eprint = {2409.16852},
 primaryClass = {astro-ph.HE},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2024arXiv240916852N},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}


@ARTICLE{2023MNRAS.tmp..259N,
       author = {{Nedora}, Vsevolod and {Dietrich}, Tim and {Shibata}, Masaru and {Pohl}, Martin and {Menegazzi}, Ludovica Crosato},
        title = "{Modeling kilonova afterglows: Effects of the thermal electron population and interaction with GRB outflows}",
      journal = "Mon. Not. Roy. Astron. Soc.,
     keywords = {neutron star mergers, stars: neutron, equation of state, gravitational waves, Astrophysics - High Energy Astrophysical Phenomena, General Relativity and Quantum Cosmology},
         year = 2023,
        month = jan,
          doi = {10.1093/mnras/stad175},
archivePrefix = {arXiv},
       eprint = {2208.01558},
 primaryClass = {astro-ph.HE},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2023MNRAS.tmp..259N},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
```
