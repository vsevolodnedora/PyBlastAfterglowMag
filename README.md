# PyBlastAfterglow

PyBlastAfterglow is a C++ code with python interface designed to model electromagnetic signatures (light curves, spectra, sky maps) of gamma-ray bursts and kilonova afterglows. The model approximates the hydrodynamics of the fluid flow by treating each fluid element as an independent element, a system comprised of forward shock, contact discontinuity, and reverse shock (a blastwave).  

At these hydrodynamic shocks, collisionless shocks form where electrons that cross shocks' front are accelerated and seed magnetic fields are amplified. The code _numerically_ follows the evolution of electron distribution in shocks' downstream, accounting for radiation and adiabatic cooling. Afterwards, the comoving synchrotron and synchrotron self-Compton (SSC) emission are computed numerically. Synchrotron-self absorption and pair-production are included in an approximate form. The observed emission is obtained via equal-time arrival surface integration. The code is modular and designed so that new physics (e.g., new radiation methods, energy injection) can be easily incorporated.  

The code also includes several formulations of the same physics for testing, comparison, and analysis. For example analytical model of electron distribution and synchrotron radiation, only forward shock dynamics, two methods to discretize a jet into lateral layers (see paper).  

_In summation:_ - PyBlastAfterglow is designed to stand in between full numeric radiation relativistic hydrodynamics codes (e.g., GAMMA, JET) and semi-analytic codes (e.g., afterglowpy, jetsimpy). It has more physics than the latter and is faster than the former. 

__The code is under active development__ 
  

## Repository structure

* `/analysis/` contains the simulation analysis  
*  `/data/` contains tabulated data used by the code (e.g., EBL table)  
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

## Acknowledgments

The code is inspired by [Gavin. P.~Lamb](https://doi.org/10.1093/mnras/stab2879) 
and [J.~Ryan](https://iopscience.iop.org/article/10.3847/1538-4357/ab93cf) afterglow models.  
Numerical microphsyics module for electron distribution evolution is inspired by [Miceli & Nava paper](https://www.mdpi.com/2075-4434/10/3/66).  

## Usage  

See examples in `/tst/` and `/notebooks/` for how to use the code. 

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
