# PyBlastAfterglow

PyBlastAfterglow is a numerical, C++ code designed to model electromagnetic signatures (light curves, spectra, sky maps) of gamma-ray bursts and kilonova afterglows. The model approximates the hydrodynamics of the fluid flow by treating each fluid element as an independent element, a system comprised of forward shock, contact discontinuity, and reverse shock. These shocks are collisionless shocks that form when ejecta propagates through the interstellar medium. At these shocks injected electrons are accelerated and seed magnetic fields are amplified. The code numerically follows the evolution of electron distribution in shocks downstream, accounting for radiation and adiabatic cooling. Afterwards, the comoving synchrotron and synchrotron self-Compton (SSC)  emission are computed numerically. Finally, the observed emission is obtained via equal-time arrival surface integration. 
The code is modular and designed so that new physics (e.g., new radiation methods, energy injection) can be easily incorporated. 
_In summation:- PyBlastAfterglow stands in between full numeric radiation relativistic hydrodynamics codes (e.g., GAMMA, JET) and semi-analytic codes (e.g., afterglowpy). It has much more physics than the latter and is faster than the former.  
__The code is under active development and is not yet released__.  

## Repository structure

The C++ source code is located in `/src/`  
The Python interface and processing tools are located in `/package/`  
The test cases are located in `/tst/` for GRB and kilonova afterglows.

## Installation  

To install C++ code base, `g++` with lhdf5 is required.  
C++ requirements are noted in `Makefile`. Then, execute  
`cd src/`  
`make main`  
or  
`cd src/`  
`make alterantive`  
depending on how h5lib was installed, it might be tricky to complile. 

Once compilation was done successfully, and `src/pba.out` was generated, install the python interface that allows easy use of the code output. 
Python requirements are noted in `package/requirements.py`. A separate conda environment is advised. Installation of the package is done via `pip` as:  
`cd package/`  
`pip uninstall --no-cache-dir PyBlastAfterglowMag & pip install .` 


## Sources

The code is inspired by [Gavin. P.~Lamb](https://doi.org/10.1093/mnras/stab2879) and [J.~Ryan](https://iopscience.iop.org/article/10.3847/1538-4357/ab93cf) afterglow models. 

## Usage  

To run the code requires a parfile. Examples of a profile can be found in `/tst/` directory. 


## Attribution  
If you use the code, please consider citing  
```
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
