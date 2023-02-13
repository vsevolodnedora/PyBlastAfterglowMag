# PyBlastAfterglowMag

Semi-analytic model for computing GRB and kilonova afterglow written in C++ with python class for postprocessing and plotting the results. 

The code is design for easy implementation of new physics and numerical methods. Specifically, the code suppoerts several formulations of a blast-wave evolution, blast-wave lateral spreading, analytic electron distribution and symchrotron emission approximants. 

The code is inspired by [Gavin. P.~Lamb](https://doi.org/10.1093/mnras/stab2879) and [J.~Ryan](https://iopscience.iop.org/article/10.3847/1538-4357/ab93cf) afterglow models. 
Overall, the code discretises ejecta into angular and velocity elements. Then, considering each of them as a part of a spherical blastwave propageting through ambient medium, it solves a system of ODEs for energy consdervation. After that, the code computes doppler-corrected synchrotron emission and integrates the observed flux via EATS integration.  

# Installation  

To install C++ code base, `g++` with lhdf5 is required.  
C++ requirements are noted in `Makefile`. Then, execute  
`cd src/`  
`make main`  
or  
`cd src/`  
`make alterantive`  
depending on how h5lib was installed, it might be tricky to complile. 

Once combilation was done successfully, and `src/pba.out` was generated, install the python interface that allows an easy use of the code output. 
Python requirements noted in `package/requirements.py`. Separate conda environment is advised. Installation of the package is done via `pip` as:  
`cd package/`  
`pip uninstall --no-cache-dir PyBlastAfterglowMag & pip install .` 


# Usage  

To run the code requires a parfile. Examples of a parfile can be found in `tst/` and `projects/` directoris. Depending on the usage, a parfile can differ a lot. For common GRB afterglow application see `projects/grbgauss_mcmc/` and `\projects/grbtophat_parallel/` . 
For all possible applications see `/tst/`.  
The jupyter notebook section with better examples is coming soon. 

# Attribution  
If you use the code, please cite  
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