'''
    List of all the parameters C++ code expects in parfile.par and their defaulr values
'''

default_parfile_main_part = dict(
    pars=dict(
        tb0 = 1.e3,     # [numeric] start of the time grid (burster frame) [s]
        tb1 = 1.e11,    # [numeric] end of the time grid (burster frame) [s]
        ntb = 1000,     # [numeric] number of grid points for ODE solver
        iout = 1,       # [numeric] keep (store) solution at 'iout'th iteration
        rtol = 1e-8,    # [numeric] relative tolerance for ODE solver
        nmax = 1000,    # [numeric] maximum number of iteration for adaptive ODE solver
        A0 = -1,        # [ISM] wind environment constant, keep =0 for uniform ISM
        s = -1,         # [ISM] wind environment slope, keep =0 for uniform ISM
        r_ej = -1,      # [ISM] radius at which first break in density profile [cm]
        r_ism = -1,     # [ISM] radius at which second break in density profile [cm]
        n_ism = 1.,     # [ISM] ism number density if it is constnat [cm^-3]
        d_l = 3.09e26,  # [source] luminocity distance to the source
        z = 0.028,      # [source] redshift of the source (used for EBL table interpolation as well)
        theta_obs = 0   # [source] observer angle with respect to the polar axis
    ),
    opts=dict(
        do_average_solution = "no", # [numeric] store only averaged ODE result for 'n_store_substeps' (see below)
        lc_freqs = "array 1e9 1e18", # [obs] frequencies to compute light curve [Hz]
        lc_times = "array logspace 3.e3 1.e10 100", # [obs] times to compute light curve [s]
        lc_use_freq_to_time = "no", # [obs] assume freq-to-time 1-1 relation (lc=len(times)*len(freqs) otherwise)
        skymap_freqs = "array 1e9 1e18", # [obs] frequencies to compute skymap [Hz]
        skymap_times = "array logspace 8.e5 2e9 50", # [obs] times to compute skymaps [s]
        integrator = "DOP853E", # [numeric] ODE solver (DOP853E uses adaptive step-size)
    )
)
default_parfile_grb_part = dict(
    pars=dict(
        # --- Forward Shock Microphysics ---
        eps_e_fs = 0.1,
        eps_b_fs = 0.01,
        eps_t_fs = 0.,
        p_fs = 2.2,
        epsilon_e_rad_fs = -1, # [FS] fraction of the Esh2 removed due to radiation per timestep (dynamics)
        gamma_max_fs = 1e7,    # [numeric] Used only if 'method_gamma_max_fs=useConst'; must be < gam2
        max_substeps_fs = 2000,# [numeric] Number of cooling substeps in electron evolution (between evol.steps)
        gam1_fs = 1.,      # [numeric] lower lim for comoving electron spectrum
        gam2_fs = 1.e8,    # [numeric] upper lim for comoving electron spectrum
        ngam_fs = 451,     # [numeric] size of the electron grid points for Chang-Cooper scheme
        freq1_fs = 1.e5,   # [numeric] lower lim for comoving synchrotron spectrum
        freq2_fs = 1.e32,  # [numeric] uppers lim for comoving synchrotron spectrum
        nfreq_fs = 601,    # [numeric] size of the freq. grid points for Chang-Cooper scheme
        # --- Reverse shock Microphsyics ---
        eps_e_rs = 0.1,
        eps_b_rs = 0.01,
        eps_t_rs = 0.,
        p_rs = 2.2,
        epsilon_e_rad_rs = -1, # [RS] fraction of the Esh2 removed due to radiation per timestep (dynamics)
        gamma_max_rs = 1e7,    # [numeric] Used only if 'method_gamma_max_fs=useConst'; must be < gam2
        max_substeps_rs = 2000,# [numeric] Number of cooling substeps in electron evolution (between evol.steps)
        gam1_rs = 1.,      # [numeric] lower lim for comoving electron spectrum
        gam2_rs = 1.e7,    # [numeric] upper lim for comoving electron spectrum
        ngam_rs = 451,     # [numeric] size of the electron grid points for Chang-Cooper scheme
        freq1_rs = 1.e5,   # [numeric] lower lim for comoving synchrotron spectrum
        freq2_rs = 1.e32,  # [numeric] uppers lim for comoving synchrotron spectrum
        nfreq_rs = 601,    # [numeric] size of the freq. grid points for Chang-Cooper scheme
        # -------------------
        n_store_substeps = 10,  # use n steps of ODE solver to average over and store (used if iout >> 1)
        tprompt = 1.e3,         # [RS] duration of the ejection (for RS initial width Delta=tprompt*c)
        a = 1,                  # [spread] if method_spread="AA", controls dtheta/dR slope
        rs_shutOff_criterion_rho = 1e-50, # [RS] criterion for rho4 when to shut down the reverse shock
        min_Gamma0_for_rs=5.,   # [RS] If initial Gamma0 of a BW (layer) < this value, use only 'fs' RHS not 'fsrs'
        mom_when_start_spread=2,# [spread] Val for \Gamma\beta below which spread is allowed
        mom0_frac_when_start_spread = 0.9, # [spread] if \Gamma\beta < frac * \Gamma_0\beta_0 spread is allowed
        rs_Gamma0_frac_no_exceed = .92, # [RS] if Gamma > frac*Gamma0; set dGammadR = 0 (prevent error acceleration)
        save_dyn_every_it = 10, # [numeric] if to save dynamics, save every it'th iteration,
        rtol_phi = 1e-6,        # [eats] relative tolerance for adaptive quadrature for EATS integration
        rtol_theta = 1e-6,      # [eats] relative tolerance for adaptive quadrature for EATS integration
        nmax_phi = 1000,        # [eats] maximum number of adaptive quadrature iterations to find solution
        nmax_theta = 1000,      # [eats] maximum number of adaptive quadrature iterations to find solution
        theta_max = 3.1415/2.,  # [eats] maximum extend of each hemispheres
        beta_min_fs = 1e-5,     # [numeric] if < betaShock; do not compute any microphsyics (FS)
        beta_min_rs = 1e-5,     # [numeric] if < betaShock; do not compute any microphsyics (RS)
        # --- Skymap adaptive calculation; resize untill 'min_sublayers' each of which has 'min_non_zero_cells'
        nsublayers = 5,         # [numeric] initial division of a theta-layer into sublayers (may give I=0 cells -> 'redo')
        frac_to_increase = 1.5, # [numeric] frac. to increase nsublayers id skymap is not resolved (see below)
        max_restarts = 12,      # [numeric] max number of increasing 'nsublayers' to resolve the jet
        min_sublayers = 3,      # [numeric] criterion, min number sublayers each of which has 'min_non_zero_cells'
        min_non_zero_cells = 5, # [numeric] criterion, min number of phi-cells with intensity > 0 for jet to be resolved
        im_max_theta = 1.5708   # [numeric] max value in intensity calculations for skymap
    ),
    opts=dict(
        run_bws = "yes",        # [task] evolve the blastwaves (if no, expected that load_dynamics=yes)
        save_dynamics = "no",   # [task] save blastwaves evolution history
        load_dynamics = "no",   # [task] load the blastwave dynamics from a file
        do_mphys_in_situ= "yes",# [task] compute electrons/syc/ssc comov.spec after each it. of BW evolution
        do_mphys_in_ppr = "no", # [task] compute electrons/syc/ssc comov.spec after full BW evolution has finisheds
        # do_ele = "yes",         # [task] compute comoving electron spectrum (numerically or analytically)
        # do_spec = "no",         # [task] compute comoving synchrotron/SSC spectrum (numerically or analytically)
        save_spec = "no",       # [task] save comoving ele/synchrotron/SSC spectrum (numerically or analytically)
        do_lc = "yes",          # [task] compute & save light curves
        do_skymap = "no",       # [task] compute & save raw skymaps
        save_raw_skymap="yes",  # [task] currently only raw, unstructured images can be saved. (UNFINISHED)
        skymap_remove_mu = "no",# [task] remove 'mu' from skymap calculation
        counter_jet = "yes",    # [numeric] do include counter jet as well in LCs and Sky Maps

        do_rs = "no",               # [RS] include RS into consideration (main switch)
        do_rs_radiation="yes",      # [RS] if RS is included, compute also the radiation from RS (adds to total LC)
        bw_type = "fs",             # [numeric] type pf the blastwave RHS to use, e.g. fs - forward shock only
        init_deltaR4="no",          # [numeric] set deltaR4[0] = c*beta0*tprompt for ODE solver
        exponential_rho4="yes",     # [numeric] use exponential rho4 decay as exp(-Delta4/Delta0)
        method_collision = "none",  # [numeric] include blastwave collision (UNFINISHED)
        method_eats = "adaptive",   # [numeric] main switch for blastwave discretezation (piece-wise or adaptive)
        method_quad = "CADRE",      # [numeric] EATS quadrature method (CADRE is adaptive)
        method_comp_mode = "comovSpec", # [numeric] interpolated comoving spectra, or compute in-situe
        allow_termination = "no",   # [numerc] continue if one of the blastwaves fails (ODE solver fails)

        do_thermrad_loss = "no",    # [numeric] include thermal radiation from ejecta (UNFINISHED)
        do_eninj_inside_rhs = "no", # [numeric] magnetar-driven ejecta; (UNFINISHED)

        use_1d_id = "yes",          # [I/O] type of the initail data, if 'yes' expects 1D arrays with E,Gamma...
        fname_ejecta_id = "id.h5",  # [I/O] file name (in working_dir) with initial data
        load_r0 = "no",             # [I/O] use R0 from the file instead of computing it as R0=beta0 * tb0 * c
        fname_dyn = "dyn.h5",       # [I/O] file name (in working_dir) to save dynamics
        fname_spectrum = "spec.h5", # [I/O] file name (in working_dir) to save the comoving ele/radition spectrum
        fname_light_curve = "lc.h5",# [I/O] file name (in working_dir) to save light curve
        fname_sky_map = "skymap.h5",# [I/O] file name (in working_dir) to save raw light curve

        ebl_tbl_fpath = "../../../data/EBL/Franceschini18/table.h5", # if not "none"  -- no EBL correction

        do_nucinj = "no", # [numeric] include r-process heating in ejecta (UNFINISHED)

        method_spread = "our",              # [spread] method for lateral spreading
        method_limit_spread="MomValAndFrac",# [numeric] how to limit spreading of the blastwave
        method_dgdr = "our",                # [numeric] choice of equation for BW dynamical evolution dGamma/dR
        method_eos = "Nava13",              # [numeric] choice of EOS for the blast wave
        method_dmdr = "usingdthdr",         # [numeric] choice of equation for accreted mass dm/dr

        use_dens_prof_behind_jet_for_ejecta = "no", # [numeric] include jet in ejecta mode (UNFINISHED)

        # --- Forward Shock ---
        # method_radius_fs = "useGammaR",       # [numeric] how to ge radius for forward shock (use GammaShock or not)
        use_adiabLoss = "yes",              # [numeric] include blast wave adiabatic lossess (FS)
        # method_Gamma_fs = "useGammaShock",  # [numeric] compute GammaShock via EOS or assume = to Gamma (not used)
        method_Up_fs = "useEint2",          # [numeric] compute internal energy from Eint2 or Gamma
        method_thickness_fs = "useJoh06",   # [numeric] compute shock thickness dR, as 1/Gamma^2 or Johannesson paper
        method_vel_fs = "sameAsBW",         # [numeric] "shockVel" in EATS, compute abberation using GammaShock or Gamma
        method_ele_fs = "numeric",          # [numeric] assume analytical electron profile or evolve
        num_ele_use_adi_loss_fs="yes",      # [numeric] include adiabatic cooling term into kinetic eq. for ele. evolution
        method_ne_fs = "useNe",             # [numeric] compute emissivities using Ne or nprime
        method_nonrel_dist_fs ="use_Sironi",# [numeric] include Deep Newtonian regime for electron dist.
        method_gamma_min_fs = "useNumeric", # [numeric] how to compute gamma_min
        method_gamma_c_fs = "useTcomov",    # [numeric] how to compute gamma_c
        method_gamma_max_fs = "useB",       # [numeric] how to compute gamma_max
        method_B_fs = "useU_b" ,            # [numeric] how to compute magnetic field
        method_synchrotron_fs = "CSYN",     # [numeric] method for the synchrotron radiation
        incl_th_in_marg21_fs = "yes",       # [numeric] if Syn.method is MARG21; allows to toggle j_th and a_th
        use_ssa_fs = "no",                  # [numeric] include SSA
        method_ssc_fs = "none",             # [numeric] method for SSC
        method_pp_fs = "none",              # [numeric] method for pair-production
        method_tau_fs = "smooth",           # [numeric] method for optical depth calculation in shock

        # --- Reverse Shock ---
        # method_radius_rs = "sameAsR",     # [numeric] how to ge radius for reverse shock (use GammaShock or not)
        use_adiabLoss_rs = "yes",           # [numeric] include blast wave adiabatic lossess (RS)
        # method_Gamma_rs = "useGammaShock",  # [numeric] compute GammaShock via EOS or assume = to Gamma (not used)
        method_Up_rs = "useEint2",          # [numeric] compute internal energy from Eint2 or Gamma
        method_thickness_rs = "useJoh06",   # [numeric] compute shock thickness dR, as 1/Gamma^2 or Johannesson paper
        method_vel_rs = "sameAsBW",         # [numeric] compute shock thickness dR, as 1/Gamma^2 or Johannesson paper
        method_ele_rs = "numeric",          # [numeric] assume analytical electron profile or evolve
        num_ele_use_adi_loss_rs="yes",      # [numeric] include adiabatic cooling term into kinetic eq. for ele. evolution
        method_ne_rs = "useNe",             # [numeric] compute emissivities using Ne or nprime
        method_nonrel_dist_rs ="use_Sironi",# [numeric] include Deep Newtonian regime for electron dist.
        method_gamma_min_rs = "useNumeric", # [numeric] how to compute gamma_min
        method_gamma_c_rs = "useTcomov",    # [numeric] how to compute gamma_c
        method_gamma_max_rs = "useB",       # [numeric] how to compute gamma_max
        method_B_rs = "useU_b",             # [numeric] how to compute magnetic field
        method_synchrotron_rs = "CSYN",     # [numeric] method for the synchrotron radiation
        incl_th_in_marg21_rs = "yes",       # [numeric] if Syn.method is MARG21; allows to toggle j_th and a_th
        use_ssa_rs = "no",                  # [numeric] include SSA
        method_ssc_rs = "none",             # [numeric] method for SSC
        method_pp_rs = "none",              # [numeric] method for pair-production
        method_tau_rs = "smooth",           # [numeric] method for optical depth calculation in shock
    )
)

# NOT IMPLEMENTED
default_parfile_magnetar_part = dict(
    pars = dict(),
    opts = dict()
)

# NOT IMPLEMENTED
default_parfile_pwn_part = dict(
    pars = dict(),
    opts = dict()
)