//
// Created by vsevolod on 02/04/23.
//

#ifndef SRC_EJECTA_H
#define SRC_EJECTA_H

//#include "../utilitites/pch.h"
//#include "../utilitites/utils.h"
//#include "../utilitites/interpolators.h"
//#include "../utilitites/ode_solvers.h"
//#include "../utilitites/quadratures.h"
//#include "../utilitites/rootfinders.h"
//#include "../image.h"
//#include "../synchrotron_an.h"
//
////#include "model_magnetar.h"
//#include "../blastwave/blastwave_components.h"
#include "../blastwave/blastwave.h"
#include "../blastwave/blastwave_collision.h"
#include "ejecta_cumshell.h"
#include "ejecta_base.h"

class Ejecta : public EjectaBase {
    std::unique_ptr<Output> p_out = nullptr;
    std::unique_ptr<logger> p_log = nullptr;

public:
    Ejecta(Vector & t_arr, int loglevel) : EjectaBase(t_arr, loglevel) {
        p_log = std::make_unique<logger>(std::cout, std::cerr, loglevel, "Ejecta");
        p_out = std::make_unique<Output>(loglevel);
    }

    /// dispatch tasks
    void processEvolved(StrDbMap & main_pars, StrStrMap & main_opts){

        /// work on GRB afterglow
        if (run_bws || load_dyn){
            bool lc_freq_to_time = getBoolOpt("lc_use_freq_to_time",main_opts,AT,p_log,false,true);
            Vector lc_freqs = makeVecFromString(getStrOpt("lc_freqs",main_opts,AT,p_log,"",true),p_log);
            Vector lc_times = makeVecFromString(getStrOpt("lc_times",main_opts,AT,p_log,"",true), p_log);
            Vector skymap_freqs = makeVecFromString(getStrOpt("skymap_freqs",main_opts,AT,p_log,"",true), p_log);
            Vector skymap_times = makeVecFromString(getStrOpt("skymap_times",main_opts,AT,p_log,"",true), p_log);

            if (load_dyn)
                loadEjectaBWDynamics(working_dir,
                                     getStrOpt("fname_dyn", m_opts, AT, p_log, "", true));

//            if (do_ele)
//                evolveElectronDistributionAndRadiation();

            if (do_spec || do_lc || do_skymap)
//                setPreComputeEjectaAnalyticSynchrotronPars();
                evolveElectronDistributionAndRadiation();

            if (save_dyn)
                saveEjectaBWsDynamics(
                        working_dir,
                        getStrOpt("fname_dyn", m_opts, AT, p_log, "", true),
                        main_pars, m_pars);

            if (do_spec && save_spec)
                computeSaveEjectaSpectrum(
                        working_dir,
                        getStrOpt("fname_spectrum", m_opts, AT, p_log, "", true),
                        main_pars, m_pars);

            if (do_lc) {
                computeSaveEjectaLightCurveAnalytic_new(
                        working_dir,
                        getStrOpt("fname_light_curve", m_opts, AT, p_log, "", true),
                        getStrOpt("fname_light_curve_layers", m_opts, AT, p_log, "", true),
                        lc_times, lc_freqs, main_pars, m_pars, lc_freq_to_time);
//                    (*p_log)(LOG_INFO, AT) << "jet analytic synch. light curve finished [" << timer.checkPoint() << " s]" << "\n";
            }
            if (do_skymap)
                computeSaveEjectaSkyImagesAnalytic_new(
                        working_dir,
                        getStrOpt("fname_sky_map", m_opts, AT, p_log, "", true),
                        skymap_times, skymap_freqs, main_pars, m_pars, m_opts);
//                    (*p_log)(LOG_INFO, AT) << "jet analytic synch. sky map finished [" << timer.checkPoint() << " s]" << "\n";

        }
    }

    /// compute optical depth along the line of sight in ejecta
    bool evalOptDepthsAlongLineOfSight(double & frac, double & tau_comp, double & tau_BH, double tau_bf,
                                       double ctheta, double rpwn, double phi, double theta,
                                       double phi_obs, double theta_obs, double r_obs, size_t il,
                                       double time,  //size_t ia, size_t ib, double ta, double tb,
                                       double freq, double z, size_t info_ilpwn,
                                       void * params) {
        auto * p_pars = (struct PWNPars *) params;
        std::vector<std::string> lightpath {};
        for (size_t ish = 0; ish < nshells(); ish++){
            auto & cumshell = p_cumShells[il];
            auto & bw = cumshell->getBW(ish);
            /// skip if this shells was not initialized or evolved
            if (bw->getPars()->i_end_r==0) { status[il][ish] = 'N'; continue; }
            size_t ia=0, ib=0;
            bool res = false;
            res = bw->evalEATSindexes(ia,ib,time,z,ctheta,phi,theta_obs,obsAngle);\
            /// skip if there is no EATS surface for the time 'time'
            if (!res){ status[il][ish] = 'T'; continue; }
            Vector & ttobs = bw->getPars()->ttobs;
            /// skip if at EATS times there the BW is no loger evolved (collided)
            if ((bw->getData()[BW::Q::iEJr][ia] == 0)||(bw->getData()[BW::Q::iEJr][ib] == 0)) { status[il][ish] = 'n'; continue;  }
            double rej = interpSegLog(ia, ib, time, ttobs, bw->getData()[BW::Q::iR]);
            /// skip if at time the PWN outruns the ejecta shell (should never occure)
            if ( rpwn > rej){ status[il][ish] = 'R'; continue; }
            /// if this is the first entry, change the value to 0
            if (tau_comp < 0){tau_comp = 0.;}
            if (tau_BH < 0){tau_BH = 0.;}
            if (tau_bf < 0){tau_bf = 0.;}
            /// evaluate the optical depths at this time and this freq.
            double _tau_comp, _tau_BH, _tau_bf;
            bw->evalOpticalDepths(_tau_comp,_tau_BH,_tau_bf,ia,ib,time,freq);
            tau_comp+=_tau_comp; tau_BH+=_tau_BH; tau_bf+=_tau_bf;
            /// save for debugging
            lightpath.push_back("il="+std::to_string(il)
                                +" ish="+std::to_string(ish)
                                + " comp="+std::to_string(_tau_comp)
                                + " BH="+std::to_string(_tau_BH)
                                + " bf="+std::to_string(_tau_bf));
            status[il][ish] = '+';
        }
        /// Combine individual optical depths into fraction
        double power_Compton=0.0;
        if (tau_comp > 1.0)
            power_Compton = tau_comp*tau_comp;
        else
            power_Compton = tau_comp;

        double e_gamma = freq * 6.62606957030463E-27;// Hz -> erg  *4.1356655385381E-15*CGS::EV_TO_ERG; // TODO this has to be red-shifted!
//        double e_gamma = freq * 4.1356655385381E-15;
        frac = exp(-(tau_BH+tau_bf))
               * (exp(-(tau_comp)) + (1.0 - exp(-(tau_comp))) * pow(1.0 - PWNradiationMurase::gamma_inelas_Compton(e_gamma),power_Compton));
//        (*p_log)(LOG_INFO,AT) << "[ilpw="<<info_ilpwn<<", ctheta="<<ctheta<<", rpwn="<<rpwn<<", "<<"t="<<time<<"] "
//                              << " tau_comp="<<tau_comp<<" tau_BH="<<tau_BH<<" tau_bf="<<tau_bf<<" -> frac="<<frac<<"\n";
        if ((frac<0)||(frac>1)||!std::isfinite(frac)){
            (*p_log)(LOG_ERR,AT) << "[ilpw=" << info_ilpwn << ", ctheta=" << ctheta << ", rpwn=" << rpwn << ", " << "t=" << time << "] "
                                 <<" tau_comp="<<tau_comp<<" tau_BH="<<tau_BH<<" tau_bf="<<tau_bf<<" -> frac="<<frac<<"\n";
//            exit(1);
            return false;
        }
        return true;
    }

private:

    /// COMPUTE electrons
    void evolveElectronDistributionAndRadiation(){//(StrDbMap pars, StrStrMap opts){
        (*p_log)(LOG_INFO,AT) << "Computing Ejecta analytic electron pars...\n";

        if ((!run_bws)&&(!load_dyn)){
            (*p_log)(LOG_ERR,AT) << " ejecta BWs were not evolved. Cannot evolveElectronDistributionAndRadiation electrons (analytic) exiting...\n";
            exit(1);
        }
        auto & models = getShells();
        for (auto & model : models)
            for (auto & bw : model->getBWs()) {
                auto & pars = bw->getPars();
                /// skip calulations for not evolved blast waves
                if (pars->end_evolution && (pars->E0 < 0))
                    continue;
                bw->evolveElectronDistAndComputeRadiation();
            }

        is_ejecta_anal_ele_computed = true;
        is_ejecta_anal_synch_computed = true;
    }

#if 0
    /// COMPUTE synchrotron
    void setPreComputeEjectaAnalyticSynchrotronPars(){//(StrDbMap pars, StrStrMap opts){
        (*p_log)(LOG_INFO,AT) << "Computing Ejecta analytic synchrotron pars...\n";

        if ((!run_bws)&&(!load_dyn)){
            (*p_log)(LOG_ERR,AT) << " ejecta BWs were not evolved. Cannot setPreComputeEjectaAnalyticSynchrotronPars electrons (analytic) exiting...\n";
            exit(1);
        }

        auto & models = getShells();
        for (auto & model : models)
            for (auto & bw : model->getBWs())
                bw->computeShockRadiation();

        is_ejecta_anal_synch_computed = true;
    }
#endif

private:
    /// OUTPUT dynamics/electrons
    void saveEjectaBWsDynamics(
            std::string workingdir, std::string fname, StrDbMap & main_pars, StrDbMap & ej_pars){


        (*p_log)(LOG_INFO,AT) << "Saving Ejecta BW dynamics...\n";

        std::string fpath = workingdir + fname;
        Output::remove_file_if_existis(fname,p_log);
        H5::H5File file(fpath, H5F_ACC_TRUNC); // "/home/m/Desktop/tryout/file.h5"

        size_t nshells_ = nshells();
        size_t nlayers_ = nlayers();

        /// add data for each shell
        for (size_t ish = 0; ish < nshells_; ish++) {
            for (size_t il = 0; il < nlayers_; il++) {
                std::string group_name = "shell=" + std::to_string(ish) + " layer=" + std::to_string(il);
                H5::Group grp(file.createGroup(group_name));
                auto &bw = getShells()[il]->getBW(ish);
                for (auto & ivn : BW::VARS.at(bw->getPars()->m_type))
                    Output::addVectorToGroup(grp,bw->getData(static_cast<BW::Q>(ivn)), BW::VARNAMES[ivn]);
                StrDbMap bw_atts {
                    {"E0",bw->getPars()->E0},
                    {"Gamma0",bw->getPars()->Gamma0},
                    {"R0",bw->getPars()->R0},
                    {"Eint0",bw->getPars()->Eint0},
                    {"Rd",bw->getPars()->Rd},
                    {"t0",bw->getPars()->tb0},
                    {"theta_b0",bw->getPars()->theta_b0},
                    {"theta_a",bw->getPars()->theta_a}

                };
                Output::addStrDbMap(bw_atts, grp);
                grp.close();
            }
        }

        std::unordered_map<std::string, double> attrs{
                {"nshells", nshells_ },
                {"nlayers", nlayers_ },
                {"ntimes", t_arr.size() },
                {"ncells", ncells() }
        };

        for (auto& [key, value]: main_pars) { attrs[key] = value; }
        for (auto& [key, value]: ej_pars) { attrs[key] = value; }

        Output::addStrDbMap(attrs, file);

        file.close();

    }

    /// OUTPUT skymap
    void computeSaveEjectaSkyImagesAnalytic_new(
            std::string workingdir, std::string fname, Vector times, Vector freqs,
            StrDbMap & main_pars, StrDbMap & ej_pars, StrStrMap & ej_opts){
        if ((!run_bws)&&(!load_dyn))
            return;

        (*p_log)(LOG_INFO,AT) << "Computing and saving Ejecta sky image with analytic synchrotron...\n";

        if (!is_ejecta_anal_synch_computed){
            std::cerr  << "ejecta analytic electrons were not evolved. Cannot computeSynchrotronEmissivityAbsorptionAnalytic images (analytic) exiting...\n";
            std::cerr << AT << " \n";
            exit(1);
        }
        if (!is_ejecta_obs_pars_set){
            std::cerr<< "ejecta observer parameters are not set. Cannot computeSynchrotronEmissivityAbsorptionAnalytic image (analytic) exiting...\n";
            std::cerr << AT << " \n";
            exit(1);
        }

        size_t nshells_ = nshells();
        size_t nlayers_ = nlayers();

        /// for adaptive method: each BW is its own radial shell there, so we split
        if (id->method_eats==EjectaID2::iadaptive) {
            nshells_ = nlayers();
            nlayers_ = 1;
        }

        std::vector<ImageExtend> ims;

        VecVector other_data{ times, freqs };
        std::vector<std::string> other_names { "times", "freqs" };

        size_t itinu = 0;
        for (size_t ifreq = 0; ifreq < freqs.size(); ifreq++) {
            /// allocate space for the fluxes associated with skymaps
            Vector tota_flux(times.size(), 0.0);
            Vector xc(times.size(), 0.0);
            Vector yc(times.size(), 0.0);

            VecVector total_flux_shell( nshells_ );
            for (auto & total_flux_shel : total_flux_shell)
                total_flux_shel.resize( times.size(), 0.0 );

            /// compute skymaps
            for (size_t it = 0; it < times.size(); ++it) {
                ims.emplace_back( ImageExtend(times[it], freqs[ifreq], nshells_, nlayers_, m_loglevel) );
                if (id->method_eats == EjectaID2::STUCT_TYPE::ipiecewise) {
                    computeEjectaSkyMapPW(ims[itinu], times[it], freqs[ifreq]);
                } else if (id->method_eats == EjectaID2::STUCT_TYPE::iadaptive) {
                    computeEjectaSkyMapA(ims[itinu], times[it], freqs[ifreq],
                                         (size_t) getDoublePar("ntheta", ej_pars, AT, p_log, 100, true),
                                         (size_t) getDoublePar("nphi", ej_pars, AT, p_log, 200, true));
                }

                if (getBoolOpt("save_raw_skymap", ej_opts, AT, p_log, false, false))
                    saveRawImage(ims[itinu], itinu,workingdir, main_pars,ej_pars,p_log);
                else {
                    (*p_log)(LOG_ERR , AT) << " Cannot save final image. Only raw images can now be saved...\n";
                    exit(1);
                }

                /// TODO this words but we don't really need to do it here, as we do not save the final image...
//                ims[itinu].evalXcYc();

                /// compute 2D histogram (aka collecting photons that arrive at a given area)
                /// TODO this words but we don't really need to do it here, as we do not save the final image...
//                ims[itinu].binImage((size_t) getDoublePar("skymap_nx", ej_pars, AT, p_log, 200, true),
//                                    (size_t) getDoublePar("skymap_ny", ej_pars, AT, p_log, 200, true) );

                /// perform bilinear interpolation of the data for a regular grid
                // TODO: Because the grid the large and ustructred one need Delaney interpolation. Too complex...
//                ims[itinu].interpImage((size_t) getDoublePar("skymap_nx", ej_pars, AT, p_log, 200, true),
//                                       (size_t) getDoublePar("skymap_ny", ej_pars, AT, p_log, 200, true) );

                itinu++;

            }
            /// TODO sacing the final image is not usefull unless we can interpolate the result...
//            saveImages(ims, times, freqs, main_pars, ej_pars, workingdir+fname,p_log);


//            other_data.emplace_back( tota_flux );
//            other_names.emplace_back( "totalflux at freq="+ string_format("%.4e", freqs[ifreq]) );
//            for (size_t ish = 0; ish < nshells_; ish++){
//                other_data.emplace_back( total_flux_shell[ish] );
//                other_names.emplace_back( "totalflux at freq="+ string_format("%.4e", freqs[ifreq])
//                                          + " shell=" + string_format("%d", ish));
//            }
//
//            /// create names for the groups of data columns
//            std::vector<std::string> group_names{};
//            for (size_t ifreq = 0; ifreq < freqs.size(); ifreq++){
//                for (size_t it = 0; it < times.size(); it++){
//                    for (size_t i_sh = 0; i_sh < nshells_; ++i_sh) {
//                        group_names.emplace_back(
//                                "shell=" + std::to_string(i_sh) +
//                                " time=" +  string_format("%.4e",times[it]) +
//                                " freq=" + string_format("%.4e",freqs[ifreq]));
//                    }
//                }
//            }

            /// names of each data column in the group
//            auto in_group_names = IMG::m_names;
//
//            /// add attributes from model parameters
//            std::unordered_map<std::string,double> attrs{ {"nshells", nshells_}, {"nlayers", nlayers_} };
//            for (auto& [key, value]: main_pars) { attrs[key] = value; }
//            for (auto& [key, value]: ej_pars) { attrs[key] = value; }
//
//            p_out->VectorOfVectorsH5()
//
//            p_out->VecVectorOfVectorsAsGroupsAndVectorOfVectorsH5(workingdir+fname,out,group_names,
//                                                                  in_group_names,
//                                                                  other_data,other_names,attrs);

        }
    }
#if 0
    /// OUTPUT spectrum
    void computeSaveEjectaSpectrumOLD(
            std::string workingdir,std::string fname, std::string fname_shells_layers,
            StrDbMap & main_pars, StrDbMap & ej_pars, bool lc_freq_to_time){

        (*p_log)(LOG_INFO,AT) << "Computing and saving Ejecta spectrum with analytic synchrotron...\n";

//        size_t nshells = p_cumShells->nshells();
//        size_t m_nlayers = p_cumShells->m_nlayers();
//        size_t ncells =  (int)p_cumShells->ncells();

        if (!is_ejecta_anal_synch_computed){
            (*p_log)(LOG_INFO,AT) << " ejecta analytic electrons were not evolved. Cannot computeSynchrotronEmissivityAbsorptionAnalytic light curve (analytic) exiting...\n";
            exit(1);
        }
        if (!is_ejecta_obs_pars_set){
            (*p_log)(LOG_INFO,AT) << " ejecta observer parameters are not set. Cannot computeSynchrotronEmissivityAbsorptionAnalytic light curve (analytic) exiting...\n";
            exit(1);
        }

//        auto & tmp = getShells()[0]->getBW(0)->getSynchAnPtr();

        std::vector< // layers / shells
                std::vector< // options
                        std::vector<double>>> // freqs*times
        out {};

        /// evaluate light curve
//        auto spectrum = evalEjectaSpectrum();
        auto & spec_freqs = p_cumShells[0]->getBW(0)->getPars()->p_syn_a->m_freq_arr;
        if (spec_freqs.size()<1){
            (*p_log)(LOG_INFO,AT) << " m_freq_arr is not initialized for a BW. Cannot compute comoving spectrum \n ";
            exit(1);
        }

        /// save total_rad lightcurve
        size_t n = t_arr.size() * spec_freqs.size();
        Vector total_power (n, 0.0);
        Vector _times, _freqs;
        cast_times_freqs(t_arr,spec_freqs,_times,_freqs,lc_freq_to_time,p_log);

        std::string var = getStrOpt("spec_var_out",m_opts,AT,p_log,"em", true);
        if ((var == "em_rs" || var == "abs_rs") && (!p_cumShells[0]->getBW(0)->getPars()->do_rs)){
            (*p_log)(LOG_INFO,AT) << " cannot ahve 'spec_var_out = 'em_rs or 'abs_rs' if 'do_rs = no' \n ";
            exit(1);
        }

        for (size_t itnu = 0; itnu < n; ++itnu) {
            size_t ishil = 0;
            for (size_t ishell = 0; ishell < nshells(); ++ishell) {
                for (size_t ilayer = 0; ilayer < nlayers(); ++ilayer) {
                    if (var == "em")
                        total_power[itnu] += p_cumShells[ilayer]->getBW(ishell)->getPars()->p_syn_a->out_spectrum[itnu];
                    else if (var == "abs")
                        total_power[itnu] += p_cumShells[ilayer]->getBW(ishell)->getPars()->p_syn_a->out_specturm_ssa[itnu];
                    else if (var == "em_rs")
                        total_power[itnu] += p_cumShells[ilayer]->getBW(ishell)->getPars()->p_syn_a_rs->out_spectrum[itnu];
                    else if (var == "abs_rs")
                        total_power[itnu] += p_cumShells[ilayer]->getBW(ishell)->getPars()->p_syn_a_rs->out_specturm_ssa[itnu];
                    else{
                        (*p_log)(LOG_INFO,AT) << " spec_var_out is not recognized. Possible options: "
                                              << " em "<< " abs "<<" em_rs " <<" abs_rs "<<" Givem="<<var<<"\n";
                        exit(1);
                    }
//                    auto & spectrum = p_cumShells[ilayer]->getBW(ishell)->getPars()->out_spectrum;
//                    if (spectrum.size()<1){
//                        (*p_log)(LOG_INFO,AT) << " spectrum is not initialized for a BW \n ";
//                        exit(1);
//                    }
//                    total_power[itnu] += spectrum[itnu];
                    ishil ++;
                }
            }
        }
        std::vector<std::string> other_names { "times", "freqs", "total_power" };
        VecVector out_data {_times, _freqs, total_power};

        std::unordered_map<std::string,double> attrs{ {"nshells", nshells()}, {"nlayers", nlayers()} };
        for (auto& [key, value]: main_pars) { attrs[key] = value; }
        for (auto& [key, value]: ej_pars) { attrs[key] = value; }
        p_out->VectorOfVectorsH5(out_data, other_names, workingdir+fname,  attrs);


        /// save light curve for each shell and layer
        if (fname_shells_layers == "none")
            return;
        std::vector<std::string> group_names;
        VecVector total_fluxes_shell_layer(nshells()*nlayers());
        size_t ii = 0;
        for (size_t ishell = 0; ishell < nshells(); ++ishell) {
            for (size_t ilayer = 0; ilayer < nlayers(); ++ilayer) {
                group_names.emplace_back("shell=" + std::to_string(ishell) + " layer=" + std::to_string(ilayer));
                total_fluxes_shell_layer[ii].resize(n,0.);
                auto & spectrum = p_cumShells[ilayer]->getBW(ishell)->getPars()->p_syn_a->out_spectrum;
                for (size_t ifnu = 0; ifnu < n; ifnu++) {

                    if (var == "em")
                        total_fluxes_shell_layer[ii][ifnu]= p_cumShells[ilayer]->getBW(ishell)->getPars()->p_syn_a->out_spectrum[ifnu];
                    else if (var == "abs")
                        total_fluxes_shell_layer[ii][ifnu]= p_cumShells[ilayer]->getBW(ishell)->getPars()->p_syn_a->out_specturm_ssa[ifnu];
                    else if (var == "em_rs")
                        total_fluxes_shell_layer[ii][ifnu]= p_cumShells[ilayer]->getBW(ishell)->getPars()->p_syn_a_rs->out_spectrum[ifnu];
                    else if (var == "abs_rs")
                        total_fluxes_shell_layer[ii][ifnu]= p_cumShells[ilayer]->getBW(ishell)->getPars()->p_syn_a_rs->out_specturm_ssa[ifnu];
                    else{
                        (*p_log)(LOG_INFO,AT) << " spec_var_out is not recognized. Possible options: "
                                              << " em "<< " abs "<<" em_rs " <<" abs_rs "<<" Givem="<<var<<"\n";
                        exit(1);
                    }

//                    total_fluxes_shell_layer[ii][ifnu] = spectrum[ifnu];
                }
                ii++;
            }
        }
        total_fluxes_shell_layer.emplace_back(_times);
        total_fluxes_shell_layer.emplace_back(_freqs);
        total_fluxes_shell_layer.emplace_back(total_power);

        group_names.emplace_back("times");
        group_names.emplace_back("freqs");
        group_names.emplace_back("total_power");
        p_out->VectorOfVectorsH5(total_fluxes_shell_layer, group_names, workingdir+fname,  attrs);
    }
#endif
    /// OUTPUT spectrum
    void computeSaveEjectaSpectrum(
            std::string workingdir,std::string fname, StrDbMap & main_pars, StrDbMap & ej_pars){

        (*p_log)(LOG_INFO,AT) << "Computing and saving ejecta spectra...\n";

        if (!is_ejecta_anal_synch_computed){
            (*p_log)(LOG_INFO,AT) << " ejecta electrons were not evolved. "
                                                " Cannot compute light curve exiting...\n";
            exit(1);
        }

        std::string fpath = workingdir+fname;
        Output::remove_file_if_existis(fname,p_log);
        H5::H5File file(fpath, H5F_ACC_TRUNC); // "/home/m/Desktop/tryout/file.h5"

        size_t nshells_ = nshells();
        size_t nlayers_ = nlayers();

        /// for each shell and layer save spectra
        for (size_t ish = 0; ish < nshells(); ++ish) {
            for (size_t il = 0; il < nlayers(); ++il) {
                std::string group_name = "shell=" + std::to_string(ish) + " layer=" + std::to_string(il);
                H5::Group grp(file.createGroup(group_name));
                auto &bw = getShells()[il]->getBW(ish);
                /// save electron spectrum if electrons are numerically evolved
                if (bw->getPars()->p_syn_a->m_eleMethod!=METHODS_SHOCK_ELE::iShockEleAnalyt)
                    Output::addVectorToGroup(grp, bw->getPars()->p_syn_a->ele.f_all,
                                             "n_ele_fs");
                /// save synchrotron spectra for forward shock
                Output::addVectorToGroup(grp, bw->getPars()->p_syn_a->syn.j_all,
                                         "synch_fs");
                Output::addVectorToGroup(grp, bw->getPars()->p_syn_a->syn.a_all,
                                         "ssa_fs");
                /// save SSC spectrum
                if (bw->getPars()->p_syn_a->m_methods_ssc!=METHOD_SSC::inoSSC)
                    Output::addVectorToGroup(grp, bw->getPars()->p_syn_a->ssc.j_all,
                                             "ssc_fs");

                /// save spectra of the reverse shock
                if (bw->getPars()->do_rs){
                    /// save electron spectrum if electrons are numerically evolved
                    if (bw->getPars()->p_syn_a->m_eleMethod!=METHODS_SHOCK_ELE::iShockEleAnalyt)
                        Output::addVectorToGroup(grp, bw->getPars()->p_syn_a_rs->ele.f_all,
                                                 "n_ele_rs");
                    Output::addVectorToGroup(grp, bw->getPars()->p_syn_a_rs->syn.j_all,
                                             "synch_fs_rs");
                    Output::addVectorToGroup(grp, bw->getPars()->p_syn_a_rs->syn.a_all,
                                             "ssa_fs_rs");
                    /// save SSC spectrum
                    if (bw->getPars()->p_syn_a_rs->m_methods_ssc!=METHOD_SSC::inoSSC)
                        Output::addVectorToGroup(grp, bw->getPars()->p_syn_a_rs->ssc.j_all,
                                                 "ssc_rs");
                }
                grp.close();
            }
        }

        /// add time and frequency. Note: specta are 1D and loop is for t in time, { for freq in freqs { } }
        auto &bw0 = getShells()[0]->getBW(0);
        size_t nn = bw0->get_tburst().size() * bw0->getPars()->p_syn_a->total_rad.numbins;

        /// make vectors for time and freq with the same structure as emissivity and absorption
        Vector _times_for_freq(nn, 0.);
        Vector _freqs(nn, 0.);
        size_t ii=0;
        for (size_t it = 0; it < bw0->get_tburst().size(); it++) {
            for (size_t ifreq = 0; ifreq < bw0->getPars()->p_syn_a->total_rad.numbins; ifreq++){
                _times_for_freq[ii]=bw0->get_tburst()[it];
                _freqs[ii]=bw0->getPars()->p_syn_a->total_rad.e[ifreq];
                ii++;
            }
        }
        Output::addVector(_times_for_freq,"times_freqs",file);
        Output::addVector(_freqs,"freqs",file);


        /// make vectors for time and freq with the same structure as emissivity and absorption
        if (bw0->getPars()->p_syn_a->m_eleMethod!=METHODS_SHOCK_ELE::iShockEleAnalyt){
            Vector _times_freq_gam(bw0->get_tburst().size() * bw0->getPars()->p_syn_a->ele.numbins, 0.);
            Vector _gams (bw0->get_tburst().size() * bw0->getPars()->p_syn_a->ele.numbins, 0.);
            ii=0;
            for (size_t it = 0; it < bw0->get_tburst().size(); it++) {
                for(size_t igam = 0; igam < bw0->getPars()->p_syn_a->ele.numbins; igam++) {
                    _times_freq_gam[ii]=bw0->get_tburst()[it];
                    _gams[ii] = bw0->getPars()->p_syn_a->ele.e[igam];
                    ii++;
                }
            }
            Output::addVector(_times_freq_gam,"times_gams",file);
            Output::addVector(_gams,"gams",file);
        }


        /// add attributes
        std::unordered_map<std::string, double> attrs{
                {"nshells", nshells_ },
                {"nlayers", nlayers_ },
                {"ntimes", t_arr.size() },
                {"ncells", ncells() }
        };

        /// add attributes read from parfile
        for (auto& [key, value]: main_pars) { attrs[key] = value; }
        for (auto& [key, value]: ej_pars) { attrs[key] = value; }
        Output::addStrDbMap(attrs, file);

        file.close();
    }

    /// OUTPUT light curves
    void computeSaveEjectaLightCurveAnalytic_new(
            std::string workingdir, std::string fname, std::string fname_shells_layers,
            Vector lc_times, Vector lc_freqs, StrDbMap & main_pars, StrDbMap & ej_pars, bool lc_freq_to_time){

        Vector _times, _freqs;
        cast_times_freqs(lc_times,lc_freqs,_times,_freqs,lc_freq_to_time,p_log);
        (*p_log)(LOG_INFO,AT) << "Computing and saving Ejecta light curve with analytic synchrotron...\n";

        VecVector out{}; // [ish*il][it*inu]
        out.resize(nshells() * nlayers());
        for (auto & arr : out)
            arr.resize(_times.size(), 0.);

        evalEjectaLightCurves(out, _times, _freqs);

        /// Collect total_rad flux at a given time/freq from all shells/layers
        Vector total_fluxes (_times.size(), 0.0);
        for (size_t itnu = 0; itnu < _times.size(); ++itnu) {
            size_t ii = 0;
            for (size_t ish = 0; ish < nshells(); ish++) {
                for (size_t il = 0; il < nlayers(); il++) {
                    total_fluxes[itnu] += out[ii][itnu];
                    ii++;
                }
            }
        }

        std::vector<std::string> other_names { "times", "freqs", "total_fluxes" };
        VecVector out_data {_times, _freqs, total_fluxes};

        /// save the total_rad fluxes (sparse output)
        std::unordered_map<std::string,double> attrs{ {"nshells", nshells()}, {"nlayers", nlayers()} };
        for (auto& [key, value]: main_pars) { attrs[key] = value; }
        for (auto& [key, value]: ej_pars) { attrs[key] = value; }
        p_out->VectorOfVectorsH5(out_data, other_names, workingdir+fname,  attrs);

        /// save light curve for each shell and layer (dense output)
        if (fname_shells_layers == "none")
            return;
        std::vector<std::string> group_names;
//        VecVector total_fluxes_shell_layer(nshells()*nlayers());
        size_t ii = 0;
        for (size_t ishell = 0; ishell < nshells(); ++ishell)
            for (size_t ilayer = 0; ilayer < nlayers(); ++ilayer)
                group_names.emplace_back("shell=" + std::to_string(ishell) + " layer=" + std::to_string(ilayer));


        out.emplace_back(_times);
        out.emplace_back(_freqs);
        out.emplace_back(total_fluxes);

        group_names.emplace_back("times");
        group_names.emplace_back("freqs");
        group_names.emplace_back("total_fluxes");
        p_out->VectorOfVectorsH5(out, group_names, workingdir+fname,  attrs);
    }

    /// INPUT dynamics
    void loadEjectaBWDynamics(std::string workingdir, std::string fname){
        if (!std::experimental::filesystem::exists(workingdir+fname)) {
            (*p_log)(LOG_ERR, AT) << "File not found. " + workingdir + fname << "\n";
            exit(1);
        }

        Exception::dontPrint();
        H5std_string FILE_NAME(workingdir+fname);
        H5File file(FILE_NAME, H5F_ACC_RDONLY);
        size_t nshells_ = (size_t)getDoubleAttr(file,"nshells");
        size_t nlayers_ = (size_t)getDoubleAttr(file, "nlayers");
        size_t ntimes_ = (size_t)getDoubleAttr(file, "ntimes");
        if (nshells_ != nshells()){
            (*p_log)(LOG_ERR,AT) << "Wring attribute: nshells_="<<nshells_<<" expected nshells="<<nshells()<<"\n";
        }
        if (nlayers_ != nlayers()){
            (*p_log)(LOG_ERR,AT) << "Wring attribute: nlayers_="<<nlayers_<<" expected nlayers_="<<nlayers()<<"\n";
        }

        auto & models = getShells();
        for (size_t il = 0; il < nlayers(); il++) {
            size_t n_empty_shells = 0;
            for (size_t ish = 0; ish < nshells(); ish++) {
                auto &bw = models[il]->getBW(ish);
                for (auto & iv : BW::VARS.at(bw->getPars()->m_type)) {
                    std::string key = "shell=" + std::to_string(ish) + " layer=" + std::to_string(il);
                    if (!pathExists(file.getId(), key))
                        continue;
                    auto grp = file.openGroup(key);
                    auto & vec = bw->getData(static_cast<BW::Q>(iv));
                    if (!vec.empty()){
                        (*p_log)(LOG_ERR,AT) << " container is not isEmpty\n";
                    }

                    Output::load1DDataFromGroup(BW::VARNAMES[iv],vec,grp);

                    if ( bw->getData(static_cast<BW::Q>(iv)).empty() ){
                        (*p_log)(LOG_ERR,AT) << " Failed loading v_n="<< key << " group="<<key<<"\n";
                        exit(1);
                    }
                    grp.close();
                }
                if (bw->getData(BW::iR)[0] == 0){
//                    (*p_log)(LOG_WARN,AT) << "Loaded not evolved shell [il="<<il<<", "<<"ish="<<ish<<"] \n";
                    n_empty_shells+=1;
                }
                bw->checkEvolutionEnd();
                bw->getPars()->nr = bw->getData()[BW::iR].size();
            }
            (*p_log)(LOG_INFO,AT) << "Loaded [il="<<il<<"] empty shells "<<n_empty_shells<<"\n";
        }
        file.close();
//        if ( p_cumShells[0]->getBW(0)->getData()[BW::iR][0] == 0 ){
//            std::cout << p_cumShells[0]->getBW(0)->getData()[BW::iR] << "\n";
//            std::cout << " faild" << std::endl;
//            exit(1);
//        }

#if 0
        LoadH5 ldata;
        ldata.setFileName(workingdir+fname);
        ldata.setVarName("nshells");
        double nshells = ldata.getDoubleAttr("nshells");
        double m_nlayers = ldata.getDoubleAttr("m_nlayers");
        auto & models = getShells();
        for (size_t ish = 0; ish < nshells-1; ish++){
            for (size_t il = 0; il < m_nlayers-1; il++){
                auto & bw = models[il]->getBW(ish);
                for (size_t ivar = 0; ivar < BW::m_vnames.size(); ivar++) {
                    std::string key = "shell=" + std::to_string(ish)
                                    + " layer=" + std::to_string(il)
                                    + " key=" + BW::m_vnames[ivar];
                    ldata.setVarName(BW::m_vnames[ivar]);
                    bw->getData().emplace_back( std::move( ldata.getDataVDouble() ) );
                    (*p_log)(LOG_INFO,AT) << "Reading: "<<key<<"\n";
//                    bw->getData(static_cast<BW::Q>(ivar))
//                        = std::move( ldata.getDataVDouble() );
                }

//                auto & bw = models[il]->getBW(ish);
//                std::string
//                bw->getData()[]
            }
        }
#endif
        (*p_log)(LOG_INFO,AT)<<" dynamics loaded successfully\n";

    }

#if 0
    void updateEjectaObsPars(StrDbMap pars) {

        auto & models = getShells();

//        size_t nshells = nshells();
//        size_t m_nlayers = p_cumShells->m_nlayers();
//        size_t ncells =  (int)p_cumShells->ncells();

        (*p_log)(LOG_ERR,AT) << "Updating Ejecta observer pars...\n";
        size_t ii = 0;
        for (size_t ishell = 0; ishell < nshells(); ishell++) {
//            auto &struc = ejectaStructs.structs[ishell];
//            size_t n_layers_ej = struc.m_nlayers;//(p_pars->ej_method_eats == LatStruct::i_pw) ? struc.nlayers_pw : struc.nlayers_a ;
            for (size_t ilayer = 0; ilayer < nlayers(); ilayer++) {
                auto & model = getShells()[ilayer]->getBW(ishell);//ejectaModels[ishell][ilayer];
                model->getFsEATS()->updateObsPars(pars);
                ii++;
            }
        }
//        p_pars->is_ejecta_obs_pars_set = true;
    }
    bool evalOptDepthsAlongLineOfSight(double & frac, double ctheta, double r,
                                       double phi, double theta,
                                       double phi_obs, double theta_obs, double r_obs,
                                       double mu, double time, double freq, size_t ilpwn, void * params){

        // Convert spherical coordinates to Cartesian coordinates
        double x1 = r * std::sin(phi) * std::cos(ctheta);
        double y1 = r * std::sin(phi) * std::sin(ctheta);
        double z1 = r * std::cos(phi);

        // do the same for observer
        double x3 = r_obs * std::sin(phi_obs) * std::cos(theta_obs);
        double y3 = r_obs * std::sin(phi_obs) * std::sin(theta_obs);
        double z3 = r_obs * std::cos(phi_obs);

        /// Calculate the direction vector of the line between the two points
        double dx = x3 - x1;
        double dy = y3 - y1;
        double dz = z3 - z1;

        /// iterate over all layers and shells and find for each layer a shell that lies on the line of sight
        double tau_comp=0., tau_BH=0., tau_bf=0.;
        bool found = false;
        double r_ej_max = 0;
        size_t tot_nonzero_layers = 0;
        size_t tot_nonzero_shells = 0;
        std::vector<std::vector<char>> status(nlayers());
        for (auto & arr : status)
            arr.resize(nshells(), '-');
        Vector tobs_ej_max (nlayers(), 0.);
        /// ---------------------------------------------------
        for (size_t il = 0; il < nlayers(); il++){
            status.resizeEachImage(nshells());
            size_t nonzero_shells = 0;
            bool found_il = false;
            auto & cumshell = p_cumShells[il];
//            Vector cphis = EjectaID2::getCphiGridPW( il );
            if (cumshell->getPars()->n_active_shells==1){
                (*p_log)(LOG_ERR,AT)<<" not implemented\n";
                exit(1);
            }
            for (size_t ish = 0; ish < cumshell->getPars()->n_active_shells-1; ish++){
                auto & bw = cumshell->getBW(ish);
                size_t idx = ish;//cumshell->getIdx()[ish]; // TODO assume sorted shells (after evolution)
                auto & bw_next = cumshell->getBW(ish+1);
                size_t idx_next = ish+1;//cumshell->getIdx()[ish+1];// TODO assume sorted shells (after evolution)

                if ((bw->getPars()->i_end_r==0)||(bw_next->getPars()->i_end_r==0)) {
//                    (*p_log)(LOG_WARN,AT) << "[il="<<il<<" ish="<<ish<<"]"
//                        << " Skipping as bw->getPars()->i_end_r="<<bw->getPars()->i_end_r
//                        <<" and bw_next->getPars()->i_end_r="<<bw_next->getPars()->i_end_r
//                        <<"\n";
                    status[il][ish] = 'n';
                    continue;
                }
                bw->getFsEATS()->parsPars(time, freq,
                                          bw->getPars()->theta_c_l, bw->getPars()->theta_c_h,
                                          0., M_PI, obsAngle);
                bw->getFsEATS()->check_pars();

                // get BW (a cell) properties
                double cphi = 0. ; // We don't care about phi diretion due to symmetry
                double ctheta_cell = bw->getPars()->ctheta0;//m_data[BW::Q::ictheta][0]; //cthetas[0];
#if 1
                size_t ia=0, ib=0;
                bool is_in_time = bw->getFsEATS()->evalEATSindexes(ia,ib,time,theta_obs, ctheta_cell,cphi,obsAngle);
                Vector & ttobs = bw->getFsEATS()->getTobs();
                if (!is_in_time){
//                    (*p_log)(LOG_WARN,AT) << "[il="<<il<<" ish="<<ish<<"]" << " Skipping as tobs="<<time
//                        <<" while ttobs is in["<<ttobs[0]<<", "<<ttobs[bw->getPars()->i_end_r-1]<<"] \n";
                    status[il][ish] = 'T';
                    continue;
                }
#endif
                /// interpolate the exact radial position of the blast that corresponds to the req. obs time
                double r_cell = interpSegLog(ia, ib, time, ttobs, bw->getData(BW::Q::iEJr));
                if ( r_cell > r_ej_max ) r_ej_max = r_cell;
                double rho_ej_cell = interpSegLog(ia, ib, time, ttobs, bw->getData(BW::Q::iEJrho));
                double delta_ej_cell = interpSegLog(ia, ib, time, ttobs, bw->getData(BW::Q::iEJdelta));
                if ((rho_ej_cell<=0.)||(!std::isfinite(rho_ej_cell))||(delta_ej_cell<=0)||(!std::isfinite(delta_ej_cell))){
                    (*p_log)(LOG_ERR,AT) << "[il="<<il<<" ish="<<ish<<"]"<<" error in opt depth along line of sight\n";
                    exit(1);
                }
                if ((r >= r_cell)){
//                    (*p_log)(LOG_WARN,AT) << "[il="<<il<<" ish="<<ish<<"]" << " Skipping as r_pwn="<<r
//                     <<" > r_ej_cell="<<r_cell<<" Overall, r_ej_max="<<bw->getData(BW::Q::iEJr)[bw->getPars()->i_end_r-1]<<"\n";
                    status[il][ish] = 'R';
                    continue;
//                    exit(1);
                }

                double e_gamma = freq*4.1356655385381E-15*CGS::EV_TO_ERG;
                double mu_e = bw->getPars()->mu_e;
                double Z_eff = bw->getPars()->Z_eff;
                int opacitymode = bw->getPars()->opacitymode;
                double albd_fac = bw->getPars()->albd_fac;



                bw_next->getFsEATS()->parsPars(time, freq,
                                          bw_next->getPars()->theta_c_l, bw_next->getPars()->theta_c_h,
                                          0., M_PI, obsAngle);
#if 1
                bw_next->getFsEATS()->check_pars();
                is_in_time = bw_next->getFsEATS()->evalEATSindexes(ia,ib,time,theta_obs, ctheta_cell,cphi,obsAngle);
                ttobs = bw_next->getFsEATS()->getTobs();
                if (!is_in_time){
//                    (*p_log)(LOG_WARN,AT) << "[il="<<il<<" ish="<<ish<<"]" << " Skipping as tobs="<<time
//                                         <<" while ttobs is in["<<ttobs[0]<<", "<<ttobs[bw_next->getPars()->i_end_r-1]<<"] \n";
                    status[il][ish] = 't';
                    continue;
                }
#endif
                double r_cell_next =interpSegLog(ia,ib, time, ttobs, bw_next->getData(BW::Q::iEJr));

                if ((r >= r_cell_next)){
//                    (*p_log)(LOG_WARN,AT)  << "[il="<<il<<" ish="<<ish<<"]" << " r > r_\n";
                    status[il][ish] = 'r';
                    continue;
//                    exit(1);
                }

                // Calculate the intersection point of the line with the middle sphere
                double a = dx*dx + dy*dy + dz*dz;
                double b = 2. * (x1*dx + y1*dy + z1*dz);
                double c = x1*x1 + y1*y1 + z1*z1 - r_cell*r_cell;
                double disc = b*b - 4.*a*c;
                double t1 = (-b - std::sqrt(disc)) / (2.*a);
                double t2 = (-b + std::sqrt(disc)) / (2.*a);
                double x = x1 + t2*dx;
                double y = y1 + t2*dy;
                double z = z1 + t2*dz;

                double r_ = std::sqrt(x*x + y*y + z*z);
                double theta_ = std::atan(y/x);
                double phi_ = std::acos(z / r);


                // Calculate the intersection point of the line with the middle sphere
                double a_next = dx*dx + dy*dy + dz*dz;
                double b_next = 2. * (x1*dx + y1*dy + z1*dz);
                double c_next = x1*x1 + y1*y1 + z1*z1 - r_cell*r_cell;
                double disc_next = b_next*b_next - 4.*a_next*c;
                double t1_next = (-b_next - std::sqrt(disc_next)) / (2.*a_next);
                double t2_next = (-b_next + std::sqrt(disc_next)) / (2.*a_next);
                double x_next = x1 + t2_next*dx;
                double y_next = y1 + t2_next*dy;
                double z_next = z1 + t2_next*dz;

                double r_next = std::sqrt(x_next*x_next + y_next*y_next + z_next*z_next);
                double theta_next = std::atan(y_next/x_next);
                double phi_next = std::acos(z_next / r_cell_next);

                if (((theta_ > bw->getPars()->theta_c_l) && (theta_ <=bw->getPars()->theta_c_h)) &&
                    ((theta_next > bw_next->getPars()->theta_c_l) && (theta_next <=bw_next->getPars()->theta_c_h))){
                    /// --- optical depth due to compton scattering
                    double Kcomp = PWNradiationMurase::sigma_kn(e_gamma)/mu_e/CGS::M_PRO;
                    double tau_comp_ = rho_ej_cell*delta_ej_cell*Kcomp;
                    tau_comp+=tau_comp_;
                    /// optical depth of BH pair production
                    double KBH = (1.0+Z_eff)*PWNradiationMurase::sigma_BH_p(e_gamma)/mu_e/CGS::M_PRO;
                    double tau_BH_ = rho_ej_cell*delta_ej_cell*KBH;
                    tau_BH+=tau_BH_;
                    /// The photoelectric absorption at high energies is taken into account, using the bound–free opacity
                    double Kbf = (1.0-albd_fac)*PWNradiationMurase::kappa_bf(e_gamma, Z_eff, opacitymode);
                    double tau_bf_ = rho_ej_cell*delta_ej_cell*Kbf;
                    tau_bf+=tau_bf_;
//                    (*p_log)(LOG_INFO,AT) << "[il="<<il<<" ish="<<ish<<"]"
//                        << " tau_comp="<<tau_comp_<<" tau_BH="<<tau_BH_<<" tau_bf="<<tau_bf_
//                        << " | case 1 | " << "theta_="<<theta_<<" is in ["
//                        << bw->getPars()->theta_c_l<< ", "<<bw->getPars()->theta_c_h<<"] "
//                        << " and theta_next="<<theta_next<<" is in ["<<bw_next->getPars()->theta_c_l
//                        << ", " << bw_next->getPars()->theta_c_h <<"] \n";
                    status[il][ish] = '1';
//                    double tau_abs = (1.0+PWNradiationMurase::gamma_inelas_Compton(e_gamma))*(tau_BH+tau_bf);
//                    double tau_eff = sqrt((tau_abs+tau_comp_)*tau_abs);
                    found_il= true;
                }

                if (((theta_ > bw->getPars()->theta_c_l) && (theta_ <=bw->getPars()->theta_c_h)) &&
                    ((theta_next < bw_next->getPars()->theta_c_l) || (theta_next > bw_next->getPars()->theta_c_h))){
                    /// --- Bottom entrance only // TODO compute how much light travelled within this cell (not just 0.5!)
                    /// --- optical depth due to compton scattering
                    double Kcomp = PWNradiationMurase::sigma_kn(e_gamma)/mu_e/CGS::M_PRO;
                    double tau_comp_ = rho_ej_cell*delta_ej_cell*Kcomp;
                    tau_comp+=tau_comp_;
                    /// optical depth of BH pair production
                    double KBH = (1.0+Z_eff)*PWNradiationMurase::sigma_BH_p(e_gamma)/mu_e/CGS::M_PRO;
                    double tau_BH_ = rho_ej_cell*delta_ej_cell*KBH;
                    tau_BH+=tau_BH_;
                    /// The photoelectric absorption at high energies is taken into account, using the bound–free opacity
                    double Kbf = (1.0-albd_fac)*PWNradiationMurase::kappa_bf(e_gamma, Z_eff, opacitymode);
                    double tau_bf_ = rho_ej_cell*delta_ej_cell*Kbf;
                    tau_bf+=tau_bf_;
//                    (*p_log)(LOG_INFO,AT) << "[il="<<il<<" ish="<<ish<<"]"
//                                          << " tau_comp="<<tau_comp_<<" tau_BH="<<tau_BH_<<" tau_bf="<<tau_bf_
//                                          << " | case 2 | " << "theta_="<<theta_<<" is in ["
//                                          << bw->getPars()->theta_c_l<< ", "<<bw->getPars()->theta_c_h<<"] "
//                                          << " but theta_next="<<theta_next<<" is NOT in ["<<bw_next->getPars()->theta_c_l
//                                          << ", " << bw_next->getPars()->theta_c_h <<"] \n";
                    status[il][ish] = '2';
                    found_il= true;
                }

                if (((theta_ < bw->getPars()->theta_c_l) || (theta_ > bw->getPars()->theta_c_h)) &&
                    ((theta_next > bw_next->getPars()->theta_c_l) && (theta_next <= bw_next->getPars()->theta_c_h))){
                    /// --- Top exit only // TODO compute how much light travelled within this cell (not just 0.5!)
                    /// --- optical depth due to compton scattering
                    double Kcomp = PWNradiationMurase::sigma_kn(e_gamma)/mu_e/CGS::M_PRO;
                    double tau_comp_ = rho_ej_cell*delta_ej_cell*Kcomp;
                    tau_comp+=tau_comp_;
                    /// optical depth of BH pair production
                    double KBH = (1.0+Z_eff)*PWNradiationMurase::sigma_BH_p(e_gamma)/mu_e/CGS::M_PRO;
                    double tau_BH_ = rho_ej_cell*delta_ej_cell*KBH;
                    tau_BH+=tau_BH_;
                    /// The photoelectric absorption at high energies is taken into account, using the bound–free opacity
                    double Kbf = (1.0-albd_fac)*PWNradiationMurase::kappa_bf(e_gamma, Z_eff, opacitymode);
                    double tau_bf_ = rho_ej_cell*delta_ej_cell*Kbf;
                    tau_bf+=tau_bf_;
//                    (*p_log)(LOG_INFO,AT) << "[il="<<il<<" ish="<<ish<<"]"
//                                          << " tau_comp="<<tau_comp_<<" tau_BH="<<tau_BH_<<" tau_bf="<<tau_bf_
//                                          << " | case 3 | " << "theta_="<<theta_<<" is NOT in ["
//                                          << bw->getPars()->theta_c_l<< ", "<<bw->getPars()->theta_c_h<<"] "
//                                          << " but theta_next="<<theta_next<<" is in ["<<bw_next->getPars()->theta_c_l
//                                          << ", " << bw_next->getPars()->theta_c_h <<"] \n";
                    status[il][ish] = '3';
                    found_il= true;
                }

                if (((theta_ < bw->getPars()->theta_c_l) || (theta_ > bw->getPars()->theta_c_h)) &&
                    ((theta_next < bw_next->getPars()->theta_c_l) || (theta_next > bw_next->getPars()->theta_c_h))){
                    /// --- No intersection
//                    (*p_log)(LOG_INFO,AT) << "[il="<<il<<" ish="<<ish<<"]"
//                        << " tau_comp="<<0<<" tau_BH="<<0<<" tau_bf="<<0
//                        << " | case 4 |"<< " theta_="<<theta_<<" is NOT in ["
//                        << bw->getPars()->theta_c_l<< ", "<<bw->getPars()->theta_c_h<<"] "
//                        << " and theta_next="<<theta_next<<" is NOT in ["<<bw_next->getPars()->theta_c_l
//                        << ", " << bw_next->getPars()->theta_c_h <<"] \n";
                    status[il][ish] = '4';

                }
                if (found_il){
                    tot_nonzero_shells+=1;
                    nonzero_shells+=1;
                }
//                (*p_log)(LOG_INFO,AT) << "[il="<<il<<" ish="<<ish<<"]" << " tau_comp="<<tau_comp<<" tau_BH="<<tau_BH<<" tau_bf="<<tau_bf<<"\n";
//                int __x = 1;
            }
            if(nonzero_shells>0)
                tot_nonzero_layers+=1;
            if (found_il)
                found = true;
        }

        auto & stream = std::cout;
        stream << "------ t="<<time<<", nu="<<freq<<" | PWN: il="<<ilpwn<<" r="<<r<<" ctheta="<<ctheta<<" ---\n";
        for (size_t il = 0; il < nlayers(); il++){
            stream << "il="<<il<<"| ";
            for (size_t ish = 0; ish < nshells(); ish++) {
                if (ish == nshells()-1)
                    stream << status[il][ish] << " | \n";
                else
                    stream << status[il][ish] << " | ";
            }
        }
        stream << "---------------------------------------------------------------------------------------------"
                  "------\n";

        if (r > r_ej_max){
            (*p_log)(LOG_ERR,AT) << "[ilpw="<<ilpwn<<", ctheta="<<ctheta<<", r="<<r<<" > "<<"r_ej_max"<<r_ej_max<<"] "<<"\n";
            int _x = 1;
        }
        if (!found){
            (*p_log)(LOG_INFO,AT) << "[ilpw="<<ilpwn<<", ctheta="<<ctheta<<", r="<<r<<", "<<"t="<<time<<"] "
                <<" not found layer/shell in which optical depth can be computed"<<"\n";
            return false;
        }


        /// Combine individual optical depths into fraction
        double power_Compton=0.0;
        if (tau_comp > 1.0)
            power_Compton = tau_comp*tau_comp;
        else
            power_Compton = tau_comp;

        double e_gamma = freq*4.1356655385381E-15*CGS::EV_TO_ERG; // TODO this has to be red-shifted!
        frac = exp(-(tau_BH+tau_bf))
               * (exp(-(tau_comp)) + (1.0-exp(-(tau_comp)))
                                         * pow(1.0-PWNradiationMurase::gamma_inelas_Compton(e_gamma),power_Compton));
        (*p_log)(LOG_INFO,AT) << "[ilpw="<<ilpwn<<", ctheta="<<ctheta<<", r="<<r<<", "<<"t="<<time<<"] "
                                         << " tau_comp="<<tau_comp<<" tau_BH="<<tau_BH<<" tau_bf="<<tau_bf<<" -> frac="<<frac<<"\n";
        if ((frac<0)||(frac>1)||!std::isfinite(frac)){
            (*p_log)(LOG_ERR,AT) << "[ilpw="<<ilpwn<<", ctheta="<<ctheta<<", r="<<r<<", "<<"t="<<time<<"] "
                <<" tau_comp="<<tau_comp<<" tau_BH="<<tau_BH<<" tau_bf="<<tau_bf<<" -> frac="<<frac<<"\n";
//            exit(1);
            return false;
        }

        return true;
    }
#endif

#if 0 // uses t_obs
    bool evalOptDepthsAlongLineOfSight(double & frac, double & tau_comp, double & tau_BH, double tau_bf,
                                       double ctheta, double r, double phi, double theta,
                                       double phi_obs, double theta_obs, double r_obs,
                                       double time, double freq, size_t info_ilpwn, void * params) {

        // Convert spherical coordinates to Cartesian coordinates
        double x1 = r * std::sin(phi) * std::cos(ctheta);
        double y1 = r * std::sin(phi) * std::sin(ctheta);
        double z1 = r * std::cos(phi);

        // do the same for observer
        double x3 = r_obs * std::sin(phi_obs) * std::cos(theta_obs);
        double y3 = r_obs * std::sin(phi_obs) * std::sin(theta_obs);
        double z3 = r_obs * std::cos(phi_obs);

        /// Calculate the direction vector of the line between the two points
        double dx = x3 - x1;
        double dy = y3 - y1;
        double dz = z3 - z1;

        /// iterate over all layers and shells and find for each layer a shell that lies on the line of sight
//        double tau_comp=0., tau_BH=0., tau_bf=0.;
        bool found = false;
        double r_ej_max = 0;
        size_t tot_nonzero_layers = 0;
        size_t tot_nonzero_shells = 0;

        /// ---------------------------------------------------
        for (size_t il = 0; il < nlayers(); il++){
//            status.resizeEachImage(nshells());
            size_t nonzero_shells = 0;
            bool found_il = false;
            auto & cumshell = p_cumShells[il];
//            Vector cphis = EjectaID2::getCphiGridPW( il );
//            if (cumshell->getPars()->n_active_shells==1){
//                (*p_log)(LOG_ERR,AT)<<" not implemented\n";
//                exit(1);
//            }
            for (size_t ish = 0; ish < nshells(); ish++){
                auto & bw = cumshell->getBW(ish);
                size_t idx = ish;//cumshell->getIdx()[ish]; // TODO assume sorted shells (after evolution)
//                auto & bw_next = cumshell->getBW(ish+1);
                size_t idx_next = ish+1;//cumshell->getIdx()[ish+1];// TODO assume sorted shells (after evolution)
                if (bw->getPars()->i_end_r==0) {
//                    (*p_log)(LOG_WARN,AT) << "[il="<<il<<" ish="<<ish<<"]"
//                        << " Skipping as bw->getPars()->i_end_r="<<bw->getPars()->i_end_r
//                        <<" and bw_next->getPars()->i_end_r="<<bw_next->getPars()->i_end_r
//                        <<"\n";
                    status[il][ish] = 'N';
                    continue;
                }
//                if (bw_next->getPars()->i_end_r==0) {
////                    (*p_log)(LOG_WARN,AT) << "[il="<<il<<" ish="<<ish<<"]"
////                        << " Skipping as bw->getPars()->i_end_r="<<bw->getPars()->i_end_r
////                        <<" and bw_next->getPars()->i_end_r="<<bw_next->getPars()->i_end_r
////                        <<"\n";
//                    status[il][ish] = 'n';
//                    continue;
//                }
                bw->getFsEATS()->parsPars(time, freq,
                                          bw->getPars()->theta_c_l, bw->getPars()->theta_c_h,
                                          0., M_PI, obsAngle);
                bw->getFsEATS()->check_pars();

                // get BW (a cell) properties
                double cphi = 0. ; // We don't care about phi diretion due to symmetry
                double ctheta_cell = bw->getPars()->ctheta0;//m_data[BW::Q::ictheta][0]; //cthetas[0];

                size_t ia=0, ib=0;
                bool is_in_time = bw->getFsEATS()->evalEATSindexes(ia,ib,time,theta_obs, ctheta_cell,cphi,obsAngle);
                Vector & ttobs = bw->getFsEATS()->getTobs();
                if (!is_in_time){
//                    (*p_log)(LOG_WARN,AT) << "[il="<<il<<" ish="<<ish<<"]" << " Skipping as tobs="<<time
//                        <<" while ttobs is in["<<ttobs[0]<<", "<<ttobs[bw->getPars()->i_end_r-1]<<"] \n";
                    status[il][ish] = 'T';
                    continue;
                }

                /// interpolate the exact radial position of the blast that corresponds to the req. obs time
                double r_cell = interpSegLog(ia, ib, time, ttobs, bw->getData(BW::Q::iEJr));
                if ( r_cell > r_ej_max ) r_ej_max = r_cell;
                double rho_ej_cell = interpSegLog(ia, ib, time, ttobs, bw->getData(BW::Q::iEJrho));
                double delta_ej_cell = interpSegLog(ia, ib, time, ttobs, bw->getData(BW::Q::iEJdelta));
                if ((rho_ej_cell<=0.)||(!std::isfinite(rho_ej_cell))||(delta_ej_cell<=0)||(!std::isfinite(delta_ej_cell))){
                    (*p_log)(LOG_ERR,AT) << "[il="<<il<<" ish="<<ish<<"]"<<" error in opt depth along line of sight\n";
                    exit(1);
                }
                if ((r >= r_cell)){
//                    (*p_log)(LOG_WARN,AT) << "[il="<<il<<" ish="<<ish<<"]" << " Skipping as r_pwn="<<r
//                     <<" > r_ej_cell="<<r_cell<<" Overall, r_ej_max="<<bw->getData(BW::Q::iEJr)[bw->getPars()->i_end_r-1]<<"\n";
                    status[il][ish] = 'R';
                    continue;
//                    exit(1);
                }

                double e_gamma = freq*4.1356655385381E-15*CGS::EV_TO_ERG;
                double mu_e = bw->getPars()->mu_e;
                double Z_eff = bw->getPars()->Z_eff;
                int opacitymode = bw->getPars()->opacitymode;
                double albd_fac = bw->getPars()->albd_fac;


                // Calculate the intersection point of the line with the middle sphere
                double a = dx*dx + dy*dy + dz*dz;
                double b = 2. * (x1*dx + y1*dy + z1*dz);
                double c = x1*x1 + y1*y1 + z1*z1 - r_cell*r_cell;
                double disc = b*b - 4.*a*c;
                double t1 = (-b - std::sqrt(disc)) / (2.*a);
                double t2 = (-b + std::sqrt(disc)) / (2.*a);
                double x = x1 + t2*dx;
                double y = y1 + t2*dy;
                double z = z1 + t2*dz;

                double r_ = std::sqrt(x*x + y*y + z*z);
                double theta_ = std::atan(y/x);
                double phi_ = std::acos(z / r);

                if (((theta_ > bw->getPars()->theta_c_l) && (theta_ <=bw->getPars()->theta_c_h))){
                    /// --- optical depth due to compton scattering
                    double Kcomp = PWNradiationMurase::sigma_kn(e_gamma)/mu_e/CGS::M_PRO;
                    double tau_comp_ = rho_ej_cell*delta_ej_cell*Kcomp;
                    tau_comp+=tau_comp_;
                    /// optical depth of BH pair production
                    double KBH = (1.0+Z_eff)*PWNradiationMurase::sigma_BH_p(e_gamma)/mu_e/CGS::M_PRO;
                    double tau_BH_ = rho_ej_cell*delta_ej_cell*KBH;
                    tau_BH+=tau_BH_;
                    /// The photoelectric absorption at high energies is taken into account, using the bound–free opacity
                    double Kbf = (1.0-albd_fac)*PWNradiationMurase::kappa_bf(e_gamma, Z_eff, opacitymode);
                    double tau_bf_ = rho_ej_cell*delta_ej_cell*Kbf;
                    tau_bf+=tau_bf_;
#if 0
                    (*p_log)(LOG_INFO,AT) << "[il="<<il<<" ish="<<ish<<"]"
                        << " rho="<<rho_ej_cell<<" delta="<<delta_ej_cell
                        << " tau_comp="<<tau_comp_<<" tau_BH="<<tau_BH_<<" tau_bf="<<tau_bf_
                        << " | case 1 | " << "theta_="<<theta_<<" is in ["
                        << bw->getPars()->theta_c_l<< ", "<<bw->getPars()->theta_c_h<<"] "
                        <<" is in ["<<bw_next->getPars()->theta_c_l << ", " << bw_next->getPars()->theta_c_h <<"] \n";
#endif
                    status[il][ish] = '1';
//                    double tau_abs = (1.0+PWNradiationMurase::gamma_inelas_Compton(e_gamma))*(tau_BH+tau_bf);
//                    double tau_eff = sqrt((tau_abs+tau_comp_)*tau_abs);
                    found_il= true;
                }
                if (((theta_ < bw->getPars()->theta_c_l) || (theta_ > bw->getPars()->theta_c_h))){
                    /// --- No intersection
//                    (*p_log)(LOG_INFO,AT) << "[il="<<il<<" ish="<<ish<<"]"
//                        << " tau_comp="<<0<<" tau_BH="<<0<<" tau_bf="<<0
//                        << " | case 4 |"<< " theta_="<<theta_<<" is NOT in ["
//                        << bw->getPars()->theta_c_l<< ", "<<bw->getPars()->theta_c_h<<"] "
//                        << " and theta_next="<<theta_next<<" is NOT in ["<<bw_next->getPars()->theta_c_l
//                        << ", " << bw_next->getPars()->theta_c_h <<"] \n";
                    status[il][ish] = '4';

                }
                if (found_il){
                    tot_nonzero_shells+=1;
                    nonzero_shells+=1;
                }
//                (*p_log)(LOG_INFO,AT) << "[il="<<il<<" ish="<<ish<<"]" << " tau_comp="<<tau_comp<<" tau_BH="<<tau_BH<<" tau_bf="<<tau_bf<<"\n";
//                int __x = 1;
            }
            if(nonzero_shells>0)
                tot_nonzero_layers+=1;
            if (found_il)
                found = true;
        }

#if 0
        auto & stream = std::cout;
        stream << "------ t="<<time<<", nu="<<freq<<" | PWN: il="<<info_ilpwn<<" r="<<r<<" ctheta="<<ctheta<<" ---\n";
        for (size_t il = 0; il < nlayers(); il++){
            stream << "il="<<il<<"| ";
            for (size_t ish = 0; ish < nshells(); ish++) {
                if (ish == nshells()-1)
                    stream << status[il][ish] << " | \n";
                else
                    stream << status[il][ish] << " | ";
            }
        }
        stream << "Tot non-zero layers = "<<tot_nonzero_layers<< " tot_non_layer_shells = "<<tot_nonzero_shells << "\n";
        stream << "---------------------------------------------------------------------------------------------"
                  "------\n";
#endif
        if (r > r_ej_max){
            (*p_log)(LOG_ERR,AT) << "[ilpw=" << info_ilpwn << ", ctheta=" << ctheta << ", r=" << r << " > " << "r_ej_max" << r_ej_max << "] " << "\n";
            int _x = 1;
        }
        if (!found){
            (*p_log)(LOG_INFO,AT) << "[ilpw=" << info_ilpwn << ", ctheta=" << ctheta << ", r=" << r << ", " << "t=" << time << "] "
                                  <<" not found layer/shell in which optical depth can be computed"<<"\n";
            return false;
        }


        /// Combine individual optical depths into fraction
        double power_Compton=0.0;
        if (tau_comp > 1.0)
            power_Compton = tau_comp*tau_comp;
        else
            power_Compton = tau_comp;

        double e_gamma = freq*4.1356655385381E-15*CGS::EV_TO_ERG; // TODO this has to be red-shifted!

        frac = exp(-(tau_BH+tau_bf))
               * (exp(-(tau_comp)) + (1.0-exp(-(tau_comp)))
                                     * pow(1.0-PWNradiationMurase::gamma_inelas_Compton(e_gamma),power_Compton));
//        (*p_log)(LOG_INFO,AT) << "[ilpw="<<info_ilpwn<<", ctheta="<<ctheta<<", r="<<r<<", "<<"t="<<time<<"] "
//                              << " tau_comp="<<tau_comp<<" tau_BH="<<tau_BH<<" tau_bf="<<tau_bf<<" -> frac="<<frac<<"\n";
        if ((frac<0)||(frac>1)||!std::isfinite(frac)){
            (*p_log)(LOG_ERR,AT) << "[ilpw=" << info_ilpwn << ", ctheta=" << ctheta << ", r=" << r << ", " << "t=" << time << "] "
                                 <<" tau_comp="<<tau_comp<<" tau_BH="<<tau_BH<<" tau_bf="<<tau_bf<<" -> frac="<<frac<<"\n";
//            exit(1);
            return false;
        }

        return true;
    }
#endif

};

#endif //SRC_EJECTA_H
