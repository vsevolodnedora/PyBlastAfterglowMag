//
// Created by vsevolod on 13/10/23.
//

#ifndef SRC_BASE_H
#define SRC_BASE_H


class ShockMicrophysics{
    VecVector & mD;
    std::unique_ptr<logger> p_log = nullptr;
public:

    ShockMicrophysics(VecVector & m_data, int loglevel) : mD(m_data){
        p_log = std::make_unique<logger>(std::cout, std::cerr, loglevel, "ShockMicrophysics");
    }

};




# if 0
class ShockMicrophysics{
    struct EatsPars{
        double eps_e=-1., eps_b=-1., eps_t=-1., p=-1., ksi_n=-1.;
    };
public:
    std::unique_ptr<EatsPars> p_pars = nullptr;
    std::unique_ptr<logger> p_log = nullptr;
    ShockMicrophysics(int loglevel){
        p_log = std::make_unique<logger>(std::cout, std::cerr, loglevel, "Base");
        p_pars = std::make_unique<EatsPars>();
    }

    void setPars(StrDbMap & pars, StrStrMap & opts, bool is_rs){
        // set parameters
        std::string fs_or_rs;
        if (is_rs)
            fs_or_rs += "_rs";

        p_pars->ksi_n = getDoublePar("ksi_n" + fs_or_rs, pars, AT, p_log, 1., false);//pars.at("ksi_n");
        p_pars->eps_e = getDoublePar("eps_e" + fs_or_rs, pars, AT, p_log, -1, true);//pars.at("eps_e");
        p_pars->eps_b = getDoublePar("eps_b" + fs_or_rs, pars, AT, p_log, -1, true);//pars.at("eps_b");
        p_pars->eps_t = getDoublePar("eps_t" + fs_or_rs, pars, AT, p_log, 0., true);//pars.at("eps_t");
        p_pars->p = getDoublePar("p" + fs_or_rs, pars, AT, p_log, -1, true);//pars.at("p");

    }
    void evalEleDist(double eprime, double Gamma, double Gamam_sh, double t_e, double n_prime){
        double p = p_pars->p;
        double eps_b = p_pars->eps_b;
        double eps_e = p_pars->eps_e;
        double eps_t = p_pars->eps_t;
        double ksi_n = p_pars->ksi_n;
    }
private:
    void checkPars(double nprime, double eprime){
        // check
        if (p_pars->eps_e <= 0){
            (*p_log)(LOG_ERR,AT) << " eps_e is not set (eps_e="<<p_pars->eps_e<<")\n";
            exit(1);
        }
        if (p_pars->eps_b <= 0){
            (*p_log)(LOG_ERR,AT) << " eps_b is not set (eps_e="<<p_pars->eps_b<<")\n";
            exit(1);
        }
        if (p_pars->eps_t < 0){
            (*p_log)(LOG_ERR,AT) << " eps_e is not set (eps_e="<<p_pars->eps_t<<") "
                                            << "(even through it is only needed in Marg21 model)\n";
            exit(1);
        }
        if (p_pars->p <= 0){
            (*p_log)(LOG_ERR,AT) << " eps_e is not set (eps_e="<<p_pars->p<<")\n";
            exit(1);
        }
        if (p_pars->ksi_n <= 0){
            (*p_log)(LOG_ERR,AT)<< " ksi_n is not set (ksi_n="<<p_pars->ksi_n<<")\n";
            exit(1);
        }
    }

};
#endif
#endif //SRC_BASE_H
