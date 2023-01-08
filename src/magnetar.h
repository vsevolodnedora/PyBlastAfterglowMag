//
// Created by vsevolod on 08/01/23.
//

#ifndef SRC_MAGNETAR_H
#define SRC_MAGNETAR_H

#include "pch.h"
#include "logger.h"

class Magnetar{
    Array & m_tb_arr;
    VecArray m_data{}; // container for the solution of the evolution
    std::unique_ptr<logger> p_log;
public:
    /// RHS pars
    const static int neq = 0;
    std::vector<std::string> vars {  };
    size_t getNeq() { return neq; }
    enum Q_SOL {  };

    /// All variables
    enum Q {
        // -- dynamics ---
        iR, iRsh, irho, idrhodr, iGammaCBM, iGammaREL, idGammaCBMdr, idGammaRELdGamma, iPcbm, idPCBMdrho, iMCBM, iCSCBM,
        iGamma, ibeta, imom, iEint2, iU_p, itheta, ictheta, iErad2, iEsh2, iEad2, iM2,
        itcomov, itburst, itt, ithickness, iadi, irho2, iGammaFsh,
        // --- electrons  ---
        igm, igc, igM, iB, iTheta, iz_cool, ix, inprime, iacc_frac,
        // --- observables ---
        imu,
        // ---
        ijl, ir_dist
    };
    std::vector<std::string> m_vnames{
            ""
    };
    static constexpr size_t NVALS = 1; // number of variables to save
    inline Array & operator[](unsigned ll){ return m_data[ll]; }
    inline double & operator()(size_t ivn, size_t ir){ return m_data[ivn][ir]; }

    Magnetar( Array & t_grid, int loglevel ) : m_tb_arr(t_grid) {
        ///
    }

    Array & getTbGrid() {return m_tb_arr;}
    Array getTbGrid(size_t every_it) {
        if ((every_it == 1)||(every_it==0)) return m_tb_arr;
        Vector tmp{};
        for (size_t it = 0; it < m_tb_arr.size(); it = it + every_it){
            tmp.push_back(m_tb_arr[it]);
        }
        Array tmp2 (tmp.data(), tmp.size());
        return std::move(tmp2);
    }

    void setPars(StrDbMap & pars, StrStrMap & opts){

    }

    void setInitConditions( double * arr, size_t i ) {

    }

    void insertSolution( const double * sol, size_t it, size_t i ){

    }

    void evaluateRhs( double * out_Y, size_t i, double x, double const * Y ){

    }

    void addOtherVariables( size_t it ){

    }

    bool isSolutionOk( double * sol, size_t i ){
        return true;
    }



    void applyUnits( double * sol, size_t i ){

    }

};

#endif //SRC_MAGNETAR_H
