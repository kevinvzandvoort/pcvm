/*
    Model implementation for the diamond shaped SIS-type pneumococcal model
*/
#pragma once

int n_comps_prevalence = 4;
int n_comps_incidence = 2;

class VaccinationGroup : public BaseVaccinationGroup {
public:
    Compartment Sus, VT, NVT, B;
    arma::rowvec sus, vt, nvt, b, dsus, dvt, dnvt, db, inf_s_vt, inf_s_nvt, inf_vt_b, inf_nvt_b,
    clear_vt_s, clear_nvt_s, clear_b_vt, clear_b_nvt, incidence;

    // Constructor
    VaccinationGroup(int n_agrp, Rcpp::List vac_parms, arma::rowvec arate, arma::rowvec arate_corr) 
        :   BaseVaccinationGroup(n_agrp, vac_parms),
            //Note that last value in Compartment constructor is proportion of newborns
            // born into this compartment
            Sus(n_agrp, arate, arate_corr, 1.0),
            VT(n_agrp, arate, arate_corr, 0.0),
            NVT(n_agrp, arate, arate_corr, 0.0),
            B(n_agrp, arate, arate_corr, 0.0)
    {
            //Rcpp::Rcout << "DEBUG: Create VaccinationGroup" << std::endl;
            compartments = {&Sus, &VT, &NVT, &B}; // Add pointers to compartment objects
            sus = vt = nvt = b = dsus = dvt = dnvt = db = inf_s_vt = inf_s_nvt = inf_vt_b = inf_nvt_b =
            clear_vt_s = clear_nvt_s = clear_b_vt = clear_b_nvt = arma::rowvec(n_agrp, arma::fill::zeros);
            incidence = arma::rowvec(n_agrp*n_comps_incidence, arma::fill::zeros);
    }

    //Deconstructor
    ~VaccinationGroup(){  }

    void calculateDerivs(arma::rowvec &foi_vt, arma::rowvec &foi_nvt, double &comp, arma::rowvec &crate_vt, arma::rowvec &crate_nvt) {
        sus = Sus.getUpdatedValue();
        vt = VT.getUpdatedValue();
        nvt = NVT.getUpdatedValue();
        b = B.getUpdatedValue();
        
        //Infections from S to VT and NVT (note FOI values already account for travel between populations)
        inf_s_vt =  foi_vt  % (1.0 - vac_eff) % sus;
        inf_s_nvt = foi_nvt % sus;
    
        //Superinfections from VT and NVT to B
        inf_vt_b =  comp * foi_nvt % vt;
        inf_nvt_b = comp * foi_vt  % (1.0 - vac_eff) % nvt;
    
        //Clearance of carriage from VT and NVT to S
        clear_vt_s =  crate_vt  % vt;
        clear_nvt_s = crate_nvt % nvt;
    
        //Clearance of superinfection from B to VT and NVT
        clear_b_vt =  crate_nvt % b;
        clear_b_nvt = crate_vt  % b;
    
        //All rates are then combined to update each ODE for each compartment
        dsus = -(inf_s_vt+inf_s_nvt) +clear_vt_s+clear_nvt_s;
        Sus.updateDerivs(dsus);
        
        dvt = +inf_s_vt -inf_vt_b +clear_b_vt -clear_vt_s;
        VT.updateDerivs(dvt);
        
        dnvt = +inf_s_nvt -inf_nvt_b +clear_b_nvt -clear_nvt_s;
        NVT.updateDerivs(dnvt);
        
        db = +(inf_vt_b+inf_nvt_b) -(clear_b_vt+clear_b_nvt); 
        B.updateDerivs(db);
    }

    arma::rowvec& getIncidence(){
        incidence = arma::rowvec(n_agrp*n_comps_incidence, arma::fill::zeros);
        incidence.subvec(n_agrp*0, (n_agrp*1 - 1)) = get_inf_s_vt();
        incidence.subvec(n_agrp*1, (n_agrp*2 - 1)) = get_inf_s_nvt();
        incidence.subvec(n_agrp*2, (n_agrp*3 - 1)) = get_inf_vt_b();
        incidence.subvec(n_agrp*3, (n_agrp*4 - 1)) = get_inf_nvt_b();

        return incidence;
    }

    arma::rowvec& get_inf_s_vt(){
        return inf_s_vt;
    }
    arma::rowvec& get_inf_s_nvt(){ 
        return inf_s_nvt;
    }
    arma::rowvec& get_inf_vt_b(){
        return inf_vt_b;
    }
    arma::rowvec& get_inf_nvt_b(){
        return inf_nvt_b;
    }
};