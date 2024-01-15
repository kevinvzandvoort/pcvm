/*
    Model implementation for the diamond shaped SIS-type pneumococcal model
*/
#pragma once

int n_comps_prevalence = 6;
int n_comps_incidence = 2;

class VaccinationGroup : public BaseVaccinationGroup {
public:
    Compartment Sus, VT, NVT, B, VT2, NVT2;
    arma::rowvec sus, vt, nvt, vt2, b, nvt2, dsus, dvt, dnvt, dvt2, db, dnvt2, inf_s_vt, inf_s_nvt, inf_vt_vt2, inf_vt_b, inf_nvt_b, inf_nvt_nvt2,
    clear_vt_s, clear_nvt_s, clear_vt2_vt, clear_b_vt, clear_b_nvt, clear_nvt2_nvt, incidence;

    // Constructor
    VaccinationGroup(int n_agrp, Rcpp::List vac_parms, arma::rowvec arate, arma::rowvec arate_corr) 
        :   BaseVaccinationGroup(n_agrp, vac_parms),
            //Note that last value in Compartment constructor is proportion of newborns
            // born into this compartment
            Sus(n_agrp, arate, arate_corr, 1.0),
            VT(n_agrp, arate, arate_corr, 0.0),
            NVT(n_agrp, arate, arate_corr, 0.0),
            VT2(n_agrp, arate, arate_corr, 0.0),
            B(n_agrp, arate, arate_corr, 0.0),
            NVT2(n_agrp, arate, arate_corr, 0.0)
    {
            //Rcpp::Rcout << "DEBUG: Create VaccinationGroup" << std::endl;
            compartments = {&Sus, &VT, &NVT, &VT2, &B, &NVT2}; // Add pointers to compartment objects
            sus = vt = nvt = vt2 = b = nvt2 = dsus = dvt = dnvt = dvt2 = db = dnvt2 = inf_s_vt = inf_s_nvt = inf_vt_vt2 = inf_vt_b = inf_nvt_b = inf_nvt_nvt2 =
            clear_vt_s = clear_nvt_s = clear_vt2_vt = clear_b_vt = clear_b_nvt = clear_nvt2_nvt = arma::rowvec(n_agrp, arma::fill::zeros);
            incidence = arma::rowvec(n_agrp*n_comps_incidence, arma::fill::zeros);
    }

    //Deconstructor
    ~VaccinationGroup(){  }

    void calculateDerivs(arma::rowvec &foi_vt, arma::rowvec &foi_nvt, double &comp, arma::rowvec &crate_vt, arma::rowvec &crate_nvt) {
        sus = Sus.getUpdatedValue();
        vt = VT.getUpdatedValue();
        nvt = NVT.getUpdatedValue();
        vt2 = VT2.getUpdatedValue();
        b = B.getUpdatedValue();
        nvt2 = NVT2.getUpdatedValue();
        
        //Infections from S to VT and NVT (note FOI values already account for travel between populations)
        inf_s_vt =  foi_vt  % (1.0 - vac_eff) % sus;
        inf_s_nvt = foi_nvt % sus;
    
        //Superinfections
        inf_vt_vt2 =  comp * foi_vt % (1.0 - vac_eff) % vt;
        inf_vt_b =  comp * foi_nvt % vt;
        inf_nvt_b = comp * foi_vt  % (1.0 - vac_eff) % nvt;
        inf_nvt_nvt2 =  comp * foi_nvt % nvt;
    
        //Clearance of carriage from VT and NVT to S
        clear_vt_s =  crate_vt  % vt;
        clear_nvt_s = crate_nvt % nvt;
    
        //Clearance of superinfection
        clear_vt2_vt = 2.0 * crate_vt % vt2;
        clear_b_vt = crate_nvt % b;
        clear_b_nvt = crate_vt  % b;
        clear_nvt2_nvt = 2.0 * crate_nvt % nvt2;
    
        //All rates are then combined to update each ODE for each compartment
        dsus = -(inf_s_vt+inf_s_nvt) +clear_vt_s+clear_nvt_s;
        Sus.updateDerivs(dsus);
        
        dvt = +inf_s_vt -inf_vt_vt2 -inf_vt_b +clear_vt2_vt +clear_b_vt -clear_vt_s;
        VT.updateDerivs(dvt);
        
        dnvt = +inf_s_nvt -inf_nvt_nvt2 -inf_nvt_b +clear_nvt2_nvt +clear_b_nvt -clear_nvt_s;
        NVT.updateDerivs(dnvt);
        
        dvt2 = +inf_vt_vt2 -clear_vt2_vt; 
        VT2.updateDerivs(dvt2);

        db = +(inf_vt_b+inf_nvt_b) -(clear_b_vt+clear_b_nvt); 
        B.updateDerivs(db);

        dnvt2 = +inf_nvt_nvt2 -clear_nvt2_nvt; 
        NVT2.updateDerivs(dnvt2);
    }

    arma::rowvec& getIncidence(){
        incidence = arma::rowvec(n_agrp*n_comps_incidence, arma::fill::zeros);
        incidence.subvec(n_agrp*0, (n_agrp*1 - 1)) = get_inf_s_vt() +get_inf_vt_vt2() +get_inf_nvt_b();
        incidence.subvec(n_agrp*1, (n_agrp*2 - 1)) = get_inf_s_nvt() +get_inf_nvt_nvt2() +get_inf_vt_b();

        return incidence;
    }

    arma::rowvec& get_inf_s_vt(){
        return inf_s_vt;
    }
    arma::rowvec& get_inf_s_nvt(){ 
        return inf_s_nvt;
    }
    arma::rowvec& get_inf_vt_vt2(){
        return inf_vt_vt2;
    }
    arma::rowvec& get_inf_vt_b(){
        return inf_vt_b;
    }
    arma::rowvec& get_inf_nvt_b(){
        return inf_nvt_b;
    }
    arma::rowvec& get_inf_nvt_nvt2(){
        return inf_vt_vt2;
    }
};

class Population : public BasePopulation {
public:
    arma::mat beta_vt, beta_nvt;
    arma::rowvec crate_vt, crate_nvt, prev_vt, prev_nvt, foi_vt, foi_nvt;
    double comp;
    arma::rowvec adj_acq_vt_on, adj_acq_nvt_on, adj_acq_vt_off, adj_acq_nvt_off;
    double adj_acq_start, adj_acq_stop;

    // Constructor
    Population(int n_agrp, Rcpp::List parms, int c) 
        :   BasePopulation(n_agrp, parms, c),
            comp(parms["comp"]),
            crate_vt(Rcpp::as<arma::rowvec>(parms["clearVT"])),
            crate_nvt(Rcpp::as<arma::rowvec>(parms["clearNVT"]))
    {
            //Rcpp::Rcout << "DEBUG: Create Population" << std::endl;
            Rcpp::List trial_arm = Rcpp::as<Rcpp::List>(parms["trial_arms"])[c];
            Rcpp::List cluster_parameters = Rcpp::as<Rcpp::List>(trial_arm["parameters"]);
            beta_vt = Rcpp::as<arma::mat>(cluster_parameters["betaVT"]);
            beta_nvt = Rcpp::as<arma::mat>(cluster_parameters["betaNVT"]);
            adj_acq_vt_off = arma::rowvec(n_agrp, arma::fill::ones);
            adj_acq_nvt_off = arma::rowvec(n_agrp, arma::fill::ones);
            adj_acq_vt_on = Rcpp::as<arma::rowvec>(cluster_parameters["adjust_acq_vt"]);
            adj_acq_nvt_on = Rcpp::as<arma::rowvec>(cluster_parameters["adjust_acq_nvt"]);
            adj_acq_start = cluster_parameters["adjust_acq_start"];
            adj_acq_stop = cluster_parameters["adjust_acq_stop"];
    }

    ~Population(){  }

    arma::rowvec& getPrevVT(){
        prev_vt = arma::rowvec(n_agrp, arma::fill::zeros);
        for(int t=0; t < n_vstrat; t++){
            prev_vt += vac_strata[t]->VT.getUpdatedValue() + 2.0 * vac_strata[t]->VT2.getUpdatedValue() + vac_strata[t]->B.getUpdatedValue();
        }
        return(prev_vt);
    }

    arma::rowvec& getPrevNVT(){
        prev_nvt = arma::rowvec(n_agrp, arma::fill::zeros);
        for(int t=0; t < n_vstrat; t++){
            prev_nvt += vac_strata[t]->NVT.getUpdatedValue() + 2.0 * vac_strata[t]->NVT2.getUpdatedValue() + vac_strata[t]->B.getUpdatedValue();
        }
        return(prev_nvt);
    }

    void calculateDerivs(std::vector<std::unique_ptr<Population>> & populations, int n_pops, int p, double & time){
        //calculate weighted prevalence across all clusters
        arma::rowvec prev_vt_calc = arma::rowvec(n_agrp, arma::fill::zeros);
        arma::rowvec prev_nvt_calc = arma::rowvec(n_agrp, arma::fill::zeros);
    
        for(int j=0; j<n_pops; j++){
          prev_vt_calc += populations[j]->getPrevVT() * travel(j, p);
          prev_nvt_calc += populations[j]->getPrevNVT() * travel(j, p);
        }
    
        //Prevalence is then used to calculate force of infection
        if(time >= adj_acq_start && time < adj_acq_stop){
            foi_vt = trans( beta_vt * trans(prev_vt_calc % adj_acq_vt_on) );
            foi_nvt = trans( beta_nvt * trans(prev_nvt_calc % adj_acq_nvt_on) );  
        } else {
            foi_vt = trans( beta_vt * trans(prev_vt_calc % adj_acq_vt_off) );
            foi_nvt = trans( beta_nvt * trans(prev_nvt_calc % adj_acq_nvt_off) );  
        }
    
        #ifdef MP_ENABLED
        //#pragma omp parallel for
        #endif
        for(int t=0; t < n_vstrat; t++){
            //We calculate the derivatives for this specific stratum through the calculateDerivs() function.
            vac_strata[t]->calculateDerivs(foi_vt, foi_nvt, comp, crate_vt, crate_nvt);
        }
    }
};