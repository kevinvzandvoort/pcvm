/*
    Model implementation for a diamond-shaped Streptococcus Pneumoniae model
    Author: Kevin van Zandvoort, London School of Hygiene and Tropical Medicine
*/

// Standard library includes
#include <vector>
#include <memory>
#ifdef MP_ENABLED
    #include <omp.h>
#endif

// External library includes
#include <R.h>
#include <RcppArmadillo.h>
//The RcppArmadillo attribute sets some macros that enable some external libraries and speed up Armadillo
// [[Rcpp::depends(RcppArmadillo)]]

#include <thread>
#include <chrono>

// Internal project includes
// - Note, there are are additional includes at the end of the file
#include "./MetaVax/Compartment.h"
#include "./MetaVax/BaseVaccinationGroup.h"
#include "./MetaVax/BasePopulation.h"

// Set number of compartments
int n_comps_prevalence = 4;
int n_comps_incidence = 2;

/* Define the VaccinationGroup class
- the following functions need to exist:
  constructor that sets the different epidemiological compartments
  - should correspond to the number provided earlier (n_comps_prevalence)
  - should pass n_agrp and vac_parms to the BaseVaccinationGroup constructor
  calculateDerivs with the calculations for the transitions between epidemiological compartments
  getIncidence to return any incidence values
  - should correspond to number provided earlier (n_comps_incidence)
*/ 
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
        incidence.subvec(n_agrp*0, (n_agrp*1 - 1)) = get_inf_s_vt() + get_inf_nvt_b();
        incidence.subvec(n_agrp*1, (n_agrp*2 - 1)) = get_inf_s_nvt() + get_inf_vt_b();

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

/* Define the Population class
- the following functions need to exist:
  constructor that sets the different epidemiological parameters
  - should pass n_agrp, parms, and c (index of population)to the BasePopulation constructor
  calculateDerivs that calculates the FOI and calls the calculateDerivs function in each VaccinationGroup
*/ 
class Population : public BasePopulation {
  public:
    double comp;
    arma::mat beta_vt, beta_nvt;
    arma::rowvec crate_vt, crate_nvt, prev_vt, prev_nvt, foi_vt, foi_nvt;
    arma::rowvec adj_acq_on, adj_acq_off;
    double adj_acq_start, adj_acq_stop;

    // Constructor
    Population(int n_agrp, Rcpp::List parms, int c) 
        :   BasePopulation(n_agrp, parms, c),
            comp(Rcpp::as<Rcpp::List>(parms["global_settings"])["comp"]),
            crate_vt(Rcpp::as<arma::rowvec>(Rcpp::as<Rcpp::List>(parms["global_settings"])["clearVT"])),
            crate_nvt(Rcpp::as<arma::rowvec>(Rcpp::as<Rcpp::List>(parms["global_settings"])["clearNVT"]))
    {
            //Rcpp::Rcout << "DEBUG: Create Population" << std::endl;
            Rcpp::List trial_arm = Rcpp::as<Rcpp::List>(parms["populations"])[c];
            Rcpp::List cluster_parameters = Rcpp::as<Rcpp::List>(trial_arm["parameters"]);
            beta_vt = Rcpp::as<arma::mat>(cluster_parameters["betaVT"]);
            beta_nvt = Rcpp::as<arma::mat>(cluster_parameters["betaNVT"]);
            adj_acq_off = arma::rowvec(n_agrp, arma::fill::ones);
            adj_acq_on = Rcpp::as<arma::rowvec>(cluster_parameters["adjust_acq"]);
            adj_acq_start = cluster_parameters["adjust_acq_start"];
            adj_acq_stop = cluster_parameters["adjust_acq_stop"];
    }

    ~Population(){  }

    arma::rowvec& getPrevVT(){
        prev_vt = arma::rowvec(n_agrp, arma::fill::zeros);
        for(int t=0; t < n_vstrat; t++){
            prev_vt += vac_strata[t]->VT.getUpdatedValue() + vac_strata[t]->B.getUpdatedValue();
        }
        return(prev_vt);
    }

    arma::rowvec& getPrevNVT(){
        prev_nvt = arma::rowvec(n_agrp, arma::fill::zeros);
        for(int t=0; t < n_vstrat; t++){
            prev_nvt += vac_strata[t]->NVT.getUpdatedValue() + vac_strata[t]->B.getUpdatedValue();
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
            foi_vt = trans( beta_vt * trans(prev_vt_calc % adj_acq_on) );
            foi_nvt = trans( beta_nvt * trans(prev_nvt_calc % adj_acq_on) );  
        } else {
            foi_vt = trans( beta_vt * trans(prev_vt_calc % adj_acq_off) );
            foi_nvt = trans( beta_nvt * trans(prev_nvt_calc % adj_acq_off) );  
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

// Additional internal project includes
// - Note, these need to be last. Other includes are at the top of the file
#include "./MetaVax/BasePopulation.cpp"
#include "./MetaVax/metapopulation.cpp"