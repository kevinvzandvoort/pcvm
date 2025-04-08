/*
 Model implementation for a cube-shaped Streptococcus Pneumoniae model
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

// Internal project includes
// - Note, there are are additional includes at the end of the file
#include "./MetaVax/Compartment.h"
#include "./MetaVax/BaseVaccinationGroup.h"
#include "./MetaVax/BasePopulation.h"

// Set number of compartments
int n_comps_prevalence = 8;
int n_comps_incidence = 3;

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
  arma::rowvec vac_eff_ot;
  Compartment Sus, VT, OT, NVT, VT_OT, VT_NVT, OT_NVT, VT_OT_NVT;
  arma::rowvec sus, vt, ot, nvt, vt_ot, vt_nvt, ot_nvt, vt_ot_nvt,
  dsus, dvt, dot, dnvt, dvt_ot, dvt_nvt, dot_nvt, dvt_ot_nvt,
  inf_s2vt, inf_s2ot, inf_s2nvt, inf_vt2vt_ot, inf_vt2vt_nvt, inf_ot2vt_ot, inf_ot2ot_nvt,
  inf_nvt2vt_nvt, inf_nvt2ot_nvt, inf_vt_ot2vt_ot_nvt, inf_vt_nvt2vt_ot_nvt, inf_ot_nvt2vt_ot_nvt,
  clear_vt2s, clear_ot2s, clear_nvt2s, clear_vt_ot2vt, clear_vt_nvt2vt, clear_vt_ot2ot, clear_ot_nvt2ot,
  clear_vt_nvt2nvt, clear_ot_nvt2nvt, clear_vt_ot_nvt2vt_ot, clear_vt_ot_nvt2vt_nvt, clear_vt_ot_nvt2ot_nvt,
  incidence;
  // Constructor
  VaccinationGroup(int n_agrp, Rcpp::List vac_parms, arma::rowvec arate, arma::rowvec arate_corr) 
    :   BaseVaccinationGroup(n_agrp, vac_parms), vac_eff_ot(Rcpp::as<arma::rowvec>(vac_parms["efficacy_transmission_ot"])),
        //Note that last value in Compartment constructor is proportion of newborns
        // born into this compartment
        Sus(n_agrp, arate, arate_corr, 1.0),
        VT(n_agrp, arate, arate_corr, 0.0),
        OT(n_agrp, arate, arate_corr, 0.0),
        NVT(n_agrp, arate, arate_corr, 0.0),
        VT_OT(n_agrp, arate, arate_corr, 0.0),
        VT_NVT(n_agrp, arate, arate_corr, 0.0),
        OT_NVT(n_agrp, arate, arate_corr, 0.0),
        VT_OT_NVT(n_agrp, arate, arate_corr, 0.0)
  {
    //Rcpp::Rcout << "DEBUG: Create VaccinationGroup" << std::endl;
    compartments = {&Sus, &VT, &OT, &NVT, &VT_OT, &VT_NVT, &OT_NVT, &VT_OT_NVT}; // Add pointers to compartment objects
    sus = vt = ot = nvt = vt_ot = vt_nvt = ot_nvt = dsus = dvt = dot = dnvt = dvt_ot= dvt_nvt = dot_nvt = dvt_ot_nvt =
      inf_s2vt = inf_s2ot = inf_s2nvt = inf_vt2vt_ot = inf_vt2vt_nvt = inf_ot2vt_ot = inf_ot2ot_nvt =
      inf_nvt2vt_nvt = inf_nvt2ot_nvt = inf_vt_ot2vt_ot_nvt = inf_vt_nvt2vt_ot_nvt = inf_ot_nvt2vt_ot_nvt =
      clear_vt2s = clear_ot2s = clear_nvt2s = clear_vt_ot2vt = clear_vt_nvt2vt = clear_vt_ot2ot = clear_ot_nvt2ot =
      clear_vt_nvt2nvt = clear_ot_nvt2nvt = clear_vt_ot_nvt2vt_ot = clear_vt_ot_nvt2vt_nvt = clear_vt_ot_nvt2ot_nvt = arma::rowvec(n_agrp, arma::fill::zeros);
    incidence = arma::rowvec(n_agrp*n_comps_incidence, arma::fill::zeros);
  }
  
  //Deconstructor
  ~VaccinationGroup(){  }
  
  void calculateDerivs(arma::rowvec &foi_vt, arma::rowvec &foi_ot, arma::rowvec &foi_nvt, double &comp, arma::rowvec &crate_vt, arma::rowvec &crate_ot, arma::rowvec &crate_nvt) {
    sus = Sus.getUpdatedValue();
    vt = VT.getUpdatedValue();
    ot = OT.getUpdatedValue();
    nvt = NVT.getUpdatedValue();
    vt_ot = VT_OT.getUpdatedValue();
    vt_nvt = VT_NVT.getUpdatedValue();
    ot_nvt = OT_NVT.getUpdatedValue();
    vt_ot_nvt = VT_OT_NVT.getUpdatedValue();
    
    //Infections from S to VT, OT and NVT (note FOI values already account for travel between populations)
    inf_s2vt =  foi_vt  % (1.0 - vac_eff) % sus;
    inf_s2ot =  foi_ot  % (1.0 - vac_eff_ot) % sus;
    inf_s2nvt = foi_nvt % sus;
    
    //Superinfections
    inf_vt2vt_ot =  comp * foi_ot % (1.0 - vac_eff_ot) % vt;
    inf_vt2vt_nvt = comp * foi_nvt % vt;
    inf_ot2vt_ot =  comp * foi_vt % (1.0 - vac_eff) % ot;
    inf_ot2ot_nvt = comp * foi_nvt % ot;
    inf_nvt2vt_nvt = comp * foi_vt % (1.0 - vac_eff) % nvt;
    inf_nvt2ot_nvt = comp * foi_ot % (1.0 - vac_eff_ot) % nvt;
    
    //assume comp scales by number of serotypes
    inf_vt_ot2vt_ot_nvt = comp * comp * foi_nvt % vt_ot;
    inf_vt_nvt2vt_ot_nvt = comp * comp * foi_ot % (1.0 - vac_eff_ot) % vt_nvt;
    inf_ot_nvt2vt_ot_nvt = comp * comp * foi_vt % (1.0 - vac_eff) % ot_nvt;
    
    //Clearance of carriage from VT, OT, and NVT to S
    clear_vt2s =  crate_vt  % vt;
    clear_ot2s =  crate_ot  % ot;
    clear_nvt2s = crate_nvt % nvt;
    
    //Clearance of superinfections to VT, OT and NVT
    clear_vt_ot2vt = crate_ot % vt_ot;
    clear_vt_ot2ot = crate_vt % vt_ot;
    clear_vt_nvt2vt = crate_nvt % vt_nvt;
    clear_vt_nvt2nvt = crate_vt % vt_nvt;
    clear_ot_nvt2ot = crate_nvt % ot_nvt;
    clear_ot_nvt2nvt = crate_ot % ot_nvt;
    
    clear_vt_ot_nvt2vt_ot = crate_nvt % vt_ot_nvt;
    clear_vt_ot_nvt2vt_nvt = crate_ot % vt_ot_nvt;
    clear_vt_ot_nvt2ot_nvt = crate_vt % vt_ot_nvt;
    
    //All rates are then combined to update each ODE for each compartment
    dsus = -(inf_s2vt+inf_s2ot+inf_s2nvt) +(clear_vt2s+clear_ot2s+clear_nvt2s);
    Sus.updateDerivs(dsus);
    
    dvt = +inf_s2vt -(inf_vt2vt_ot+inf_vt2vt_nvt) +(clear_vt_ot2vt+clear_vt_nvt2vt) -clear_vt2s;
    VT.updateDerivs(dvt);
    
    dot = +inf_s2ot -(inf_ot2vt_ot+inf_ot2ot_nvt) +(clear_vt_ot2ot+clear_ot_nvt2ot) -clear_ot2s;
    OT.updateDerivs(dot);
    
    dnvt = +inf_s2nvt -(inf_nvt2vt_nvt+inf_nvt2ot_nvt) +(clear_vt_nvt2nvt+clear_ot_nvt2nvt) -clear_nvt2s;
    NVT.updateDerivs(dnvt);
    
    dvt_ot = +(inf_vt2vt_ot+inf_ot2vt_ot) -inf_vt_ot2vt_ot_nvt +clear_vt_ot_nvt2vt_ot -(clear_vt_ot2vt+clear_vt_ot2ot); 
    VT_OT.updateDerivs(dvt_ot);
    
    dvt_nvt = +(inf_vt2vt_nvt+inf_nvt2vt_nvt) -inf_vt_nvt2vt_ot_nvt +clear_vt_ot_nvt2vt_nvt -(clear_vt_nvt2vt+clear_vt_nvt2nvt); 
    VT_NVT.updateDerivs(dvt_nvt);
    
    dot_nvt = +(inf_ot2ot_nvt+inf_nvt2ot_nvt) -inf_ot_nvt2vt_ot_nvt +clear_vt_ot_nvt2ot_nvt -(clear_ot_nvt2ot+clear_ot_nvt2nvt); 
    OT_NVT.updateDerivs(dot_nvt);
    
    dvt_ot_nvt = +(inf_vt_ot2vt_ot_nvt+inf_vt_nvt2vt_ot_nvt+inf_ot_nvt2vt_ot_nvt) -(clear_vt_ot_nvt2vt_ot+clear_vt_ot_nvt2vt_nvt+clear_vt_ot_nvt2ot_nvt); 
    VT_OT_NVT.updateDerivs(dvt_ot_nvt);
  }
  
  arma::rowvec& getIncidence(){
    incidence = arma::rowvec(n_agrp*n_comps_incidence, arma::fill::zeros);
    incidence.subvec(n_agrp*0, (n_agrp*1 - 1)) = get_inf_s2vt() + get_inf_ot2vt_ot() + get_inf_nvt2vt_nvt() + get_inf_ot_nvt2vt_ot_nvt();
    incidence.subvec(n_agrp*1, (n_agrp*2 - 1)) = get_inf_s2ot() + get_inf_vt2vt_ot() + get_inf_nvt2ot_nvt() + get_inf_vt_nvt2vt_ot_nvt();
    incidence.subvec(n_agrp*2, (n_agrp*3 - 1)) = get_inf_s2nvt() + get_inf_vt2vt_nvt() + get_inf_ot2ot_nvt() + get_inf_vt_nvt2vt_ot_nvt();
    
    return incidence;
  }
  
  arma::rowvec& get_inf_s2vt(){
    return inf_s2vt;
  }
  arma::rowvec& get_inf_s2ot(){
    return inf_s2ot;
  }
  arma::rowvec& get_inf_s2nvt(){ 
    return inf_s2nvt;
  }
  arma::rowvec& get_inf_vt2vt_ot(){ 
    return inf_vt2vt_ot;
  }
  arma::rowvec& get_inf_vt2vt_nvt(){ 
    return inf_vt2vt_nvt;
  }
  arma::rowvec& get_inf_ot2vt_ot(){ 
    return inf_ot2vt_ot;
  }
  arma::rowvec& get_inf_ot2ot_nvt(){ 
    return inf_ot2ot_nvt;
  }
  arma::rowvec& get_inf_nvt2vt_nvt(){ 
    return inf_nvt2vt_nvt;
  }
  arma::rowvec& get_inf_nvt2ot_nvt(){ 
    return inf_nvt2ot_nvt;
  }
  arma::rowvec& get_inf_vt_ot2vt_ot_nvt(){ 
    return inf_vt_ot2vt_ot_nvt;
  }
  arma::rowvec& get_inf_vt_nvt2vt_ot_nvt(){ 
    return inf_vt_nvt2vt_ot_nvt;
  }
  arma::rowvec& get_inf_ot_nvt2vt_ot_nvt(){ 
    return inf_ot_nvt2vt_ot_nvt;
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
  arma::mat beta_vt, beta_ot, beta_nvt;
  arma::rowvec crate_vt, crate_ot, crate_nvt, prev_vt, prev_ot, prev_nvt, foi_vt, foi_ot, foi_nvt;
  arma::rowvec adj_acq_on, adj_acq_off;
  double adj_acq_start, adj_acq_stop;
  
  // Constructor
  Population(int n_agrp, Rcpp::List parms, int c) 
    :   BasePopulation(n_agrp, parms, c),
        comp(Rcpp::as<Rcpp::List>(parms["global_settings"])["comp"]),
        crate_vt(Rcpp::as<arma::rowvec>(Rcpp::as<Rcpp::List>(parms["global_settings"])["clearVT"])),
        crate_ot(Rcpp::as<arma::rowvec>(Rcpp::as<Rcpp::List>(parms["global_settings"])["clearOT"])),
        crate_nvt(Rcpp::as<arma::rowvec>(Rcpp::as<Rcpp::List>(parms["global_settings"])["clearNVT"]))
  {
    //Rcpp::Rcout << "DEBUG: Create Population" << std::endl;
    Rcpp::List trial_arm = Rcpp::as<Rcpp::List>(parms["populations"])[c];
    Rcpp::List cluster_parameters = Rcpp::as<Rcpp::List>(trial_arm["parameters"]);
    beta_vt = Rcpp::as<arma::mat>(cluster_parameters["betaVT"]);
    beta_ot = Rcpp::as<arma::mat>(cluster_parameters["betaOT"]);
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
      prev_vt += vac_strata[t]->VT.getUpdatedValue() + vac_strata[t]->VT_OT.getUpdatedValue() + vac_strata[t]->VT_NVT.getUpdatedValue() + vac_strata[t]->VT_OT_NVT.getUpdatedValue();
    }
    return(prev_vt);
  }
  
  arma::rowvec& getPrevOT(){
    prev_ot = arma::rowvec(n_agrp, arma::fill::zeros);
    for(int t=0; t < n_vstrat; t++){
      prev_ot += vac_strata[t]->OT.getUpdatedValue() + vac_strata[t]->VT_OT.getUpdatedValue() + vac_strata[t]->OT_NVT.getUpdatedValue() + vac_strata[t]->VT_OT_NVT.getUpdatedValue();
    }
    return(prev_ot);
  }
  
  arma::rowvec& getPrevNVT(){
    prev_nvt = arma::rowvec(n_agrp, arma::fill::zeros);
    for(int t=0; t < n_vstrat; t++){
      prev_nvt += vac_strata[t]->NVT.getUpdatedValue() + vac_strata[t]->VT_NVT.getUpdatedValue() + vac_strata[t]->OT_NVT.getUpdatedValue() + vac_strata[t]->VT_OT_NVT.getUpdatedValue();
    }
    return(prev_nvt);
  }
  
  void calculateDerivs(std::vector<std::unique_ptr<Population>> & populations, int n_pops, int p, double & time){
    //calculate weighted prevalence across all clusters
    arma::rowvec prev_vt_calc = arma::rowvec(n_agrp, arma::fill::zeros);
    arma::rowvec prev_ot_calc = arma::rowvec(n_agrp, arma::fill::zeros);
    arma::rowvec prev_nvt_calc = arma::rowvec(n_agrp, arma::fill::zeros);
    
    for(int j=0; j<n_pops; j++){
      prev_vt_calc += populations[j]->getPrevVT() * travel(j, p);
      prev_ot_calc += populations[j]->getPrevOT() * travel(j, p);
      prev_nvt_calc += populations[j]->getPrevNVT() * travel(j, p);
    }
    
    //Prevalence is then used to calculate force of infection
    if(time >= adj_acq_start && time < adj_acq_stop){
      foi_vt = trans( beta_vt * trans(prev_vt_calc % adj_acq_on) );
      foi_ot = trans( beta_ot * trans(prev_ot_calc % adj_acq_on) );
      foi_nvt = trans( beta_nvt * trans(prev_nvt_calc % adj_acq_on) );  
    } else {
      foi_vt = trans( beta_vt * trans(prev_vt_calc % adj_acq_off) );
      foi_ot = trans( beta_ot * trans(prev_ot_calc % adj_acq_off) );
      foi_nvt = trans( beta_nvt * trans(prev_nvt_calc % adj_acq_off) );  
    }
    
#ifdef MP_ENABLED
    //#pragma omp parallel for
#endif
    for(int t=0; t < n_vstrat; t++){
      //We calculate the derivatives for this specific stratum through the calculateDerivs() function.
      vac_strata[t]->calculateDerivs(foi_vt, foi_ot, foi_nvt, comp, crate_vt, crate_ot, crate_nvt);
    }
  }
};

// Additional internal project includes
// - Note, these need to be last. Other includes are at the top of the file
#include "./MetaVax/BasePopulation.cpp"
#include "./MetaVax/metapopulation.cpp"