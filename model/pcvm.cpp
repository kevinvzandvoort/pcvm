/*
 * Pneumococcal transmission model for vaccine trials with multiple vaccination strategies
 *  - groups pneumococcal serotypes as vaccine type or non-vaccine type
 *  - transmission compartments within any stratum are Susceptible (S), Vaccine type carrier (VT), Non-Vaccine type
 *    carrier (NVT), and Both Vaccine type and Non-Vaccine type carrier (B)
 *    
 *    Classes:
  *    - TransComp: Transmission Compartments with states and related functions within one single stratum, should be
 *        used within a Cluster object
 *     - Cluster: Represents one individual arm in the trial. Multiple TransComp objects can be nested within one
 *        cluster object to represent different vaccine strata (e.g. those unvaccinated, who received 1 dose, received
 *        two doses, etc.).
 *    
 *    Functions:
 *     - get_deSolve_gparms_Rcpp: Retrieves parms argument passed to deSolve in R
 *     - initmod: Sets up model before the first iteration of the ODE solver
 *     - derivs: Ran by deSolve at each iteration (timestep) of the model. Updates compartment sizes with new states,
 *        and calls functions to calculate ODEs.
 */

//#include <memory>
#include <vector>
#include <R.h>
#include <RcppArmadillo.h>

//tmp for debugging
#include <chrono>
#include <thread>
//Rcpp::Rcout << "DEBUG message" << std::endl;
//std::this_thread::sleep_for(std::chrono::milliseconds(500));

//The RcppArmadillo attribute sets some macros that enable some external libraries and speed up Armadillo
// [[Rcpp::depends(RcppArmadillo)]]

//Function to return the parms list passed to deSolve as a SEXP object
SEXP attribute_hidden get_deSolve_gparms_Rcpp() {
  static SEXP(*fun)() = NULL;
  if (fun == NULL)
    fun = (SEXP(*)()) R_GetCCallable("deSolve","get_deSolve_gparms");
  return fun();
}

//Ensure C linkage so functions keep their names and can be identified by deSolve in the DLL/SO
extern "C" {
  void derivs(int *neq, double *t, double *y, double *ydot, double *yout, int *ip);
  void initmod(void (* odeparms)(int *, double *));
  void cleanUp();
}

//A TransComp object represents the age-stratified transmission compartments for one individual vaccine stratum
// (e.g. those unvaccinated, or vaccinated with a single dose) within an individual cluster or arm in the trial.
// It is nested within a Cluster object, and holds the states, parameters, and functions to calculate the ODEs for any
// stratum (diamond-shaped type SIS model).
class TransComp{
private:
  arma::rowvec Sus, VT, NVT, B, N, VT_B, NVT_B, VT_NVT_B;
  arma::rowvec dSus, dVT, dNVT, dB;
  arma::rowvec inf_s_vt, inf_s_nvt, inf_vt_b, inf_nvt_b;
  arma::rowvec clear_vt_s, clear_nvt_s, clear_b_vt, clear_b_nvt;
  arma::rowvec wane_s, wane_vt, wane_nvt, wane_b;
  arma::rowvec vac_out_s, vac_out_vt, vac_out_nvt, vac_out_b;
  arma::rowvec age_out_s, age_out_vt, age_out_nvt, age_out_b;
  arma::rowvec age_in_s, age_in_vt, age_in_nvt, age_in_b;
  arma::rowvec age_in_s_vac, age_in_vt_vac, age_in_nvt_vac, age_in_b_vac;
  arma::rowvec age_in_s_unvac, age_in_vt_unvac, age_in_nvt_unvac, age_in_b_unvac;
  arma::rowvec mgr_out_s_vac, mgr_out_vt_vac, mgr_out_nvt_vac, mgr_out_b_vac;
  arma::rowvec mgr_out_s_unvac, mgr_out_vt_unvac, mgr_out_nvt_unvac, mgr_out_b_unvac;
  arma::rowvec vac_cov, vac_eff, vac_waning;
  int n_agrp;
  
public:
  //Constructor
  TransComp(int &n_agrp, Rcpp::List &vac_parms)
    : n_agrp(n_agrp), vac_cov(Rcpp::as<arma::rowvec>(vac_parms["coverage"])),
      vac_eff(Rcpp::as<arma::rowvec>(vac_parms["efficacy"])), vac_waning(Rcpp::as<arma::rowvec>(vac_parms["waning"]))
  {
    //Preallocate memory for efficiency
    mgr_out_s_unvac = mgr_out_vt_unvac = mgr_out_nvt_unvac = mgr_out_b_unvac = 
      mgr_out_s_vac = mgr_out_vt_vac = mgr_out_nvt_vac = mgr_out_b_vac = wane_s = wane_vt = wane_nvt = wane_b = vac_out_s = vac_out_vt = vac_out_nvt = vac_out_b = age_in_b =
      age_in_nvt = age_in_vt = age_in_s = age_out_b = age_out_nvt = age_out_vt = age_out_s = inf_nvt_b =
      inf_vt_b = inf_s_vt = inf_s_nvt = clear_vt_s = clear_nvt_s = clear_b_vt = clear_b_nvt = 
      dSus = dVT = dNVT = dB = Sus = VT = NVT = B = VT_B = NVT_B =
      VT_NVT_B = arma::rowvec(n_agrp, arma::fill::zeros);
  }
  
  //Deconstructor
  ~TransComp(){  }
  
  //This function assigns the current state of each compartment from values calculated in deSolve. Values are stored
  // in one long array. *y points to the address of the first value.
  void setState(double *y, int start){
    for(int a=0; a<n_agrp; a++){
      Sus(a) = y[start + n_agrp * 0 + a];
      VT(a) = y[start + n_agrp * 1 + a];
      NVT(a) = y[start + n_agrp * 2 + a];
      B(a) = y[start + n_agrp * 3 + a];
    }
    
    //Also calculate the total number of people in all compartments, by agegroup
    N = Sus + VT + NVT + B;
    //And the total number that carry combinations of any serotype, by agegroup
    VT_B = VT + B;
    NVT_B = NVT + B;
    VT_NVT_B = VT + NVT + B;
    
    //Set to 0 at start of new iteration
    dSus = dVT = dNVT = dB = arma::rowvec(n_agrp, arma::fill::zeros);
  }
  
  //This function updates the population (processes ageing and migration)
  void updateDemographics(arma::rowvec &arate, arma::rowvec &arate_corr, double N_all_ageing, arma::rowvec &mrate, arma::rowvec &vac_in_s, arma::rowvec &vac_in_vt, arma::rowvec &vac_in_nvt, arma::rowvec &vac_in_b){
    // Process ageing and vaccination
    age_out_s = arate % Sus;
    //Correct ageing rates to account for age groups of different sizes (arate_corr)
    age_in_s.subvec(1, (n_agrp-1)) = age_out_s.subvec(0, (n_agrp-2)) % arate_corr.subvec(1, (n_agrp-1));
    //Birthrate is equal to deathrate (those that age in the last agegroup in all vaccine arms in this this cluster).
    //Births only occur in unvaccinated stratum, so value passed to N_all_ageing is 0 for all other strata
    age_in_s(0) = arate(n_agrp-1) * arate_corr(0) * N_all_ageing;
    //Total number that are not vaccinated and remain in this arm
    age_in_s_unvac = age_in_s % (1.0 - vac_cov);
    //Total number that are effectively vaccinated (and move to next vaccine stratum in a cluster)
    age_in_s_vac = age_in_s % (vac_cov);
    
    age_out_vt = arate % VT;
    age_in_vt.subvec(1, (n_agrp-1)) = age_out_vt.subvec(0, (n_agrp-2)) % arate_corr.subvec(1, (n_agrp-1));
    age_in_vt_unvac = age_in_vt % (1.0 - vac_cov);
    age_in_vt_vac = age_in_vt % (vac_cov);
    
    age_out_nvt = arate % NVT;
    age_in_nvt.subvec(1, (n_agrp-1)) = age_out_nvt.subvec(0, (n_agrp-2)) % arate_corr.subvec(1, (n_agrp-1));
    age_in_nvt_unvac = age_in_nvt % (1.0 - vac_cov);
    age_in_nvt_vac = age_in_nvt % (vac_cov);
    
    age_out_b = arate % B;
    age_in_b.subvec(1, (n_agrp-1)) = age_out_b.subvec(0, (n_agrp-2)) % arate_corr.subvec(1, (n_agrp-1));
    age_in_b_unvac = age_in_b % (1.0 - vac_cov);
    age_in_b_vac = age_in_b % (vac_cov);
    
    // Unvaccinated/Vaccinated individuals who migrate out
    mgr_out_s_unvac = (Sus - age_out_s + age_in_s_unvac) % mrate;
    mgr_out_s_vac = age_in_s_vac % mrate;
    mgr_out_vt_unvac = (VT - age_out_vt + age_in_vt_unvac) % mrate;
    mgr_out_vt_vac = age_in_vt_vac % mrate;
    mgr_out_nvt_unvac = (NVT - age_out_nvt + age_in_nvt_unvac) % mrate;
    mgr_out_nvt_vac = age_in_nvt_vac % mrate;
    mgr_out_b_unvac = (B - age_out_b + age_in_b_unvac) % mrate;
    mgr_out_b_vac = age_in_b_vac % mrate; 
    
    // vaccinated individuals who do not migrate out (assigned to next arm in this cluster)
    age_in_s_vac = age_in_s_vac % (1.0 - mrate);
    age_in_vt_vac = age_in_vt_vac % (1.0 - mrate);
    age_in_nvt_vac = age_in_nvt_vac % (1.0 - mrate);
    age_in_b_vac = age_in_b_vac % (1.0 - mrate);
    
    //The change in states with ageing, in-migration, and in-vaccination, but without in-migration  
    //Those who are vaccinated will move to next arm
    
    dSus += Sus % (-arate -mrate +arate % mrate) + (1.0 - mrate) % age_in_s_unvac +vac_in_s;
    dVT += VT % (-arate -mrate +arate % mrate) + (1.0 - mrate) % age_in_vt_unvac +vac_in_vt;
    dNVT += NVT % (-arate -mrate +arate % mrate) + (1.0 - mrate) % age_in_nvt_unvac +vac_in_nvt;
    dB += B % (-arate -mrate +arate % mrate) + (1.0 - mrate) % age_in_b_unvac +vac_in_b;
    
  }
  
  //This function actually calculates the ODEs
  // Nb all variables are rowvectors (one cell per agegroup)
  void calculateDerivs(
      arma::rowvec &foi_vt, arma::rowvec &foi_nvt, double &comp, arma::rowvec &crate_vt, arma::rowvec &crate_nvt
  ){
    
    //Process changes to ageing, migration, and vaccination
    Sus += dSus;
    VT += dVT;
    NVT += dNVT;
    B += dB;
    
    //Infections from S to VT and NVT (note FOI values already account for travel between clusters)
    inf_s_vt =  foi_vt  % (1.0 - vac_eff) % Sus;
    inf_s_nvt = foi_nvt % Sus;
    
    //Superinfections from VT and NVT to B
    inf_vt_b =  comp * foi_nvt % VT;
    inf_nvt_b = comp * foi_vt  % (1.0 - vac_eff) % NVT;
    
    //Clearance of carriage from VT and NVT to S
    clear_vt_s =  crate_vt  % VT;
    clear_nvt_s = crate_nvt % NVT;
    
    //Clearance of superinfection from B to VT and NVT
    clear_b_vt =  crate_nvt % B;
    clear_b_nvt = crate_vt  % B;
    
    //Total number whose vaccine protection wanes (these will move back to the unvaccinated stratum in this cluster)
    wane_s = vac_waning % Sus;
    wane_vt = vac_waning % VT;
    wane_nvt = vac_waning % NVT;
    wane_b = vac_waning % B;
    
    //All rates are then combined to calculate to add to each ODE. Note, rates for those with waning vaccine protection flowing
    // back into the first stratum (unvaccinated) is added to these ODEs by the updateDerivs() function
    dSus += -(inf_s_vt+inf_s_nvt) +(clear_vt_s+clear_nvt_s) -wane_s;
    dVT += +inf_s_vt -inf_vt_b +clear_b_vt -clear_vt_s -wane_vt;
    dNVT += +inf_s_nvt -inf_nvt_b +clear_b_nvt -clear_nvt_s -wane_nvt;
    dB += +(inf_vt_b+inf_nvt_b) -(clear_b_vt+clear_b_nvt) -wane_b;
  }
  
  //This function updates the ODEs calculated by the calculateDerivs() function by adding the rates for those with 
  // waning vaccine protection flowing back into the first stratum (unvaccinated)
  void updateDerivs(arma::rowvec &delta_Sus, arma::rowvec &delta_VT, arma::rowvec &delta_NVT, arma::rowvec &delta_B){
    dSus += delta_Sus;
    dVT += delta_VT;
    dNVT += delta_NVT;
    dB += delta_B;
  }
  
  //Getter functions to return private variables. Note all return values are passed by reference!
  arma::rowvec& get_inf_s_vt(){ return inf_s_vt; }
  arma::rowvec& get_inf_s_nvt(){ return inf_s_nvt; }
  arma::rowvec& get_inf_vt_b(){ return inf_vt_b; }
  arma::rowvec& get_inf_nvt_b(){ return inf_nvt_b; }
  arma::rowvec& getVT_B(){ return VT_B; }
  arma::rowvec& getNVT_B(){ return NVT_B; }
  arma::rowvec& getN(){ return N; }
  arma::rowvec& getdSus(){ return dSus; }
  arma::rowvec& getdVT(){ return dVT; }
  arma::rowvec& getdNVT(){ return dNVT; }
  arma::rowvec& getdB(){ return dB; }
  arma::rowvec& get_vac_out_s(){ return age_in_s_vac; }
  arma::rowvec& get_vac_out_vt(){ return age_in_vt_vac; }
  arma::rowvec& get_vac_out_nvt(){ return age_in_nvt_vac; }
  arma::rowvec& get_vac_out_b(){ return age_in_b_vac; }
  arma::rowvec& get_wane_s(){ return wane_s; }
  arma::rowvec& get_wane_vt(){ return wane_vt; }
  arma::rowvec& get_wane_nvt(){ return wane_nvt; }
  arma::rowvec& get_wane_b(){ return wane_b; }
  arma::rowvec& get_migr_s_unvac(){ return mgr_out_s_unvac; }
  arma::rowvec& get_migr_s_vac(){ return mgr_out_s_vac; }
  arma::rowvec& get_migr_vt_unvac(){ return mgr_out_vt_unvac; }
  arma::rowvec& get_migr_vt_vac(){ return mgr_out_vt_vac; }
  arma::rowvec& get_migr_nvt_unvac(){ return mgr_out_nvt_unvac; }
  arma::rowvec& get_migr_nvt_vac(){ return mgr_out_nvt_vac; }
  arma::rowvec& get_migr_b_unvac(){ return mgr_out_b_unvac; }
  arma::rowvec& get_migr_b_vac(){ return mgr_out_b_vac; }
  
  //Getter function for the total number of carriers (with any serotype) by agegroup
  arma::rowvec& getCarriers(){ return VT_NVT_B; }
};

//A Cluster object represents an individual cluster or arm in the trial. It will hold multiple vaccinated strata
class Cluster{
private:
  arma::mat beta_vt, beta_nvt;
  arma::rowvec crate_vt, crate_nvt;
  arma::rowvec arate, arate_corr;
  int n_agrp;
  double comp;
  std::vector<TransComp> vac_strata;
  arma::rowvec carriers;
  arma::rowvec deqs_all, deqs;
  arma::rowvec incidence_all, incidence;
  arma::rowvec adj_acq_vt_on, adj_acq_nvt_on;
  arma::rowvec adj_acq_vt_off, adj_acq_nvt_off;
  double adj_acq_start, adj_acq_stop;
  arma::rowvec population_size, N, VT_B, NVT_B, prev_vt, prev_nvt, foi_vt, foi_nvt;
  arma::rowvec population_change;
  double N_ageing;
  arma::mat travel, migration;
  arma::rowvec vac_out_s, vac_out_vt, vac_out_nvt, vac_out_b;
  arma::rowvec mgr_out_s, mgr_out_vt, mgr_out_nvt, mgr_out_b;
  arma::rowvec mgr_out_s_cl, mgr_out_vt_cl, mgr_out_nvt_cl, mgr_out_b_cl;
  arma::rowvec wane_s, wane_vt, wane_nvt, wane_b;
  std::vector<arma::rowvec> mrates;
  arma::rowvec mrates_total;
public:
  //Constructor
  Cluster(int n_agrp, Rcpp::List parms, int c)
    : n_agrp(n_agrp), comp(parms["comp"]),
      crate_vt(Rcpp::as<arma::rowvec>(parms["clearVT"])), crate_nvt(Rcpp::as<arma::rowvec>(parms["clearNVT"])),
      arate(Rcpp::as<arma::rowvec>(parms["ageout"])), travel(Rcpp::as<arma::mat>(parms["travel"])),
      migration(Rcpp::as<arma::mat>(parms["migration"]))
  {
    //We take the total number of required vaccinated strata from the trial_arms element in the parms list passed to
    // deSolve. Total number of arms will be the number provided + 1 for the unvaccinated
    Rcpp::List trial_arm = Rcpp::as<Rcpp::List>(parms["trial_arms"])[c];
    Rcpp::List cluster_parameters = Rcpp::as<Rcpp::List>(trial_arm["parameters"]);
    Rcpp::List vstrata = Rcpp::as<Rcpp::List>(trial_arm["arms"]);
    
    //set population size
    population_size = Rcpp::as<arma::rowvec>(cluster_parameters["N"]);
    population_change = arma::rowvec(n_agrp, arma::fill::zeros);
    
    //set contact matrix
    beta_vt = Rcpp::as<arma::mat>(cluster_parameters["betaVT"]);
    beta_nvt = Rcpp::as<arma::mat>(cluster_parameters["betaNVT"]);
    
    //Set cluster acquisition parameters
    adj_acq_vt_off = arma::rowvec(n_agrp, arma::fill::ones);
    adj_acq_nvt_off = arma::rowvec(n_agrp, arma::fill::ones);
    adj_acq_vt_on = Rcpp::as<arma::rowvec>(cluster_parameters["adjust_acq_vt"]);
    adj_acq_nvt_on = Rcpp::as<arma::rowvec>(cluster_parameters["adjust_acq_nvt"]);
    adj_acq_start = cluster_parameters["adjust_acq_start"];
    adj_acq_stop = cluster_parameters["adjust_acq_stop"];
    
    //Need to correct age rates to account for strata of different size
    arate_corr = arma::rowvec(n_agrp, arma::fill::zeros);
    arate_corr(0) = arate(0)/arate(n_agrp-1);
    arate_corr.subvec(1, (n_agrp-1)) = arate.subvec(1, (n_agrp-1))/arate.subvec(0, (n_agrp-2));
    
    //Vaccinated strata are stored in a vector with TransComp objects
    vac_strata.reserve(vstrata.size()+1);
    
    //We add one stratum for unvaccinated people by default in all arms
    Rcpp::List vac_parms_unvac = Rcpp::List::create(
      Rcpp::Named("coverage", Rcpp::wrap(arma::rowvec(n_agrp, arma::fill::zeros))),
      Rcpp::Named("efficacy", Rcpp::wrap(arma::rowvec(n_agrp, arma::fill::zeros))),
      Rcpp::Named("waning", Rcpp::wrap(arma::rowvec(n_agrp, arma::fill::zeros)))
    );
    
    //We assign each stratum the coverage of the subsequent dose, e.g. the coverage set in the first stratum
    // (unvaccinated) is that of the 1st vaccine dose.
    if(vstrata.size() > 0){
      arma::rowvec coverage_next_dose = Rcpp::as<Rcpp::List>(vstrata[0])["coverage"];
      vac_parms_unvac["coverage"] = Rcpp::wrap(coverage_next_dose);
    }
    
    //Adjust mrates so population size will remain fixed
    mrates_total = arma::rowvec(n_agrp, arma::fill::zeros);
    
    vac_strata.emplace_back(TransComp(n_agrp, vac_parms_unvac));
    
    //Add subsequent strata for vaccinated transmission compartments
    for(int t=0; t < vstrata.size(); t++){
      Rcpp::List vac_parms = vstrata[t];
      //assign coverage for next dose
      if(t < (vstrata.size()-1)){
        arma::rowvec coverage_next_dose = Rcpp::as<Rcpp::List>(vstrata[t+1])["coverage"];
        vac_parms["coverage"] = Rcpp::wrap(coverage_next_dose);
      } else {
        vac_parms["coverage"] = Rcpp::wrap(arma::rowvec(n_agrp, arma::fill::zeros));
      }
      vac_strata.emplace_back(TransComp(n_agrp, vac_parms));
    }
    
    //We preallocate memory for some variables for efficiency
    deqs_all = incidence_all = arma::rowvec(4 * n_agrp *vac_strata.size(), arma::fill::zeros);
    deqs = incidence = arma::rowvec(4 * n_agrp, arma::fill::zeros);
    wane_s = wane_vt = wane_nvt = wane_b = prev_vt = prev_nvt = foi_vt = foi_nvt = vac_out_s = vac_out_vt =
      vac_out_nvt = vac_out_b = mgr_out_s = mgr_out_vt = mgr_out_nvt = mgr_out_b = mgr_out_s_cl = mgr_out_vt_cl = mgr_out_nvt_cl = mgr_out_b_cl = arma::rowvec(n_agrp, arma::fill::zeros);
    N_ageing = 0.0;
    
  }
  
  //Deconstructor, empties the vector with the vaccine strata
  ~Cluster(){
    vac_strata.clear();
  }
  
  void setMigrationRates(int n_clus, int c, std::vector<std::unique_ptr<Cluster>> & clusters){
    for(int j=0; j < n_clus; j++){
      arma::rowvec mrate_by_age;
      if(migration(j, c) < 0.0){
        //Make sure total number of migrants is consistent
        mrate_by_age = migration(c, j) * clusters[j]->getPopSize()/population_size;
      } else {
        mrate_by_age = migration(j, c) * arma::rowvec(n_agrp, arma::fill::ones);
      }
      mrates.emplace_back(mrate_by_age);
      mrates_total += mrate_by_age;
    }
  }
  
  //This function assigns the current state of each compartment from values calculated in deSolve. Values are stored
  // in one long array. *y points to the address of the first value.
  int setState(double *y, int start){
    for(int t=0; t < vac_strata.size(); t++){
      vac_strata[t].setState(y, start + n_agrp * 4 * t);
    }
    
    //also make sure incidence is empty
    incidence_all = arma::rowvec(4 * n_agrp *vac_strata.size(), arma::fill::zeros);
    
    return vac_strata.size() * 4 * n_agrp;
  }
  
  //This function returns the total number of carriers by age, aggregated over all vaccine strata in this arm
  arma::rowvec& getCarriers(){
    carriers = arma::rowvec(n_agrp, arma::fill::zeros);
    for(int t=0; t < vac_strata.size(); t++){
      carriers += vac_strata[t].getCarriers();
    }
    return carriers;
  }
  
  void calculatePrev(){
    N = VT_B = NVT_B = arma::rowvec(n_agrp, arma::fill::zeros);
    
    //We need to know VT (including B) and NVT (including B) prevalence in all vaccine strata in this arm
    for(int t=0; t < vac_strata.size(); t++){
      N += vac_strata[t].getN();
      VT_B += vac_strata[t].getVT_B();
      NVT_B += vac_strata[t].getNVT_B();
    }
    //We are already modelling proportions, so shouldn't divide by N
    //prev_vt = (VT_B);
    //prev_nvt = (NVT_B);
    prev_vt = (VT_B);
    prev_nvt = (NVT_B);
  }
  
  arma::rowvec& getPrevVT(){ return prev_vt; }
  arma::rowvec& getPrevNVT(){ return prev_nvt; }
  
  //This function updates the demographics for each arm in the cluster (process ageing, migration, and vaccination)
  void updateDemographics(std::vector<std::unique_ptr<Cluster>> & clusters, int n_clus, int c){
    N = arma::rowvec(n_agrp, arma::fill::zeros);
    for(int t=0; t < vac_strata.size(); t++){
      N += vac_strata[t].getN();
    }
    
    //Loop through all vaccine strata to calculate the derivatives at this timestep
    for(int t=0; t < vac_strata.size(); t++){

      if(t == 0){
        //No vaccinated individuals move into the first stratum (unvaccinated).
        vac_out_s = vac_out_vt = vac_out_nvt = vac_out_b = arma::rowvec(n_agrp, arma::fill::zeros);
        mgr_out_s = mgr_out_vt = mgr_out_nvt = mgr_out_b = arma::rowvec(n_agrp, arma::fill::zeros);
        //Newborns are all born into the first stratum (unvaccinated), and calculated as a proportion of all those who
        // age in all strata (N_ageing).
        N_ageing = N(n_agrp -1);
      } else {
        //Vaccination happens as people age. Effectively vaccinated people move from the previous stratum to the next.
        vac_out_s = vac_strata[t-1].get_vac_out_s();
        vac_out_vt = vac_strata[t-1].get_vac_out_vt();
        vac_out_nvt = vac_strata[t-1].get_vac_out_nvt();
        vac_out_b = vac_strata[t-1].get_vac_out_b();
        //No newborns enter the vaccinated strata.
        N_ageing = 0.0;
        //Migrate who are vaccinated
        mgr_out_s = vac_strata[t-1].get_migr_s_vac();
        mgr_out_vt = vac_strata[t-1].get_migr_vt_vac();
        mgr_out_nvt = vac_strata[t-1].get_migr_nvt_vac();
        mgr_out_b = vac_strata[t-1].get_migr_b_vac();
      }

      //We update the derivatives for this specific stratum through the calculateDerivs() function.
      vac_strata[t].updateDemographics(arate, arate_corr, N_ageing, mrates_total, vac_out_s, vac_out_vt, vac_out_nvt, vac_out_b);
      
      //Also add those migrating but not vaccinated (remain in this arm in other cluster)
      mgr_out_s += vac_strata[t].get_migr_s_unvac();
      mgr_out_vt += vac_strata[t].get_migr_vt_unvac();
      mgr_out_nvt += vac_strata[t].get_migr_nvt_unvac();
      mgr_out_b += vac_strata[t].get_migr_b_unvac();
      
      for(int j=0; j<n_clus; j++){
        
        if(c == j || arma::accu(mrates_total) == 0 || arma::accu(mrates[j]) == 0){ 
          continue;
        }
        
        //Calculate proportion that goes to this cluster
        mgr_out_s_cl = mgr_out_s % (mrates[j] / mrates_total);
        mgr_out_vt_cl = mgr_out_vt % (mrates[j] / mrates_total);
        mgr_out_nvt_cl = mgr_out_nvt % (mrates[j] / mrates_total);
        mgr_out_b_cl = mgr_out_b % (mrates[j] / mrates_total);
        
        //Adjust as population size may differ
        mgr_out_s_cl = mgr_out_s_cl / mrates[j];
        mgr_out_vt_cl = mgr_out_vt_cl / mrates[j];
        mgr_out_nvt_cl = mgr_out_nvt_cl / mrates[j];
        mgr_out_b_cl = mgr_out_b_cl / mrates[j];
        
        clusters[j]->addMigrants(c, t, mgr_out_s_cl, mgr_out_vt_cl, mgr_out_nvt_cl, mgr_out_b_cl);
      }
      
    }
  }
  
  void addMigrants(int origin, int t, arma::rowvec &delta_Sus, arma::rowvec &delta_VT, arma::rowvec &delta_NVT, arma::rowvec &delta_B){
    
    //rates are multiplied by mrate for this cluster, to keep population sizes constant
    arma::rowvec delta_Sus_adj = delta_Sus % mrates[origin];
    arma::rowvec delta_VT_adj = delta_VT % mrates[origin];
    arma::rowvec delta_NVT_adj = delta_NVT % mrates[origin];
    arma::rowvec delta_B_adj = delta_B % mrates[origin];
    
    updateDerivsArm(t, delta_Sus_adj, delta_VT_adj, delta_NVT_adj, delta_B_adj);
  }
  
  void updateDerivsArm(int t, arma::rowvec &delta_Sus, arma::rowvec &delta_VT, arma::rowvec &delta_NVT, arma::rowvec &delta_B){
    vac_strata[t].updateDerivs(delta_Sus, delta_VT, delta_NVT, delta_B);
  }
  
  //This function calculates the derivatives in all vaccine strata within this arm
  void calculateDerivs(std::vector<std::unique_ptr<Cluster>> & clusters, int n_clus, int c, double & time){
    //calculate weighted prevalence across all clusters
    arma::rowvec prev_vt_calc = arma::rowvec(n_agrp, arma::fill::zeros);
    arma::rowvec prev_nvt_calc = arma::rowvec(n_agrp, arma::fill::zeros);
    
    for(int j=0; j<n_clus; j++){
      prev_vt_calc += clusters[j]->getPrevVT() * travel(j, c);
      prev_nvt_calc += clusters[j]->getPrevNVT() * travel(j, c);
    }
    
    //Prevalence is then used to calculate force of infection
    if(time >= adj_acq_start && time < adj_acq_stop){
      foi_vt = trans( beta_vt * trans(prev_vt_calc % adj_acq_vt_on) );
      foi_nvt = trans( beta_nvt * trans(prev_nvt_calc % adj_acq_nvt_on) );  
    } else {
      foi_vt = trans( beta_vt * trans(prev_vt_calc % adj_acq_vt_off) );
      foi_nvt = trans( beta_nvt * trans(prev_nvt_calc % adj_acq_nvt_off) );  
    }
    //foi_vt = trans( beta_vt * trans(prev_vt_calc % (time >= adj_acq_start && time < adj_acq_stop ? adj_acq_vt_on : adj_acq_vt_off)) );
    //foi_nvt = trans( beta_nvt * trans(prev_nvt_calc % (time >= adj_acq_start && time < adj_acq_stop ? adj_acq_nvt_on : adj_acq_nvt_off)) );
    
    //Loop through all vaccine strata to calculate the derivatives at this timestep
    wane_s = wane_vt = wane_nvt = wane_b = arma::rowvec(n_agrp, arma::fill::zeros);
    for(int t=0; t < vac_strata.size(); t++){
      //We calculate the derivatives for this specific stratum through the calculateDerivs() function.
      vac_strata[t].calculateDerivs(foi_vt, foi_nvt, comp, crate_vt, crate_nvt);
      //We add all those whose vaccine efficacy wanes across all vaccinated strata
      wane_s += vac_strata[t].get_wane_s();
      wane_vt += vac_strata[t].get_wane_vt();
      wane_nvt += vac_strata[t].get_wane_nvt();
      wane_b += vac_strata[t].get_wane_b();
    }
    //All individuals whose vaccine efficacy has waned are added to the first stratum (unvaccinated), through the
    // updateDerivs() function
    vac_strata[0].updateDerivs(wane_s, wane_vt, wane_nvt, wane_b);
  }
  
  //This function gets the derivatives from all strata, and returns them as one long vector.
  arma::rowvec& getDerivs(){
    for(int t=0; t < vac_strata.size(); t++){
      deqs.subvec(0 + n_agrp*0, (n_agrp+n_agrp*0 - 1)) = vac_strata[t].getdSus();
      deqs.subvec(0 + n_agrp*1, (n_agrp+n_agrp*1 - 1)) = vac_strata[t].getdVT();
      deqs.subvec(0 + n_agrp*2, (n_agrp+n_agrp*2 - 1)) = vac_strata[t].getdNVT();
      deqs.subvec(0 + n_agrp*3, (n_agrp+n_agrp*3 - 1)) = vac_strata[t].getdB();
      deqs_all.subvec(0 + 4*n_agrp*t, (4*n_agrp + 4*n_agrp*t - 1)) = deqs;
    }
    return deqs_all;
  }
  
  //This function gets the incidence from all strata, and returns them as one long vector.
  arma::rowvec& getIncidence(){
    for(int t=0; t < vac_strata.size(); t++){
      incidence.subvec(0 + n_agrp*0, (n_agrp+n_agrp*0 - 1)) = vac_strata[t].get_inf_s_vt();
      incidence.subvec(0 + n_agrp*1, (n_agrp+n_agrp*1 - 1)) = vac_strata[t].get_inf_s_nvt();
      incidence.subvec(0 + n_agrp*2, (n_agrp+n_agrp*2 - 1)) = vac_strata[t].get_inf_vt_b();
      incidence.subvec(0 + n_agrp*3, (n_agrp+n_agrp*3 - 1)) = vac_strata[t].get_inf_nvt_b();
      //Calculate cumulative incidence
      incidence_all.subvec(0 + 4*n_agrp*t, (4*n_agrp + 4*n_agrp*t - 1)) = incidence;//+= incidence;
    }
    return incidence_all;
  }
  
  //Return the current population size for this cluster
  arma::rowvec& getPopSize(){ return population_size; }
};

//Keep global parameters in memory so they can be accessed in initMod() and derivs()
//std::vector<Cluster*> global_clusters;
std::vector<std::unique_ptr<Cluster>> clusters;
int n_clus, n_agrp;

//This function sets the model up, and stores the parameter values in memory. It is only called once when setting up
// the model
void initmod(void (* odeparms)(int *, double *)) {
  //Rcpp::Rcout << "DEBUG: initializing model" << std::endl;
  //We get the parms argument passed to deSolve as SEXP object
  SEXP sparms = get_deSolve_gparms_Rcpp();
  
  try {
    //Parse parameters passed to deSolve as Rcpp::List
    Rcpp::List parms = Rcpp::clone(Rcpp::as<Rcpp::List>(sparms));
    
    //Define the number of trial arms/clusters and agegroups from parameter list passed to deSolve
    n_clus = Rcpp::as<Rcpp::List>(parms["trial_arms"]).size();
    if(n_clus == 0) n_clus = 1;
    n_agrp = parms["n_agrp"];
    
    //We can't do garbage collection at end of model run in deSolve, so we do it if the same DLL/SO is still loaded
    // and deSolve is ran again
    //for(int c = 0; c < global_clusters.size(); c++){
    //  delete global_clusters[c];
    //}
    //global_clusters.clear();
    clusters.clear();
    
    //We create a new Cluster object for every arm in the trial, and store a pointer to it in the clusters vector
    for(int c = 0; c < n_clus; c++){
      //Cluster* cluster = new Cluster(n_agrp, parms, c);
      //global_clusters.emplace_back(cluster);
      clusters.emplace_back(std::make_unique<Cluster>(n_agrp, parms, c));
    }
    
    //Now all cluster objects exist, process the migration rates so these are consistent with the population size
    for(int c = 0; c < n_clus; c++){
      //Cluster* cluster = new Cluster(n_agrp, parms, c);
      //global_clusters.emplace_back(cluster);
      clusters[c]->setMigrationRates(n_clus, c, clusters);
    }
  } catch(std::exception& __ex__){
    forward_exception_to_r(__ex__);
  } catch(...){
    ::Rf_error( "c++ exception (unknown reason)" );
  }
}

void cleanUp() {
  for(int c = 0; c < clusters.size(); c++){
    //delete clusters[c];
    clusters[c].reset();
  }
  clusters.clear();
  
}

//This function is called by deSolve in every iteration of the integrator
void derivs(int *neq, double *t, double *y, double *ydot, double *yout, int *ip) {
  /* Something strange going on here, I need to copy the global_clusters variable to a local variable for this to work
   *  otherwise the model 1) runs very slow and 2) returns the wrong results.
   * Another thing that makes this work is by printing something (i.e. Rcpp::Rcout or std::cout) at every step, such
   *  as the time-step.
   */
  double time = t[0];
  
  //TODO: check if this is necessary
  if (ip[0] < 1) error("nout should be at least 1");
  //if(time < 1){
  //  Rcpp::Rcout << "DEBUG: time: " << time << "; ip[0]: " << ip[0] << std::endl;
  //}
  
  //We loop through every cluster in the model, and update the state in each compartment. State is passed from deSolve
  // by the y array. The states are ordered as: cluster > vaccine_arm > compartment (S, VT, NVT, B) > agegroup.
  //The setState() function returns the total number of compartments updated by that cluster (which is
  // the number of vaccine strata in that cluster * 4 infection compartments * the number of agegroups)
  int start = 0;
  for(int c = 0; c < n_clus; c++){
    start += clusters[c]->setState(y, start);
  }
  
  //Process demographic changes (ageing and migration) and vaccinations
  for(int c = 0; c < n_clus; c++){
    clusters[c]->updateDemographics(clusters, n_clus, c); 
    clusters[c]->calculatePrev(); 
  }
  
  //We calculate the total number of carriers with any serotype, aggregated over all agegroups and vaccine arms
  double n_carr = 0;
  for(int c = 0; c < n_clus; c++){
    n_carr += arma::accu(clusters[c]->getCarriers());
  }
  
  //Only continue with the model if the total number of carriers is larger than 2
  if(n_carr > 0){
    //We calculate the ODEs within all clusters, which is done in the calculateDerivs() function
    for(int c = 0; c < n_clus; c++){
      clusters[c]->calculateDerivs(clusters, n_clus, c, time);
    }
    
    //We return the model output to deSolve by copying the model output in the ydot array. The states are ordered
    // as: cluster > vaccine_arm > compartment (S, VT, NVT, B) > agegroup.
    int i = 0;
    for(int c = 0; c < n_clus; c++){
      arma::rowvec deqs = clusters[c]->getDerivs();
      for(int d = 0; d < deqs.size(); d++){
        ydot[i] = deqs(d);
        i += 1;
      }
    }
  } else {
    //If there are only two or less carriers, the values of the ODEs are set to 0 in all compartments
    for(int i=0; i<*neq; i++){
      ydot[i] = 0.0;
    }
  }
  
  //Additional return value to deSolve, not used. TODO: check if this is necessary
  if (ip[0] > 1){
    int i = 0;
    for(int c = 0; c < n_clus; c++){
      arma::rowvec incidence = clusters[c]->getIncidence();
      for(int d = 0; d < incidence.size(); d++){
        yout[i] = incidence(d);
        i += 1;
      }
    }
    //Population size never changes in current model, this is not needed (for now)
    /*for(int c = 0; c < n_clus; c++){
      arma::rowvec total_people = clusters[c]->getPopSize();
      for(int d = 0; d < total_people.size(); d++){
        yout[i] = total_people(d);
        i += 1;
      }
    }*/
    //Rcpp::Rcout << "Finished" << std::endl;
  }
//yout[0] = 1.0;
  
  //for(int c = 0; c < n_clus; c++){
  //  delete clusters[c];
  //}
  //clusters.clear();
}