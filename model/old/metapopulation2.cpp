//#include <memory>
#include <vector>
#include <R.h>
#include <RcppArmadillo.h>

//tmp for debugging
#include <chrono>
#include <thread>
#include <iostream>

#include "Compartment.h"
#include "BaseVaccinationGroup.h"
#include "BasePopulation.h"

#include "epimodel_streppneumo_neutral.h"

#include "BasePopulation.cpp"

#ifdef MP_ENABLED
#include <omp.h>
#endif

//The RcppArmadillo attribute sets some macros that enable some external libraries and speed up Armadillo
// [[Rcpp::depends(RcppArmadillo)]]

//Function to return the parms list passed to deSolve as a SEXP object
SEXP attribute_hidden get_deSolve_gparms_Rcpp() {
  static SEXP(*fun)() = NULL;
  if (fun == NULL)
    fun = (SEXP(*)()) R_GetCCallable("deSolve","get_deSolve_gparms");
  return fun();
}

//Function to return the parms list passed to rootSolve as a SEXP object
//Nb. This requires the adapted rootSolve version
SEXP attribute_hidden get_rootSolve_gparms_Rcpp() {
  static SEXP(*fun)() = NULL;
  if (fun == NULL)
    fun = (SEXP(*)()) R_GetCCallable("rootSolve","get_rootSolve_gparms");
  return fun();
}

//Ensure C linkage so functions keep their names and can be identified by deSolve in the DLL/SO
extern "C" {
  void derivs(int *neq, double *t, double *y, double *ydot, double *yout, int *ip);
  void initmod(void (* odeparms)(int *, double *));
  void rt_initmod(void (* odeparms)(int *, double *));
  void cleanUp();
  void vaccineCampaignEvent(int *n, double *t, double *y);
}

//Keep global parameters in memory so they can be accessed in initMod() and derivs()
//std::vector<Cluster*> global_clusters;
std::vector<std::unique_ptr<Population>> populations;
int n_pops, n_agrp;

void vaccineCampaignEvent(int *n, double *t, double *y) {
  double time = t[0];
  
  //first set state of all compartments
  int start = 0;
  for(int p = 0; p < n_pops; p++){
    start += populations[p]->setState(y, start, time);
  }
  
  //now process vaccination
  start = 0;
  for(int p = 0; p < n_pops; p++){
    start += populations[p]->setStateVaccineCampaign(y, start, time);
  }
}

bool debug = false;

//This function sets the model up, and stores the parameter values in memory. It is only called once when setting up
// the model
void initmod(void (* odeparms)(int *, double *)) {
//std::cerr << "TEST IN initmod" << std::endl;
  //We get the parms argument passed to deSolve as SEXP object
  SEXP sparms = get_deSolve_gparms_Rcpp();
  
  try {
    //Parse parameters passed to deSolve as Rcpp::List
    Rcpp::List parms = Rcpp::clone(Rcpp::as<Rcpp::List>(sparms));
    
    //Define the number of trial arms/clusters and agegroups from parameter list passed to deSolve
    n_pops = Rcpp::as<Rcpp::List>(parms["trial_arms"]).size();
    if(n_pops == 0) n_pops = 1;
    n_agrp = parms["n_agrp"];
    
    //We can't do garbage collection at end of model run in deSolve, so we do it if the same DLL/SO is still loaded
    // and deSolve is ran again
    populations.clear();
    
    //We create a new Cluster object for every arm in the trial, and store a pointer to it in the clusters vector
    for(int p = 0; p < n_pops; p++){
      //Cluster* cluster = new Cluster(n_agrp, parms, c);
      //global_clusters.emplace_back(cluster);
      populations.emplace_back(std::make_unique<Population>(n_agrp, parms, p));
    }
    
    //Now all cluster objects exist, process the migration rates so these are consistent with the population size
    for(int p = 0; p < n_pops; p++){
      populations[p]->setMigrationRates(n_pops, p, populations);
    }
  } catch(std::exception& __ex__){
    forward_exception_to_r(__ex__);
  } catch(...){
    ::Rf_error( "c++ exception (unknown reason)" );
  }
}

//This function sets the model up, and stores the parameter values in memory. It is only called once when setting up
// the model
void rt_initmod(void (* odeparms)(int *, double *)) {

  //We get the parms argument passed to deSolve as SEXP object
  SEXP sparms = get_rootSolve_gparms_Rcpp();
  
  try {
    //Parse parameters passed to deSolve as Rcpp::List
    Rcpp::List parms = Rcpp::clone(Rcpp::as<Rcpp::List>(sparms));
    
    //Define the number of trial arms/clusters and agegroups from parameter list passed to deSolve
    n_pops = Rcpp::as<Rcpp::List>(parms["trial_arms"]).size();
    if(n_pops == 0) n_pops = 1;
    n_agrp = parms["n_agrp"];
    
    //We can't do garbage collection at end of model run in deSolve, so we do it if the same DLL/SO is still loaded
    // and deSolve is ran again
    populations.clear();
    
    //We create a new Cluster object for every arm in the trial, and store a pointer to it in the clusters vector
    for(int p = 0; p < n_pops; p++){
      populations.emplace_back(std::make_unique<Population>(n_agrp, parms, p));
    }
    
    //for(int p = 0; p < n_pops; p++){
    //  for(int t = 0; t < populations[p]->n_vstrat; t++){
    //    for(int c = 0; c < n_comps_prevalence; c++){
    //      Rcpp::Rcout << "DEBUG A pop: " << p << "; t: " << t << "; c: " << c << "; val: " << populations[p]->vac_strata[t]->compartments[c]->n_agrp;
    //    }
    //  }
    //}

    //Now all cluster objects exist, process the migration rates so these are consistent with the population size
    for(int p = 0; p < n_pops; p++){
      populations[p]->setMigrationRates(n_pops, p, populations);
    }

    //for(int p = 0; p < n_pops; p++){
    //  for(int t = 0; t < populations[p]->n_vstrat; t++){
    //    for(int c = 0; c < n_comps_prevalence; c++){
    //      Rcpp::Rcout << "DEBUG B pop: " << p << "; t: " << t << "; c: " << c << "; val: " << populations[p]->vac_strata[t]->compartments[c]->n_agrp;
    //    }
    //  }
    //}
  } catch(std::exception& __ex__){
    forward_exception_to_r(__ex__);
  } catch(...){
    ::Rf_error( "c++ exception (unknown reason)" );
  }
}

void cleanUp() {
  for(int p = 0; p < n_pops; p++){
    populations[p].reset();
  }
  populations.clear();
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
  //  Rcpp::Rcout << "DEBUG: time: " << time << std::endl;
  //  std::this_thread::sleep_for(std::chrono::milliseconds(5));
  //}
  
  //for(int p = 0; p < n_pops; p++){
  //    for(int t = 0; t < populations[p]->n_vstrat; t++){
  //      for(int c = 0; c < n_comps_prevalence; c++){
  //        Rcpp::Rcout << "DEBUG C pop: " << p << "; t: " << t << "; c: " << c << "; val: " << populations[p]->vac_strata[t]->compartments[c]->n_agrp;
  //      }
  //    }
  //  }

  //We loop through every cluster in the model, and update the state in each compartment. State is passed from deSolve
  // by the y array. The states are ordered as: cluster > vaccine_arm > compartment (S, VT, NVT, B) > agegroup.
  //The setState() function returns the total number of compartments updated by that cluster (which is
  // the number of vaccine strata in that cluster * 4 infection compartments * the number of agegroups)
  //We also update any time-varying parameters in the model
  //Rcpp::Rcout << "DEBUG: Set States" << std::endl;
  //std::this_thread::sleep_for(std::chrono::milliseconds(5));
  int start = 0;
  for(int p = 0; p < n_pops; p++){
    start += populations[p]->setState(y, start, time);
  }
  
  //Process demographic changes (ageing and migration) and vaccinations
  //can run in parallel if no migration in the model
  //std::cout << "DEBUG: Update Demographics" << std::endl;
  //std::this_thread::sleep_for(std::chrono::milliseconds(500));
  for(int p = 0; p < n_pops; p++){
    //std::cout << "DEBUG: Update Demographics, p: " << p << std::endl;
    //std::this_thread::sleep_for(std::chrono::milliseconds(50));
    populations[p]->updateDemographics(populations, n_pops, p); 
  }
  
  //We calculate the ODEs within all clusters, which is done in the calculateDerivs() function
  //if(n_pops > 1){
  //  Rcpp::Rcout << "DEBUG: Calculate Derivs" << std::endl;
  //  std::this_thread::sleep_for(std::chrono::milliseconds(5));
  //}
  #ifdef MP_ENABLED
  //#pragma omp parallel for
  #endif
  for(int p = 0; p < n_pops; p++){
    populations[p]->calculateDerivs(populations, n_pops, p, time);
  }
    
  //We return the model output to deSolve by copying the model output in the ydot array. The states are ordered
  // as: cluster > vaccine_arm > compartment (S, VT, NVT, B) > agegroup.
  int i = 0;
  //Rcpp::Rcout << "DEBUG: Return to deSolve" << std::endl;
  //std::this_thread::sleep_for(std::chrono::milliseconds(5));
  for(int p = 0; p < n_pops; p++){
    arma::rowvec deqs = populations[p]->getDerivs();
    for(int d = 0; d < deqs.size(); d++){
      ydot[i] = deqs(d);
      i += 1;
    }
  }
  
  //Additional return value to deSolve if incidence is returned
  if (ip[0] > 1){
    int i = 0;
    for(int p = 0; p < n_pops; p++){
      arma::rowvec incidence = populations[p]->getIncidence();
      for(int d = 0; d < incidence.size(); d++){
        yout[i] = incidence(d);
        i += 1;
      }
    }
  }
}