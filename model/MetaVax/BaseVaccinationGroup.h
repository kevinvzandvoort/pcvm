#pragma once

//A VaccinationGroup object represents the age-stratified transmission compartments for one individual vaccine stratum
// (e.g. those unvaccinated, or vaccinated with a single dose) within an individual population in the trial.
// It is nested within a Population object, and holds the states, parameters, and functions to calculate the ODEs for this specific stratum.
//base class for VaccinationGroup
class BaseVaccinationGroup{
private:
  arma::rowvec vac_cov_r, vac_cov_c, vac_waning, vac_cov_r_to, vac_cov_c_to, N;
  int vac_cov_r_index, vac_cov_c_index;
  Rcpp::List vac_cov_r_values, vac_cov_c_values;
  bool vac_cov_r_change_final, vac_cov_c_change_final, vac_cov_c_implemented;
  double vac_cov_r_change_time, vac_cov_c_change_time;
public:
  int n_agrp;
  arma::rowvec vac_eff;
  std::vector<Compartment*> compartments;
  //Constructor
  BaseVaccinationGroup(int &n_agrp, Rcpp::List &vac_parms)
    : vac_waning(Rcpp::as<arma::rowvec>(vac_parms["waning"])), vac_cov_r_index(0), vac_cov_c_index(0),
      vac_cov_r_values(Rcpp::as<Rcpp::List>(vac_parms["coverage_r"])), vac_cov_c_values(Rcpp::as<Rcpp::List>(vac_parms["coverage_c"])),
      vac_cov_c_implemented(false), n_agrp(n_agrp), vac_eff(Rcpp::as<arma::rowvec>(vac_parms["efficacy_transmission"]))
  {
    vac_cov_r_change_time = (vac_cov_r_values.size() == 0 ? 999999 : Rcpp::as<Rcpp::List>(vac_cov_r_values[vac_cov_r_index])["time"]);
    if(vac_cov_r_change_time > 0.0){
      //if first vaccination time is in the future, initialize vaccination coverage at 0
      vac_cov_r = arma::rowvec(n_agrp, arma::fill::zeros);
      vac_cov_r_to = arma::rowvec(n_agrp, arma::fill::zeros);
    } else {
      vac_cov_r = Rcpp::as<arma::rowvec>(Rcpp::as<Rcpp::List>(vac_cov_r_values[vac_cov_r_index])["value"]);
      vac_cov_r_to = Rcpp::as<arma::rowvec>(Rcpp::as<Rcpp::List>(vac_cov_r_values[vac_cov_r_index])["coverage_to"]);
    }
    vac_cov_r_change_final = vac_cov_r_index == (vac_cov_r_values.size() - 1); //check if this is the final value to be updated
    //vac_cov_r_change_time = (vac_cov_r_change_final ? 999999 : Rcpp::as<Rcpp::List>(vac_cov_r_values[vac_cov_r_index+1])["time"]);
    
    vac_cov_c = Rcpp::as<arma::rowvec>(Rcpp::as<Rcpp::List>(vac_cov_c_values[vac_cov_c_index])["value"]);
    vac_cov_c_to = Rcpp::as<arma::rowvec>(Rcpp::as<Rcpp::List>(vac_cov_c_values[vac_cov_c_index])["coverage_to"]);
    vac_cov_c_change_final = vac_cov_c_index == (vac_cov_c_values.size() - 1); //check if this is the final value to be updated
    vac_cov_c_change_time = (vac_cov_c_values.size() == 0 ? 999999 : Rcpp::as<Rcpp::List>(vac_cov_c_values[vac_cov_c_index])["time"]);

    N = arma::rowvec(n_agrp, arma::fill::zeros);
    //Compartments are stored in a vector with Compartment objects, to be created when class is extended
  }
  
  //Deconstructor
  ~BaseVaccinationGroup(){  }
  
  arma::rowvec getCoverageVaccineCampaign(double & time){
    //no campaign coverage unless we are doing a campaign
    vac_cov_c = arma::rowvec(n_agrp, arma::fill::zeros);
    
    if(time >= vac_cov_c_change_time){
      
      //Rcpp::Rcout << "DEBUG - Doing a campaign; TIME: " << time << std::endl;
      vac_cov_c = Rcpp::as<arma::rowvec>(Rcpp::as<Rcpp::List>(vac_cov_c_values[vac_cov_c_index])["value"]);
      vac_cov_c_to = Rcpp::as<arma::rowvec>(Rcpp::as<Rcpp::List>(vac_cov_c_values[vac_cov_c_index])["coverage_to"]);
      vac_cov_c_change_final = vac_cov_c_index == (vac_cov_c_values.size()-1); //check if this is the final value to be updated
      
      vac_cov_c_index++;
      if(!vac_cov_c_change_final){
        vac_cov_c_change_time = Rcpp::as<Rcpp::List>(vac_cov_c_values[vac_cov_c_index])["time"]; //check at which time the value changes next 
      } else {
        vac_cov_c_change_time = 999999;
      }
    }

    return vac_cov_c;
  }
  
  void updateParams(double & time, bool & solver_difference){
    //update any time parameters on the transcomp level
    if(time >= vac_cov_r_change_time){
      
      vac_cov_r = Rcpp::as<arma::rowvec>(Rcpp::as<Rcpp::List>(vac_cov_r_values[vac_cov_r_index])["value"]);
      vac_cov_r_to = Rcpp::as<arma::rowvec>(Rcpp::as<Rcpp::List>(vac_cov_r_values[vac_cov_r_index])["coverage_to"]);
      vac_cov_r_change_final = vac_cov_r_index == (vac_cov_r_values.size()-1); //check if this is the final value to be updated
      //if(!vac_cov_r_change_final) vac_cov_r_change_time = Rcpp::as<Rcpp::List>(vac_cov_r_values[vac_cov_r_index+1])["time"]; //check at which time the value changes next
      
      vac_cov_r_index++;
      if(!vac_cov_r_change_final){
        vac_cov_r_change_time = Rcpp::as<Rcpp::List>(vac_cov_r_values[vac_cov_r_index])["time"]; //check at which time the value changes next 
      } else {
        vac_cov_r_change_time = 999999;
      }
    }
    
    //can do the same for catch-up coverage if using difference equations, otherwise needs to be implemented using events
    //if(solver_difference){
    //  if(!vac_cov_c_change_final && time >= vac_cov_c_change_time){
    //    Rcpp::Rcout << "DEBUG: at time " << time << " update vac_cov_c " << vac_cov_c_index << " to " << vac_cov_c_index+1 << "; change time " << vac_cov_c_change_time << std::endl;
    //    std::this_thread::sleep_for(std::chrono::milliseconds(50));
    //    vac_cov_c_index++;
    //    vac_cov_c = Rcpp::as<arma::rowvec>(Rcpp::as<Rcpp::List>(vac_cov_c_values[vac_cov_c_index])["value"]);
    //    vac_cov_c_to = Rcpp::as<arma::rowvec>(Rcpp::as<Rcpp::List>(vac_cov_c_values[vac_cov_c_index])["coverage_to"]);
    //    vac_cov_c_change_final = vac_cov_c_index == (vac_cov_c_values.size()-1); //check if this is the final value to be updated
    //    if(!vac_cov_c_change_final) vac_cov_c_change_time = Rcpp::as<Rcpp::List>(vac_cov_c_values[vac_cov_c_index+1])["time"]; //check at which time the value changes next
    //  } 
    //}
  }
  
  arma::rowvec& getN(){
    N = arma::rowvec(n_agrp, arma::fill::zeros);
    for(int c = 0; c < (int) compartments.size(); c++){
      N += compartments[c]->getValue();
    }
    return N;
  }
  
  arma::rowvec& get_vac_cov_r_to(){ return vac_cov_r_to; }
  arma::rowvec& get_vac_cov_c_to(){ return vac_cov_c_to; } 
  arma::rowvec& get_vac_waning(){ return vac_waning; }
  arma::rowvec& get_vac_cov_r(){ return vac_cov_r; }
};