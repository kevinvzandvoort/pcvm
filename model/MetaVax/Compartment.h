#pragma once

class Compartment{
    public:
        int n_agrp;
        double prop_born;
        arma::rowvec value, delta, delta2;
        arma::rowvec arate, arate_corr;
        arma::rowvec age_none, age_out, age_in;
        arma::rowvec vac_c_out;
        arma::rowvec a0_v0_w1_m0, a0_v0_w1_m1, a0_v0_w0_m1, a1_v0_w0_m0, a1_v0_w1_m0, a1_v0_w1_m1, a1_v0_w0_m1, a1_v1_w0_m0, a1_v1_w0_m1;
    public:
    //Constructor
    Compartment(int n_agrp, arma::rowvec &arate, arma::rowvec &arate_corr, double prop_born)
        : n_agrp(n_agrp), prop_born(prop_born), arate(arate), arate_corr(arate_corr)
    {
        value = delta = delta2 = age_none = age_out = age_in = a0_v0_w1_m0 = a0_v0_w1_m1 = a0_v0_w0_m1 = a1_v0_w0_m0 = a1_v0_w1_m0 =
        a1_v0_w1_m1 = a1_v0_w0_m1 = a1_v1_w0_m0 = a1_v1_w0_m1 = vac_c_out = arma::rowvec(n_agrp, arma::fill::zeros);
    }
  
    //Deconstructor
    ~Compartment(){ }

  //This function assigns the current state of this compartment from the values calculated in deSolve. Values are stored
  // in one long array. *y points to the address of the first value.
  void setState(double *y, int c, int start){
    for(int a=0; a<n_agrp; a++){
      value(a) = y[start + n_agrp * c + a];
    }

    //Set to 0 at start of new iteration
    delta = delta2 = arma::rowvec(n_agrp, arma::fill::zeros);
  }

    //double check this
  void setStateVaccineCampaign(arma::rowvec &vac_cov_c, arma::rowvec &vac_in, double *y, int start, int c){
    //vaccinees moving out of this stratum
    vac_c_out = value % vac_cov_c;
    
    //update states in this TransComp
    for(int a=0; a<n_agrp; a++){
      //if(a == 6) Rcpp::Rcout << "DEBUG: Compartment::setStateVaccineCampaign start: " << start << "; c: " << c << "; vac_in[6]" << vac_in[6] << "; vac_cov_c[6]" << vac_cov_c[6] << "; vac_c_out[6]" << vac_c_out[6] << std::endl;
      y[start + n_agrp * c + a] = y[start + n_agrp * c + a] -vac_c_out(a) +vac_in(a);
    }
  }

  arma::rowvec getVacCOut(){
    return vac_c_out;
  }

  arma::rowvec getVacROutNoMigrate(){
    return a1_v1_w0_m0;
  }
  arma::rowvec getVacROutWhoMigrate(){
    return a1_v1_w0_m1;
  }
  arma::rowvec get_vac_none_wane_none_migr_out(){
    return a0_v0_w0_m1 + a1_v0_w0_m1;
  }
  arma::rowvec get_vac_none_wane_out_migr_out(){
    return a1_v0_w1_m1 + a0_v0_w1_m1;
  }
  arma::rowvec get_wane_out_migr_none(){
    return a1_v0_w1_m0 + a0_v0_w1_m0;
  }

  //This function updates the population (processes ageing and migration)
  void updateDemographics(double &N_all_ageing, arma::rowvec &vac_cov_r,  arma::rowvec &wrate, arma::rowvec &mrate, arma::rowvec &vac_in){
    //calculate ageing
    //those that do not age
    age_none = (1.0 - arate) % value;
    //those that do age (out)
    age_out = arate % value;
    //those that age in
    age_in.subvec(1, (n_agrp-1)) = age_out.subvec(0, (n_agrp-2)) % arate_corr.subvec(1, (n_agrp-1));
    //proportion of newborns that flow to this compartment (often all assumed to be susceptible)
    age_in(0) = arate(n_agrp-1) * arate_corr(0) * prop_born * N_all_ageing;
    
    //do not age, do not migrate, do wane
    a0_v0_w1_m0 = age_none % wrate % (1.0 - mrate);
    //do not age, do wane and migrate
    a0_v0_w1_m1 = age_none % wrate % mrate;
    //do not age, do not wane, do migrate
    a0_v0_w0_m1 = age_none % (1.0 - wrate) % mrate;
    //do age, are not vaccinated, do not wane or migrate
    a1_v0_w0_m0 = age_in % (1.0 - vac_cov_r) % (1.0 - wrate) % (1.0 - mrate);
    //do age, are not vaccinated, do not migrate, do wane
    a1_v0_w1_m0 = age_in % (1.0 - vac_cov_r) % wrate % (1.0 - mrate);
    //do age, are not vaccinated, do wane and migrate
    a1_v0_w1_m1 = age_in % (1.0 - vac_cov_r) % wrate % mrate;
    //do age, are not vaccinated, do not wane, do migrate
    a1_v0_w0_m1 = age_in % (1.0 - vac_cov_r) % (1.0 - wrate) % mrate;
    //do age, are vaccinated, do not migrate
    a1_v1_w0_m0 = age_in % vac_cov_r % (1.0 - mrate);
    //do age, are vaccinated, do migrate
    a1_v1_w0_m1 = age_in % vac_cov_r % mrate;

    //The change in states with in- and out-ageing, in- and out-vaccination (routine), out-migration, and out-waning, but without in-migration and in-waning
    delta += -(age_out +a0_v0_w1_m0 +a0_v0_w1_m1 +a0_v0_w0_m1) +(a1_v0_w0_m0 +vac_in);
  }

  void updateDerivs(arma::rowvec &additional_value){
    delta += additional_value;
  }
      
  void updateDerivs2(arma::rowvec &additional_value){
    delta2 += additional_value;
  }    

  arma::rowvec getDelta(){
    //Rcpp::Rcout << "DEBUG: Compartment::getDelta " << std::endl;
    return(delta + delta2);
  }

  arma::rowvec getValue(){
    return(value);
  }

  //returns value after processing any already calculated delta
  arma::rowvec getUpdatedValue(){
    return(value + delta + delta2);
  }
      
  //returns value after processing any already calculated delta
  arma::rowvec getUpdatedValueOnlyDemo(){
    return(value + delta);
  }
};