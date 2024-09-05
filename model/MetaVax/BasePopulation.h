#pragma once

//forward declarations
class Population;
class VaccinationGroup;

//A BasePopulation object represents an individual cluster or arm in the trial. It will hold multiple vaccinated strata
class BasePopulation{
public:
  int n_agrp;
  arma::rowvec arate, arate_corr;
  int n_vstrat;
  std::vector<std::unique_ptr<VaccinationGroup>> vac_strata; 
  arma::rowvec deqs_all, deqs;
  arma::rowvec incidence_all, incidence;
  arma::rowvec population_size, population_change;
  arma::rowvec N;
  double N_ageing;
  arma::mat travel, migration;
  arma::rowvec vac_out, mgr_out, mgr_out_cl, mgr_out_wane, mgr_out_wane_cl, wane_out_nomigr;
  std::vector<arma::rowvec> mrates;
  arma::rowvec mrates_total;
 
  BasePopulation(int n_agrp, Rcpp::List parms, int p);
  ~BasePopulation();
  void setMigrationRates(int n_pops, int p, std::vector<std::unique_ptr<Population>> & populations);
  void updateDemographics(std::vector<std::unique_ptr<Population>> & populations, int n_pops, int p);
  void updateParams(double & time);
  int setState(double *y, int start, double & time, bool & solver_difference);
  int setStateVaccineCampaign(double *y, int start, double & time);
  void addMigrants(int origin, int t, int c, arma::rowvec &delta);
  arma::rowvec& getDerivs();
  arma::rowvec& getIncidence();
  arma::rowvec& getPopSize();
};