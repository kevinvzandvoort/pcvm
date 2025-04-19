#pragma once

BasePopulation::BasePopulation(int n_agrp, Rcpp::List parms, int p)
    : n_agrp(n_agrp),
      arate(Rcpp::as<arma::rowvec>(Rcpp::as<Rcpp::List>(parms["global_settings"])["ageout"])),
      travel(Rcpp::as<arma::mat>(Rcpp::as<Rcpp::List>(parms["global_settings"])["travel"])),
      migration(Rcpp::as<arma::mat>(Rcpp::as<Rcpp::List>(parms["global_settings"])["migration"]))
  {
    //We take the total number of required vaccinated strata from the trial_arms element in the parms list passed to
    // deSolve.
    Rcpp::List trial_arm = Rcpp::as<Rcpp::List>(parms["populations"])[p];
    Rcpp::List cluster_parameters = Rcpp::as<Rcpp::List>(trial_arm["parameters"]);
    Rcpp::List vstrata = Rcpp::as<Rcpp::List>(trial_arm["vaccination_groups"]);
    
    //set population size
    population_size = Rcpp::as<arma::rowvec>(cluster_parameters["N"]);
    population_change = arma::rowvec(n_agrp, arma::fill::zeros);
    
    //Need to correct age rates to account for strata of different size
    arate_corr = arma::rowvec(n_agrp, arma::fill::zeros);
    arate_corr(0) = arate(0)/arate(n_agrp-1);
    arate_corr.subvec(1, (n_agrp-1)) = arate.subvec(1, (n_agrp-1))/arate.subvec(0, (n_agrp-2));
    
    //Vaccinated strata are stored in a vector with TransComp objects
    vac_strata.reserve(vstrata.size());
    n_vstrat = vstrata.size();
  
    //list with empty vaccine coverage
    Rcpp::List coverage_null_inner = Rcpp::List::create(
      Rcpp::Named("value", Rcpp::wrap(arma::rowvec(n_agrp, arma::fill::zeros))),
      Rcpp::Named("time", 0.0),
      Rcpp::Named("coverage_to", Rcpp::wrap(arma::rowvec(n_agrp, arma::fill::zeros)))
    );
    Rcpp::List coverage_null = Rcpp::List::create(Rcpp::clone(coverage_null_inner));
    
    //Adjust mrates so population size will remain fixed
    mrates_total = arma::rowvec(n_agrp, arma::fill::zeros);
    
    //Add subsequent strata for vaccinated transmission compartments
    //Rcpp::Rcout << "BasePopulation:: Q - create vstrata " << std::endl;
    for(int t=0; t < vstrata.size(); t++){
      Rcpp::List vac_parms = vstrata[t];
      
      //unvaccinated stratum has no waning or VE
      if(t == 0){
        vac_parms["waning"] = Rcpp::wrap(arma::rowvec(n_agrp, arma::fill::zeros));
        vac_parms["efficacy_transmission"] = Rcpp::wrap(arma::rowvec(n_agrp, arma::fill::zeros));
      }
      
      //final stratum has no coverage
      if(t == (vstrata.size()-1)){
        vac_parms["coverage_c"] = Rcpp::clone(coverage_null);
        vac_parms["coverage_r"] = Rcpp::clone(coverage_null);
      }
      
      vac_strata.emplace_back(std::make_unique<VaccinationGroup>(n_agrp, vac_parms, arate, arate_corr)); //object is already destroyed here...
    }
    
    //We preallocate memory for some variables for efficiency
    deqs_all = arma::rowvec(n_comps_prevalence * n_agrp * n_vstrat, arma::fill::zeros);
    incidence_all = arma::rowvec(n_comps_incidence * n_agrp * n_vstrat, arma::fill::zeros);
    deqs = arma::rowvec(n_comps_prevalence * n_agrp, arma::fill::zeros);
    incidence = arma::rowvec(n_comps_incidence * n_agrp, arma::fill::zeros);
    wane_out_nomigr = vac_out = mgr_out = mgr_out_cl = mgr_out_wane = mgr_out_wane_cl = arma::rowvec(n_agrp, arma::fill::zeros);
    N_ageing = 0.0;
}

BasePopulation::~BasePopulation(){
    //Rcpp::Rcout << "Destroy cluster" << std::endl;
    //Rcpp::Rcout << "Deconstructing BasePopulation" << std::endl;
    vac_strata.clear();
  }

void BasePopulation::setMigrationRates(int n_pops, int p, std::vector<std::unique_ptr<Population>> & populations){
    for(int j=0; j < n_pops; j++){
      arma::rowvec mrate_by_age;
      if(migration(j, p) < 0.0){
        //Make sure total number of migrants is consistent
        mrate_by_age = migration(p, j) * populations[j]->getPopSize()/population_size;
      } else {
        mrate_by_age = migration(j, p) * arma::rowvec(n_agrp, arma::fill::ones);
      }
      mrates.emplace_back(mrate_by_age);
      mrates_total += mrate_by_age;
      //Rcpp::Rcout << "DEBUG BasePopulation::setMigrationRates pop p: " << p << "; j: " << j << "; mrate_by_age: " << mrate_by_age << std::endl;
    }
  }

  //This function updates the demographics for each arm in the cluster (process ageing, migration, and vaccination)
  void BasePopulation::updateDemographics(std::vector<std::unique_ptr<Population>> & populations, int n_pops, int p){
    N = arma::rowvec(n_agrp, arma::fill::zeros);
    for(int t=0; t < n_vstrat; t++){
      N += vac_strata[t]->getN();
    }
    
    //Loop through all vaccine strata to calculate the derivatives at this timestep
    for(int t=0; t < n_vstrat; t++){
      arma::rowvec vac_cov_r = vac_strata[t]->get_vac_cov_r();
      arma::rowvec wrate = vac_strata[t]->get_vac_waning();
      
      for(int c=0; c < n_comps_prevalence; c++){
        //initially, no people move into this stratum
        vac_out = arma::rowvec(n_agrp, arma::fill::zeros);
        mgr_out = arma::rowvec(n_agrp, arma::fill::zeros);

        if(t == 0){
          //Newborns are all born into the first stratum (unvaccinated), and calculated as a proportion of all those who
          // age in all strata (N_ageing).
          N_ageing = N(n_agrp -1);

          //No vaccinated individuals move into the first stratum (unvaccinated).
        } else {
          //No newborns enter the vaccinated strata.
          N_ageing = 0.0;
          
          //Effectively vaccinated people move from previous strata to the current stratum (if they do not migrate).
          for(int v = 0; v < t; v++){
            
            //check if people move into this arm
            arma::umat cov_this_arm = vac_strata[v]->get_vac_cov_r_to() == arma::rowvec(n_agrp, arma::fill::value(t));
            if(accu(cov_this_arm) > 0){
              vac_out += vac_strata[v]->compartments[c]->getVacROutNoMigrate() % cov_this_arm;
              mgr_out += vac_strata[v]->compartments[c]->getVacROutWhoMigrate() % cov_this_arm;
            }
          }
        }
        
        //We update the derivatives for this specific stratum through the calculateDerivs() function.
        vac_strata[t]->compartments[c]->updateDemographics(N_ageing, vac_cov_r, wrate, mrates_total, vac_out);

        //get those who wane, but do not migrate
        wane_out_nomigr = vac_strata[t]->compartments[c]->get_wane_out_migr_none();
        vac_strata[0]->compartments[c]->updateDerivs(wane_out_nomigr);

        //Also add those migrating but not vaccinated (remain in this arm in other cluster)
        mgr_out += vac_strata[t]->compartments[c]->get_vac_none_wane_none_migr_out();
        mgr_out_wane = vac_strata[t]->compartments[c]->get_vac_none_wane_out_migr_out();

        //Process migration
        for(int j=0; j<n_pops; j++){
          //no migration to the same population
          if(p == j || arma::accu(mrates_total) == 0 || arma::accu(mrates[j]) == 0){ 
            continue;
          }
          //Calculate proportion that goes to this cluster
          //-could probably just divide by mrates_total
          mgr_out_cl = mgr_out % (mrates[j] / mrates_total) / mrates[j];
          populations[j]->addMigrants(p, t, c, mgr_out_cl);
          mgr_out_wane_cl = mgr_out_wane % (mrates[j] / mrates_total);
          populations[j]->addMigrants(p, 0, c, mgr_out_wane_cl); //those that wane always go to first stratum
        }
      }
    }
  }

  void BasePopulation::updateParams(double & time){
    //update any time parameters on the cluster level
  }

  int BasePopulation::setState(double *y, int start, double & time, bool & solver_difference){
    //Rcpp::Rcout << "DEBUG: Population::setState" << std::endl;
    //std::this_thread::sleep_for(std::chrono::milliseconds(5));
    updateParams(time);
    
    for(int t=0; t < n_vstrat; t++){
      for(int c=0; c < n_comps_prevalence; c++){
        //Rcpp::Rcout << "DEBUG: Population::setState; pre - n_agrp: " << vac_strata[t]->compartments[c]->n_agrp << std::endl;
        //Rcpp::Rcout << "DEBUG: Population::setState t: " << t << "; c: " << c << "; start: " << start << "; n_agrp: " << n_agrp << "; n_comps_prevalence" << n_comps_prevalence << "; pass: " << start + n_agrp * n_comps_prevalence * t << std::endl;
        //std::this_thread::sleep_for(std::chrono::milliseconds(50));
        vac_strata[t]->compartments[c]->setState(y, c, start + n_agrp * n_comps_prevalence * t);
      }
      
      //Rcpp::Rcout << "DEBUG: Population::setState t: " << t << " - updateParams" << std::endl;
      //std::this_thread::sleep_for(std::chrono::milliseconds(5));
      vac_strata[t]->updateParams(time, solver_difference);
    }
    
    //Rcpp::Rcout << "DEBUG: Population::setState incidence" << std::endl;
    //std::this_thread::sleep_for(std::chrono::milliseconds(5));
    //also make sure incidence is empty
    incidence_all = arma::rowvec(n_comps_incidence * n_agrp * n_vstrat, arma::fill::zeros);
    
    //Rcpp::Rcout << "DEBUG: Population::setState finished and return: " << vac_strata.size() * n_comps_prevalence * n_agrp << std::endl;
    //std::this_thread::sleep_for(std::chrono::milliseconds(5));
    return n_vstrat * n_comps_prevalence * n_agrp;
  }

  int BasePopulation::setStateVaccineCampaign(double *y, int start, double & time){
    arma::rowvec vac_in;
    
    for(int t=0; t < n_vstrat; t++){
      arma::rowvec vac_cov_c = vac_strata[t]->getCoverageVaccineCampaign(time);
      
      for(int c=0; c < n_comps_prevalence; c++){
        vac_in = arma::rowvec(n_agrp, arma::fill::zeros);
        
        if(t > 0){
          //Effectively vaccinated people move from previous strata to the current stratum
          for(int v = 0; v < t; v++){
            //check if people move into this arm
            arma::umat cov_this_arm = vac_strata[v]->get_vac_cov_c_to() == arma::rowvec(n_agrp, arma::fill::value(t));
            if(accu(cov_this_arm) > 0){
              vac_in += vac_strata[v]->compartments[c]->getVacCOut() % cov_this_arm;
            }
          }
        }
      
        vac_strata[t]->compartments[c]->setStateVaccineCampaign(vac_cov_c, vac_in, y, start + n_agrp * n_comps_prevalence * t, c);
      }
    }
    
    //return new start value
    return n_vstrat * n_comps_prevalence * n_agrp;
  }

  void BasePopulation::addMigrants(int origin, int t, int c, arma::rowvec &delta){
    //rates are multiplied by mrate for this cluster, to keep population sizes constant
    arma::rowvec delta_adj = delta % mrates[origin];
    //Rcpp::Rcout << "DEBUG BasePopulation::addMigrants origin: " << origin << "; t: " << t << "; c: " << c << "; delta[10]: " << delta[10] << "; mrate[origin][10]: " << mrates[origin][10] << "; delta_adj: " << delta_adj[10] << std::endl;
    vac_strata[t]->compartments[c]->updateDerivs(delta_adj);
  }

  arma::rowvec& BasePopulation::getDerivs(){
    //Rcpp::Rcout << "DEBUG: BasePopulation::getDerivs " << std::endl;
    for(int t=0; t < n_vstrat; t++){
      for(int c=0; c < n_comps_prevalence; c++){
        deqs.subvec(0 + n_agrp*c, (n_agrp+n_agrp*c - 1)) = vac_strata[t]->compartments[c]->getDelta();  
      }
      //Rcpp::Rcout << "DEBUG: BasePopulation::getDerivs - finish loop c" << std::endl;
      deqs_all.subvec(0 + n_comps_prevalence*n_agrp*t, (n_comps_prevalence*n_agrp + n_comps_prevalence*n_agrp*t - 1)) = deqs;
    }
    //Rcpp::Rcout << "DEBUG: BasePopulation::getDerivs - finish loop t" << std::endl;
    return deqs_all;
  }

  arma::rowvec& BasePopulation::getIncidence(){
    for(int t=0; t < n_vstrat; t++){
      incidence_all.subvec(0 + n_comps_incidence*n_agrp*t, (n_comps_incidence*n_agrp + n_comps_incidence*n_agrp*t - 1)) = vac_strata[t]->getIncidence();
    }
    
    return incidence_all;
  }

  arma::rowvec& BasePopulation::getPopSize(){ return population_size; }