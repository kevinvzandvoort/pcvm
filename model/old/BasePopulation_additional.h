#pragma once

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
    //std::cout << "DEBUG: BasePopulation::updateDemographics, p: " << p << std::endl;
    //std::this_thread::sleep_for(std::chrono::milliseconds(50));

    N = arma::rowvec(n_agrp, arma::fill::zeros);
    for(int t=0; t < n_vstrat; t++){
      N += vac_strata[t]->getN();
    }
    //std::cout << "DEBUG: N: " << N << std::endl;
    //std::this_thread::sleep_for(std::chrono::milliseconds(50));
    
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