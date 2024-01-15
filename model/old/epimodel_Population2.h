#pragma once

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