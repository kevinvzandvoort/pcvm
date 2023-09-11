#' Adjust a case carrier ratio based on malnutrition prevalence
adjustCaseCarrierRatio = function(ccr, kilifi_maln_prev = 0.2, maln_RR = 2){
  base_risk = ccr/(1 - kilifi_maln_prev + kilifi_maln_prev*maln_RR)
  return(base_risk)
}

#' Calculate the multinomial log likelihood  
multinomial_log_ll = function(data){
  if(nrow(data[modelled <= 0]) > 0){
    warning("At least one modelled value was <0, rejecting the proposal.")
    return(-Inf)
  } else
    data[, LL := observed * log(modelled)] %>%
    .[, .(LL = sum(LL)), by = c("age")] %>%
    .[, LL] %>%
    return
}

#' tryCatch to wrap around the deSolve functions
#' - there may be runs where deSolve returns an error, in which case we do not want the MCMC to fail
tryCatchWE = function(expr){
  status = 0
  wHandler = function(w) status <<- 1
  eHandler = function(e) status <<- 2
  value = list(value = tryCatch(expr, error = eHandler, warning = wHandler), status = status)
  return(value)
}

#' Should work to compile on all platforms
#' Mostly based on inline::cxxfunction
#' Compiled model is slightly slower compared to manual compilation
#' - TODO: figure out why
#' - system('g++ -c -o ./model/build/Fit_by_arm_quick_v5.o ./model/Fit_by_arm_quick_v5.cpp -L/usr/lib/R/lib -lR -std=gnu++11 -shared -L/usr/lib/R/lib -Wl,-Bsymbolic-functions,-z,relro -I"/usr/share/R/include" -I"/home/lsh1604011/workspace/Dosing_schedule_and_fitting" -I"/home/lsh1604011/R/x86_64-pc-linux-gnu-library/4.0/Rcpp/include" -I"/home/lsh1604011/R/x86_64-pc-linux-gnu-library/4.0/RcppArmadillo/include" -fPIC -DNDEBUG -fopenmp -g -O3 -fdebug-prefix-map=/build/r-base-8T8CYO/r-base-4.0.3=. -fstack-protector-strong -Wformat -Werror=format-security -Wdate-time -D_FORTIFY_SOURCE=2 -g')
#' - system('g++ -o ./model/build/Fit_by_arm_quick_v5.so ./model/build/Fit_by_arm_quick_v5.o -L/usr/lib/R/lib -lR -std=gnu++11 -shared -L/usr/lib/R/lib -Wl,-Bsymbolic-functions,-z,relro -I"/usr/share/R/include" -I"/home/lsh1604011/workspace/espicc_model" -I"/home/lsh1604011/R/x86_64-pc-linux-gnu-library/4.0/Rcpp/include" -fPIC -fopenmp -llapack -lblas -lgfortran -lm -lquadmath -L/usr/lib/R/lib -lR')
compileModel = function(cpp_file, shrd_lib_loc){
  if(!file.exists(cpp_file))
    stop("Cpp file does not exist")
  shrd_lib_name = gsub(".cpp", .Platform$dynlib.ext, basename(cpp_file))
  if(file.exists(sprintf("%s/%s", shrd_lib_loc, shrd_lib_name)))
    file.remove(sprintf("%s/%s", shrd_lib_loc, shrd_lib_name))
  
  cmd = paste0(R.home(component = "bin"), "/R")
  
  #' Need to include these directories  
  paths = sapply(c("Rcpp", "RcppArmadillo"), find.package)
  flag = paste(paste0("-I\"", paths, "/include\""), collapse = " ")
  
  #' Need to set additional compiler flags
  do.call(Sys.setenv, inline::getPlugin("RcppArmadillo")$env)
  Sys.setenv(CLINK_CPPFLAGS = flag)
  Sys.setenv(PKG_CXXFLAGS="-std=c++14")
  
  #' Compile model
  system2(cmd, args = paste(" CMD SHLIB -o", sprintf("%s/%s", shrd_lib_loc, shrd_lib_name), cpp_file))
  if(!file.exists(sprintf("%s/%s", shrd_lib_loc, shrd_lib_name)))
    stop("Something went wrong")
  
  #' Remove object file
  if(file.exists(gsub(".cpp", ".o", cpp_file)))
    file.remove(gsub(".cpp", ".o", cpp_file))
}

reshapeModelOutput = function(model_output, param.list){
  if(PCVM_VERSION == 2){
    compartments = c("S", "VT", "NVT", "VT2", "B", "NVT2")
  } else {
    compartments = c("S", "VT", "NVT", "B")
  }
  N_agegroups = param.list$n_agrp
  cluster_arms = param.list$trial_arms %>%
    lapply(function(x){
      c("unvaccinated", names(x[["arms"]])) %>% lapply(function(arm_name){
        data.table(
          dose = factor(arm_name, c("unvaccinated", names(x[["arms"]]))) %>% rep(length(compartments) * N_agegroups),
          compartment = factor(compartments, compartments) %>% rep(each = N_agegroups),
          age = 1:N_agegroups %>% rep(length(compartments)))
      }) %>% rbindlist})
  for(cl in 1:length(param.list$trial_arms)){
    clname = names(param.list$trial_arms)[cl]
    cluster_arms[[cl]][, cluster := factor(clname, names(param.list$trial_arms))]
  }
  cluster_arms = cluster_arms %>% rbindlist
  cluster_arms = cluster_arms[, c("cluster", "dose", "compartment", "age")] %>% setorder(cluster, dose, compartment, age)
  cluster_arms[, compartment_index := as.character(1:.N)]
  
  out_prevalence = model_output[, c(colnames(model_output) %>% subset(!grepl("incidence", .) & !grepl("N_", .))), with=FALSE] %>%
    melt(id.vars = "time", variable.name = "compartment_index", value.name = "prevalence")
  
  if(any(colnames(model_output) %>% grepl("incidence", .))){
    out_incidence = model_output[, c("time", colnames(model_output) %>% subset(grepl("incidence", .))), with=FALSE] %>%
      melt(id.vars = "time", variable.name = "compartment_index", value.name = "incidence") %>%
      .[, compartment_index := compartment_index %>% gsub("incidence_", "", .)]
    
    out = out_prevalence %>% merge(out_incidence, by=c("compartment_index", "time"))
    out = out %>% melt(measure.vars = c("prevalence", "incidence"))
  } else {
    out = out_prevalence %>% melt(measure.vars = c("prevalence"))
  }
  
  out = out %>% merge(cluster_arms, by = c("compartment_index"), all = TRUE)
  out %>% setorder(variable, cluster, dose, compartment, age, time)
  if(PCVM_VERSION == 2){
    out[variable == "incidence", compartment := switch(as.character(compartment), "S" = "S_VT", "VT" = "S_NVT", "NVT" = "VT_VT2", "VT2" = "VT_B", "B" = "NVT_B", "NVT2" = "NVT_NVT2"), by = c("compartment")]
  } else {
    out[variable == "incidence", compartment := switch(as.character(compartment), "S" = "S_VT", "VT" = "S_NVT", "NVT" = "VT_B", "B" = "NVT_B"), by = c("compartment")]
  }
  
  out = out[, c("variable", "cluster", "dose", "compartment", "age", "time", "value")]
  
  if(any(colnames(model_output) %>% grepl("N_", .))){
    out_population = model_output[, c("time", colnames(model_output) %>% subset(grepl("N_", .))), with=FALSE] %>%
      melt(id.vars = "time", variable.name = "cluster_age_index", value.name = "prevalence") %>%
      .[, cluster_age_index := cluster_age_index %>% gsub("N_", "", .) %>% as.numeric()] %>%
      merge(data.table(age = 1:N_agegroups %>% rep(length(param.list$trial_arms)),
                       cluster = names(param.list$trial_arms) %>% rep(each = N_agegroups),
                       cluster_age_index = seq_len(N_agegroups * length(param.list$trial_arms))),
            by = "cluster_age_index") %>% .[, c("compartment", "variable", "value") := .("N", "prevalence", prevalence)] %>%
      .[, c("variable", "cluster", "compartment", "age", "time", "value")]
    
    out = out %>% rbind(out_population, fill=TRUE)
  }
  
  return(out)
}

statePropModelOutput = function(model_output, param.list, catchup_ages = 0:14){
  if(length(model_output[, unique(time)]) != 1) stop("Table model_output needs to be for a single timestep")
  if(length(param.list$trial_arms) == 0) stop("There are no vaccine arms specified")
  n_catchup_arms = sum(sapply(param.list$trial_arms, function(x) "catchup" %in% names(x)))
  apply_catchup = TRUE
  if(n_catchup_arms == 0){
    apply_catchup = FALSE
    warning("There is no catchup dose specified in the trial_arms argument. No catchup campaign will be implemented and variable catchup_ages will be ignored")
  } else if(n_catchup_arms > 0 & n_catchup_arms < length(param.list$trial_arms)){
    stop("Some trial_arms have catchup campaign specified and some have not. States will not be correctly calculated for arms where no catchup campaign is specified")
  }
  
  n_agrp = param.list$n_agrp
  n_clus = length(param.list$trial_arms)
  state_prop = numeric(sum(sapply(param.list$trial_arms, length) + 1) * ifelse(PCVM_VERSION == 2, 6, 4) * n_agrp)
  i=1
  for(cl in 1:n_clus)
    for(d in c("unvaccinated", names(param.list$trial_arms[[cl]]))){
      state_prop[i:(i + ifelse(PCVM_VERSION == 2, 6, 4) * n_agrp - 1)] = switch(
        d,
        "unvaccinated" = model_output[order(compartment, age), .(value = value * ifelse(apply_catchup, !age %in% catchup_ages, 1)), by=seq_len(model_output[, .N])][, value],
        "catchup" = model_output[order(compartment, age), .(value = value * ifelse(apply_catchup, age %in% catchup_ages, 1)), by=seq_len(model_output[, .N])][, value],
        rep(0, model_output[, .N])
      )
      i = i + ifelse(PCVM_VERSION == 2, 6, 4) * n_agrp
    }
  return(state_prop)
}

eqStatesVaccinate = function(model_output, param.list){
  if(length(model_output[, unique(time)]) != 1) stop("Table model_output needs to be for a single timestep")
  
  model_output %>% setorder(cluster, dose, compartment, age)
  
  compartments = c("S", "VT", "NVT", "B")
  cluster_arms = param.list$trial_arms %>%
    lapply(function(x){
      c("unvaccinated", names(x[["arms"]])) %>% lapply(function(arm_name){
        data.table(
          dose = factor(arm_name, c("unvaccinated", names(x[["arms"]]))) %>% rep(length(compartments) * age_groups_model[, .N]),
          compartment = factor(compartments, compartments) %>% rep(each = age_groups_model[, .N]),
          age = 1:age_groups_model[, .N] %>% rep(length(compartments)),
          value = switch(arm_name,
                         "unvaccinated" = rep(1, length(compartments) * age_groups_model[, .N]) - {
                           if(length(x[["arms"]]) == 0) 0 else
                             x[["arms"]] %>% sapply(function(x) x[["catchup_coverage"]]) %>%
                             rowSums %>% rep(length(compartments))},
                         x[["arms"]][[arm_name]][["catchup_coverage"]]))
      }) %>% rbindlist})
  for(cl in 1:length(param.list$trial_arms)){
    clname = names(param.list$trial_arms)[cl]
    cluster_arms[[cl]][, cluster := factor(clname, names(param.list$trial_arms))]
  }
  cluster_arms = cluster_arms %>% rbindlist
  cluster_arms = cluster_arms[, c("cluster", "dose", "compartment", "age", "value")] %>% setorder(cluster, dose, compartment, age)
  model_output[, c("cluster", "dose", "compartment") := .(factor(cluster, levels(cluster_arms$cluster)),
                                                          factor(dose, levels(cluster_arms$dose)),
                                                          factor(compartment, levels(cluster_arms$compartment)))]
  cluster_arms = cluster_arms %>% merge(model_output[,-"dose"], by = c("cluster", "compartment", "age")) %>%
    .[, state := value.x * value.y] %>% .[, -c("value.x", "value.y")]
  
  cluster_arms = cluster_arms[, c("cluster", "dose", "compartment", "age", "state")] %>% setorder(cluster, dose, compartment, age)
  return(cluster_arms)
}

adaptiveMCMC = function(
  calculateLL, model_params, updateParameters, model_params_prior, model_params_prior_covmat = NULL,
  mcmc_steps = 500000, mcmc_adapt_size_start = 1000, mcmc_adapt_size_cooling = 0.999, mcmc_adapt_shape_start = 1000,
  mcmc_adapt_shape_stop = 1000, mcmc_scaling_factor_max = 50,
  output_file = sprintf("%s/MCMC_with6A_1.txt", OUTPUT_FOLDER),
  model_output_folder = sprintf("%s/%s", OUTPUT_FOLDER, "model_output")
){
  dir.create(model_output_folder)
  
  model_params_current = model_params_prior[, c("variable", "mean")] %>% (function(x){z=x[, mean]; names(z)=x[, variable]; return(z)})
  model_params = updateParameters(
    model_params_current %>% t %>% as.data.table %>% melt(id.vars=character(0)),
    model_params)
  
  #Placeholders for return values
  mcmc_adapt_shape_start_at = NA_integer_
  mcmc_adapt_shape_stop_at = NA_integer_
  
  mcmc_out_settings = list(
    list(param = "mcmc_steps", value = mcmc_steps),
    list(param = "mcmc_adapt_size_start", value = mcmc_adapt_size_start),
    list(param = "mcmc_adapt_size_cooling", value = mcmc_adapt_size_cooling),
    list(param = "mcmc_adapt_shape_start", value = mcmc_adapt_shape_start),
    list(param = "mcmc_adapt_shape_start_at", value = NA_integer_),
    list(param = "mcmc_adapt_shape_stop", value = mcmc_adapt_shape_stop),
    list(param = "mcmc_adapt_shape_stop_at", value = NA_integer_),
    list(param = "mcmc_scaling_factor_max", value = mcmc_scaling_factor_max)) %>%
    lapply(as.data.table) %>% rbindlist
  fwrite(mcmc_out_settings, sprintf("%s/mcmc_out_settings.csv", OUTPUT_FOLDER))
  
  updateCovMat = function(model_params_covmat_empirical, model_params_mean, model_params_current, i){
    model_params_residual = model_params_current - model_params_mean
    
    model_params_covmat_empirical = (
      model_params_covmat_empirical * (i - 1) +
        (i - 1)/i * model_params_residual %*% t(model_params_residual)
    )/i
    
    model_params_mean = model_params_mean + model_params_residual/i
    
    return(list(model_params_covmat_empirical = model_params_covmat_empirical, model_params_mean = model_params_mean))
  }
  
  mcmc_params = matrix(NA, nrow = mcmc_steps, ncol = nrow(model_params_prior) + 3)
  colnames(mcmc_params) = c(model_params_prior[, variable], c("logprior", "loglikelihood", "logposterior"))
  
  mcmc_adapting_size = FALSE
  mcmc_adapting_shape = FALSE
  mcmc_adapt_shape_i = 0
  mcmc_acceptance_rate = 0
  mcmc_scaling_factor = 1
  proposed_accepted = 0
  
  model_params_covmat_proposed = model_params_prior_covmat
  if(is.null(model_params_covmat_proposed)){
    model_params_covmat_proposed = matrix(
      diag(model_params_prior[, sd]^2, nrow = nrow(model_params_prior)), nrow = nrow(model_params_prior),
      dimnames = list(model_params_prior[, variable], model_params_prior[, variable]))
  }
  
  model_params_covmat_proposed_init = model_params_covmat_proposed
  model_params_current = model_params_prior[, c("variable", "mean")] %>% (function(x){z=x[, mean]; names(z)=x[, variable]; return(z)})
  model_params_mean = model_params_current
  
  model_params = updateParameters(
    model_params_current %>% t %>% as.data.table %>% melt(id.vars=character(0)),
    model_params)
  
  #' Log-likelihood initial parameters
  model_output = calculateLL(model_params)
  
  priorLikelihood = function(theta, param_priors = priors, combine=TRUE){
    pl = numeric(length(theta))
    for(i in 1:length(theta)){
      pl[i] = param_priors$prior[[i]](theta[i])
    }
    if(combine) return(sum(pl)) else return(pl)
  }
  
  log_likelihood_current = model_output$log_ll_total
  log_prior_current = priorLikelihood(model_params_current)
  log_posterior_current = log_likelihood_current + log_prior_current
  
  #' Initial empirical covariance matrix is 0 everywhere
  model_params_covmat_empirical = model_params_covmat_proposed
  model_params_covmat_empirical[, ] = 0
  
  z = 0
  for(i in 1:mcmc_steps){
    
    if(i >= mcmc_adapt_size_start && mcmc_adapt_shape_i <= mcmc_adapt_shape_stop){
      if(!mcmc_adapting_shape && mcmc_adapt_shape_i == 0){
        if(!mcmc_adapting_size){
          mcmc_adapting_size = TRUE
          message(sprintf("%s (%s%%) - %s - Start adapting size of covariance matrix",
                          i, round(100*i/mcmc_steps, 1), Sys.time())) 
        }
        #' At every iteration, update initial covariance matrix with scaling factor, aim for acceptance rate of 0.234
        mcmc_scaling_multiplier = exp(mcmc_adapt_size_cooling^(i - mcmc_adapt_size_start) * (mcmc_acceptance_rate - 0.234))
        mcmc_scaling_factor = mcmc_scaling_factor * mcmc_scaling_multiplier
        mcmc_scaling_factor = min(c(mcmc_scaling_factor, mcmc_scaling_factor_max))
        model_params_covmat_proposed_new = mcmc_scaling_factor^2 * model_params_covmat_proposed_init
        if(!(any(diag(model_params_covmat_proposed_new) < .Machine$double.eps))){
          model_params_covmat_proposed = model_params_covmat_proposed_new
        }
      }
      
      if(mcmc_acceptance_rate*i >= mcmc_adapt_shape_start){
        if(mcmc_adapting_size){
          mcmc_adapting_size = FALSE
          mcmc_adapting_shape = TRUE
          
          mcmc_adapt_shape_start_at = i
          mcmc_out_settings[param == "mcmc_adapt_shape_start_at", value := mcmc_adapt_shape_start_at]
          fwrite(mcmc_out_settings, sprintf("%s/mcmc_out_settings.csv", OUTPUT_FOLDER))
          
          message(sprintf("%s (%s%%) - %s - Stop adapting size of covariance matrix",
                          i, round(100*i/mcmc_steps, 1), Sys.time()))
          message(sprintf("%s (%s%%) - %s - Start adapting shape of covariance matrix",
                          i, round(100*i/mcmc_steps, 1), Sys.time()))
          
          #WHY update again? Keep old scaling fator
          mcmc_scaling_factor = 2.38/sqrt(nrow(model_params_prior))
        }
        
        if(mcmc_adapting_shape){
          if(mcmc_adapt_shape_i < mcmc_adapt_shape_stop){
            #' Scale empirical covariance matrix
            model_params_covmat_proposed = mcmc_scaling_factor^2 * model_params_covmat_empirical
            mcmc_adapt_shape_i = mcmc_adapt_shape_i + mcmc_proposed_accepted  
          } else {
            mcmc_adapting_shape = FALSE
            message(sprintf("%s (%s%%) - %s - Stop adapting shape of covariance matrix",
                            i, round(100*i/mcmc_steps, 1), Sys.time()))
            
            mcmc_adapt_shape_stop_at = i
            mcmc_out_settings[param == "mcmc_adapt_shape_stop_at", value := mcmc_adapt_shape_stop_at]
            fwrite(mcmc_out_settings, sprintf("%s/mcmc_out_settings.csv", OUTPUT_FOLDER))
          } 
        }
      }
    }
    
    #' Ensure none of the variances are 0
    if(any(diag(model_params_covmat_proposed) < .Machine$double.eps)){
      print(model_params_covmat_proposed)
      stop("Non-positive definite covariance matrix")
    }
    
    #' Sample new set of parameters
    model_params_proposed = tmvtnorm::rtmvnorm(
      n = 1, mean = model_params_current, sigma = model_params_covmat_proposed, lower = model_params_prior[, min],
      upper = model_params_prior[, max]) %>% as.numeric %>% (function(x){names(x)=priors[, variable]; return(x)})
    
    model_params = updateParameters(
      model_params_proposed %>% t %>% as.data.table %>% melt(id.vars=character(0)),
      model_params)
    
    #' Run the model and calculate the Log-Likelihood
    model_output = calculateLL(model_params)
    log_likelihood_proposed = model_output$log_ll_total
    log_prior_proposed = priorLikelihood(model_params_proposed)
    log_posterior_proposed = log_likelihood_proposed + log_prior_proposed
    
    #' Compute the Metropolis-Hastings ratio (log-scale)
    log_acceptance_ratio = log_posterior_proposed - log_posterior_current
    
    #' Adjust the MH ratio to account for truncated mvnorm prior
    log_acceptance_ratio = log_acceptance_ratio + tmvtnorm::dtmvnorm(
      x = model_params_current, mean = model_params_proposed, sigma = model_params_covmat_proposed, lower = model_params_prior[, min],
      upper = model_params_prior[, max], log = TRUE)
    
    #' Check if proposal is accepted
    mcmc_proposed_accepted = log(runif(1)) < log_acceptance_ratio
    if(mcmc_proposed_accepted){
      model_params_current = model_params_proposed
      log_likelihood_current = log_likelihood_proposed
      log_prior_current = log_prior_proposed
      log_posterior_current = log_posterior_proposed
      
      #' Save model output
      res_postvacc_all = model_output$save_data
      #res_postvacc_all = res_postvacc_all %>% as.data.table %>% .[, -"output"] %>% reshapeModelOutput(model_params$params_vac)
      #res_postvacc_all = res_postvacc_all %>%
      #  merge(age_groups[, c("id", "name", "age_group_data_trial", "age_group_data_2006")], by.x="age", by.y="id") %>%
      #  merge(popsize_model, by="name") %>%
      #  .[, value := value * N] %>% .[!is.na(age_group_data_trial), .(value = sum(value)), by=c("cluster", "compartment", "age_group_data_trial", "time")] %>%
      #  .[, N := sum(value), by=c("cluster", "age_group_data_trial", "time")] %>% .[, prev := value / N] %>% .[, -c("N", "value")] %>%
      #  dcast(...~compartment, value.var="prev") %>% .[, VT := VT + B] %>% .[, -"B"] %>% melt(measure.vars = c("S", "VT", "NVT"), variable.name = "compartment", value.name = "modelled")
      
      qs::qsave(res_postvacc_all, sprintf("%s/out_%s.qs", model_output_folder, i))
    }
    
    mcmc_params[i, ] = c(model_params_current, log_prior_current, log_likelihood_current, log_posterior_current)
    
    #' Update acceptance rate
    if(i == 1)
      mcmc_acceptance_rate = mcmc_proposed_accepted
    else
      mcmc_acceptance_rate = mcmc_acceptance_rate + (mcmc_proposed_accepted - mcmc_acceptance_rate)/i
    
    #' Update empirical covariance matrix
    if(mcmc_adapt_shape_i <= mcmc_adapt_shape_stop){
      tmp = updateCovMat(model_params_covmat_empirical, model_params_mean, model_params_current, i)
      model_params_covmat_empirical = tmp$model_params_covmat_empirical
      model_params_mean = tmp$model_params_mean
    }
    
    if( (i %% 1000) == 0 )
      write.table(mcmc_params, output_file)
    
    if( (i %% 10) == 0 ){
      #message(sprintf("%s/%s (%s%%) - %s - LLp1: %s; LL: %s; a: %s (%s%%); s: %s",
      #                i, mcmc_steps, round(100*i/mcmc_steps, 1), Sys.time(),
      #                round(log_likelihood_proposed, 1), round(log_likelihood_current, 1),
      #                round(mcmc_acceptance_rate*i, 1), round(100*mcmc_acceptance_rate,1),
      #                round(mcmc_scaling_factor, 5)))
      
      if((z %% 15) == 0){
        message(paste0(sprintf("%15.15s", c("Time", "Iteration", "LLProp/Curr", "LPrProp/Curr", "LPoProp/Curr", "Accepted", "Adapt",
                                            colnames(mcmc_params)[-ncol(mcmc_params)])), collapse = " "))
      }
      message(paste0(c(format(Sys.time(), format="%H:%M:%S"), sprintf("%s (%s%%)", i, round(100*i/mcmc_steps,1)),
                       sprintf("%s/%s", round(log_likelihood_proposed, 1), round(log_likelihood_current, 1)),
                       sprintf("%s/%s", round(log_prior_proposed, 1), round(log_prior_current, 1)),
                       sprintf("%s/%s", round(log_posterior_proposed, 1), round(log_posterior_current, 1)),
                       sprintf("%s (%s%%)", round(mcmc_acceptance_rate*i, 1), round(100*mcmc_acceptance_rate,1)),
                       round(mcmc_scaling_factor, 5), sprintf("%1.4e", (mcmc_params[i, -ncol(mcmc_params)]))) %>%
                       sprintf("%15.15s", .), collapse = " "))
      z = z + 1
    }
  }
  
  return(list(mcmc_params = mcmc_params, acceptance_rate = mcmc_acceptance_rate, covmat = model_params_covmat_proposed,
              times = c(adapt_size_start_at = mcmc_adapt_shape_start, adapt_shape_start_at = mcmc_adapt_shape_start_at,
                        adapt_shape_stop_at = mcmc_adapt_shape_stop_at)))
}

#' Process the contact matrix to have the same age groups as the model, and return as a matrix with
#' contactor age-groups in columns and contactee age-groups in rows
adjustContactMatrixAgeGroups = function(age_groups_model, contact_matrix_data, contact_matrix_data_agegroups, population_size_model){
  contact_matrix = age_groups_model %>%
    combineAgeBreaks(contact_matrix_data_agegroups[, -"name"] %>%
                       cbind(contact_matrix_data %>% dcast(contactee_age_group ~ contactor_age_group)),
                     method = "sum", value.var = contact_matrix_data_agegroups[, name]) %>%
    .[, -c("from", "to")] %>%
    .[, contactee_age_group := name] %>% .[, -"name"] %>%
    melt(id.vars = "contactee_age_group", variable.name = "contactor_age_group") %>%
    dcast(contactor_age_group ~ contactee_age_group) %>%
    cbind(contact_matrix_data_agegroups[, -"name"]) %>%
    combineAgeBreaks(x = age_groups_model, y = ., method = "mean", value.var = age_groups_model[, name]) %>%
    .[, -c("from", "to", "name")] %>%
    as.matrix(rownames = age_groups_model[, name]) %>%
    t()
  
  #'Make symmetrical for this population
  #' multiply columns with total number of contactors by age
  #' calculate average total number of contacts
  #' divide columns by total number of contactors by age
  total_matrix = contact_matrix %*% diag(population_size_model[, value])
  total_matrix = (total_matrix + t(total_matrix))/2
  contact_matrix = total_matrix %*% diag(1/population_size_model[, value])
  colnames(contact_matrix) = age_groups_model[, name]
  
  return(contact_matrix)
}

#' Get dominant eigenvalue
dominantEigenValue = function(contact_matrix){
  return(Re(eigen(contact_matrix, only.values = TRUE)[["values"]][1]))
}

sampleCaseCarrierRatio = function(age_groups_model, case_carrier_data){
  i = sample(case_carrier_data[, iter], 1)
  
  case_carrier_ratio_model = age_groups_model %>%
    combineAgeBreaks(case_carrier_data[iter == i] %>% dcast(from+to+name ~ st),
                     value.var = c("NVT", "VT"))
  
  return(case_carrier_ratio_model)
}

#' Updates rates by timestep
#check if any > 1
#for each arm
adjustForTimeStep = function(value, MODEL_TIMESTEP.=MODEL_TIMESTEP){
  MODEL_TIMESTEP = MODEL_TIMESTEP.
  multTimestep = function(x, checkMatrix = FALSE, popsize = NULL){
    x = x * MODEL_TIMESTEP
    if(!is.matrix(x) & any(x >= 1)) stop("Rate per timestep > 1, use smaller timestep.")
    
    #' Check if contact rate > 1
    #if(is.matrix(x) & checkMatrix){
    #  if(!is.null(popsize)){
    #    if(any(t(x) %*% diag(1/popsize) > 1)) stop("Rate per timestep > 1, use smaller timestep.")
    #  } else {
    #    if(any(x > 1)) stop("Rate per timestep > 1, use smaller timestep.") 
    #  }
    #}
          
    return(x)
  }
  
  #' If adjusting the contact matrix only, this line is used
  if(is.matrix(value)) return(multTimestep(value))
  
  model_params = value
  model_params$params_unvac$clearVT %<>% multTimestep()
  model_params$params_unvac$clearNVT %<>% multTimestep()
  model_params$params_unvac$ageout %<>% multTimestep()
  model_params$params_unvac$migration[which(model_params$params_unvac$migration != -1)] %<>% multTimestep()
  model_params$params_unvac$trial_arms %<>% lapply(function(x){
    x$parameters$betaVT %<>% multTimestep(TRUE, x$parameters$N)
    x$parameters$betaNVT %<>% multTimestep(TRUE, x$parameters$N)
    x$arms %<>% lapply(function(z){
      z$waning %<>% multTimestep()
      return(z)})
    return(x)})
  
  model_params$params_vac$clearVT %<>% multTimestep()
  model_params$params_vac$clearNVT %<>% multTimestep()
  model_params$params_vac$ageout %<>% multTimestep()
  model_params$params_vac$migration[which(model_params$params_vac$migration != -1)] %<>% multTimestep()
  model_params$params_vac$trial_arms %<>% lapply(function(x){
    x$parameters$betaVT %<>% multTimestep(TRUE, x$parameters$N)
    x$parameters$betaNVT %<>% multTimestep(TRUE, x$parameters$N)
    x$arms %<>% lapply(function(z){
      z$waning %<>% multTimestep()
      return(z)})
    return(x)})
  
  return(model_params)
}

getVaccineCoverage = function(age_groups_model, age_breaks, coverage, return_table = FALSE){
  if(length(coverage) == 1 & length(age_breaks) > 2){
    warning("Applying same coverage to all age breaks provided. Consider widening the age-break provided to a single group.")
  } else if(length(age_breaks) != 1 & (length(age_breaks) - 1 != length(coverage))){
    stop("Length of coverage does not match age-breaks, and more than one coverage is specified.")
  }
  
  if(length(coverage) == 1 & length(age_breaks) == 1){
    coverage_table = copy(age_groups_model) %>% .[, value := 0]
    if(coverage_table[ageeq(from, age_breaks), .N] == 0 & coverage != 0){
      stop(sprintf("Age break %s %s does not exit in age_groups_model", age_breaks, attr(age_breaks, "units")$numerator))
    }
    coverage_table[ageeq(from, age_breaks), value := coverage]
  } else {
    coverage_by_age = setAgeBreaks(age_breaks) %>% .[, value := 0]
    
    for(a in 1:(length(age_breaks)-1)){
      if(length(coverage) == 1){
        coverage_by_age[ageeq(from, age_breaks[a]), value := coverage]
      } else {
        coverage_by_age[ageeq(from, age_breaks[a]), value := coverage[a]]  
      }
    }
    
    coverage_table = age_groups_model %>% combineAgeBreaks(coverage_by_age)
  }
  
  if(return_table) return(coverage_by_age)
  else return(coverage_table[, value])
}

#' Create expected contact matrix for a given population size
#' - assuming contacts are made completely at random
#' - and assuming everyone makes exactly one contact with every other individual
createExpectedMatrix = function(popsize){
    expected = matrix(rep(popsize/sum(popsize), each=length(popsize)), length(popsize))
    return(expected)
}

#' Create prior for BayesianTools
createBTPrior = function(priors){
  BayesianTools::createPrior(
    density = function(par){
      seq_len(length(par)) %>%
        sapply(function(p, par, priors) (priors[p, density][[1]])(par[p]),
               par=par, priors = priors) %>% sum},
    sampler = function(n=1){
      values = seq_len(priors[, .N]) %>%
        sapply(function(p, n, priors) (priors[p, sampler[[1]]])(n),
               n=n, priors=priors)
      if(n == 1) values = t(values)
      colnames(values) = priors[, variable]
      return(values)},
    lower = priors[, min],
    upper = priors[, max],
    best = NULL)
}
