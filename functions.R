#' tryCatch to wrap around the deSolve functions
#' - there may be runs where deSolve returns an error, in which case we do not want the MCMC to fail
tryCatchWE = function(expr){
  status = 0
  wHandler = function(w) {status <<- 1; return(paste0(w, collapse = "; "))}
  eHandler = function(e) {status <<- 1; return(paste0(e, collapse = "; "))}
  value = list(value = tryCatch(expr, error = eHandler, warning = wHandler), status = status)
  return(value)
}

#' Should work to compile on all platforms
#' Mostly based on inline::cxxfunction
#' Compiled model is slightly slower compared to manual compilation
#' - TODO: figure out why
#' - system('g++ -c -o ./model/build/Fit_by_arm_quick_v5.o ./model/Fit_by_arm_quick_v5.cpp -L/usr/lib/R/lib -lR -std=gnu++11 -shared -L/usr/lib/R/lib -Wl,-Bsymbolic-functions,-z,relro -I"/usr/share/R/include" -I"/home/lsh1604011/workspace/Dosing_schedule_and_fitting" -I"/home/lsh1604011/R/x86_64-pc-linux-gnu-library/4.0/Rcpp/include" -I"/home/lsh1604011/R/x86_64-pc-linux-gnu-library/4.0/RcppArmadillo/include" -fPIC -DNDEBUG -fopenmp -g -O3 -fdebug-prefix-map=/build/r-base-8T8CYO/r-base-4.0.3=. -fstack-protector-strong -Wformat -Werror=format-security -Wdate-time -D_FORTIFY_SOURCE=2 -g')
#' - system('g++ -o ./model/build/Fit_by_arm_quick_v5.so ./model/build/Fit_by_arm_quick_v5.o -L/usr/lib/R/lib -lR -std=gnu++11 -shared -L/usr/lib/R/lib -Wl,-Bsymbolic-functions,-z,relro -I"/usr/share/R/include" -I"/home/lsh1604011/workspace/espicc_model" -I"/home/lsh1604011/R/x86_64-pc-linux-gnu-library/4.0/Rcpp/include" -fPIC -fopenmp -llapack -lblas -lgfortran -lm -lquadmath -L/usr/lib/R/lib -lR')
compileModel = function(cpp_file, shrd_lib_loc, shrd_lib_name = ""){
  if(!file.exists(cpp_file))
    stop("Cpp file does not exist")
  if(shrd_lib_name == "") shrd_lib_name = gsub(".cpp", .Platform$dynlib.ext, basename(cpp_file))
  if(file.exists(sprintf("%s/%s", shrd_lib_loc, shrd_lib_name)))
    file.remove(sprintf("%s/%s", shrd_lib_loc, shrd_lib_name))
  
  cmd = paste0(R.home(component = "bin"), "/R")
  
  #' Need to include these directories  
  paths = sapply(c("Rcpp", "RcppArmadillo"), find.package)
  flag = paste(paste0("-I\"", paths, "/include\""), collapse = " ")
  
  #' Need to set additional compiler flags
  do.call(Sys.setenv, inline::getPlugin("RcppArmadillo")$env)
  #change this to make omp optional
  #could add -g flag to ease debugging
  Sys.setenv(CLINK_CPPFLAGS = paste0(c(flag, "-lomp -DMP_ENABLED -fopenmp"), collapse = " "))
  Sys.setenv(PKG_CXXFLAGS="-std=c++14")
  
  #' Compile model
  system2(cmd, args = paste(" CMD SHLIB -o", sprintf("%s/%s", shrd_lib_loc, shrd_lib_name), cpp_file))
  if(!file.exists(sprintf("%s/%s", shrd_lib_loc, shrd_lib_name)))
    stop("Something went wrong")
  
  #' Remove object file
  if(file.exists(gsub(".cpp", ".o", cpp_file)))
    file.remove(gsub(".cpp", ".o", cpp_file))
}

reshapeModelOutput2 = function(result, model_params){
  result = result %>% as.data.table()
  incidence = !"output" %in% colnames(result)
  if(!incidence) result = result[, -"output"]
  result = result %>% melt(id.vars="time", variable.name = "id")
  
  #' create data table to match to molten data
  columns_prevalence = model_params$trial_arms %>%
    (function(populations){
      lapply(names(populations), function(pop, populations){
        v_strata = populations[[pop]][["arms"]]
        data.table(
          population = pop %>% rep(length(v_strata) * length(compartments_prevalence) * age_groups_model[, .N]),
          vaccination_group = names(v_strata) %>% rep(each = length(compartments_prevalence) * age_groups_model[, .N]),
          compartment = compartments_prevalence %>% rep(each=age_groups_model[, .N]) %>% rep(length(v_strata)),
          age = seq_len(age_groups_model[, .N]) %>% rep(length(compartments_prevalence)) %>% rep(length(v_strata))
        )}, populations)}) %>% rbindlist()
  columns_prevalence[, c("id", "outcome") := .(as.character(.I), "prevalence")]
  columns_match = columns_prevalence
  
  if(incidence){
    columns_incidence = model_params$trial_arms %>%
      (function(populations){
        lapply(names(populations), function(pop, populations){
          v_strata = populations[[pop]][["arms"]]
          data.table(
            population = pop %>% rep(length(v_strata) * length(compartments_incidence) * age_groups_model[, .N]),
            vaccination_group = names(v_strata) %>% rep(each = length(compartments_incidence) * age_groups_model[, .N]),
            compartment = compartments_incidence %>% rep(each=age_groups_model[, .N]) %>% rep(length(v_strata)),
            age = seq_len(age_groups_model[, .N]) %>% rep(length(compartments_incidence)) %>% rep(length(v_strata))
          )}, populations)}) %>% rbindlist()
    columns_incidence[, c("id", "outcome") := .(paste0("inc_", .I), "incidence")]
    columns_match = rbind(columns_match, columns_incidence)
  }
  
  columns_match[, population := factor(population, names(model_params$trial_arms))]
  columns_match[, vaccination_group := factor(vaccination_group, model_params$trial_arms %>% sapply(function(p) names(p[["arms"]])) %>% unlist() %>% as.vector() %>% unique())]
  if(incidence){
    columns_match[, outcome := factor(outcome, c("prevalence", "incidence"))]
    columns_match[, compartment := factor(compartment, c(compartments_prevalence, compartments_incidence))]
  } else {
    columns_match[, outcome := factor(outcome, c("prevalence"))]
    columns_match[, compartment := factor(compartment, compartments_prevalence)]
  }
  
  columns_match = columns_match %>% merge(age_groups_model[, c("age", "name")] %>% setNames(c("age", "age_group")), by="age")
  
  result = result %>% merge(columns_match, by="id") %>%
    .[, c("outcome", "population", "vaccination_group", "compartment", "age_group", "age", "time", "value")]
  
  setorder(result, outcome, population, vaccination_group, compartment, age_group, age, time)
  
  return(result)
}

eqStatesVaccinate2 = function(model_output, model_params, pop_unvacc = NULL){
  if(length(model_output[, unique(time)]) != 1) stop("Table model_output needs to be for a single timestep")
  if(is.null(pop_unvacc)){
    if(length(model_output[, unique(population)]) == 1){
      pop_unvacc = rep(model_output[, unique(population)], length(model_params$trial_arms))
    } else if(all(names(model_params$trial_arms) %in% model_output[, unique(population)])){
      pop_unvacc = names(model_params$trial_arms)
    } else {
      stop("eqStatesVaccinate2: method not yet implemented")
    }
  }
  names(pop_unvacc) = names(model_params$trial_arms)
  
  model_output = model_output[outcome == "prevalence"]
  model_output %>% setorder(population, vaccination_group, compartment, age)
  
  model_input = model_params$trial_arms %>%
    (function(populations){
      lapply(names(populations), function(pop, populations){
        v_strata = names(populations[[pop]][["arms"]])
        v_strata %>% lapply(function(vstrat, pop){
          copy(model_output[population == pop_unvacc[[pop]]]) %>%
            .[, c("population", "vaccination_group", "value") := .(pop, vstrat, ifelse(vstrat == "unvaccinated", value, 0)), by=c("compartment", "age")] %>% .[]
        }, pop) %>% rbindlist()}, populations) %>% rbindlist})
  
  model_input[, population := factor(population, names(model_params$trial_arms))]
  model_input[, vaccination_group := factor(vaccination_group, model_params$trial_arms %>% sapply(function(p) names(p[["arms"]])) %>% unlist() %>% as.vector() %>% unique())]
  model_input[, outcome := factor(outcome, c("prevalence"))]
  model_input[, compartment := factor(compartment, compartments_prevalence)]
  
  setorder(model_input, outcome, population, vaccination_group, compartment, age_group, age, time)
  
  return(model_input %>% .[])
}

#' Process the contact matrix to have the same age groups as the model, and return as a matrix with
#' contactor age-groups in columns and contactee age-groups in rows
adjustContactMatrixAgeGroups = function(age_groups_model, contact_matrix_data, contact_matrix_data_agegroups, population_size_model){
  contact_matrix = age_groups_model %>% .[, -"age"] %>%
    combineAgeBreaks(contact_matrix_data_agegroups[, -"name"] %>%
                       cbind(contact_matrix_data %>% dcast(contactee_age_group ~ contactor_age_group)),
                     method = "sum", value.var = contact_matrix_data_agegroups[, name]) %>%
    .[, -c("from", "to")] %>%
    .[, contactee_age_group := name] %>% .[, -c("name")] %>%
    melt(id.vars = "contactee_age_group", variable.name = "contactor_age_group") %>%
    dcast(contactor_age_group ~ contactee_age_group) %>%
    cbind(contact_matrix_data_agegroups[, -"name"]) %>%
    combineAgeBreaks(x = age_groups_model[, -"age"], y = ., method = "mean", value.var = age_groups_model[, name]) %>%
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

#sampleCaseCarrierRatio = function(age_groups_model, case_carrier_data){
#  i = sample(case_carrier_data[, iter], 1)
#  
#  case_carrier_ratio_model = age_groups_model %>%
#    combineAgeBreaks(case_carrier_data[iter == i] %>% dcast(from+to+name ~ st),
#                     value.var = c("NVT", "VT"))
#  
#  return(case_carrier_ratio_model)
#}

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

createBetaPriorBT = function(name, min = 0, max = 1, plotmin = NULL, plotmax = NULL, shape1, shape2){
  range = max - min
  if(is.null(plotmin)) plotmin = min
  if(is.null(plotmax)) plotmax = max
  data.table(variable = name, min = min, max = max, plotmin = plotmin, plotmax = plotmax,
             density = function(x, uselog=TRUE) dbeta((x - min)/range, shape1 = shape1, shape2 = shape2, log=uselog),
             sampler = function(n) rbeta(n, shape1 = shape1, shape2 = shape2) * range + min)
}

createUnifPriorBT = function(name, min = 0, max = 1, plotmin = NULL, plotmax = NULL){
  if(is.null(plotmin)) plotmin = min
  if(is.null(plotmax)) plotmax = max
  data.table(variable = name, min = min, max = max, plotmin = plotmin, plotmax = plotmax,
             density = function(x, uselog=TRUE) dunif(x, min = min, max = max, log=uselog),
             sampler = function(n) runif(n, min, max))
}

createLogNormPriorBT = function(name, min = 0, max = 1, plotmin = NULL, plotmax = NULL, meanlog, sdlog, flippedx = FALSE){
  if(is.null(plotmin)) plotmin = min
  if(is.null(plotmax)) plotmax = max
  data.table(variable = name, min = min, max = max, plotmin = plotmin, plotmax = plotmax,
             density = function(x, uselog=TRUE) dlnorm(ifelse(flippedx, 1 - x, x), meanlog = meanlog, sdlog = sdlog),
             sampler = function(n){
               in_range = FALSE
               while(!in_range){
                 val = rlnorm(n, meanlog = meanlog, sdlog = sdlog)
                 if(flippedx) val = 1 - val
                 in_range = (val >= min & val <= max)
               }
               return(val)
             })
}

uniqueSharedObject = function(){
  #' copy compiled model
  main_file_path = sprintf("./model/build/%s%s", MODEL_NAME, .Platform$dynlib.ext)
  new_file_name = sprintf("%s_%s", MODEL_NAME, Sys.getpid())
  new_file_path = sprintf("./model/build/%s%s", new_file_name, .Platform$dynlib.ext)
  
  if(!file.exists(new_file_path)){
    if(Sys.info()["sysname"] == "Windows"){
      #' copy compiled model
      file.copy(main_file_path, new_file_path)
    } else {
      #' recompile model
      #' - need to first copy cpp file so the .o file will be unique
      #' - need to update this, first compile .o file, then only need to be linked in unique .so
      file.copy(sprintf("%s/model/%s.cpp", PCVM_FOLDER, MODEL_NAME), sprintf("%s/model/%s.cpp", PCVM_FOLDER, new_file_name))
      compileModel(sprintf("%s/model/%s.cpp", PCVM_FOLDER, new_file_name), "./model/build/",
                   sprintf("%s%s", new_file_name, .Platform$dynlib.ext))
      file.remove(sprintf("%s/model/%s.cpp", PCVM_FOLDER, new_file_name))
    }
  }
  if(!is.loaded("derivs", new_file_name)) dyn.load(new_file_path)
  if(!is.loaded("derivs", new_file_name)) stop("MetaVax is not loaded")
  
  return(new_file_name)
}

#' Function that will be used in MCMC algorithm
runModel = function(initial_state, model_params, steady_state = FALSE, times = c(0, 1), hmin = 0, hmax = NULL, rtol = 1e-06, atol = 1e-06, incidence = FALSE, parallel = FALSE){
  nout_incidence = model_params$trial_arms %>% sapply(function(x) length(x[["arms"]])) %>% sum() * age_groups_model[, .N] * length(compartments_incidence)
  
  if(steady_state){
    #if(incidence) warning("currently not running incidence in steady state")
    result = runsteady(
      y=initial_state, func = "derivs",
      initpar = model_params, dllname = {if(parallel) uniqueSharedObject() else MODEL_NAME},
      nout = ifelse(incidence, nout_incidence, 1), outnames = {if(incidence) paste0("inc_", seq_len(nout_incidence)) else "output"},
      initfunc = "rt_initmod", jactype = "fullint",
      hmin = hmin, hmax = hmax, atol = atol, rtol = rtol) %>% tryCatchWE()
  } else {
    campaign_times = sapply(model_params$trial_arms,
                            function(cluster){
                              sapply(cluster$arms,
                                     function(arm) sapply(arm$coverage_c, "[[", "time")) %>%
                                unlist %>% unique %>% sort}) %>% unlist() %>% sort() %>% unique()
    if(length(campaign_times) == 0){
      result = lsode(
        y=initial_state, times=times, func = "derivs",
        parms = model_params, dllname = {if(parallel) uniqueSharedObject() else MODEL_NAME},
        nout = ifelse(incidence, nout_incidence, 1), outnames = {if(incidence) paste0("inc_", seq_len(nout_incidence)) else "output"},
        initfunc = "initmod", jactype = "fullint",
        hmin = hmin, hmax = hmax, atol = atol, rtol = rtol) %>% tryCatchWE
    } else {
      events = list(func="vaccineCampaignEvent",
                    time=campaign_times)
      times = sort(unique(c(times, campaign_times)))
      result = lsode(
        y=initial_state, times=times, func = "derivs",
        parms = model_params, dllname = {if(parallel) uniqueSharedObject() else MODEL_NAME},
        nout = ifelse(incidence, nout_incidence, 1), outnames = {if(incidence) paste0("inc_", seq_len(nout_incidence)) else "output"},
        initfunc = "initmod", jactype = "fullint",
        events = events, hmin = hmin, hmax = hmax, atol = atol, rtol = rtol) %>% tryCatchWE  
    }
  }
  
  if(result$status == 0){
    if(steady_state){
      if(incidence){
        result$value = data.table(time = attr(result$value, "time")) %>%
          cbind(data.table(i = c(1:length(result$value$y), length(result$value$y) + 1:length(result$value$var)), val = c(result$value$y, result$value$var))) %>%
          dcast(time~i, value.var="val") %>%
          setNames(c("time", 1:length(result$value$y), paste0("inc_", seq_len(nout_incidence)))) %>%
          reshapeModelOutput2(model_params)  
      } else {
        result$value = data.table(time = attr(result$value, "time")) %>%
          cbind(data.table(i = c(1:length(result$value$y)), val = result$value$y)) %>%
          dcast(time~i, value.var="val") %>%
          .[, output := 1] %>%
          reshapeModelOutput2(model_params)
      }
      
    } else{
      result$value = result$value %>%
        reshapeModelOutput2(model_params) 
    }
  }
  
  return(result)
}

#' coverage_to can be provided as the name of another vaccination_group to improve readability, but needs to be provided as an index to
#' the model. This function matches names provided in coverage_to, to the correct index 
renameCoverageTo = function(model_populations){
  lapply(model_populations, function(population){
    arm_names = names(population$arms)
    population$arms = lapply(arm_names, function(name, arms){
      arm = arms[[name]]
      arm$coverage_r = arm$coverage_r %>% lapply(function(x, arms, name){
        if(!is.null(x$coverage_to)){
          if(is.numeric(x$coverage_to)){
            message("coverage_to is numeric, assuming correct indices are already provided")
          } else {
            x$coverage_to = x$coverage_to %>% sapply(function(z, arms){
              #cpp index starts at 0
              which(arms == z) - 1}, arms)
          }
        } else {
          #move to the next vaccination_group, if not provided
          x$coverage_to = rep(which(arms == name), age_groups_model[, .N])
        }
        
        return(x)
      }, names(arms), name)
      
      arm$coverage_c = arm$coverage_c %>% lapply(function(x, arms, name){
        if(!is.null(x$coverage_to)){
          if(is.numeric(x$coverage_to)){
            message("coverage_to is numeric, assuming correct indices are already provided")
          } else {
            x$coverage_to = x$coverage_to %>% sapply(function(z, arms){
              #cpp index starts at 0
              which(arms == z) - 1}, arms)
          }
        } else {
          x$coverage_to = rep(which(arms == name), age_groups_model[, .N])
        }
        
        return(x)
      }, names(arms), name)
      
      return(arm)
    }, population$arms)
    
    names(population$arms) = arm_names
    
    return(population)
  }) 
}

#' alternative summary for BayesianTools
altSummary = function(out){
  
  ##update Parameters
  #for(j in 1:length(out$chain)){
  #  x = out$chain[[j]]  
  #  x = apply(x[, 1:priors[, .N]], 1, function(z){
  #    names(z) = priors$variable
  #    return(as.vector(updateDepParameters(z)))
  #  }) %>% t()
  #  out$chain[[j]][, 1:priors[, .N]] = x
  #}
  
  #try(DInf <- DIC(sampler), silent = TRUE)
  sampler <- out
  MAPvals <- round(MAP(sampler)$parametersMAP, 3)
  psf <- FALSE
  mcmcsampler <- sampler$settings$sampler
  runtime <- sampler$settings$runtime[3]
  correlations <- round(cor(getSample(sampler)), 3)
  chain <- getSample(sampler, parametersOnly = T, coda = T)
  if ("mcmc.list" %in% class(chain)) {
    psf <- TRUE
    nrChain <- length(chain)
    nrIter <- nrow(chain[[1]])
    conv <- tryCatchWE(round(coda::gelman.diag(chain)$mpsrf, 3))# ifelse(chain$setup$numPars > 1, round(coda::gelman.diag(chain)$mpsrf, 3), round(coda::gelman.diag(chain)$mpsrf, 3)$psrf[1])
    if(conv$status == 0) conv = conv$value else conv = 9999
    npar <- sampler$setup$numPars
    lowerq <- upperq <- numeric(npar)
    medi <- numeric(npar)
    parnames <- colnames(chain[[1]])
    for (i in 1:npar) {
      if (nchar(parnames[i]) > 8) 
        parnames[i] <- paste(substring(parnames[i], 1, 
                                       6), "...", sep = "")
    }
    for (i in 1:npar) {
      tmp <- unlist(chain[, i])
      tmp <- quantile(tmp, probs = c(0.025, 0.5, 0.975))
      lowerq[i] <- round(tmp[1], 3)
      medi[i] <- round(tmp[2], 3)
      upperq[i] <- round(tmp[3], 3)
    }
  } else {
    nrChain <- 1
    nrIter <- nrow(chain)
    npar <- sampler$setup$numPars
    conv <- "Only one chain; convergence cannot be determined!"
    medi <- numeric(npar)
    lowerq <- upperq <- numeric(npar)
    parnames <- colnames(chain)
    for (i in 1:npar) {
      tmp <- quantile(chain[, i], probs = c(0.025, 0.5, 
                                            0.975))
      lowerq[i] <- round(tmp[1], 3)
      medi[i] <- round(tmp[2], 3)
      upperq[i] <- round(tmp[3], 3)
    }
  }
  parOutDF <- cbind(MAPvals, lowerq, medi, upperq)
  colnames(parOutDF) <- c("MAP", "2.5%", "median", "97.5%")
  if (psf == TRUE) {
    psf <- round(gelmanDiagnostics(sampler)$psrf[, 1], 3)
    parOutDF <- cbind(psf, parOutDF)
  }
  
  parnames = priors$variable
  row.names(parOutDF) <- parnames
  cat(rep("#", 25), "\n")
  cat("## MCMC chain summary ##", "\n")
  cat(rep("#", 25), "\n", "\n")
  #cat("# MCMC sampler: ", mcmcsampler, "\n")
  #cat("# Nr. Chains: ", nrChain, "\n")
  cat("# Iterations per chain: ", nrIter, "\n")
  cat("# Acceptance rate: ", 1 - ifelse(out$setup$numPars == 
                                          1 & class(chain) == "mcmc.list", round(mean(sapply(chain, 
                                                                                             coda::rejectionRate)), 3), round(mean(coda::rejectionRate(chain)), 
                                                                                                                              3)), "\n")
  cat("# Effective sample size: ", ifelse(sampler$setup$numPars == 
                                            1, round(coda::effectiveSize(chain), 0), round(mean(coda::effectiveSize(chain)), 
                                                                                           0)), "\n")
  cat("# Runtime: ", runtime, " sec.", "\n", "\n")
  cat("# Parameters (adjusted)\n")
  print(parOutDF)
  cat("\n")
  #try(cat("## DIC: ", round(DInf$DIC, 3), "\n"), silent = TRUE)
  cat("## Convergence", "\n", "Gelman Rubin multivariate psrf: ", 
      conv, "\n", "\n")
}

#' Increases or decreases the brightness of a color by a percentage of the current brightness.
#' adjust_by: A number between -1 and 1. E.g. 0.3 = 30% lighter; -0.4 = 40% darker.
#' See https://stackoverflow.com/a/54393956
adjustColBrightness = function(hex_code, adjust_by=c(-1, 0, 1)){
  if(any(adjust_by < -1) | any(adjust_by > 1)) stop("adjust_by needs to be between -1 and 1")
  
  if (nchar(hex_code) == 4) {
    hex_code = gsub("#", "", hex_code) %>%
      strsplit("") %>% .[[1]] %>% rep(each=2) %>% paste0("#", .)
  }
  
  rgb = col2rgb(hex_code)
  
  col_matrix = matrix(rep(rgb, length(adjust_by)), 3)
  adjustable_limit = col_matrix
  adjustable_limit[, which(adjust_by > 0)] = (255 - as.matrix(adjustable_limit[, which(adjust_by > 0)]))
  adjustable_limit = ceiling(adjustable_limit * matrix(rep(adjust_by, each=3), 3))
  col_matrix = col_matrix + adjustable_limit
  
  apply(col_matrix, MARGIN = 2, FUN = function(x) rgb(x[1]/255, x[2]/255, x[3]/255, 1))
}

# combine multiple BayesianTools output files from the same folder
combineOutFiles = function(OUTPUT_FOLDER){
  out_files = list.files(OUTPUT_FOLDER) %>%
    subset(grepl(pattern = "out(\\w)+.RDS", x = .))
  chains = sapply(strsplit(out_files, "_"), "[[", 2) %>% unique %>% as.numeric()
  
  # read max iteration for each chain
  out = lapply(chains, function(chain){
    #get max iteration for each chain
    max_iter = out_files %>%
      subset(grepl(pattern = sprintf("out_%s", chain), x=.)) %>%
      sapply(function(x){
        gsub("^[^0-9]*[0-9]+[^0-9]*", "", x) %>%
          gsub("^([0-9]+).*", "\\1", .) %>%
          as.numeric
      }) %>%
      max()
    
    out = readRDS(sprintf("%s/out_%s_%s.RDS", OUTPUT_FOLDER, chain, max_iter))  
    return(out)})
  
  #process mcmclist to use in analyses
  if(length(out) == 1) {
    out = out[[1]]
  } else {
    #make sure all chains are of the same length
    min_it = min(sapply(out, function(x) min(dim(x$chain[[1]])[1])))
    for(i in 1:length(out)){
      out[[i]]$chain = out[[i]]$chain %>% lapply(function(x){
        x = x[1:min_it, ]
        class(x) = "mcmc"
        return(x)})
      
      class(out[[i]]$chain) = "mcmc.list"
    }
    
    class(out) = c("mcmcSamplerList", "bayesianOutput")
  }
}

createTracePlot = function(posterior){
  posterior %>%
    .[, -c("LP", "LL", "LPr")] %>%
    .[, i := 1:.N, by="chain"] %>%
    melt(id.vars=c("chain", "i")) %>%
    ggplot(aes(x=i, y=value, colour=as.factor(chain)))+
    facet_grid(factor(variable, priors$variable)~., scales = "free")+
    geom_line(alpha=0.75)+
    theme_bw()+
    labs(x="Iteration", y="Value", colour="Chain")+
    guides(colour = guide_legend(override.aes = list(size = 5))) 
}

createPriorPosteriorPlot = function(priors, posterior){
  posterior_long = posterior %>%
    .[, -c("LP", "LL", "LPr")] %>%
    .[, i := 1:.N, by="chain"] %>%
    melt(id.vars=c("chain", "i")) %>%
    .[, type := "posterior"]
  
  prior_evaluated = lapply(1:nrow(priors), function(i){
    eval_at = seq(from=priors[i, plotmin] + 1e-6, to=priors[i, plotmax], length=1000)
    data.table(variable = priors[i, variable],
               value = eval_at,
               density = sapply(eval_at, function(x) (priors[i, density][[1]])(x, uselog = FALSE)))}) %>%
    rbindlist()
  prior_evaluated[, scaled_density := density/max(density), by="variable"] %>% .[, type := "prior"]
  prior_evaluated[variable %in% c("VE_waning_full", "VE_waning_minor"), value := value+2]
  
  scales = lapply(1:priors[, .N], function(i, priors) {
    scale_x_continuous(limits = c(priors[i, plotmin], priors[i, plotmax]))
  }, priors)
  
  posterior_long %>%
    ggplot(aes(x=value, fill=type, colour=type))+
    facet_wrap(facets = ~ factor(variable, priors$variable), scales = "free")+
    geom_area(data = prior_evaluated, aes(y=density), alpha=0.75)+
    geom_density(alpha=0.75)+
    geom_vline(data = posterior_long %>%
                 .[, .(med=median(value), l95=quantile(value, 0.025), u95=quantile(value, 0.975)), by=c("variable")],
               aes(xintercept=med))+
    geom_vline(data = posterior_long %>%
                 .[, .(med=median(value), l95=quantile(value, 0.025), u95=quantile(value, 0.975)), by=c("variable")],
               aes(xintercept=l95), linetype=2)+
    geom_vline(data = posterior_long %>%
                 .[, .(med=median(value), l95=quantile(value, 0.025), u95=quantile(value, 0.975)), by=c("variable")],
               aes(xintercept=u95), linetype=2)+
    theme_bw()+
    labs(x="Value", y="Density")+
    scale_fill_manual(values=c("prior" = "#DDDDDD", "posterior" = "#777777"))+
    scale_colour_manual(values=c("prior" = "#DDDDDD", "posterior" = "#777777"))+
    ggh4x::facetted_pos_scales(x = scales)
}

parallelModelRuns = function(singleRun, iterations = 10, cores = parallel::detectCores()){
  #' setup parallel environment
  if (cores > parallel::detectCores()) stop("More cores specified than available on this machine")
  cl = parallel::makeCluster(cores)
  
  #' load environments
  packages = (.packages())
  tmpdlls = getLoadedDLLs()
  dlls = vector(mode = "character", length = length(tmpdlls))
  counter = 0
  for (i in tmpdlls) {
    counter = counter + 1
    dlls[counter] = i[[2]]
  }
  objects = ls(envir = .GlobalEnv)
  objects = c(objects)
  packageFun = function(packages = NULL, dlls = NULL) {
    if (!is.null(packages)) {
      for (i in packages) library(i, character.only = TRUE)
    }
    if (!is.null(dlls)) {
      for (i in dlls) try(dyn.load(i), silent = T)
    }
  }
  
  parallel::clusterCall(cl, packageFun, packages, dlls)
  parallel::clusterExport(cl, varlist = objects)
  
  #' get modelled estimates for each iteration (may take a while)
  posterior_prevalence = parallel::parLapplyLB(cl, seq_len(iterations), function(i){
    if(i %% floor(iterations/10) == 0) message(sprintf("%s/%s (%s%%)", i, iterations, round(i/iterations * 100)))
    set.seed(i)
    
    #' sample values from the posterior
    sample_posterior = unlist(posterior[sample(1:.N, 1), ])
    
    #' Run model
    model_run = singleRun(sample_posterior)
    
    model_run[, run := i]
    return(model_run)}) %>% rbindlist
  
  #' stop parallel environment
  parallel::stopCluster(cl = cl)
  
  # process output
  colnames(posterior_prevalence)[2] = "age_group"
  
  return(posterior_prevalence)
}