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
  
  compartments_prevalence = model_params$global_settings$compartments_prevalence
  compartments_incidence = model_params$global_settings$compartments_incidence
  
  #' create data table to match to molten data
  columns_prevalence = model_params$populations %>%
    (function(populations){
      lapply(names(populations), function(pop, populations){
        v_strata = populations[[pop]][["vaccination_groups"]]
        data.table(
          population = pop %>% rep(length(v_strata) * length(compartments_prevalence) * age_groups_model[, .N]),
          vaccination_group = names(v_strata) %>% rep(each = length(compartments_prevalence) * age_groups_model[, .N]),
          compartment = compartments_prevalence %>% rep(each=age_groups_model[, .N]) %>% rep(length(v_strata)),
          age = seq_len(age_groups_model[, .N]) %>% rep(length(compartments_prevalence)) %>% rep(length(v_strata))
        )}, populations)}) %>% rbindlist()
  columns_prevalence[, c("id", "outcome") := .(as.character(.I), "prevalence")]
  columns_match = columns_prevalence
  
  if(incidence){
    columns_incidence = model_params$populations %>%
      (function(populations){
        lapply(names(populations), function(pop, populations){
          v_strata = populations[[pop]][["vaccination_groups"]]
          data.table(
            population = pop %>% rep(length(v_strata) * length(compartments_incidence) * age_groups_model[, .N]),
            vaccination_group = names(v_strata) %>% rep(each = length(compartments_incidence) * age_groups_model[, .N]),
            compartment = compartments_incidence %>% rep(each=age_groups_model[, .N]) %>% rep(length(v_strata)),
            age = seq_len(age_groups_model[, .N]) %>% rep(length(compartments_incidence)) %>% rep(length(v_strata))
          )}, populations)}) %>% rbindlist()
    columns_incidence[, c("id", "outcome") := .(paste0("inc_", .I), "incidence")]
    columns_match = rbind(columns_match, columns_incidence)
  }
  
  columns_match[, population := factor(population, names(model_params$populations))]
  columns_match[, vaccination_group := factor(vaccination_group, model_params$populations %>% sapply(function(p) names(p[["vaccination_groups"]])) %>% unlist() %>% as.vector() %>% unique())]
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
      pop_unvacc = rep(model_output[, unique(population)], length(model_params$populations))
    } else if(all(names(model_params$populations) %in% model_output[, unique(population)])){
      pop_unvacc = names(model_params$populations)
    } else {
      stop("eqStatesVaccinate2: method not yet implemented")
    }
  }
  names(pop_unvacc) = names(model_params$populations)
  
  model_output = model_output[outcome == "prevalence"]
  model_output %>% setorder(population, vaccination_group, compartment, age)
  
  model_input = model_params$populations %>%
    (function(populations){
      lapply(names(populations), function(pop, populations){
        v_strata = names(populations[[pop]][["vaccination_groups"]])
        v_strata %>% lapply(function(vstrat, pop){
          copy(model_output[population == pop_unvacc[[pop]]]) %>%
            .[, c("population", "vaccination_group", "value") := .(pop, vstrat, ifelse(vstrat == "unvaccinated", value, 0)), by=c("compartment", "age")] %>% .[]
        }, pop) %>% rbindlist()}, populations) %>% rbindlist})
  
  model_input[, population := factor(population, names(model_params$populations))]
  model_input[, vaccination_group := factor(vaccination_group, model_params$populations %>% sapply(function(p) names(p[["vaccination_groups"]])) %>% unlist() %>% as.vector() %>% unique())]
  model_input[, outcome := factor(outcome, c("prevalence"))]
  model_input[, compartment := factor(compartment, model_params$global_settings$compartments_prevalence)]
  
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
  model_params$params_unvac$populations %<>% lapply(function(x){
    x$parameters$betaVT %<>% multTimestep(TRUE, x$parameters$N)
    x$parameters$betaNVT %<>% multTimestep(TRUE, x$parameters$N)
    x$vaccination_groups %<>% lapply(function(z){
      z$waning %<>% multTimestep()
      return(z)})
    return(x)})
  
  model_params$params_vac$clearVT %<>% multTimestep()
  model_params$params_vac$clearNVT %<>% multTimestep()
  model_params$params_vac$ageout %<>% multTimestep()
  model_params$params_vac$migration[which(model_params$params_vac$migration != -1)] %<>% multTimestep()
  model_params$params_vac$populations %<>% lapply(function(x){
    x$parameters$betaVT %<>% multTimestep(TRUE, x$parameters$N)
    x$parameters$betaNVT %<>% multTimestep(TRUE, x$parameters$N)
    x$vaccination_groups %<>% lapply(function(z){
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
  data.table(variable = name, min = min, max = max,
             density = function(x, uselog=TRUE) dbeta((x - min)/range, shape1 = shape1, shape2 = shape2, log=uselog),
             sampler = function(n) rbeta(n, shape1 = shape1, shape2 = shape2) * range + min)
}

createUnifPriorBT = function(name, min = 0, max = 1, plotmin = NULL, plotmax = NULL){
  if(is.null(plotmin)) plotmin = min
  if(is.null(plotmax)) plotmax = max
  data.table(variable = name, min = min, max = max,
             density = function(x, uselog=TRUE) dunif(x, min = min, max = max, log=uselog),
             sampler = function(n) runif(n, min, max))
}

createLogNormPriorBT = function(name, min = 0, max = 1, plotmin = NULL, plotmax = NULL, meanlog, sdlog, flippedx = FALSE){
  if(is.null(plotmin)) plotmin = min
  if(is.null(plotmax)) plotmax = max
  data.table(variable = name, min = min, max = max,
             density = function(x, uselog=TRUE) dlnorm(ifelse(flippedx, 1 - x, x), meanlog = meanlog, sdlog = sdlog, log = uselog),
             sampler = function(n){
               #in_range = FALSE
               #while(!in_range){
              #   val = rlnorm(n, meanlog = meanlog, sdlog = sdlog)
              #   if(flippedx) val = 1 - val
              #   in_range = (val >= min & val <= max)
              # }
              # return(val)
               rlnorm(n, meanlog = meanlog, sdlog = sdlog)
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
      file.copy(sprintf("%s/model/%s.cpp", METAVAX_FOLDER, MODEL_NAME), sprintf("%s/model/%s.cpp", METAVAX_FOLDER, new_file_name))
      compileModel(sprintf("%s/model/%s.cpp", METAVAX_FOLDER, new_file_name), "./model/build/",
                   sprintf("%s%s", new_file_name, .Platform$dynlib.ext))
      file.remove(sprintf("%s/model/%s.cpp", METAVAX_FOLDER, new_file_name))
    }
  }
  if(!is.loaded("derivs", new_file_name)) dyn.load(new_file_path)
  if(!is.loaded("derivs", new_file_name)) stop("MetaVax is not loaded")
  
  return(new_file_name)
}

#' Function that will be used in MCMC algorithm
runModel = function(initial_state, model_params, steady_state = FALSE, times = c(0, 1), hmin = 0, hmax = NULL, rtol = 1e-06, atol = 1e-06, incidence = FALSE, parallel = FALSE){
  difference_equations = model_params$global_settings$model_solver_difference
  if(difference_equations & !"solver_difference_delta_t" %in% names(model_params$global_settings)) stop("solver_difference_delta_t not defined")
  
  nout_incidence = model_params$populations %>%
    sapply(function(x) length(x[["vaccination_groups"]])) %>%
    sum() * model_params$global_settings$age_groups_model[, .N] * length(model_params$global_settings$compartments_incidence)
  
  if(steady_state){
    #always use ODE model to calculate steady state
    model_params$global_settings$solver_difference = FALSE
    result = runsteady(
      y = initial_state, func = "derivs",
      initpar = model_params, dllname = {if(parallel) uniqueSharedObject() else MODEL_NAME},
      nout = ifelse(incidence, nout_incidence, 1), outnames = {if(incidence) paste0("inc_", seq_len(nout_incidence)) else "output"},
      initfunc = "rt_initmod", jactype = "fullint",
      hmin = hmin, hmax = hmax, atol = atol, rtol = rtol) %>% tryCatchWE()
  } else {
    model_params$global_settings$solver_difference = model_params$global_settings$model_solver_type == "DIFF"
    
    campaign_times = sapply(model_params$populations,
                            function(p){
                              sapply(p$vaccination_groups,
                                     function(vaccination_group) sapply(vaccination_group$coverage_c, "[[", "time")) %>%
                                unlist %>% unique %>% sort}) %>% unlist() %>% sort() %>% unique()
    if(length(campaign_times) == 0){
      if(difference_equations){
        result = ode(
          y=initial_state, times=times, func = "derivs",
          parms = model_params, dllname = {if(parallel) uniqueSharedObject() else MODEL_NAME},
          nout = ifelse(incidence, nout_incidence, 1), outnames = {if(incidence) paste0("inc_", seq_len(nout_incidence)) else "output"},
          initfunc = "initmod", method = "iteration") %>% tryCatchWE  
      } else {
        result = lsode(
          y=initial_state, times=times, func = "derivs",
          parms = model_params, dllname = {if(parallel) uniqueSharedObject() else MODEL_NAME},
          nout = ifelse(incidence, nout_incidence, 1), outnames = {if(incidence) paste0("inc_", seq_len(nout_incidence)) else "output"},
          initfunc = "initmod", jactype = "fullint",
          hmin = hmin, hmax = hmax, atol = atol, rtol = rtol) %>% tryCatchWE
      }
    } else {
      events = list(func="vaccineCampaignEvent",
                    time=campaign_times)
      times = sort(unique(c(times, campaign_times)))
      if(difference_equations){
        result = ode(
          y=initial_state, times=times, func = "derivs",
          parms = model_params, dllname = {if(parallel) uniqueSharedObject() else MODEL_NAME},
          nout = ifelse(incidence, nout_incidence, 1), outnames = {if(incidence) paste0("inc_", seq_len(nout_incidence)) else "output"},
          initfunc = "initmod",
          events = events, method = "iteration") %>% tryCatchWE  
      } else {
        result = lsode(
          y=initial_state, times=times, func = "derivs",
          parms = model_params, dllname = {if(parallel) uniqueSharedObject() else MODEL_NAME},
          nout = ifelse(incidence, nout_incidence, 1), outnames = {if(incidence) paste0("inc_", seq_len(nout_incidence)) else "output"},
          initfunc = "initmod", jactype = "fullint",
          events = events, hmin = hmin, hmax = hmax, atol = atol, rtol = rtol) %>% tryCatchWE   
      }
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
    arm_names = names(population$vaccination_groups)
    population$vaccination_groups = lapply(arm_names, function(name, arms){
      arm = arms[[name]]
      arm$coverage_r = arm$coverage_r %>% lapply(function(x, arms, name){
        if(!is.null(x$coverage_to)){
          if(is.numeric(x$coverage_to)){
            #message("coverage_to is numeric, assuming correct indices are already provided")
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
            #message("coverage_to is numeric, assuming correct indices are already provided")
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
    }, population$vaccination_groups)
    
    names(population$vaccination_groups) = arm_names
    
    return(population)
  }) 
}

renameCoverageTo2 = function(model_params){
  model_params$populations = lapply(model_params$populations, function(population){
    arm_names = names(population$vaccination_groups)
    population$vaccination_groups = lapply(arm_names, function(name, arms){
      arm = arms[[name]]
      arm$coverage_r = arm$coverage_r %>% lapply(function(x, arms, name){
        if(!is.null(x$coverage_to)){
          if(is.numeric(x$coverage_to)){
            #message("coverage_to is numeric, assuming correct indices are already provided")
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
            #message("coverage_to is numeric, assuming correct indices are already provided")
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
    }, population$vaccination_groups)
    
    names(population$vaccination_groups) = arm_names
    
    return(population)
  })
  
  return(model_params)
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
  
  #' add for last iteration
  chain_last = getSample(sampler, parametersOnly = T, coda = T, start = max(0, nrow(chain[[1]]) - sampler$settings$iterations/length(chain)))
  medi_last <- numeric(npar)
  lowerq_last <- upperq_last <- numeric(npar)
  if ("mcmc.list" %in% class(chain)) {
    for (i in 1:npar) {
      tmp <- unlist(chain_last[, i])
      tmp <- quantile(tmp, probs = c(0.025, 0.5, 0.975))
      lowerq_last[i] <- round(tmp[1], 3)
      medi_last[i] <- round(tmp[2], 3)
      upperq_last[i] <- round(tmp[3], 3)
    }
  } else {
    for (i in 1:npar) {
      tmp <- quantile(chain_last[, i],
                      probs = c(0.025, 0.5, 0.975))
      lowerq_last[i] <- round(tmp[1], 3)
      medi_last[i] <- round(tmp[2], 3)
      upperq_last[i] <- round(tmp[3], 3)
    }
  }
  
  divider = seq_len(length(lowerq))
  parOutDF <- cbind(MAPvals, lowerq, medi, upperq, divider, lowerq_last, medi_last, upperq_last)
  colnames(parOutDF) <- c("MAP", "2.5%", "median", "97.5%", "", "2.5%", "median", "97.5%")
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
combineOutFiles = function(OUTPUT_FOLDER, align_start = FALSE){
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
        if(align_start){
          x = x[1:min_it, ]  
        } else {
          x = x[seq(to = nrow(x), length.out = min_it), ]
        }
        class(x) = "mcmc"
        return(x)})
      
      class(out[[i]]$chain) = "mcmc.list"
    }
    
    class(out) = c("mcmcSamplerList", "bayesianOutput")
  }
  
  return(out)
}

createTracePlot = function(posterior){
  #' Let's also assess the cross correlation table for fitted parameters
  posterior = as.matrix(posterior)
  posterior[, 1:priors[, .N]] = apply(posterior[, 1:priors[, .N]], 1, function(z){
    names(z) = priors$variable
    return(as.vector(updateDepParameters(z)))
  }) %>% t()
  posterior = as.data.table(posterior)
  
  posterior %>%
    .[, i := 1:.N, by="chain"] %>%
    melt(id.vars=c("chain", "i")) %>%
    ggplot(aes(x=i, y=value, colour=as.factor(chain)))+
    facet_wrap(factor(variable, priors$variable)~., scales = "free")+
    geom_line(alpha=0.75)+
    theme_bw()+
    labs(x="Iteration", y="Value", colour="Chain")+
    guides(colour = guide_legend(override.aes = list(size = 5))) 
}

createPriorPosteriorPlot = function(priors, posterior){
  #' Let's also assess the cross correlation table for fitted parameters
  posterior = as.matrix(posterior)
  posterior[, 1:priors[, .N]] = apply(posterior[, 1:priors[, .N]], 1, function(z){
    names(z) = priors$variable
    return(as.vector(updateDepParameters(z)))
  }) %>% t()
  posterior = as.data.table(posterior)
  
  posterior_long = posterior %>%
    .[, i := 1:.N, by="chain"] %>%
    melt(id.vars=c("chain", "i")) %>%
    .[, type := "posterior"]
  
  priors = priors %>% cbind(lapply(priors[, sampler], function(s){
    samples = s(1000)
    plotmin = quantile(samples, 0.25)
    plotmax = quantile(samples, 0.75)
    data.table(plotmin = plotmin,
               plotmax = plotmax)
  }) %>% rbindlist())
  
  priors = priors %>% merge(posterior_long[, .(postmin = quantile(value, 0.001), postmax = quantile(value, 0.999)), by = "variable"],
                   by = "variable")
  priors[, c("plotmin", "plotmax") := .(min(plotmin, postmin), max(plotmax, postmax)), by = "variable"]
  
  prior_evaluated = lapply(1:nrow(priors), function(i){
    eval_at = seq(from=priors[i, plotmin], to=priors[i, plotmax], length=1000)
    data.table(variable = priors[i, variable],
               value = eval_at,
               density = sapply(eval_at, function(x) (priors[i, density][[1]])(x, uselog = FALSE)))}) %>%
    rbindlist()
  prior_evaluated[, scaled_density := density/max(density), by="variable"] %>% .[, type := "prior"]
  prior_evaluated = prior_evaluated %>%
    merge(posterior_long[, .(maxdensity = max(hist(value, breaks = 20, plot = FALSE)$density)), by=variable]) %>%
    .[, scaled_density := scaled_density * maxdensity]
  
  scales = lapply(1:priors[, .N], function(i, priors) {
    scale_x_continuous(limits = c(priors[i, plotmin], priors[i, plotmax]))
  }, priors)
  
  posterior_long %>%
    ggplot(aes(x=value, fill=type, colour=type))+
    facet_wrap(facets = ~ factor(variable, priors$variable), scales = "free")+
    geom_area(data = prior_evaluated, aes(y=scaled_density), alpha=0.75)+
    geom_histogram(alpha=0.50, bins = 20, aes(y = ..density..))+
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
    labs(x="Value", y="Density (scaled)")+
    scale_fill_manual(values=c("prior" = "#DDDDDD", "posterior" = "#777777"))+
    scale_colour_manual(values=c("prior" = "#DDDDDD", "posterior" = "#777777"))+
    ggh4x::facetted_pos_scales(x = scales)+
    theme(axis.text.x = element_text(angle=45, hjust=1),
          axis.text.y = element_text(colour = "#00000000"),
          axis.ticks.y = element_blank())
}

parallelModelRuns = function(singleRun, parameter_values, cores = parallel::detectCores()){
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
  posterior_prevalence = parallel::parLapplyLB(cl, seq_len(nrow(parameter_values)), function(i){
    set.seed(i)
    
    #' sample values from the posterior
    sample_posterior = unlist(parameter_values[i, ])
    
    #' Run model
    model_run = singleRun(sample_posterior)
    
    model_run[, run := i]
    return(model_run)}) %>% rbindlist
  
  #' stop parallel environment
  parallel::stopCluster(cl = cl)
  
  return(posterior_prevalence)
}

#' functions to create model_params
createParameters = function(model_params,
                            age_groups_model = setAgeBreaks(0),
                            populations = c("pop1"),
                            vaccination_groups = c("unvaccinated"),
                            model_start_date = Sys.Date()){
  #' set start date
  model_params[["global_settings"]][["model_start_date"]] = model_start_date
  
  #' set age group related values
  model_params[["global_settings"]][["age_groups_model"]] = age_groups_model
  model_params[["global_settings"]][["n_agrp"]] = age_groups_model[, .N]
  model_params[["global_settings"]][["ageout"]] = age_groups_model %>%
    .[, .(duration = (to - from) %>% set_units("days"))] %>%
    .[, 1/as.numeric(duration)]
  
  #' set populations and vaccination strata
  model_params[["populations"]] = lapply(populations, function(population){
    list(parameters = list(),
         vaccination_groups = lapply(vaccination_groups, function(vaccination_group) list()) %>% setNames(vaccination_groups))
  }) %>% setNames(populations)
  
  return(model_params)
}

setParameter = function(model_params,
                        key,
                        value,
                        level = c("global", "population", "vaccination_group")[1],
                        population = "all",
                        vaccination_group = "all",
                        strict = TRUE){
  #' helper function to check population level
  validatePopulation = function(model_params, population){
    validate_result = population %in% names(model_params$populations)
    
    if(!validate_result)
      if(strict) stop(sprintf("Population %s does not exist", population))
      else warning(sprintf("Population %s does not exist", population))
    
    return(validate_result)
  }
  
  #' helper function to check vaccination stratum level
  validateVaccinationGroup = function(model_params, population, vaccination_group){
    validate_result = vaccination_group %in% names(model_params$populations[[population]]$vaccination_groups)
    
    if(!validate_result)
      if(strict) stop(sprintf("Vaccination group %s does not exist in population %s", vaccination_group, population))
      else warning(sprintf("Vaccination group %s does not exist in population %s", vaccination_group, population))
    
    return(validate_result)
  }
  
  #' check that correct level is specified
  if(!level %in% c("global", "population", "vaccination_group"))
    stop("Wrong level")
  
  #' overwrite with population names if needed
  if(length(population) == 1) if(population == "all") population = names(model_params$populations)
  
  #' set key with value
  if(level == "global"){
    model_params$global_settings[[key]] = value
  } else if(level == "population"){
    for(p in population){
      if(!validatePopulation(model_params, p)) next
      model_params$populations[[p]]$parameters[[key]] = value
    }
  } else if(level == "vaccination_group"){
    for(p in population){
      if(!validatePopulation(model_params, p)) next
      
      #' overwrite with vaccination stratum names if needed
      if(length(vaccination_group) == 1) if(vaccination_group == "all") vaccination_group = names(model_params$populations[[p]]$vaccination_groups)
      
      for(v in vaccination_group){
        if(!validateVaccinationGroup(model_params, p, v)) next
        
        model_params$populations[[p]]$vaccination_groups[[v]][[key]] = value
      }
    }
  }
  
  return(model_params)
}

createBTAdditionalPosterior = function(additional_posteriors){
  additional_posteriors = list(
    sampler = function(n = 1){
      result = additional_posteriors$setup[, `sampler`] %>%
        sapply(function(s, n){ s(n) }, n = n)
      
      if(n == 1){
        names(result) = (additional_posteriors$setup[, `variable`])
      } else {
        colnames(result) = additional_posteriors$setup[, `variable`]
      }
      
      return(result)
    },
    setup = additional_posteriors)
  
  return(additional_posteriors)
}

getParameter = function(model_params,
                        key,
                        level = c("global", "population", "vaccination_group")[1],
                        population = "all",
                        vaccination_group = "all"){
  #' helper function to check population level
  validatePopulation = function(model_params, population){
    if(!population %in% names(model_params$populations))
      stop(sprintf("Population %s does not exist"))
  }
  
  #' helper function to check vaccination stratum level
  validateVaccinationGroup = function(model_params, population, vaccination_group){
    if(!vaccination_group %in% names(model_params$populations[[population]]$vaccination_groups))
      stop(sprintf("Vaccination stratum %s does not exist in population %s", vaccination_group, population))
  }
  
  #' check that correct level is specified
  if(!level %in% c("global", "population", "vaccination_group"))
    stop("Wrong level")
  
  #' overwrite with population names if needed
  if(length(population) == 1) if(population == "all") population = names(model_params$populations)
  
  #' set key with value
  if(level == "global"){
    result = data.table(var1 = model_params$global_settings[[key]]) %>% setNames(key)
  } else if(level == "population"){
    result = lapply(population, function(p){
      validatePopulation(model_params, p)
      
      data.table(population = p,
                 var1 = model_params$populations[[p]]$parameters[[key]]) %>%
        .[, age := .I] %>% setNames(c("population", key, "age"))
    }) %>% rbindlist()
  } else if(level == "vaccination_group"){
    result = lapply(population, function(p){
      validatePopulation(model_params, p)
      
      #' overwrite with vaccination stratum names if needed
      if(length(vaccination_group) == 1) if(vaccination_group == "all") vaccination_group = names(model_params$populations[[p]]$vaccination_groups)
      
      lapply(vaccination_group, function(v){
        validateVaccinationGroup(model_params, p, v)
        
        data.table(population = p,
                   vaccination_group = v,
                   var1 = model_params$populations[[p]]$vaccination_groups[[v]][[key]]) %>%
          .[, age := .I] %>% setNames(c("population", "vaccination_group", key, "age"))
      }) %>% rbindlist()
    }) %>% rbindlist()
  }
  
  return(result)
}

checkModelOutput = function(modelled_result, .eps = .Machine$double.eps){
  #' use eps when comparing to ignore very small floating point errors
  if(modelled_result[outcome == "prevalence", any(value < (0 - .eps)) | any(value > (1 + .eps))]){
    warning(sprintf("Some modelled values <0 (min: %s) or >1 (max: %s).",
                    modelled_result[, min(value)], modelled_result[, max(value)]))
    return(FALSE)
  }
  
  return(TRUE)
}

aggregateModelOutput = function(modelled_result, model_params, aggregate_agegroups,
                                by_vaccination_group = TRUE, by_population = TRUE, by_compartment = TRUE, by_time = TRUE, additional_by = NULL){
  #' specify columns to group results by
  #' - note outcome and age_group are always included
  by_cols = c("population", "vaccination_group", "outcome", "age_group", "compartment", "time", additional_by) %>% unique()
  by_cols_N = c("population", "age_group")
  if(!by_population){
    by_cols = by_cols %>% subset(. != "population")
    by_cols_N = by_cols_N %>% subset(. != "population")
  }
  if(!by_vaccination_group) by_cols = by_cols %>% subset(. != "vaccination_group")
  if(!by_compartment) by_cols = by_cols %>% subset(. != "compartment")
  if(!by_time) by_cols = by_cols %>% subset(. != "time")
  
  #' see what age groups match
  matching_age_groups = model_params$global_settings$age_groups_model %>%
    matchingAgeBreaks(aggregate_agegroups)
  
  #' sum N by age group
  N = getParameter(model_params, level = "population", population = unique(modelled_result$population), key = "N")
  N_aggregated = N %>%
    merge(model_params$global_settings$age_groups_model[, c("age", "name")], by = "age") %>%
    merge(matching_age_groups, by.x = "name", by.y = "name.x") %>%
    .[, age_group := name.y] %>%
    .[, .(N = sum(N)), by = by_cols_N]
  
  #' sum values by age group  
  result_aggregated = modelled_result %>%
    merge(N) %>% .[, value := value * N] %>% .[, -"N"] %>%
    merge(matching_age_groups, by.x = "age_group", by.y = "name.x") %>%
    .[, age_group := name.y] %>%
    .[, .(value = sum(value)), by = by_cols]
  
  #' revert back to proportion for prevalence, keep incidence as absolute number
  result_aggregated = result_aggregated %>% merge(N_aggregated) %>%
    #' when aggregating time, want the average prevalence over all timesteps
    .[, N := ifelse(by_time, N, N * modelled_result[, length(unique(time))]), by = by_cols] %>%
    .[, value := ifelse(outcome == "prevalence", value/N, value)] %>%
    .[, -"N"]
  
  return(result_aggregated)
}

testLL = function(ll, prior){
  i = 0
  x = prior$sampler()
  y = ll(x)
  
  while(is.infinite(y)){
    i = i + 1
    message(sprintf("i: %s", i))
    x = prior$sampler()
    y = ll(x)
  }
  
  message(sprintf("Found value after %s iterations; ll was %s", i, y))
}

setOutputFolder = function(working_directory, subdirectory){
  OUTPUT_FOLDER = sprintf("%s/output", working_directory)
  if(!dir.exists(OUTPUT_FOLDER)) dir.create(OUTPUT_FOLDER)
  
  OUTPUT_FOLDER = sprintf("%s/%s", OUTPUT_FOLDER, subdirectory)
  if(!dir.exists(OUTPUT_FOLDER)) dir.create(OUTPUT_FOLDER)
  
  return(OUTPUT_FOLDER)
}

fitBT = function(bayesianSetup, settings, output_folder, chain, i){
  if(i == 1){
    burned_in = FALSE
  } else {
    burned_in = TRUE
    out = readRDS(sprintf("%s/out_%s_%s.RDS", output_folder, chain, i))
    
    #out$setup$likelihood$cl = bayesianSetup$likelihood$cl
    
    message(sprintf("chain %s - continue with i %s", chain, i))
    altSummary(out)
    message("-----")
    
    i = i+1
  }
  
  while(!file.exists("./stop_sampler")){
    tstart = proc.time()
    
    if(!burned_in){
      out = runMCMC(bayesianSetup = bayesianSetup, settings = settings)
    } else {
      out = runMCMC(bayesianSetup = out)
    }
    
    tend = proc.time()
    elapsed = tend-tstart
    
    if(!burned_in){
      burned_in = TRUE
      
      out$settings$burnin = 0
      message(sprintf("%s - Finished with %s, continue and start saving output... Resetting i to t=1", Sys.time(), i))
    }
    
    altSummary(out)
    message(sprintf("%s - Finished with %s, continue... Last run took %s minutes", Sys.time(), i, round(elapsed[3]/60, 2)))
    
    saveRDS(out, sprintf("%s/out_%s_%s.RDS", output_folder, chain, i))
    if(file.exists(sprintf("%s/out_%s_%s.RDS", output_folder, chain, i - 1))){
      file.remove(sprintf("%s/out_%s_%s.RDS", OUTPUT_FOLDER, chain, i - 1))
    }
    
    i = i + 1
  }
}

samplePosterior = function(out, start = 0, thin = 0, variable_names){
  if("mcmcSamplerList" %in% class(out)){
    posterior = lapply(seq_along(out),
                       function(x, out) getSample(out[[x]], start = start, thin = 10, coda = FALSE) %>%
                         as.data.table %>% .[, chain := x] %>% return, out) %>% rbindlist  
  } else {
    posterior = getSample(out, start = start, thin = 10, coda = FALSE) %>%
      as.data.table %>% .[, chain := 1]
  }
  
  
  colnames(posterior) = c(variable_names, "chain")  
  
  return(posterior)
}

showPriorsTable = function(priors_table){
  lapply(seq_len(priors_table[, .N]), function(i){
    samples = priors_table[i, sampler][[1]](1e6)
    data.table(variable = priors_table[i, variable],
               sampler = priors_table[i, sampler] %>%
                 as.character() %>%
                 gsub("function (n) \nr", "", x = ., fixed = TRUE) %>%
                 gsub("(n, ", "(", x = ., fixed = TRUE),
               med = median(samples),
               low95 = quantile(samples, 0.025),
               high95 = quantile(samples, 0.975))
  }) %>% rbindlist() 
}
