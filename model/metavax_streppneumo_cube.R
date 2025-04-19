MODEL_NAME = "metavax_streppneumo_cube" #which metavax version should be used?
source(sprintf("%s/index.R", METAVAX_FOLDER))

model_params = list(
  global_settings = list(
    model_name = MODEL_NAME,
    compartments_prevalence = c("S", "VT", "OT", "NVT", "VT_OT", "VT_NVT", "OT_NVT", "VT_OT_NVT"),
    compartments_incidence = c("iVT", "iOT", "iNVT"),
    model_solver_type = "DIFF",
    model_solver_difference = TRUE,
    solver_difference_delta_t = 1))

#' this function can be used to process all incidence results
processIncidence = function(modelled_result, model_params){
  #' get ccr_values, accounting for vaccine efficacy against disease
  ccr_values = getParameter(model_params, level = "population", population = modelled_result[, unique(population)], key = "ccrVT") %>%
    .[, compartment := "iVT"] %>%
    setNames(c("population", "ccr", "age", "compartment")) %>%
    rbind(getParameter(model_params, level = "population", population = modelled_result[, unique(population)], key = "ccrNVT") %>%
            .[, compartment := "iNVT"] %>%
            setNames(c("population", "ccr", "age", "compartment"))) %>%
    rbind(getParameter(model_params, level = "population", population = modelled_result[, unique(population)], key = "ccrNVT") %>%
            .[, compartment := "iOT"] %>%
            setNames(c("population", "ccr", "age", "compartment"))) %>%
    merge(getParameter(model_params, level = "vaccination_group", population = modelled_result[, unique(population)], key = "efficacy_disease"),
          by = c("population", "age"), allow.cartesian = TRUE) %>%
    merge(getParameter(model_params, level = "vaccination_group", population = modelled_result[, unique(population)], key = "efficacy_disease_OT"),
          by = c("population", "age"), allow.cartesian = TRUE) %>%
    .[, efficacy_disease := ifelse(compartment == "iVT", efficacy_disease,
                                   ifelse(compartment == "iOT", efficacy_disease_OT, 0))] %>%
    .[, ccr := ccr * (1 - efficacy_disease)] %>%
    .[, -"efficacy_disease"]
  
  modelled_result = modelled_result[outcome == "incidence"] %>%
    merge(ccr_values, by = c("population", "vaccination_group", "compartment", "age")) %>%
    .[, value := value * ccr] %>% .[, -"ccr"] %>%
    rbind(modelled_result[outcome == "prevalence"])
  
  #' reorder columns and rows
  modelled_result = modelled_result[, c("outcome", "population", "vaccination_group", "compartment", "age_group", "age", "time", "value")]
  setorder(modelled_result, outcome, population, vaccination_group, compartment, age_group, time)
  
  return(modelled_result)
}
