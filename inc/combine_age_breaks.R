library(units)
library(data.table)

setAgeBreaks = function(age_breaks, minage = set_units(0, "years"), maxage = set_units(120, "years")){
  isAnInteger = function(value, tol = 1/1e10){
    sapply(value, function(x, tol) (abs(as.numeric(x) - round(as.numeric(x))) < tol), tol) 
  }
  age_types = c(hours = "h", days = "d", weeks = "w", months = "m", years = "y")
  
  age_breaks_years = age_breaks %>% set_units("years") %>% sort %>% unique
  age_breaks_years = age_breaks_years %>% subset(. >= minage)
  if(age_breaks_years[1] > set_units(minage, "years")) age_breaks_years = c(minage, age_breaks_years)
  age_breaks_years = age_breaks_years %>% subset(. < maxage)
  
  #' Check if there are any groups of this type
  age_breaks_types = character(length(age_breaks_years))
  age_breaks_unresolved = rep(TRUE, length(age_breaks_years))
  for(i in rev(seq_along(age_types))){
    if(sum(age_breaks_unresolved) == 0) break
    unit = names(age_types)[i]
    label = age_types[i]
    
    matching_age_groups = which(age_breaks_years[age_breaks_unresolved] %>% set_units(unit, mode = "standard") %>%
                                  as.numeric %>% isAnInteger)
    age_breaks_types[age_breaks_unresolved][matching_age_groups] = label
    age_breaks_unresolved = age_breaks_types == "" 
  }
  
  if(any(age_breaks_unresolved))
    stop(sprintf("Could not resolve type of all age-breaks provided: age-breaks %s",
                 which(age_breaks_unresolved) %>% paste0(collapse = ", ")))
  
  #' Ensure logical order
  start = 1
  for(i in seq_along(age_types)){
    label = age_types[i]
    
    if(any(age_breaks_types == label)){
      age_breaks_types[start:max(which(age_breaks_types == label))] = label
      start = max(which(age_breaks_types == label)) + 1
    }
  }
  
  prev_i = NULL
  for(i in seq_along(age_types)){
    label = age_types[i]
    
    if(any(age_breaks_types == label)){
      if(is.null(prev_i)){
        age_breaks_types[1:max(which(age_breaks_types == label))] = label
      } else {
        age_breaks_types[(max(which(age_breaks_types == age_types[prev_i])) + 1):
                           max(which(age_breaks_types == label))] = label
      } 
      
      prev_i = i
    }
  }
  
  age_group_names = character(length(age_breaks_years))
  for(i in seq_along(age_types)){
    unit = names(age_types)[i]
    label = age_types[i]
    
    if(!any(age_breaks_types == label)) next
    
    is_seq = age_breaks_years[age_breaks_types == label] %>%
      set_units(unit, mode = "standard") %>% diff %>% round %>% (function(x) as.integer(x) > 1)
    if(i < length(age_types)){
      is_seq_last = age_breaks_years[which(age_breaks_types == label) %>% rev %>% .[1] %>%
                                       (function(x) x + c(0, 1))] %>%
        set_units(unit, mode = "standard") %>% diff() %>% round %>% (function(x) as.integer(x) > 1)
      is_seq = c(is_seq, is_seq_last) 
    } else {
      is_seq_last = FALSE
    }
    
    if(any(!is_seq)){
      no_seqs = paste0(age_breaks_years[age_breaks_types == label][which(!is_seq)] %>%
                         set_units(unit, mode = "standard"), label)
      age_group_names[age_breaks_types[-length(age_breaks_types)] == label][which(!is_seq)] = no_seqs
    }
    if(any(is_seq)){
      seqs = paste0(age_breaks_years[age_breaks_types == label][which(is_seq)] %>%
                      set_units(unit, mode = "standard"), "-",
                    (age_breaks_years[age_breaks_types == label][1 + which(is_seq)] %>%
                       set_units(unit, mode = "standard") - set_units(1, unit, mode = "standard")), label)
      if(is_seq_last){
        #' Really only makes sense with other labelling
        last_name = paste0(age_breaks_years[age_breaks_types == label][max(which(is_seq))] %>%
                             set_units(unit, mode = "standard"), label, "- <",
                           age_breaks_years[which(age_breaks_types == label) %>% rev %>% .[1] + 1] %>%
                             set_units(age_types %>%
                                         subset(. == age_breaks_types[which(age_breaks_types == label) %>%
                                                                        rev %>% .[1] + 1]) %>%
                                         names, mode = "standard"),
                           age_breaks_types[which(age_breaks_types == label) %>% rev %>% .[1] + 1])
        
        warning(sprintf("Age-group '%s' spans a range between different types, check to make sure this is correct",
                        last_name))
        seqs[length(seqs)] = last_name
      }
      age_group_names[age_breaks_types[-length(age_breaks_types)] == label][which(is_seq)] = seqs
    }
  }
  if(length(age_group_names) > 1){
    age_group_names[1] = paste0("<",
                                age_breaks_years[2] %>%
                                  set_units(age_types %>% subset(. == age_breaks_types[2]) %>% names, mode = "standard"),
                                age_breaks_types[2])
    age_group_names[length(age_group_names)] = paste0(age_breaks_years %>% rev %>% .[1] %>%
                                                        set_units(age_types %>% subset(. == age_breaks_types %>% rev %>%
                                                                                         .[1]) %>% names,
                                                                  mode = "standard"), "+",
                                                      age_breaks_types %>% rev %>% .[1])
  } else {
    age_group_names[1] = paste0("<", maxage, "y")
  }
  
  age_groups = data.table(
    name = age_group_names,
    from = age_breaks_years,
    to = c(age_breaks_years[-1], maxage)
  )
  
  age_groups[, name := factor(name, age_group_names)]
  
  return(age_groups)
}

#by = c("group1", "group2")
by = NULL
method = c("sum", "mean")[2]
measure.vars = c("v1", "v2")

combineAgeBreaks = function(x, y, measure.vars, method = c("sum", "mean")[2], by = NULL){
  #' calculate overlapping age groups - this is done once for all groups, but would break if different age groups are present for different groups
  y_ages = unique(y[, c("name", "from", "to")])
  if(method == "sum"){
    #' Generate overlap matrix
    overlap = pmax(pmin(matrix(rep(x[, to], y_ages[, .N]), ncol=y_ages[, .N]),
              matrix(rep(y_ages[, to], x[, .N]), nrow=x[, .N], byrow = TRUE)) - 
           pmax(matrix(rep(x[, from], y_ages[, .N]), ncol=y_ages[, .N]),
                matrix(rep(y_ages[, from], x[, .N]), nrow=x[, .N], byrow = TRUE)), 0)/
      matrix(rep(y_ages[, to - from], x[, .N]), nrow=x[, .N], byrow = TRUE)
    
    #' Nb, this is the same as this
    #overlap = matrix(nrow = x[, .N], ncol = y_ages[, .N])
    #for(i in 1:x[, .N]){
    #  for(j in 1:y_ages[, .N]){
    #    overlap[i, j] = max(0, (min(y_ages[j, to], x[i, to]) - max(x[i, from], y_ages[j, from])))/y_ages[j, to-from]
    #  }
    #}
    
    #' And the same as this
    #overlap = matrix(nrow = x[, .N], ncol = y_ages[, .N])
    #for(i in 1:x[, .N]){
    #  overlap[i, ] = as.numeric(pmax(set_units(0, "year"), pmin(x[i, to], y_ages[, to]) - pmax(x[i, from], y_ages[, from]))/y_ages[, to-from])
    #}
  } else if(method == "mean"){
    #' This is the same, but with x and y swapped
    overlap = t(pmax(t(pmin(matrix(rep(x[, to], y_ages[, .N]), ncol=y_ages[, .N]),
                  matrix(rep(y_ages[, to], x[, .N]), nrow=x[, .N], byrow = TRUE))) - 
             t(pmax(matrix(rep(x[, from], y_ages[, .N]), ncol=y_ages[, .N]),
                    matrix(rep(y_ages[, from], x[, .N]), nrow=x[, .N], byrow = TRUE))), 0)/
        matrix(rep(x[, to - from], y_ages[, .N]), nrow=y_ages[, .N], byrow = TRUE))
  } else {
    stop(sprintf("Method %s is not implemented", method))
  }
  
  y[, x %>% cbind(overlap %*% as.matrix(.SD[, measure.vars, with=FALSE])), by=by]
}

x = setAgeBreaks(c(c(1:12) %>% set_units("month"), c(1:5, 10, 15, 20, 50, 70) %>% set_units("year")))
y = setAgeBreaks(sample(1:100, 4, replace=FALSE)) %>% .[, c("v1", "v2") := .(rpois(5, 3), rpois(5, 9))]
combineAgeBreaks(x, y, c("v2"), "mean")
combineAgeBreaks(x, y, c("v1", "v2"), "sum")

y = setAgeBreaks(sample(1:100, 4, replace=FALSE)) %>%
  rbind(., ., ., .) %>%
  .[, c("v1", "v2", "group1", "group2") := .(rpois(5*2*2, 3), rpois(5*2*2, 9), rep(rep(c("A", "B"), each=5), 2), rep(rep(c("Q", "R"), each=10)))]
combineAgeBreaks(x, y, c("v1", "v2"), "sum", by=c("group1", "group2"))
