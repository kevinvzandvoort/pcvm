library(units)
library(data.table)
library(magrittr)

#' This function creates a data table with age groups, according to the provided age breaks
setAgeBreaks = function(age_breaks, minage = set_units(0, "years"), maxage = set_units(120, "years")){
  #' Helper function to check whether a value is an Integer for a given tolerance level
  #' - alternative to the is.integer function of base R, which does not check for a specific tolerance
  isAnInteger = function(value, tol = 1/1e10){
    sapply(value, function(x, tol) (abs(as.numeric(x) - round(as.numeric(x))) < tol), tol) 
  }
  age_types = c(hours = "h", days = "d", weeks = "w", months = "m", years = "y")
  
  #' Convert all age breaks to years, and order all unique values
  age_breaks_years = age_breaks %>% set_units("years") %>% sort %>% unique
  
  #' Exclude any age breaks lower than the minimum population age (default 0 years)
  age_breaks_years = age_breaks_years %>% subset(. >= minage)
  
  #' Add the minimum age if not included in the age breaks
  if(age_breaks_years[1] > set_units(minage, "years")) age_breaks_years = c(minage, age_breaks_years)
  
  #' Exclude any age breaks over the maximum age (default 120 years, will be added automatically if not present)
  age_breaks_years = age_breaks_years %>% subset(. < maxage)
  
  #' Vector to store type of age break
  age_breaks_types = character(length(age_breaks_years))
  age_breaks_unresolved = rep(TRUE, length(age_breaks_years))
  
  #' loop through all age types, and assign to each age breaks
  #' - e.g. can it be represented as an integer when represented as a year, month, week, or day?
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
  
  #' Ensure logical order of age
  #' - make sure they are always ordered as days, weeks, months, years
  start = 1
  for(i in seq_along(age_types)){
    label = age_types[i]
    
    if(any(age_breaks_types == label)){
      age_breaks_types[start:max(which(age_breaks_types == label))] = label
      start = max(which(age_breaks_types == label)) + 1
    }
  }
  
  #' assign correct label type to each age break
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
  
  #' create user-friendly name for each age group
  age_group_names = character(length(age_breaks_years))
  for(i in 1:length(age_breaks_years)){
    from = age_breaks_years[i]
    from_unit = age_breaks_types[i]
    
    if(i < length(age_breaks_years)){
      to = age_breaks_years[i + 1]
      to_unit = age_breaks_types[i + 1]
    } else {
      to = maxage
      to_unit = "y"
    }
    
    age_group_names[i] = sprintf("[%s%s, %s%s)", 
                                 set_units(from, names(age_types)[which(age_types == from_unit)], mode = "standard"),
                                 from_unit,
                                 set_units(to, names(age_types)[which(age_types == to_unit)], mode = "standard"),
                                 to_unit)
    
  }
  
  #' Construct data table with age groups
  #' - format can be easily read by other functions
  age_groups = data.table(
    name = age_group_names,
    from = age_breaks_years,
    to = c(age_breaks_years[-1], maxage)
  )
  
  #' make name a factor for easy sorting
  age_groups[, name := factor(name, age_group_names)]
  
  return(age_groups[])
}

#' combine two age breaks tables
#' - useful to distribute data to different age groups
#' - TODO: rewrite in cpp
combineAgeBreaks = function(x, y, value.var="value", method = c("sum", "mean")[2], by = NULL){
  #' calculate overlapping age groups - this is done once for all groups, but would break if different age groups are present for different groups
  y_ages = unique(y[, c("from", "to")])
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
  
  y[, x %>% cbind(overlap %*% as.matrix(.SD[, value.var, with=FALSE])), by=by]
}

matchingAgeBreaks = function(x, y){
  onesec = set_units(1, "second")
  fivesec = set_units(5, "second")
  
  #' function to check whether one age group is lower than the other
  #' - tolerance is five seconds
  lt = function(a, b){
    (a - b) < fivesec
  }
  #' function to check whether one age group is greater than the other
  #' - tolerance is five seconds
  gt = function(a, b){
    (a - b) > fivesec
  }
  
  x = copy(x) %>% setorder(from, to)
  y = copy(y) %>% setorder(from, to)
  
  x[, row := 1:.N]
  x[, name.x := name]
  y[, row := 1:.N]
  
  for(a in 1:nrow(x)){
    #' get the age groups in y that overlap with the age group in x
    ages_that_overlap = y[(from + fivesec) <= x[a, to] & (to - fivesec) >= x[a, from]]
    
    if(ages_that_overlap[, .N] > 1) stop("More than 1 row in y to match row in x")
    
    #' get duration of each age group
    dur_x = x[a, to - from]
    dur_y = ages_that_overlap[, to - from]
    
    if(dur_y >= dur_x){
      x[a, name.y := ages_that_overlap[, name]]  
    } else {
      stop("Age-break in y is smaller than age break in x")
    }
  }
  
  return(x[, c("name.x", "name.y")])
}

#' Fill age gaps if there are gaps between age groups in a dataset
fillAgeGaps = function(data, min_age = set_units(0, "years"), max_age = set_units(120, "years")){
  if(data[1, from] > min_age){
    
    #' create age-data for row
    row_age_data = setAgeBreaks(min_age, minage = min_age, maxage = data[2, from])
    
    #' set all measurement values to NA
    row_data = copy(data[1, -c("name", "from", "to")])
    for(j in 1:ncol(row_data)) row_data[, j] = NA
    
    #' create new row
    row_data = cbind(row_age_data, row_data)
    
    #' update dataset
    data = rbind(row_data, data)
  }
  
  if(data[.N, to] < max_age){
    
    #' create age-data for row
    row_age_data = setAgeBreaks(data[.N, to], minage = data[.N, to], maxage = max_age)
    
    #' set all measurement values to NA
    row_data = copy(data[1, -c("name", "from", "to")])
    for(j in 1:ncol(row_data)) row_data[, j] = NA
    
    #' create new row
    row_data = cbind(row_age_data, row_data)
    
    #' update dataset
    data = rbind(data, row_data)
  }
  
  #' N rows updated during the for-loop
  #' - run this several times until the dataset is completed
  updated = TRUE
  while(updated){
    updated = FALSE
    for(i in 2:nrow(data)){
      if(data[i, from] != data[i-1, to]){
        updated = TRUE
        
        #' create age-data for row
        row_age_data = setAgeBreaks(data[i-1, to], minage = data[i-1, to], maxage = data[i, from])
        
        #' set all measurement values to NA
        row_data = copy(data[1, -c("name", "from", "to")])
        for(j in 1:ncol(row_data)) row_data[, j] = NA
        
        #' create new row
        row_data = cbind(row_age_data, row_data)
        
        #' update dataset
        data = rbind(data[1:(i-1)], row_data, data[i:.N])
      }
    }
  }
  
  
  return(data)
}

#' function to check whether one age group is lower than the other
#' - tolerance is five seconds
agelt = function(a, b){
  fivesec = set_units(5, "second")
  (a - b) < fivesec
}
#' function to check whether one age group is greater than the other
#' - tolerance is five seconds
agegt = function(a, b){
  fivesec = set_units(5, "second")
  (a - b) > fivesec
}

ageeq = function(a, b){
  fivesec = set_units(5, "second")
  (b >= (a - fivesec)) & (b <= (a + fivesec))
}


