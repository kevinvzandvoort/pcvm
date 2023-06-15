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
combineAgeBreaks = function(x, y, method = c("mean", "sum")[1], value.var = "value"){
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
  
  if(any(value.var %in% colnames(x))){
    for(val in value.var[which(value.var %in% colnames(x))]){
      x[, paste0(val,".x") := get(val)]
      x = x[, -val, with=FALSE]
      
      y[, paste0(val,".y") := get(val)]
      y = y[, -val, with=FALSE]
      
      value.var[which(value.var == val)] = paste0(val, ".y") 
    }
  }
  
  #' for each age group in x
  for(a in 1:nrow(x)){
    #' get the age groups in y that overlap with the age group in x
    ages_that_overlap = y[(from + fivesec) <= x[a, to] & (to - fivesec) >= x[a, from]]
    
    #' get duration of each age group
    dur_x = x[a, to - from]
    dur_y = ages_that_overlap[, to - from]
    
    if(identical(dur_x, dur_y)){
      #' all are equal, value can be copied
      for(val in value.var){
        x[a, val] = ages_that_overlap[, get(val)] 
      }
    } else if(sum(dur_y) == dur_x){
      #' sum is equal, value can be summed/averaged
      if(method == "mean"){
        for(val in value.var){
          x[a, val] = weighted.mean(ages_that_overlap[, get(val)], w = as.numeric(dur_y))
        }
      } else if(method == "sum"){
        for(val in value.var){
          x[a, val] = sum(ages_that_overlap[, get(val)])
        }
      } else {
        stop(sprintf("Method %s not yet implemented", method))
      }
    } else if(sum(dur_y) > dur_x){
      #' can only use part of one age-group
      #' age groups that are not first or last will always fully contribute
      value_rel = rep(1, ages_that_overlap[, .N])
      value_rel[c(1, ages_that_overlap[, .N])] = NA
      if(ages_that_overlap[1, from] < x[a, from]){
        #' first y age group starts before x age group
        value_rel[1] = (ages_that_overlap[1, to - from] - (x[a, from] - ages_that_overlap[1, from]))/ages_that_overlap[1, to - from]
      } else {
        value_rel[1] = 1
      }
      if(ages_that_overlap[.N, to] > x[a, to]) {
        value_rel[ages_that_overlap[, .N]] = (ages_that_overlap[.N, to - from] - (ages_that_overlap[.N, to] - x[a, to]))/ages_that_overlap[.N, to - from]
      } else if(ages_that_overlap[, .N] > 1) {
        #' if not larger than 1, already processed in the previous step
        value_rel[ages_that_overlap[, .N]] = 1
      }
      
      if(method == "mean"){
        for(val in value.var){
          x[a, val] = weighted.mean(ages_that_overlap[, get(val)], w = value_rel)
        }
      } else if(method == "sum"){
        for(val in value.var){
          x[a, val] = sum(ages_that_overlap[, get(val)] * value_rel)
        }
      } else {
        stop(sprintf("Method %s not yet implemented", method))
      }
    } else if(sum(dur_y) < dur_x){
      #' don't have enough data to fill age group
      stop(sprintf("Not enough data to fill age group %s. Do you need to fillAgeGaps on y?", x[a, name]))
    } else {
      stop(sprintf("Unexpected error when filling value for age group %s", x[a, name]))
    }
  }
  
  return(x)
}

x = c(0, 2, 4) %>% set_units("year") %>% setAgeBreaks() %>% .[, value := 1:.N] %>% .[]
y = c(0:10) %>% set_units("year") %>% setAgeBreaks() %>% .[, value := 1:.N] %>% .[]
combineAgeBreaks(x, y, value.var = "value")
combineAgeBreaks(x, y, value.var = "value", method = "sum")
combineAgeBreaks(y, x, value.var = "value")
combineAgeBreaks(y, x, value.var = "value", method = "sum")

#agelt = function(a, b){
#  onesec = set_units(1, "second")
#  return((a - b) < -onesec)
#}
#agegt = function(a, b){
#  onesec = set_units(1, "second")
#  return((a - b) > onesec)
#}

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


