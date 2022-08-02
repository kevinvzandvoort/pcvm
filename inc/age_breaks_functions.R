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
        ##last_name = paste0(age_breaks_years[age_breaks_types == label][max(which(is_seq))] %>%
        ##                     set_units(unit, mode = "standard"), label, "-",
        ##                            (age_breaks_years[which(age_breaks_types == label) %>% rev %>% .[1] + 1] %>%
        ##                     set_units(unit, mode = "standard") - set_units(1, unit, mode = "standard")),
        ##                   age_breaks_types[which(age_breaks_types == label) %>% rev %>% .[1] + 1])
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

#' Test age breaks
#' with b_length = a_start
#' with b_length = 2*a_start
#' with b_length = a_start/2
#
#a_start = 4
#a_length = 2
#a_end = a_start+a_length
#
#b_length = a_length*2
#b_start = a_start - b_length - 1
#out = list()
#i = 1
#
#tol = 0.05
#
#gt = function(a, b){
#  (a - b) > tol
#}
#lt = function(a, b){
#  (a - b) < tol
#}
#
#setPlotOption = function(data){
#  b_start = data$start
#  b_end = data$end
#  #if(gt(a_start, (b_end-tol*2))){
#  #  msg = "B fully before A"
#  #} else if(lt(a_end-tol*2, b_start)){
#  #  msg = "B fully after A"
#  #} else if(lt(a_start, b_end-tol*2)) {
#  if(!gt(a_start, b_end-tol*2) & !lt(a_end-tol*2, b_start) & lt(a_start, b_end-tol*2)){
#    msg = "B partly/fully in A"
#  } else {
#    msg = "Not Yet Implemented"
#  }
#  data[, message := msg]
#  return(data)
#}
#while(b_start <= a_end){
#  out[[i]] = data.table(y = i, start = b_start, end = b_start + b_length)
#  out[[i]] = setPlotOption(out[[i]])
#  b_start = b_start + 1
#  i = i+1
#}
#out[[i]] = data.table(y = i, start = b_start, end = b_start + b_length)
#out[[i]] = setPlotOption(out[[i]])
#ggplot(data = NULL, aes(y=y, yend=y, x=start, xend=end))+
#  geom_segment(data = data.table(y = 0, start = a_start, end = a_end), colour="#000000", size = 1)+
#  geom_point(data = data.table(y = 0, start = a_start, end = a_end), colour="#000000", fill="#000000", size = 2, shape = 21)+
#  geom_point(data = data.table(y = 0, start = a_start, end = a_end), aes(x=end), colour="#000000", fill="#FFFFFF", size = 2, shape = 21)+
#  geom_segment(data = rbindlist(out), size = 1, aes(colour=message))+
#  geom_point(data = rbindlist(out), size = 2, shape = 21, aes(colour=message, fill=message))+
#  geom_point(data = rbindlist(out), aes(x=end, colour=message), fill="#FFFFFF", size = 2, shape = 21)

combineAgeBreaks = function(target, additional, method = c("mean", "sum")[1], value.var = "value", additional_group = NULL){
  onesec = set_units(1, "second")
  
  lt = function(a, b){
    (a - b) < onesec
  }
  gt = function(a, b){
    (a - b) > onesec
  }
  
  #  if(any(class(target[, from], target[, to], additional[, from], additional[, to]) != "units")) stop("Need to input units")
  
  if(!is.null(additional_group)){
    groupings = additional[, additional_group, with=FALSE] %>% unique()
    ngroups = nrow(groupings)
  } else {
    ngroups = 1
  }
  
  output = list()
  
  for(g in 1:ngroups){
    z_target = copy(target) %>% setorder(from, to)
    if(!is.null(additional_group)){
      z_additional = groupings[g] %>% merge(additional) %>% setorder(from, to)
    } else {
      z_additional = copy(additional) %>% setorder(from, to)  
    }
    
    z_target[, row := 1:.N]
    z_additional[, row := 1:.N]
    
    first_prev = NULL
    last_prev = NULL
    dur_prev = NULL
    for(a in 1:nrow(z_target)){
      #first = first(z_additional[!gt(from, z_target[a, from]) & !lt(to - onesec, z_target[a, from]), row])
      #last = last(z_additional[!gt(from, z_target[a, to] - onesec) & !lt(to - onesec, (z_target[a, to] - onesec)), row])
      ages_that_overlap = z_additional[!gt(z_target[a, from], to-onesec*2) & !lt(z_target[a, to]-onesec*2, from) & lt(z_target[a, from], to-onesec*2)]
      first = first(ages_that_overlap[, row])
      last = last(ages_that_overlap[, row])
      dur = z_target[a, to - from]
      
      if(!is.null(first_prev) & !is.null(last_prev) & !is.null(dur_prev)){
        if(first == first_prev & last == last_prev & dur == dur_prev){
          #can skip this iteration as values will be the same
          for(val in value.var){
            z_target[a, val] = z_target[a-1, get(val)]
          }
          next
        }
      }
      
      value_rows = first:last
      value_duration = z_additional[value_rows, to - from]
      value_duration[1] = min(z_additional[value_rows[1], to], z_target[a, to]) - z_target[a, from]
      if(length(value_rows) > 1) value_duration[length(value_rows)] = z_target[a, to] - z_additional[value_rows[length(value_rows)], from]
      
      if(sum(value_duration) != sum(z_target[a, to - from]))
        stop("Something went wrong - durations not equal")
      if(any(as.numeric(value_duration/z_additional[value_rows, to - from]) > 1 | as.numeric(value_duration/z_additional[value_rows, to - from]) < 0))
        stop("Something went wrong - duration greater than 1")
      
      for(val in value.var){
        if(method == "mean")
          z_target[a, val] = weighted.mean(z_additional[value_rows, get(val)], as.numeric(value_duration)) 
        else if(method == "sum")
          z_target[a, val] = sum(as.numeric(value_duration/z_additional[value_rows, to - from]) * z_additional[value_rows, get(val)])
      }
      
      first_prev = first
      last_prev = last
      dur_prev = dur
    }
    
    if(!is.null(additional_group)){
      output[[g]] = z_target[, -"row"] %>% cbind(z_additional[, additional_group, with=FALSE] %>% unique)
    } else {
      output[[g]] = z_target[, -"row"]  
    }
  }
  
  return(rbindlist(output))
}

agelt = function(a, b){
  onesec = set_units(1, "second")
  return((a - b) < -onesec)
}
agegt = function(a, b){
  onesec = set_units(1, "second")
  return((a - b) > onesec)
}

matchingAgeBreaks = function(target, additional){
  onesec = set_units(1, "second")
  lt = function(a, b){
    (a - b) < onesec
  }
  gt = function(a, b){
    (a - b) > onesec
  }
  
  z_target = copy(target) %>% setorder(from, to)
  z_additional = copy(additional) %>% setorder(from, to)  
    
  z_target[, row := 1:.N]
  z_target[, name.x := name]
  z_additional[, row := 1:.N]
    
  for(a in 1:nrow(z_target)){
    ages_that_overlap = z_additional[!gt(z_target[a, from], to-onesec*2) & !lt(z_target[a, to]-onesec*2, from) & lt(z_target[a, from], to-onesec*2)]
    first = first(ages_that_overlap[, row])
    last = last(ages_that_overlap[, row])
    #first = first(z_additional[from <= z_target[a, from] & (to - onesec) >= z_target[a, from], row])
    #last = last(z_additional[from <= (z_target[a, to] - onesec) & (to - onesec) >= (z_target[a, to] - onesec), row])
      
    if(first != last) stop("No unique age-breaks in y to match with x")
    z_target[a, name.y := z_additional[first, name]]
  }
  
  return(z_target[, c("name.x", "name.y")])
}

combineAgeBreaks2 = function(target, additional, method = c("mean", "sum")[1], value.var = "value", additional_group = NULL){
  onesec = set_units(1, "second")
  
  lt = function(a, b){
    (a - b) < onesec
  }
  gt = function(a, b){
    (a - b) > onesec
  }
  
  additional_age_groups = additional %>% .[, c("name", "from", "to")] %>% unique()
  z_target = copy(target) %>% setorder(from, to)
  z_target[, row := 1:.N]
  z_additional_age_groups = copy(additional_age_groups) %>% setorder(from, to) %>% .[, row := 1:.N]
  
  output = list()
  for(a in 1:nrow(z_target)){
    ages_that_overlap = z_additional_age_groups[!gt(z_target[a, from], to-onesec*2) & !lt(z_target[a, to]-onesec*2, from) & lt(z_target[a, from], to-onesec*2)]
    first = first(ages_that_overlap[, row])
    last = last(ages_that_overlap[, row])
    dur = z_target[a, to - from]
    
    value_rows = first:last
    value_duration = z_additional_age_groups[value_rows, to - from]
    value_duration[1] = min(z_additional_age_groups[value_rows[1], to], z_target[a, to]) - z_target[a, from]
    if(length(value_rows) > 1) value_duration[length(value_rows)] = z_target[a, to] - z_additional_age_groups[value_rows[length(value_rows)], from]
    
    if(sum(value_duration) != sum(z_target[a, to - from]))
      stop("Something went wrong - durations not equal")
    if(any(as.numeric(value_duration/z_additional_age_groups[value_rows, to - from]) > 1 | as.numeric(value_duration/z_additional_age_groups[value_rows, to - from]) < 0))
      stop("Something went wrong - duration greater than 1")
    
    if(method == "mean"){
      ages_that_overlap[, weight := value_duration]  
    } else if(method == "sum"){
      ages_that_overlap[, weight := as.numeric(value_duration/z_additional_age_groups[value_rows, to - from])]  
    }
    
    z_additional = copy(additional)
    z_additional = z_additional %>% merge(ages_that_overlap[, c("name", "weight")], by="name")
    
    z_additional = value.var %>% lapply(function(val){
      if(method == "mean")
        tmp = z_additional[, .(val = weighted.mean(get(val), weight)), by=additional_group]
      else if(method == "sum")
        tmp = z_additional[, .(val = sum(get(val) * weight)), by=additional_group]
        #z_target[a, val] = sum(as.numeric(value_duration/z_additional[value_rows, to - from]) * z_additional[value_rows, get(val)])
      colnames(tmp)[which(colnames(tmp) == "val")] = val
      return(tmp)
      }) %>% (function(x){
        tmp = x[[1]]
        if(length(x) > 1){
          for(i in 2:length(x)){
            tmp %<>% merge(x[[i]], by = additional_group)
          }
        }
        return(tmp)
      })
    
    output[[length(output) + 1]] = z_target[a, ] %>% cbind(z_additional) %>% .[, -"row"]
  }
  
  return(rbindlist(output))
}
