#' This script is to be called by another script, that sets up PCVM

#' Source PCVM_functions
#' - Note that these are in addition to project-specific functions specified in the parent
source(sprintf("%s/inc/age_breaks_functions.R", PCVM_FOLDER))
source(sprintf("%s/functions.R", PCVM_FOLDER))

#' Create some folders to store compiled models in directory of parent
for(x in c("model", "model/build"))
  if(!dir.exists(x)) dir.create(x)

#' Compile model and store in directory of parent
if(!file.exists(sprintf("./model/build/pcvm%s", .Platform$dynlib.ext)) | PCVM_COMPILE)
  compileModel(sprintf("%s/model/pcvm.cpp", PCVM_FOLDER), "./model/build/")

#' Load compiled model
dyn.load(sprintf("./model/build/pcvm%s", .Platform$dynlib.ext))
if(is.loaded("derivs", "pcvm")){
  message("PCVm succesfully loaded") 
} else {
  stop ("PCVm is not loaded")
}