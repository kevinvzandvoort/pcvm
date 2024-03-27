#' Nb OLD PCVm version - make sure to checkout correct branch
#' This script is to be called by another script, that sets up PCVM

#' Source PCVM_functions
#' - Note that these are in addition to project-specific functions specified in the parent
source(sprintf("%s/inc/age_breaks_functions.R", PCVM_FOLDER))
source(sprintf("%s/functions.R", PCVM_FOLDER))

#' Create some folders to store compiled models in directory of parent
for(x in c("model", "model/build"))
  if(!dir.exists(x)) dir.create(x)

#' Compile model and store in directory of parent
if(!file.exists(sprintf("./model/build/pcvm%s%s", ifelse(PCVM_VERSION == 2, "2", ""), .Platform$dynlib.ext)) | PCVM_COMPILE)
  compileModel(sprintf("%s/model/pcvm%s.cpp", PCVM_FOLDER, ifelse(PCVM_VERSION == 2, "2", "")), "./model/build/")

#' Load compiled model
dyn.load(sprintf("./model/build/pcvm%s%s", ifelse(PCVM_VERSION == 2, "2", ""), .Platform$dynlib.ext))
if(is.loaded("derivs", ifelse(PCVM_VERSION == 2, "pcvm2", "pcvm"))){
  message("PCVm succesfully loaded") 
} else {
  stop ("PCVm is not loaded")
}
