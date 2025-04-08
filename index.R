#' This script is to be called by another script, that sets up PCVM

#' Source METAVAX functions
#' - Note that these are in addition to project-specific functions specified in the parent
source(sprintf("%s/inc/age_breaks_functions.R", METAVAX_FOLDER))
source(sprintf("%s/functions.R", METAVAX_FOLDER))

if(!exists("METAVAX_COMPILE")) METAVAX_COMPILE = FALSE

#' Create some folders to store compiled models in directory of parent
for(x in c("model", "model/build"))
  if(!dir.exists(x)) dir.create(x)

#' Compile model and store in directory of parent
if(!file.exists(sprintf("./model/build/%s%s", MODEL_NAME, .Platform$dynlib.ext)) | METAVAX_COMPILE)
  compileModel(sprintf("%s/model/%s.cpp", METAVAX_FOLDER, MODEL_NAME), "./model/build/")

#' Load compiled model
dyn.load(sprintf("./model/build/%s%s", MODEL_NAME, .Platform$dynlib.ext))
if(is.loaded("derivs", MODEL_NAME)){
  message(sprintf("Model '%s' succesfully loaded", MODEL_NAME))
} else {
  stop ("Model is not loaded")
}
