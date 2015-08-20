.onAttach <- function(lib, pkg)  {
  packageStartupMessage("This is a compilation of functions from the Soil Security Laboratory of the University of Sydney. Version ", utils::packageDescription("soilsecuritylab", field="Version"), appendLF = TRUE)
}
