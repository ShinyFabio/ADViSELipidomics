
.onAttach <- function(libname, pkgname){
  
  version = utils::packageDescription(pkgname, fields = "Version")
  
  msg = paste0(
    "
   |----------------------------------------------------------------------------|
   |                WELCOME TO ", pkgname, " (Ver ", version,")                     |
   |----------------------------------------------------------------------------|
   |             This package has been developed at IAC-CNR/ICB-CNR             |
   |         under the financial support of the Regione Campania project:       |
   |        Piattaforma tecnologica per la lotta alle patologie oncologiche     |
   |             Antitumor Drugs and Vaccines from the SEa (ADViSE)             |
   |----------------------------------------------------------------------------|
                                               * Package version ",version," loaded *"
  )
  
  
  packageStartupMessage(msg)
  
}

