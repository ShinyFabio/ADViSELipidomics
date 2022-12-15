#' check_lipidversion
#'
#' @description \code{check_lipidversion} Run a check on LipidSearch version.
#'
#' @param file_paths vector containing all the path files to be checked. 
#' @param colnames_4.2 Colnames used in the 4.2.29 LipidSearch version. By default c("LipidIon","Class","FattyAcid","Ion","ObsMz","TopRT","Area").
#' @param colnames_5.0 Colnames used in the 4.2.29 LipidSearch version. By default c("LipidID", "ClassKey", "Adduct", "ObsMz" ,"TopRT","Area").
#'
#'
#' @import dplyr
#' @import shiny
#' @importFrom readr read_delim
#' @importFrom utils compareVersion
#'
#'
#' @details Run a check on LipidSearch version. If there is a problem, returns NULL. LipidSearch older than 4.2.29 are not supported.
#' Multiple LipidSearch version returns NULL.
#'
#' @note Last change 05/07/2022


check_lipidversion = function(file_paths,
                              colnames_4.2 = c("LipidIon","Class","FattyAcid","Ion","ObsMz","TopRT","Area"), 
                              colnames_5.0 = c("LipidID", "ClassKey", "Adduct", "ObsMz" ,"TopRT","Area")
                              ){
  #### Check versions
  message("Checking LipidSearch version...")
  if(shiny::isRunning()){
    showNotification(tagList(icon("cogs"), HTML("&nbsp;Checking LipidSearch version...")), type = "default")
  }
  
  colnames_4.2 = c("LipidIon","Class","FattyAcid","Ion","ObsMz","TopRT","Area")
  colnames_5.0 = c("LipidID", "ClassKey", "Adduct", "ObsMz" ,"TopRT","Area")
  
  versions <- lapply(file_paths, function(x) {
    version = readr::read_delim(x, " ", escape_double = FALSE, trim_ws = TRUE, n_max = 1, col_names = FALSE, show_col_types = FALSE)
    version %>% dplyr::pull(ncol(version))
  }) %>% unlist() %>% unique()
  
  message("Version of LipidSearch found: ", paste0(versions, collapse = ", "))
  if(shiny::isRunning()){
    showNotification(duration = 6, type = "default",
                     tagList(icon("info"), HTML("&nbsp;Version of LipidSearch found: ", paste0(versions, collapse = ", "))))
  }
  
  if(length(versions) > 1){
    message("Data files belonging to different versions of LipidSearch.")
    if(shiny::isRunning()){
      showNotification( type = "warning", tagList(icon("exclamation-circle"), 
                                                  HTML("&nbsp;Data files belonging to different versions of LipidSearch.")))
    }
    return(NULL)
  }
  
  #check if there are version less than 4.2.29
  comparing_4.2 = lapply(versions, function(x) utils::compareVersion(x, "4.2.29")) %>% unlist()
  if(-1 %in% comparing_4.2){
    message("LipidSearch version not supported. ADViSELipidomics supports a LipidSearch version 4.2.29 or greater. 
            Please use a newer version of LipidSearch.")
    if(shiny::isRunning()){
      showNotification(tagList(icon("times-circle"), HTML("&nbsp;LipidSearch version not supported! ADViSELipidomics supports LipidSearch version 4.2.29 or greater. 
            Please use a newer version.")), type = "error")
    }
    return(NULL)
  }
  
  
 return(1) 
}
