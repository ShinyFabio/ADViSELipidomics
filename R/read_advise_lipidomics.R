#' read_advise_lipidomics
#'
#' @description \code{read_advise_lipidomics} is the function for reading and storing data
#' produced by LipidSearch, and to be analyzed with ADViSELipidomics pipeline.
#'
#' @param out List. It is the result from the \code{init_advise_lipidomics} function.
#' @param target_file Dataframe. It should be the target file already stored in the list 'out'
#' (as out$targets$targetfile_lipidomics).
#' @param data_type Character string. The data type to read. Can be "LipidSearch" or "Liquid".
#' By default "LipidSearch".
#'
#' @return res: a list with result from the initialization step, updated with the data to be
#' analyzed with the forthcoming steps.
#'
#' @import dplyr
#' @import shiny
#' @importFrom readr read_delim parse_number
#' @importFrom tools file_path_sans_ext
#' @importFrom purrr flatten set_names map
#' @importFrom stringr str_count
#' @importFrom crayon red bgYellow
#'
#' @export
#'
#' @details The reading step reads and stores samples data (created by LipidSearch software),
#' checking if all the files listed in the target file are present in the suitable folder.
#'
#' @note Last change 17/12/2021

read_advise_lipidomics <- function(out,
                                   datapath,    ###AGGIUNTO datapath
                                   data_type = "LipidSearch",   ###AGGIUNTO data_type
                                   target_file){

  message("---> READING AND STORING DATA ADVISE-LIPIDOMICS PIPELINE START <---")

  file_list <- unlist(strsplit(target_file$File_name, split = ";"))

  ###AGGIUNTA
  datapath = paste0(datapath, "/")
  showNotification(tagList(icon("info"), HTML("&nbsp;Number of files to be imported: ", length(file_list))), type = "default")
  ###FINE AGGIUNTA
  
  message(paste0("Number of files to be imported: ", length(file_list)))
  
  message(paste0("Name of files to be imported: ", "\n"))
  for (k in 1:length(file_list)){
    message(paste0(file_list[k]))
  }

  if (!all(file_list %in% list.files(datapath))){      ###MODIFICA list.files(datapath) anzichè list.files(paste0(input_path,"Data\\")
    message("At least one file is missing or reported with the wrong name!")
    
    ###AGGIUNTA
    showNotification(tagList(icon("times-circle"), HTML("&nbsp;At least one file is missing or reported with the wrong name!")), type = "error")
    ###FINE AGGIUNTA
    
  } else {
    
    ###AGGIUNTA
    showNotification(tagList(icon("cogs"), HTML("&nbsp;Reading and storing data...")), type = "default")
    percentage <- 0
    ###FINE AGGIUNTA
    
    aux_name <- paste0(datapath, file_list)         ###MODIFICA paste0(datapath, file_list) anzichè paste0(input_path,"Data\\",file_list)
    if(data_type == "LipidSearch"){          #####aggiunta if
      ###AGGIUNTA
      raw_list <- lapply(aux_name, function(x) {
        percentage <<- percentage + 1/length(file_list)*100
        incProgress(1/length(file_list), detail = paste0("Progress: ",round(percentage,0), " %"))
        temp_file = readr::read_delim(x, "\t", escape_double = FALSE, trim_ws = TRUE, skip = 5, n_max = 1000000)
        temp_file %>% dplyr::mutate(Area = if(class(temp_file$Area) == "character") {readr::parse_number(Area)}else{Area}) %>%
          dplyr::mutate(Area = replace(Area, Area <= 0 , NA))
        })
      
      ###FINE AGGIUNTA
      
      #prima era così:
      #raw_list <- lapply(aux_name, function(x) readr::read_delim(x, "\t", escape_double = FALSE, trim_ws = TRUE, skip = 5))
      
      data_list <- lapply(raw_list, function(x) x %>%
                            dplyr::select(LipidIon,Class,FattyAcid,Ion,ObsMz,TopRT,Area))
      names(raw_list) <- tools::file_path_sans_ext(file_list)
      names(data_list) <- tools::file_path_sans_ext(file_list)
      
      ####AGGIUNTA
      nas = data.frame("NAs" = unlist(lapply(raw_list, function(x) sum(is.na(x$Area))))) %>% dplyr::filter(NAs != 0)
      if(sum(nas$NAs) > 0){
        showNotification(duration = 8, tagList(icon("exclamation-circle"), 
                                 HTML("&nbsp;Some lipid areas (",sum(nas$NAs),"in total) are equal or less than 0 and will be replaced with NA. 
                                      Check the console to see where NAs are introduced.")), type = "warning")
        cat(crayon::bgYellow(crayon::red("NAs are introduced in the following samples.")))
        print(nas)
      }
      ###FINE AGGIUNTA
      
    } ### fine aggiunta if
    
    
    if(data_type == "Liquid"){   #####inizio aggiunta if
      raw_list <- lapply(aux_name, function(x) {
        percentage <<- percentage + 1/length(file_list)*100
        incProgress(1/length(file_list), detail = paste0("Progress: ",round(percentage,0), " %"))
        readr::read_delim(x, "\t", escape_double = FALSE, trim_ws = TRUE, n_max = 1000000)})
      
      names(raw_list) <- tools::file_path_sans_ext(file_list)
      
      data_list <- lapply(raw_list, function(x) x %>% 
                            dplyr::select('Common Name','Adduct','Apex RT','Intensity') %>% 
                            dplyr::rename_with(.fn = ~gsub(" ","_",.), .cols = everything())
      )
      
      #ora devo unire i positivi e i negativi
      rem_names = gsub("_negative","",names(data_list))
      rem_names = gsub("_positive","", rem_names)
      names(data_list) <- rem_names
      data_list = purrr::map(purrr::set_names(unique(names(data_list))), ~Reduce(rbind, data_list[names(data_list)==.]))
    }  ####fine aggiunta if

  }

  
  message("Checking replicates...")
  
  ###AGGIUNTA
  showNotification(tagList(icon("cogs"), HTML("&nbsp;Checking replicates...")), type = "default")
  ###FINE AGGIUNTA
  
  aux_sam <- target_file$SampleID #MODIFICA prima eraout$targets$targetfile_lipidomics$SampleID
  aux_cou <- stringr::str_count(aux_sam, "_")
  if(sum(aux_cou) != 0){
    message("Technical replicates are present")
    ###AGGIUNTA
    showNotification(tagList(icon("info"), HTML("&nbsp;Technical replicates are present.")), type = "default")
    ###FINE AGGIUNTA
    tec_rep = TRUE
  } else {
    message("Technical replicates are absent")
    ###AGGIUNTA
    showNotification(tagList(icon("info"), HTML("&nbsp;Technical replicates are absent.")), type = "default")
    ###FINE AGGIUNTA
    tec_rep = FALSE
  }
  
  
  # Updating
  aux_out <- list(raw_data = raw_list, lipid_data = data_list, replicates = tec_rep)

  message("---> READING AND STORING DATA ADVISE-LIPIDOMICS PIPELINE END <---")
  
  res <- purrr::flatten(list(out,aux_out))
  
  ###AGGIUNTA
  showNotification(tagList(icon("check"), HTML("&nbsp;Data successfully loaded!")), type = "message")
  return(res)
  ###FINE AGGIUNTA
  
  
}
