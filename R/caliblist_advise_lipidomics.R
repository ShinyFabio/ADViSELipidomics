#' caliblist_advise_lipidomics
#'
#' @description \code{caliblist_advise_lipidomics} is the function for creation of the calibration
#' data (concentration and area).
#'
#' @param out List. It is the result from the \code{filter_advise_lipidomics} function.
#' @param calibration_path Character string. It should be the folder in which calibration files
#' are located, already stored in the list 'out' (as out$folders$calibration_path).
#' @param calibration_targetfile Character string. It is the calibration target file.
#' @param info a string that tells if it's the deuterated or nonlabeled calibration
#'
#' @return A calibration matrix
#'
#' @import dplyr
#' @import shiny
#' @importFrom readxl read_xlsx
#' @importFrom readr read_delim
#' @importFrom tibble as_tibble rownames_to_column
#' @importFrom tidyr pivot_wider
#' @importFrom data.table rbindlist
#' @importFrom crayon red bgYellow
#'
#' @export
#'
#' @details The first part of the calibration step creates a list of calibration information
#' from the calibration target files. This first part should be applied both to deuterated
#' files and to non-labeled files. The creation of the final calibration matrix is performed
#' outside the function. The calibration matrix has the correspondence between concentration
#' and area for each internal standard lipid species.
#'
#' @note Last change 20/01/2022

caliblist_advise_lipidomics <- function(out,
                                        info,               ###AGGIUNTO info
                                        calibration_path,
                                        calibration_targetfile
                                        ){

  message("---> CALIBRATION LIST ADVISE-LIPIDOMICS PIPELINE START <---")


  #Reading calibration target file
  message("Reading ", info, " calibration files...")
  if(shiny::isRunning()){
    showNotification(tagList(icon("cogs"), HTML("&nbsp;Reading", "<b>", info, "</b>", "calibration files...")), type = "default")
  }
  
  
  file_list <- unlist(strsplit(calibration_targetfile$Name, split = ";"))
  
  if (!all(file_list %in% list.files(calibration_path))){
    
    message("At least one file is missing or reported with the wrong name!")
    if(shiny::isRunning()){
      showNotification(tagList(icon("times-circle"), HTML("&nbsp;At least one file is missing or reported with the wrong name! Check the console.")), 
                       type = "error", duration = 7)
    }
    
    #check the separator
    if(!all(grepl(";",file_list))){
      message("Semi-colon not detected in Name column. In case you have replicates, remember to use ';' to separate the different files.")
      if(shiny::isRunning()){
        showNotification(tagList(icon("times-circle"), HTML("&nbsp;Semi-colon not detected in Name column. In case you have replicates, remember to use ';' to separate the different files.")), 
                         type = "warning", duration = 7)
      }
    }
    
    #check the txt extension
    if(!all(grepl(".txt", file_list))){
      message("File extension '.txt' not found in Name column. Did you forget to add it? LipidSearch outputs should be a txt file.")
      if(shiny::isRunning()){
        showNotification(tagList(icon("times-circle"), HTML("&nbsp;File extension '.txt' not found in Name column. Did you forget to add it? LipidSearch outputs should be a txt file.")), 
                         type = "warning", duration = 7)
      }
    }

    filemissing = file_list[which(!file_list %in% list.files(datapath))]
    message("Number of file missing or with wrong name: ",length(filemissing))
    message("Check the following files:")
    print(filemissing)
    stop("Calibration file missing!")
  }

  #check class column
  if(!any(grepl(";",calibration_targetfile$Class))){
    message("Semi-colon not detected in Class column. Remember to use ';' if you want to use more than one lipid class for each concentration.")
    if(shiny::isRunning()){
      showNotification(tagList(icon("times-circle"), HTML("&nbsp;Semi-colon not detected in Class column. Remember to use ';' if you want to use more than one lipid class for each concentration.")), 
                       type = "warning", duration = 7)
    }
  }
  
  calibration = calibration_targetfile
  calibration$Class = strsplit(calibration$Class, split = ";")
  calibration$Name = strsplit(calibration$Name, split = ";")
  calibration_list = list()
  
  # Creating calibration list
  message("Creating calibration list...")
  if(shiny::isRunning()){
    showNotification(tagList(icon("cogs"), HTML("&nbsp;Creating calibration list...")), type = "default")
  }
  
  nas_list = list()
  rt_range <- out$targets$internal_standard %>% dplyr::select(Class,Ion,MinRt,MaxRt,InternalStandardLipidIon)

  # #check files
  # calib_files_tocheck = calibration$Name
  # if(!all(unlist(calibration$Name) %in% list.files(calibration_path))){
  #   if(shiny::isRunning()){
  #     showNotification(tagList(icon("times-circle"), HTML("&nbsp;Calibration file missing!")), type = "error")
  #   }
  #   stop("Calibration file missing!")
  # }
  
  
  ###check lipidsearch version
  colnames_4.2 = c("LipidIon","Class","FattyAcid","Ion","ObsMz","TopRT","Area")
  colnames_5.0 = c("LipidID", "ClassKey", "Adduct", "ObsMz" ,"TopRT","Area")
  files_to_check = paste0(calibration_path, unlist(calibration$Name))
  if(is.null(check_lipidversion(file_paths = unique(files_to_check), colnames_4.2 = colnames_4.2, colnames_5.0 = colnames_5.0))){
    return(NULL)
  }
  


  for (k in 1:nrow(calibration)){

    conc_list <- list()
    conc_temp <- calibration[k,]
    our_file <- unlist(conc_temp$Name)


    for(kk in our_file){
      conc_list[[kk]] <- readr::read_delim(paste0(calibration_path,kk), "\t",
                                           escape_double = FALSE, trim_ws = TRUE, skip = 5, n_max = 1000000)
      
      colnames_tocheck =  colnames(conc_list[[kk]])
      if(all(colnames_4.2 %in% colnames_tocheck)){
        conc_list[[kk]] <- conc_list[[kk]] %>% dplyr::select(dplyr::all_of(colnames_4.2))
      }else if(all(colnames_5.0 %in% colnames_tocheck)){
        conc_list[[kk]] <- conc_list[[kk]] %>% dplyr::select(dplyr::all_of(colnames_5.0))
        #add column FattyAcid
        conc_list[[kk]] = tibble::add_column(conc_list[[kk]], FattyAcid =  paste0("(",stringr::str_extract(string = conc_list[[kk]]$LipidID, pattern = "(?<=\\()(.*?)(?=\\))"),")"), .after = 2)
        colnames(conc_list[[kk]]) <- colnames_4.2
      }else{
        return(0)
      }
      
      conc_list[[kk]] <- conc_list[[kk]] %>% 
        dplyr::mutate(Area = if(inherits(conc_list[[kk]]$Area, "character")) {readr::parse_number(Area)}else{Area}) %>%
        dplyr::mutate(Area = replace(Area, Area <= 0 , NA)) %>%
        dplyr::filter(Class %in% calibration$Class[[k]]) %>%
        dplyr::left_join(rt_range, by = c("Class","Ion")) %>%
        dplyr::filter(TopRT >= MinRt & TopRT <= MaxRt) %>%
        dplyr::filter(LipidIon %in% out$targets$internal_standard$InternalStandardLipidIon) %>%
        dplyr::select(LipidIon,Area)
    }
    
    conc_list <- conc_list[sapply(conc_list, function(x) dim(x)[1]) > 0]
    
    
    ###Removing duplicated lipids (very important for LipidSearch 5.0)
    fun_f <- function(f){
      if(!tibble::is_tibble(f)){
        return(f)
      }
      dup <- f$LipidIon[duplicated(f$LipidIon)] %>% unique()
      not_dup = setdiff(f$LipidIon, dup)
      if(length(dup) != 0){
        all_dup = list()
        max_dup = list()
        for(k in 1:length(dup)){
          all_dup[[k]] <- f[f$LipidIon == dup[k],]
          max_dup[[k]] <- all_dup[[k]] %>% dplyr::slice_max(Area,n = 1,with_ties = FALSE)
        }
        max_dup_df <- as.data.frame(do.call(rbind, max_dup))
        tibble::as_tibble(rbind(max_dup_df, dplyr::filter(f, LipidIon %in% not_dup)))
      }else{
        f
      }
    }
    #Filter duplicated lipids. Only lipids with maximum area will be picked. Convert to tibble or problems later.
    conc_list = lapply(conc_list, fun_f)
    


    aux_fun <- function(a,b) {merge(a, b, by = "LipidIon", all.x = TRUE, all.y = TRUE)}
    
    if (length(conc_list) != 0){
      
      nas_list[[k]] = data.frame("NAs" = unlist(lapply(conc_list, function(x) sum(is.na(x$Area)))))

      mean_df <- base::Reduce(aux_fun, conc_list)
      mean_df <- cbind(mean_df$LipidIon,rowMeans(mean_df[,-1], na.rm = TRUE))
      colnames(mean_df) <- c("LipidIon",calibration$`Concentration (ng/mL)`[k])
      calibration_list[[k]] <- tibble::as_tibble(mean_df)
    }
  }
  
  if(length(calibration_list) == 0){
    message(paste("Calibration list for",info,"is empty"))
    if(shiny::isRunning()){
      showNotification(tagList(icon("exclamation-circle"), HTML("&nbsp;Calibration list for",info,"is empty")), type = "warning")
    }
    return(NULL)
  }
  
  nas_list = base::Filter(Negate(is.null), nas_list)
  df_nas = lapply(nas_list, function(x) tibble::rownames_to_column(x, "Sample")) %>% data.table::rbindlist() %>% dplyr::filter(NAs != 0)
  if(sum(df_nas$NAs) > 0){
    if(shiny::isRunning()){
      showNotification(duration = 8, tagList(icon("exclamation-circle"), 
                                           HTML("&nbsp;Some lipid areas (",sum(nas$NAs),"in total) are equal or less than 0 and will be replaced with NA. 
                                    Check the console to see where NAs are introduced.")), type = "warning")
    }
    cat(crayon::bgYellow(crayon::red("NAs are introduced in the following samples.")))
    print(df_nas)
  }

  
  #check if there are not supported files (returned 0)
  check_too_new = lapply(calibration_list, function(x) if(length(x) == 1) x == 0) %>% unlist()
  if(!is.null(check_too_new)){
    message("Version of LipidSearch not supported. If your version is newer than 5.0
            please open an issue on the GitHub page. Probably they changed some column names.")
    if(shiny::isRunning()){
      showNotification(duration = 8, tagList(icon("times-circle"), HTML("&nbsp;Version of LipidSearch not supported. If your version is newer than 5.0
            please open an issue on our GitHub page. Probably they changed some column names.")), type = "error")
    }
    return(NULL)
  }
  

  
  
  calibration_mat <- base::Filter(Negate(is.null), calibration_list) %>% 
    purrr::reduce(dplyr::full_join,by = "LipidIon") %>% 
    tidyr::gather(Concentration, Area, -c(LipidIon))
  calibration_mat$Concentration <- sub("*\\.[a-zA-Z]","",calibration_mat$Concentration)
  calibration_mat <- calibration_mat[complete.cases(calibration_mat),]
  calibration_mat <- tidyr::pivot_wider(calibration_mat, id_cols = LipidIon, 
                                 names_from = Concentration, values_from = Area)
  
  message("---> CALIBRATION LIST ADVISE-LIPIDOMICS PIPELINE END <---")
  
  
  
  if(shiny::isRunning()){
    showNotification(tagList(icon("check"), HTML("&nbsp;Calibration list successfully created!")), type = "message")
  }
  
  return(calibration_mat)
}
