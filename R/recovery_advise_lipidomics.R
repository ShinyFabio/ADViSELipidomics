#' recovery_advise_lipidomics
#'
#' @description \code{recovery_advise_lipidomics} is the function for the computation of the
#' recovery percentages, and their application to define the concentration matrix (i.e. feature
#' matrix)
#'
#' @param out List. It is the result from the \code{calibplot_advise_lipidomics} function.
#' @param intercept_flag Logical value. If set to TRUE the intercept will be set to zero.
#'
#' @return res: a list with results from calibration steps, updated with the recovery percentages
#' and the concentration matrix (feature matrix).
#'
#' @import shiny
#' @import dplyr
#' @importFrom purrr flatten
#'
#' @export
#'
#' @details The recovery step is the computation of the recovery percentage for each internal
#' standard lipid species, and the application of the recovery percentage on all the filter
#' lipid species (non-labeled species). The complete algorithm of the computation and the
#' application of the recovery percentage is composed of two main steps:
#' a) the recovery percentage is calculated by the formula (q intercept, m slope,
#' NominalStdConcentration is a value of concentration from the internal standard file):
#'    - (Area/m)*100/NominalStdConcentration for intercept as zero;
#'    - (Area-q)/m)*100/NominalStdConcentration for intercept as not zero;
#' b) the recovery percentage is applied to the lipid species by the formula (ConcExpPerc
#' is the recovery percentage):
#'    - (Area/m)*100.00/ConcExpPerc for intercept as zero;
#'    - (Area-q)/m)*100.00/ConcExpPerc for intercept as not zero.
#'
#'
#' @note Last change 14/01/2022

recovery_advise_lipidomics <- function(out,
                                       intercept_flag = TRUE
                                       ){

  message("---> RECOVERY PERCENTAGE ADVISE-LIPIDOMICS PIPELINE START <---")

  aux_mat <- cbind(out$calibration_matrix,out$coefficients)
  aux_mat$InternalStandardLipidIon <- rownames(out$calibration_matrix)

  message("Calculating recovery percentage...")
  
  if(shiny::isRunning()){
    showNotification(tagList(icon("cogs"), HTML("&nbsp;Calculating recovery percentage...")), type = "default")
  }

  #---Subfunction---#
  fun_g <- function(g){
    if(intercept_flag == FALSE){
      gg <- g %>%
        dplyr::left_join(aux_mat,by = c("LipidIon" = "InternalStandardLipidIon")) %>%
        dplyr::left_join(out$targets$internal_standard,by = c("LipidIon" = "InternalStandardLipidIon")) %>%
        dplyr::select(Class,Ion,LipidIon,Area,m,NominalStdConcentration) %>%
        dplyr::mutate(ConcExpPerc = (Area/m)*100/NominalStdConcentration)
      return(gg)
    } else {
      gg <- g %>%
        dplyr::left_join(aux_mat,by = c("LipidIon" = "InternalStandardLipidIon")) %>%
        dplyr::left_join(out$targets$internal_standard,by = c("LipidIon" = "InternalStandardLipidIon")) %>%
        dplyr::select(Class,Ion,LipidIon,Area,m,q,NominalStdConcentration) %>%
        dplyr::mutate(ConcExpPerc = ((Area-q)/m)*100/NominalStdConcentration)
      return(gg)
    }
  }

  aux_rec_perc <- lapply(out$lipid_data, function(x) x %>%
                           dplyr::select(LipidIon,Area) %>% fun_g %>%
                           dplyr::filter(!is.na(ConcExpPerc)))


  #check if there are both deuterated and nonlabeled
  if(sum(grepl("deuterated", names(aux_rec_perc)), na.rm = TRUE) == sum(grepl("nonlabeled", names(aux_rec_perc)), na.rm = TRUE)){
    step_for = 2
    #if one of them doesn't contain nonlabeled OR deuterated
  }else if(sum(grepl("deuterated", names(aux_rec_perc)), na.rm = TRUE) == 0 || sum(grepl("nonlabeled", names(aux_rec_perc)), na.rm = TRUE) == 0){
    step_for = 1
  }else{
    stop("Something wrong in the recovery percentage. Error 88")
  }

  
  rec_perc <- list()
  for(k in seq(1,length(aux_rec_perc),by = step_for)){
    
    #if step_for == 1 just copy aux_rec_perc to recperc without merging deuterated and nonlabeled (since exists only one)
    if(step_for == 2){
      rec_perc[[k]] <- rbind(aux_rec_perc[[k]],aux_rec_perc[[k+1]])
    }else{
      rec_perc[[k]] = aux_rec_perc[[k]]
    }
    
    aux_code <- unlist(strsplit(names(aux_rec_perc)[k],split = "_"))[c(1,3)]
    if(!is.na(aux_code[2])){
      names(rec_perc)[k] <- paste0(aux_code[1],"_",aux_code[2])
    } else {
      names(rec_perc)[k] <- aux_code[1]
    }
  }
  
  if(step_for == 2){
    rec_perc <- rec_perc[-(seq(1,length(rec_perc),by = 2)+1)]
  }
  
  rec_perc <- lapply(rec_perc,function(x) x %>%
                       dplyr::select(Class,LipidIon,ConcExpPerc) %>%
                       dplyr::rename(InternalStandardLipidIon = LipidIon) %>%
                       dplyr::mutate(ConcExpPerc = as.numeric(format(round(ConcExpPerc,2),nsmall = 2))) %>%
                       dplyr::mutate(ConcExpPerc = dplyr::case_when(ConcExpPerc > 100.00 ~ 100.00, TRUE ~ ConcExpPerc)))


  message("Correction by recovery percentage...")
  
  if(shiny::isRunning()){
    showNotification(tagList(icon("cogs"), HTML("&nbsp;Correction by recovery percentage...")), type = "default")
  }
  
  
  aux_conc_mat <- lapply(out$lipid_filtered, function(x) x %>%
                           dplyr::select(Class,LipidIon,InternalStandardLipidIon,Area) %>%
                           dplyr::left_join(out$coefficients,by = "InternalStandardLipidIon"))
  

  for (k in 1:length(aux_conc_mat)){
    aux_conc_mat[[k]] <- merge(aux_conc_mat[[k]],rec_perc[[k]], all.x = TRUE)
    aux_code <- unlist(strsplit(names(aux_conc_mat)[k],split = "_"))[c(1,3)]
    if(!is.na(aux_code[2])){
      names(aux_conc_mat)[k] <- paste0(aux_code[1],"_",aux_code[2])
    } else {
      names(aux_conc_mat)[k] <- aux_code[1]
    }
  }



  if(intercept_flag == FALSE){
    aux_conc_mat <- lapply(aux_conc_mat,function(x) {x %>%
                             dplyr::mutate(Conc = (Area/m)) %>% 
                             dplyr::left_join(dplyr::distinct(dplyr::select(out$targets$internal_standard, MinLinearity,	MaxLinearity,InternalStandardLipidIon)), by = "InternalStandardLipidIon") %>% 
                             dplyr::mutate(Conc = dplyr::case_when(Conc < MinLinearity | Conc > MaxLinearity ~ 0, TRUE ~ Conc)) %>% 
                             dplyr::mutate(Conc = replace(Conc, Conc == 0 , NA)) %>%
                             dplyr::mutate(ConcExpPerc = replace(ConcExpPerc, ConcExpPerc == 0 , NA)) %>%
                             dplyr::mutate(Conc = Conc*100.00/ConcExpPerc)
      })
    
  } else {
    aux_conc_mat <- lapply(aux_conc_mat,function(x) x %>%
                             dplyr::mutate(Conc = ((Area-q)/m)) %>% 
                             dplyr::left_join(dplyr::distinct(dplyr::select(out$targets$internal_standard, MinLinearity,	MaxLinearity,InternalStandardLipidIon)), by = "InternalStandardLipidIon") %>% 
                             dplyr::mutate(Conc = dplyr::case_when(Conc < MinLinearity | Conc > MaxLinearity ~ 0, TRUE ~ Conc)) %>% 
                             dplyr::mutate(Conc = replace(Conc, Conc == 0 , NA)) %>%
                             dplyr::mutate(ConcExpPerc = replace(ConcExpPerc, ConcExpPerc == 0 , NA)) %>%
                             dplyr::mutate(Conc = Conc*100.00/ConcExpPerc))
  }
  
  

  conc_mat <- lapply(aux_conc_mat, function(x) x %>% dplyr::select(LipidIon,Conc))
  aux_fun <- function(a,b){ base::merge(a, b, by = "LipidIon", all.x = TRUE, all.y = TRUE) }
  conc_mat <- base::Reduce(aux_fun, conc_mat)
  colnames(conc_mat) <- c("Lipids",names(aux_conc_mat))
  
  

  
  #normalization factor
  if("Norm_factor" %in% colnames(out$targets$targetfile_lipidomics)){
    if(0 %in% out$targets$targetfile_lipidomics$Norm_factor){
      message("There's a 0 inside the Norm_factor column. Check the values.")
      if(shiny::isRunning()){
        showNotification(tagList(icon("times-circle"), HTML("&nbsp;There's a 0 inside the Norm_factor column. Check the values.")), type = "error")
      }
      return(NULL)
    }
    
    if(TRUE %in% is.na(out$targets$targetfile_lipidomics$Norm_factor)){
      message("One or more NA in the Norm_factor column. Normalization will not be performed.")
      if(shiny::isRunning()){
        showNotification(tagList(icon("info"), 
                                 HTML("&nbsp;One or more NA in the Norm_factor column. Normalization will not be performed.")), type = "default")
      }
    }else{
      for(i in colnames(conc_mat[,-1])){
        norm_fact = out$targets$targetfile_lipidomics %>% dplyr::filter(SampleID == i)
        conc_mat[,i] <- conc_mat[,i] / norm_fact$Norm_factor
      }
      message("Applied normalization.")
      if(shiny::isRunning()){
        showNotification(tagList(icon("cogs"), HTML("&nbsp;Applied normalization with the formula new_conc = old_conc / norm_fact")), type = "default")
      }
    }
    
  }else{
    message("Column 'Norm_factor' not present in the target file. Normalization will not be performed.")
    if(shiny::isRunning()){
      showNotification(tagList(icon("info"), HTML("&nbsp;Column 'Norm_factor' not present in the target file. Normalization will not be performed.")), type = "default")
    }
  }

  ##table with lipid out of linearity
  list_not_lin = lapply(aux_conc_mat, function(x) x %>% dplyr::filter(is.na(Conc)) %>% dplyr::select(LipidIon, Conc) %>% 
                        dplyr::mutate(Conc = replace(Conc, is.na(Conc), "out")))
  table_not_lin <- Reduce(aux_fun, list_not_lin)
  colnames(table_not_lin) <- c("Lipids", names(list_not_lin))
  

  # Updating
  aux_out <- list(concentration_matrix = conc_mat, recovery_percentage = rec_perc, table_not_lin = table_not_lin)

  message("---> RECOVERY PERCENTAGE ADVISE-LIPIDOMICS PIPELINE END <---")

  res <- purrr::flatten(list(out,aux_out))
  
  if(shiny::isRunning()){
    showNotification(tagList(icon("check"), HTML("&nbsp;Recovery percentage for each internal standard lipid species calculated!")), type = "message")
  }
  
  return(res)

}
