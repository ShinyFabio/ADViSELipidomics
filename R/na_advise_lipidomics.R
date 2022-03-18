#' na_advise_lipidomics
#'
#' @description \code{na_advise_lipidomics} is the function for coping with NA values,
#' filtering rows and columns (replicates/samples and lipid species) for a selected percentage
#' of NA values, and imputing survived NA values with different imputation methods.
#'
#' @param out List. It is the result from the \code{recovery_advise_lipidomics} function.
#' @param na_filter_lip Numeric value. It is the upper bound percentage of accepted NA values
#' for lipid species, in decimal format. Default = 0.3.
#' @param na_filter_sam Numeric value. It is the upper bound percentage of accepted NA values
#' for replicates/samples, in decimal format. Default = 0.6.
#' @param imputation_met Character string. It is the imputation method used to cope with the
#' presence of the survived NA values after the filtering. The imputation method can be: "mean",
#' "median", "knn" or "irmi". Default = "median".
#' @param imputation_val Numeric value. It is the imputation value to substitute in the place
#' of NA values when "median" or "mean" imputation method is selected, in decimal format.
#' Default = 0.001.
#'
#' @return res: a list with results from recovery step, updated with the imputed concentration
#' matrix and the list of the filtered out samples.
#'
#' @import dplyr
#' @importFrom VIM kNN irmi
#' @importFrom purrr flatten
#' @importFrom tidyr replace_na separate
#' @importFrom tibble rownames_to_column
#'
#' @export
#'
#' @details The imputation step is composed of two substeps:
#' a) the NA values are filtered per replicates/samples and lipid species, taking into account
#' two percentage values (upper bounds) selected by the user;
#' b) the survived NA values are imputed with a method selected by the user.
#' More in details, four imputation methods are available:
#'   - "median", the NA values are imputed with the median from related replicates
#'     (for a sample values) per lipid species, and survived NA are corrected with
#'     the imputation value;
#'   - "mean", the NA values are imputed with the mean from related replicates
#'     (for a sample values) per lipid species, and survived NA are corrected with
#'     the imputation value;
#'   - "knn", k-Nearest Neighbour imputation based on a variation of the
#'     Gower Distance, applied on replicates/samples;
#'   - "irmi", Iterative Robust Model-based Imputation, where in each step of the
#'     iteration, one variable is used as a response variable, and the remaining variables
#'     serve as the regressors, applied on replicates/samples.
#' The list of the filtered samples is stored with the results.
#'
#' @note Last change 17/12/2021

na_advise_lipidomics <- function(out,
                                 na_filter_lip = 0.3,
                                 na_filter_sam = 0.6,
                                 imputation_met = "median",
                                 imputation_val = 0.001){

  message("---> NA FILTERING AND IMPUTATION ADVISE-LIPIDOMICS PIPELINE START <---")
  

  aux_conc_mat <- out$concentration_matrix[,-1]
  rownames(aux_conc_mat) <- out$concentration_matrix[,1]

  if (na_filter_lip != 0){

    dim_conc_mat <- dim(aux_conc_mat)

    message("Filtering NA by lipid...")
    all_na_lipid <- aux_conc_mat[rowSums(is.na(aux_conc_mat)) == dim_conc_mat[2],]
    message(paste0("Number of lipids totally absent: ", nrow(all_na_lipid)))
    message(paste0("Name of lipids totally absent: ", "\n"))
    message(paste0(rownames(all_na_lipid), "\n"))
    
    if(na_filter_lip != 1 && shiny::isRunning()){
      showNotification(tagList(icon("cogs"), HTML("&nbsp;Filtering NA by lipid...")), type = "default")
      showNotification(tagList(icon("info"), HTML("&nbsp;Number of lipids totally absent: ", nrow(all_na_lipid))), type = "default")
    }
    
    aux_conc_mat$na_perc_lip = apply(aux_conc_mat, 1, function(x) sum(is.na(x))/dim_conc_mat[2])
    aux_conc_mat <- aux_conc_mat %>% dplyr::filter(na_perc_lip < na_filter_lip)
    aux_conc_mat <- aux_conc_mat[,-ncol(aux_conc_mat)]

  }

  if (na_filter_sam != 0){

    aux_conc_mat <- as.data.frame(t(aux_conc_mat))
    dim_conc_mat_2 <- dim(aux_conc_mat)

    message("Filtering NA by sample...")
    all_na_sample <- aux_conc_mat[rowSums(is.na(aux_conc_mat)) == dim_conc_mat_2[2],]
    message(paste0("Number of samples totally absent: ", nrow(all_na_sample)))
    message(paste0("Name of samples totally absent: ", "\n"))
    message(paste0(rownames(all_na_sample), "\n"))

    if(na_filter_sam != 1 && shiny::isRunning()){
      showNotification(tagList(icon("cogs"), HTML("&nbsp;Filtering NA by sample...")), type = "default")
      showNotification(tagList(icon("info"), HTML("&nbsp;Number of samples totally absent: ", nrow(all_na_sample))), type = "default")
    }
    
    aux_conc_mat$na_perc_sam = apply(aux_conc_mat, 1, function(x) sum(is.na(x))/dim_conc_mat_2[2])
    sample_filt <- aux_conc_mat %>%
      dplyr::filter(na_perc_sam > na_filter_sam) %>%
      rownames()
    aux_conc_mat <- aux_conc_mat %>% dplyr::filter(na_perc_sam < na_filter_sam)
    aux_conc_mat <- aux_conc_mat[,-ncol(aux_conc_mat)]
    aux_conc_mat <- as.data.frame(t(aux_conc_mat))

  }


  if(imputation_met != "none"){
    message("Imputing NA values...")
    if(shiny::isRunning()){
      showNotification(tagList(icon("cogs"), HTML("&nbsp;Imputing NA values...")), type = "default")
    }
  }

  if(imputation_met == "knn"){
    imp_conc_mat <-  VIM::kNN(t(aux_conc_mat), imp_var = FALSE)
    imp_conc_mat <- as.data.frame(t(imp_conc_mat))
    colnames(imp_conc_mat) <- colnames(aux_conc_mat)
  }
  
  # Model-based imputation
  if(imputation_met == "irmi"){
    imp_conc_mat <-  VIM::irmi(t(aux_conc_mat), imp_var = FALSE)
    imp_conc_mat <- as.data.frame(t(imp_conc_mat))
    colnames(imp_conc_mat) <- colnames(aux_conc_mat)
  }
  
  if(imputation_met == "mean" | imputation_met == "median"){
    na_conc_mat <- as.data.frame(t(aux_conc_mat)) %>%
      tibble::rownames_to_column("Sample") %>%
      tidyr::separate(Sample, into = c("Sample", "rep"), sep = "_", fill = "right") %>%
      dplyr::select(-rep)
    na_conc_list <- split(na_conc_mat,
                          factor(na_conc_mat$Sample,levels = unique(na_conc_mat$Sample)))
    
    na_conc_list <-  lapply(na_conc_list, function(x) x %>% dplyr::select(-Sample))

    #---Subfunction---#
    fun_h <- function(h){
      if(imputation_met == "mean"){
        col_means <- lapply(h, mean, na.rm = TRUE)
        hh <- tidyr::replace_na(h, col_means)
        hh[is.na(hh)] <- imputation_val
        return(hh)
      }
      if(imputation_met == "median"){
        col_medians <- lapply(h, median, na.rm = TRUE)
        hh <- tidyr::replace_na(h, col_medians)
        hh[is.na(hh)] <- imputation_val
        return(hh)
      }
    }
    
    na_conc_list <- lapply(na_conc_list,fun_h)
    imp_conc_mat <- as.data.frame(t(dplyr::bind_rows(na_conc_list))) %>%
      dplyr::mutate_all(as.numeric)
    colnames(imp_conc_mat) <- colnames(aux_conc_mat)
  }
  
  if(imputation_met == "none"){
    imp_conc_mat <- aux_conc_mat
  }

  # Updating
  imp_conc_mat <- cbind(Lipids = rownames(imp_conc_mat), imp_conc_mat)
  imp_conc_mat = imp_conc_mat %>% dplyr::arrange(Lipids)
  
  
  if(na_filter_sam != 0){
    aux_out <- list(concentration_matrix_filt = imp_conc_mat,
                    sample_filtered = sample_filt)
  } else {
    aux_out <- list(concentration_matrix_filt = imp_conc_mat,
                    sample_filtered = "")
  }
  
  message("---> NA FILTERING AND IMPUTATION ADVISE-LIPIDOMICS PIPELINE END <---")

  res <- purrr::flatten(list(out,aux_out))
  
  if(shiny::isRunning()){
    if(imputation_met == "none" && na_filter_sam == 1 && na_filter_lip == 1){
      showNotification(tagList(icon("info"), HTML("&nbsp;No NA imputation or filtering performed.")), type = "default")
    }else{
      showNotification(tagList(icon("check"), HTML("&nbsp;NA filtering and imputation done!")), type = "message")
    }
  }
  
  return(res)

}
