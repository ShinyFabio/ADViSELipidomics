#' expdesign_advise_lipidomics
#'
#' @description \code{expdesign_advise_lipidomics} is the function
#'
#' @param out List. It is the result from the \code{sumexp_advise_lipidomics} function.
#' @param design_vars Character vector, with the names of the variables from the target file,
#' considered as conditions of comparison, related to replicates/samples. Default = "Product".
#' @param batch_vars Character vector, with the names of the variables from the target file,
#' considered as batch effects for the entire experiment. Default = "".
#' @param batch_type Character string. Different methodologies to cope with the eventual presence
#' of batch effects: "none", "remove", "fit". Default = "fit".
#' @param batch_method Character string. Different methodologies to cope with the selection of
#' batch_type = "fit". It can be "given" or "estimate". Default = "given".
#' @param rep_mean Logical value. It is for the selection of the dataset with the replicates
#' (FALSE) or the dataset with the samples (TRUE). Default = FALSE.
#' @param file_contrast Character string. It is the name of the file with the list of the
#' comparisons (contrasts) of interest. Default = ""contrast_list.txt".
#'
#' @return res: a list with results from the Summarized Experiment step, updated with the
#' design matrix, the contrast matrix and the contrast vector (comparisons).
#'
#' @import dplyr
#' @import shiny
#' @importFrom readr read_csv
#' @importFrom SummarizedExperiment colData
#' @importFrom stringr str_replace_all
#' @importFrom limma is.fullrank makeContrasts
#' @importFrom knitr kable
#' @importFrom purrr flatten
#' @import statmod
#' @importFrom sva num.sv sva
#'
#' @export
#'
#' @details This function generates the experimental design for the subsequent differential
#' analysis. At first, it provides the linear model formula with the selected design variables,
#' considering the eventual presence of batch effects. Then, after the generation of the
#' experimental design, it reads a list of comparisons of interest (contrast list) and generates
#' the contrast matrix for the next steps. In case of batch_type == "fit" and batch_method == "estimate",
#' the selected methodology creates surrogate variables to be added to the design matrix.
#' Further information on removing batch effects is in 'diffexp_advise_lipidomics" function.
#'
#' @note Last change 17/12/2021
#' 

expdesign_advise_lipidomics <- function(out,
                                        design_vars = "",
                                        batch_vars = "",
                                        batch_type = "",
                                        batch_method = "",
                                        rep_mean = FALSE,
                                        file_contrast = "contrast_list_mod.txt"){
  
  message("---> EXPERIMENTAL DESIGN ADVISE-LIPIDOMICS PIPELINE START <---")
  
  # Selection of replicates/samples
  if(rep_mean == TRUE && out$replicates == TRUE){
    se_data = out$sumexp_data_mean
  } else {
    se_data = out$sumexp_data
  }
  
  message("Design matrix creation...")
  
  if(shiny::isRunning()){
    showNotification(tagList(icon("cogs"), HTML("&nbsp;Design matrix creation...")), type = "default")
  }
  
  # Considering batch effects in the model
  if(batch_type == "none"){
    message("Batch effect variables are not present.")
    message(paste0("Design matrix with formula: ~0+",
                   paste(design_vars, collapse = ":")))
    if(shiny::isRunning()){
      showNotification(tagList(icon("info"), HTML("&nbsp;Batch effect variables are not present.")), type = "default")
      showNotification(tagList(icon("info"), HTML("&nbsp;Design matrix with formula: ~0+",paste(design_vars, collapse = ":"))), type = "default")
    }
    
    tot_vars = design_vars
  }
  
  if(batch_type == "remove"){
    message("Batch effects will be 'removed' before the creation of the contrast matrix.")
    message(paste0("Design matrix with formula: ~0+",
                   paste(design_vars, collapse = ":")))
    if(shiny::isRunning()){
      showNotification(tagList(icon("info"), HTML("&nbsp;Batch effects will be 'removed' before the creation of the contrast matrix.")), type = "default")
      showNotification(tagList(icon("info"), HTML("&nbsp;Design matrix with formula: ~0+",paste(design_vars, collapse = ":"))), type = "default")
    }
    tot_vars = design_vars
  }
  
  
  if(batch_type == "fit" & batch_method == "given"){
    message(paste0("Batch effect variables are explicitly provided: ", batch_vars))
    message(paste0("Design matrix with formula: ~0+",
                   paste(design_vars, collapse = ":"), "+", paste(batch_vars, collapse = "+")))
    if(shiny::isRunning()){
      showNotification(tagList(icon("info"), HTML("&nbsp;Batch effect variables are explicitly provided: ", paste(batch_vars, collapse = ","))), type = "default")
      showNotification(tagList(icon("info"), HTML("&nbsp;Design matrix with formula: ~0+",paste(design_vars, collapse = ":"), "+", paste(batch_vars, collapse = "+"))), type = "default")
    }
    tot_vars = c(design_vars,batch_vars)
  }
  
  if(batch_type == "fit" & batch_method == "estimate"){
    message(paste0("Batch effect variables are estimated."))
    message(paste0("Design matrix with formula: ~0+",
                   paste(design_vars, collapse = ":")), " (Batch effect added in model matrix)")
    if(shiny::isRunning()){
      showNotification(tagList(icon("info"), HTML("&nbsp;Batch effect variables are estimated.")), type = "default")
      showNotification(tagList(icon("info"), HTML("&nbsp;Design matrix with formula: ~0+",paste(design_vars, collapse = ":"), " (Batch effect added in model matrix)")), type = "default")
    }
    tot_vars = design_vars
  }
  
  n_vars <- length(tot_vars)
  aux_vars <- SummarizedExperiment::colData(se_data)[,tot_vars]
  if(n_vars != 1){
    aux_vars <- lapply(aux_vars, function(x) as.factor(x))
  } else {
    aux_vars <- list(as.factor(aux_vars))
  }
  names(aux_vars) <- tot_vars
  n_level <- sapply(aux_vars, nlevels)
  if(all(n_level > 1)){
    
    if(batch_type == "none" | batch_type == "remove"){
      formula_des <- as.formula(paste("~ 0 +", paste(design_vars, collapse = ":")))
      design <- model.matrix(formula_des, aux_vars)
      design <- design[,which(colSums(design) > 0)]
    }
    
    if(batch_type == "fit" & batch_method == "given"){
      formula_des <- as.formula(paste("~ 0 +", paste(design_vars, collapse = ":"), "+",
                                      paste(batch_vars, collapse = "+")))
      design <- model.matrix(formula_des, aux_vars)
      design <- design[,which(colSums(design) > 0)]
    }
    
    if(batch_type == "fit" & batch_method == "estimate"){
      formula_des <- as.formula(paste("~ 0 +", paste(design_vars, collapse = ":")))
      design <- model.matrix(formula_des, aux_vars)
      design <- design[,which(colSums(design) > 0)]
      
      design_null <- model.matrix(~1, data = colData(se_data))
      n_sv <- sva::num.sv(as.matrix(assay(se_data)), design, method = "leek")
      sv_obj = sva::sva(as.matrix(assay(se_data)), design, design_null, n.sv = n_sv-1)
      batch_mat <- sv_obj$sv
      colnames(batch_mat) <- paste0("Surr", seq(1, ncol(batch_mat)))
      
      design <- cbind(design,batch_mat)
    }
    
    aux_reg <- paste0(paste(tot_vars, collapse = "|"),"|-")
    colnames(design) <- gsub(aux_reg,"",colnames(design)) %>%
      stringr::str_replace_all("[:]", "_") %>%
      stringr::str_replace_all("[;]", "_and_")
  } else {
    if(shiny::isRunning()){
      showNotification(tagList(icon("times-circle"), HTML("&nbsp;At least one factor has only one level -> ",paste0(names(which(!n_level > 1)), collapse = ";"))), type = "error")
    }
    stop(paste0("At least one factor has only one level -> ",
                paste0(names(which(!n_level > 1)), collapse = ";")))
  }
  
  # Check on design matrix rank
  if(!limma::is.fullrank(design)){
    if(shiny::isRunning()){
      showNotification(tagList(icon("times-circle"), HTML("&nbsp;Design matrix has not full rank: contrast matrix cannot be created!")), type = "error")
    }
    stop("Design matrix has not full rank: contrast matrix cannot be created!")
  }
  
  # Setting the list of contrasts to be investigated
  contrast_list = file_contrast
  
  message("List of contrast")
  print(contrast_list)
  contrast_vec = unlist(strsplit(contrast_list$X1,"="))[seq(from = 2,
                                                            to = 2*length(contrast_list$X1),
                                                            by = 2)]
  contrast_nam = unlist(strsplit(contrast_list$X1,"="))[seq(from = 1,
                                                            to = 2*length(contrast_list$X1),
                                                            by = 2)]
  names(contrast_vec) <- contrast_nam
  
  message("Summary of contrast vector")
  print(knitr::kable(contrast_vec, format = "rst"))
  
  contrast_mat <- limma::makeContrasts(contrasts = contrast_vec, levels = design)
  colnames(contrast_mat) <- names(contrast_vec)
  
  message("Summary of contrast matrix")
  print(contrast_mat)
  
  # Updating
  aux_out <- list(exp_design = design, contrast_matrix = contrast_mat, contrasts = contrast_vec)
  
  message("---> EXPERIMENTAL DESIGN ADVISE-LIPIDOMICS PIPELINE END <---")
  
  res <- purrr::flatten(list(out,aux_out))
  
  return(res)
}

