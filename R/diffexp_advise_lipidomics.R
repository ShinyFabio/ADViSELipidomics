#' diffexp_advise_lipidomics
#'
#' @description \code{diffexp_advise_lipidomics} is the function for creation of the experimental
#' design.
#'
#' @description \code{diffexp_advise_lipidomics} is the function for eventually removing the
#' batch effects and applying the differential analysis per lipids among replicates/samples,
#' taking into account the comparisons of interest.
#'
#' @param out List. It is the result from the \code{expdesign_advise_lipidomics} function.
#' @param rep_mean Logical value. It is for the selection of the dataset with the replicates
#' (FALSE) or the dataset with the samples (TRUE). Default = FALSE.
#' @param rep_effect Logical value. If rep_mean == FALSE, it is possible to consider the
#' influence of the replicates presence for fitting the linear model, taking into account a
#' value of correlation among the replicates. Default = TRUE.
#' @param bs_norm Character. It is possible to normalize data between replicates
#' or samples. Possible options are "none", "scale" and "quantile" (only if it isn't a concentration matrix). 
#' Default = "none".
#' @param batch_vars Character vector, with the names of the variables from the target file,
#' considered as batch effects for the entire experiment. Default = "".
#' @param batch_type Character string. Different methodologies to cope with the eventual presence
#' of batch effects: "none", "remove", "fit". Default = "fit".
#' @param eb_robust Logical value for inner use of empirical Bayes moderation of the
#' standard errors towards a global value. More information in "limma' package
#' (TRUE for RNAseq data). Default = FALSE.
#' @param eb_trend Logical value for inner use of empirical Bayes moderation of the
#' standard errors towards a global value. More information in "limma' package
#' (TRUE for RNAseq data). Default = FALSE.
#' @param thresh Numerical value. Threshold on adjust P-value in order to consider a lipid
#' differentially expressed or not. Default = 0.05.
#' @param decide_met Character string specifying how genes and contrasts are to be combined
#' in the multiple testing scheme. Choices are "separate", "global", "hierarchical" or "nestedF".
#' More information in "limma' package. Default = "separate".
#' @param batch_method Algorithm used for the batch effect. Choices are: "limma", "combat_nonparam" and "combat_param".
#'
#' @return res: a list with results from the experimental design, updated with the results
#' from fitting the linear model and appling decisional test on the results from the model.
#' Moreover, all the tables and plots for differential expressed lipids are available.
#'
#' @return res: a list with results from the experimental design, updated with the results
#' from fitting the linear model and appling decisional test on the results from the model.
#' Moreover, all the tables and plots for differential expressed lipids are available.
#'
#' @import dplyr
#' @import shiny
#' @importFrom SummarizedExperiment colData assay
#' @import limma
#' @importFrom purrr flatten
#' @import statmod
#' @importFrom sva ComBat
#' @importFrom utils write.table head
#'
#' @export
#'
#' @details This function has many tasks. At first, it log-transforms data for subsequent analysis
#' and provides a simple method for normalizing data (scaling) per replicates or samples.
#' Then, it is possible to remove batch effects by "limma" method (at most two batch variables) or
#' by ComBat method (from "sva" package), in parametric or non-parametric modes. In particular
#' the option "combat_param_full" will run ComBat with parametric adjustment and the null model
#' with design matrix, "combat_nonparam_full" will run ComBat with nonparametric adjustment and
#' the null model with design matrix (non-parametric approach can be very time consuming).
#' After this correction, the linear model is fitted and the result are reported in form of table
#' (raw tables and filtered tables) and in form of plots (MA-plots, Volcano plots, number of
#' differential expressed lipids plots), separated for all the comparisons of interest.
#'
#' @note last change 17/12/2021
#' 
diffexp_advise_lipidomics <- function(out,
                                      rep_mean = FALSE,
                                      rep_effect = TRUE,
                                      bs_norm = "none",
                                      batch_type = "remove",
                                      batch_method = "limma",
                                      batch_vars = "",
                                      eb_robust = FALSE,
                                      eb_trend = FALSE,
                                      thresh = 0.05,
                                      decide_met = "separate"
                                      ){
  

  message("---> DIFFERENTIAL EXPRESSION ADVISE-LIPIDOMICS PIPELINE START <---")
  
  # Selection of replicates/samples
  if(rep_mean == TRUE && out$replicates == TRUE){
    se_data <- out$sumexp_data_mean
    se_data_log <- log2(assay(out$sumexp_data_mean))
    rep_effect = FALSE
  } else {
    se_data <- out$sumexp_data
    se_data_log <- log2(assay(out$sumexp_data))
  }
  
  # Normalization between replicates or samples
  if(bs_norm != "none"){
    message("Normalization between replicates or samples...")
    if(shiny::isRunning()){
      showNotification(tagList(icon("cogs"), HTML("&nbsp;Normalization between replicates or samples...")), type = "default")
    }
    se_data_log <- limma::normalizeBetweenArrays(se_data_log, method = bs_norm)
  }
  
  # Removing batch effect
  
  # Limma methods
  if(batch_type == "remove" & batch_method == "limma"){
    if (length(batch_vars) == 1){
      message(paste0("Removing one batch effect: ", batch_vars, " ..."))
      if(shiny::isRunning()){
        showNotification(tagList(icon("cogs"), HTML(paste0("&nbsp;Removing one batch effect: ", batch_vars, " ..."))), type = "default")
      }
      aux_batch <- colData(se_data)[, batch_vars]
      se_data_log <- removeBatchEffect(se_data_log, batch = aux_batch, design = out$exp_design)
    } else if (length(batch_vars) == 2){
      message(paste0("Removing two batch effect: ", batch_vars, " ..."))
      if(shiny::isRunning()){
        showNotification(tagList(icon("cogs"), HTML(paste0("&nbsp;Removing two batch effect: ", batch_vars, " ..."))), type = "default")
      }
      aux_batch_1 <- colData(se_data)[, batch_vars[1]]
      aux_batch_2 <- colData(se_data)[, batch_vars[2]]
      se_data_log <- removeBatchEffect(se_data_log, batch = aux_batch_1, batch2 = aux_batch_2,
                                       design = out$exp_design)
    } else if (length(batch_vars) > 2){
      if(shiny::isRunning()){
        showNotification(tagList(icon("times-circle"), HTML("&nbsp;Limma method for removing batch effects can cope at most with two variables.")), type = "error")
      }
      stop("Limma method for removing batch effects can cope at most with two variables.")
    }
  }
  
  # ComBat methods
  if(batch_type == "remove" & batch_method == "combat_param"){
    mod_comb <- model.matrix(~1, data = SummarizedExperiment::colData(se_data))
    aux_batch <- colData(se_data)[, batch_vars]
    se_data_log <- sva::ComBat(dat = se_data_log, batch = aux_batch, mod = mod_comb,
                               par.prior = TRUE)
  } else if(batch_type == "remove" & batch_method == "combat_nonparam"){
    mod_comb <- model.matrix(~1, data = SummarizedExperiment::colData(se_data))
    aux_batch <- colData(se_data)[, batch_vars]
    se_data_log <- sva::ComBat(dat = se_data_log, batch = aux_batch, mod = mod_comb,
                               par.prior = FALSE)
  }
  
  
  # Fitting linear model
  message("Fitting model...")
  if(shiny::isRunning()){
    showNotification(tagList(icon("cogs"), HTML("&nbsp;Fitting model...")), type = "default")
  }
  
  if( (rep_mean == TRUE) | (rep_mean == FALSE & rep_effect == FALSE) ){
    
    message("Fitting standard-method Limma model (considering samples) ...")
    if(shiny::isRunning()){
      showNotification(tagList(icon("cogs"), HTML("&nbsp;Fitting standard-method Limma model (considering samples) ...")), type = "default")
    }
    
    fit_a <- limma::lmFit(se_data_log,out$exp_design)
    fit_b <- limma::contrasts.fit(fit_a, out$contrast_matrix)
    fit_c <- limma::eBayes(fit_b, robust = eb_robust, trend = eb_trend)
    
  } else if(rep_mean == FALSE & rep_effect == TRUE) {
    
    message("Fitting replicates-method Limma model (replicates correlation) ...")
    if(shiny::isRunning()){
      showNotification(tagList(icon("cogs"), HTML("&nbsp;Fitting replicates-method Limma model (replicates correlation) ...")), type = "default")
    }
    
    tec_rep <- as.factor(sapply(strsplit(colData(se_data)$SampleID,"_"), `[`, 1))
    cor_fit <- limma::duplicateCorrelation(se_data_log, out$exp_design, block = tec_rep)
    fit_a <- limma::lmFit(se_data_log, out$exp_design, block = tec_rep, cor = cor_fit$consensus)
    fit_b <- limma::contrasts.fit(fit_a, out$contrast_matrix)
    fit_c <- limma::eBayes(fit_b, robust = eb_robust, trend = eb_trend)
    
  }
  

  output_limma <- list()
  for(k in colnames(out$contrast_matrix)){
    
    # Differential expressed lipids toptable
    output <- limma::topTable(fit_c, adjust.method = "BH", coef = k, number = Inf, sort.by = "none")
    

    output = output %>% dplyr::mutate(DE = dplyr::case_when(adj.P.Val < thresh ~ 1, TRUE ~ 0)) %>% 
      dplyr::mutate(DE = DE * sign(logFC))
    
    
    output$Class <- SummarizedExperiment::rowData(se_data)$Class
    output$Lipids <- SummarizedExperiment::rowData(se_data)$Lipids
    
    output_limma[[k]] <- output

  }
  
  # Differential expressed lipids barplot
  output_2 <- limma::decideTests(fit_c, method = decide_met, adjust.method = "BH", p.value = thresh)
  
  # Updating
  aux_out <- list(test_result = output_2, limma_result = output_limma)
  
  message("---> DIFFERENTIAL EXPRESSION ADVISE-LIPIDOMICS PIPELINE END <---")
  
  res <- purrr::flatten(list(out,aux_out))
  
  if(shiny::isRunning()){
    showNotification(tagList(icon("check"), HTML("&nbsp;Differential expression completed!")), type = "message")
  }
  
  return(res)

}