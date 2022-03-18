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
#' @param batch_alg Character string. Different methodologies to cope with the selection of
#' batch_type = "remove". It can be "limma", "combat_param_full" or "combat_nonparam_full".
#' Default = "limma".
#' @param eb_robust Logical value for inner use of empirical Bayes moderation of the
#' standard errors towards a global value. More information in "limma' package
#' (TRUE for RNAseq data). Default = FALSE.
#' @param eb_trend Logical value for inner use of empirical Bayes moderation of the
#' standard errors towards a global value. More information in "limma' package
#' (TRUE for RNAseq data). Default = FALSE.
#' @param tag_name Character string. A string to tag all the output files for the ongoing
#' analysis from this function. Default = "analysis".
#' @param thresh Numerical value. Threshold on adjust P-value in order to consider a lipid
#' differentially expressed or not. Default = 0.05.
#' @param lfc Numerical value. Threshold on logFC values in the plots. Default = 2.
#' @param col_de Character vector. It defines the correspondence between differential expressed
#' lipids (-1 down-regulated, 1 up-regulated) and colors in the plots.
#' Default = c("-1" = "green","0" = "black","1" = "red").
#' @param decide_met Character string specifying how genes and contrasts are to be combined
#' in the multiple testing scheme. Choices are "separate", "global", "hierarchical" or "nestedF".
#' More information in "limma' package. Default = "separate".
#' @param interactive Logical value. If TRUE, main plots are interactive by the use of "plotly"
#' package. Default = TRUE.
#' @param de_table_path Character string for the folder in which storing tables resulting from
#' the differential analysis. Default = "DETable".
#' @param de_plot_path Character string for the folder in which storing tables resulting from
#' the differential analysis. Default = "DEPlot".
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
                                      #tag_name = "analysis",
                                      thresh = 0.05,
                                      #lfc = 2,
                                      #col_de = c("-1" = "green","0" = "black","1" = "red"),
                                      decide_met = "separate"#,
                                      #interactive = TRUE,
                                      #de_table_path = "DETable",
                                      #de_plot_path = "DEPlot"
                                      ){
  
  # if (!dir.exists(paste0(out$folders$output_path_tables,de_table_path))){
  #   dir.create(file.path(paste0(out$folders$output_path_tables,de_table_path)))
  # }
  # 
  # if (!dir.exists(paste0(out$folders$output_path_plots,de_plot_path))){
  #   dir.create(file.path(paste0(out$folders$output_path_plots,de_plot_path)))
  # }
  
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
    ###AGGIUNTA
    showNotification(tagList(icon("cogs"), HTML("&nbsp;Normalization between replicates or samples...")), type = "default")
    ###FINE AGGIUNTA
    se_data_log <- limma::normalizeBetweenArrays(se_data_log, method = bs_norm)
  }
  
  # Removing batch effect
  
  # Limma methods
  if(batch_type == "remove" & batch_method == "limma"){
    if (length(batch_vars) == 1){
      message(paste0("Removing one batch effect: ", batch_vars, " ..."))
      ###AGGIUNTA
      showNotification(tagList(icon("cogs"), HTML(paste0("&nbsp;Removing one batch effect: ", batch_vars, " ..."))), type = "default")
      ###FINE AGGIUNTA
      aux_batch <- colData(se_data)[, batch_vars]
      se_data_log <- removeBatchEffect(se_data_log, batch = aux_batch, design = out$exp_design)
    } else if (length(batch_vars) == 2){
      message(paste0("Removing two batch effect: ", batch_vars, " ..."))
      ###AGGIUNTA
      showNotification(tagList(icon("cogs"), HTML(paste0("&nbsp;Removing two batch effect: ", batch_vars, " ..."))), type = "default")
      ###FINE AGGIUNTA
      aux_batch_1 <- colData(se_data)[, batch_vars[1]]
      aux_batch_2 <- colData(se_data)[, batch_vars[2]]
      se_data_log <- removeBatchEffect(se_data_log, batch = aux_batch_1, batch2 = aux_batch_2,
                                       design = out$exp_design)
    } else if (length(batch_vars) > 2){
      ###AGGIUNTA
      showNotification(tagList(icon("times-circle"), HTML("&nbsp;Limma method for removing batch effects can cope at most with two variables.")), type = "error")
      ###FINE AGGIUNTA
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
  ###AGGIUNTA
  showNotification(tagList(icon("cogs"), HTML("&nbsp;Fitting model...")), type = "default")
  ###FINE AGGIUNTA
  
  if( (rep_mean == TRUE) | (rep_mean == FALSE & rep_effect == FALSE) ){
    
    message("Fitting standard-method Limma model (considering samples) ...")
    ###AGGIUNTA
    showNotification(tagList(icon("cogs"), HTML("&nbsp;Fitting standard-method Limma model (considering samples) ...")), type = "default")
    ###FINE AGGIUNTA
    
    fit_a <- limma::lmFit(se_data_log,out$exp_design)
    fit_b <- limma::contrasts.fit(fit_a, out$contrast_matrix)
    fit_c <- limma::eBayes(fit_b, robust = eb_robust, trend = eb_trend)
    
  } else if(rep_mean == FALSE & rep_effect == TRUE) {
    
    message("Fitting replicates-method Limma model (replicates correlation) ...")
    ###AGGIUNTA
    showNotification(tagList(icon("cogs"), HTML("&nbsp;Fitting replicates-method Limma model (replicates correlation) ...")), type = "default")
    ###FINE AGGIUNTA
    
    tec_rep <- as.factor(sapply(strsplit(colData(se_data)$SampleID,"_"), `[`, 1))
    cor_fit <- limma::duplicateCorrelation(se_data_log, out$exp_design, block = tec_rep)
    fit_a <- limma::lmFit(se_data_log, out$exp_design, block = tec_rep, cor = cor_fit$consensus)
    fit_b <- limma::contrasts.fit(fit_a, out$contrast_matrix)
    fit_c <- limma::eBayes(fit_b, robust = eb_robust, trend = eb_trend)
    
  }
  
  ####COMMENTATO
  # message("Writing raw ebayes table...")
  # 
  # write.table(fit_c,
  #             file = paste0(out$folders$output_path_tables, de_table_path, "\\raw_ebayes_",
  #                           tag_name, "_", eb_robust, "_", eb_trend, ".txt"),
  #             quote = FALSE, sep = "\t", row.names = FALSE, dec = ".")
  # 
  # message("Writing toptables and plotting graphs...")
  
  output_limma <- list()
  for(k in colnames(out$contrast_matrix)){
    
    # Differential expressed lipids toptable
    output <- limma::topTable(fit_c, adjust.method = "BH", coef = k, number = Inf, sort.by = "none")
    
    #output$DE <- as.integer(output$adj.P.Val < thresh) * sign(output$logFC)
    
    
    output = output %>% dplyr::mutate(DE = dplyr::case_when(adj.P.Val < thresh ~ 1, TRUE ~ 0)) %>% 
      dplyr::mutate(DE = DE * sign(logFC))
    
    

    output$Class <- SummarizedExperiment::rowData(se_data)$Class
    output$Lipids <- SummarizedExperiment::rowData(se_data)$Lipids
    
    output_limma[[k]] <- output
    
    # write.table(output,
    #             file = paste0(out$folders$output_path_tables, de_table_path, "\\toptable_",
    #                           tag_name, "_", k, "_", eb_robust, "_", eb_trend, "_", thresh, ".txt"),
    #             quote = FALSE, sep = "\t", row.names = FALSE, dec = ".")
    # 
    # # Differential expressed lipids ma-plot and volcano plot
    # ma_plot <- ggplot2::ggplot(output, aes(x = AveExpr, y = logFC,
    #                                        col = as.factor(DE), label = Lipids)) +
    #   geom_point() +
    #   scale_color_manual(name = "DE Lipids",
    #                      labels = c("DW-reg", "Unchanged", "UP-reg"),
    #                      values = col_de) +
    #   geom_hline(yintercept = 0) +
    #   geom_hline(yintercept = lfc, linetype = "dotted", size = 1) +
    #   geom_hline(yintercept = -lfc, linetype = "dotted", size = 1) +
    #   ggtitle(paste("MA Plot:",k))
    # 
    # ma_plot_2 <- ggplot2::ggplot(output, aes(x = AveExpr, y = logFC,
    #                                          col = as.factor(Class), label = Lipids)) +
    #   geom_point(aes(shape = as.factor(DE))) +
    #   geom_hline(yintercept = 0) +
    #   geom_hline(yintercept = lfc, linetype = "dotted", size = 1) +
    #   geom_hline(yintercept = -lfc, linetype = "dotted", size = 1) +
    #   ggtitle(paste("MA Plot:",k)) +
    #   labs(shape = "DE", color = "Class")
    # 
    # volcano_plot <- ggplot2::ggplot(output, aes(x = logFC, y = -log10(adj.P.Val),
    #                                             col = as.factor(DE), label = Lipids)) +
    #   geom_point() +
    #   scale_color_manual(name = "DE Lipids",
    #                      labels = c("DW-reg", "Unchanged", "UP-reg"),
    #                      values = col_de) +
    #   geom_vline(xintercept = 0) +
    #   geom_vline(xintercept = lfc, linetype = "dotted", size = 1) +
    #   geom_vline(xintercept = -lfc, linetype = "dotted", size = 1) +
    #   geom_hline(yintercept = -log10(thresh), linetype = "dotted", size = 1) +
    #   ggtitle(paste("Volcano Plot:",k))
    # 
    # volcano_plot_2 <- ggplot2::ggplot(output, aes(x = logFC, y = -log10(adj.P.Val),
    #                                               col = as.factor(Class), label = Lipids)) +
    #   geom_point(aes(shape = as.factor(DE))) +
    #   geom_vline(xintercept = 0) +
    #   geom_vline(xintercept = lfc, linetype = "dotted", size = 1) +
    #   geom_vline(xintercept = -lfc, linetype = "dotted", size = 1) +
    #   geom_hline(yintercept = -log10(thresh), linetype = "dotted", size = 1) +
    #   ggtitle(paste("Volcano Plot:",k)) +
    #   labs(shape = "DE", colour = "Class")
    # 
    # ggplot2::ggsave(paste0(out$folders$output_path_plots,de_plot_path,"\\MA_plot_", k, "_",
    #                        eb_robust, "_", eb_trend, "_", thresh, "_", tag_name, ".pdf"),
    #                 plot = ma_plot, device = "pdf", width = 7, height = 7)
    # 
    # ggplot2::ggsave(paste0(out$folders$output_path_plots,de_plot_path,"\\MA_plot_2_", k, "_",
    #                        eb_robust, "_", eb_trend, "_", thresh, "_", tag_name, ".pdf"),
    #                 plot = ma_plot_2, device = "pdf", width = 7, height = 7)
    # 
    # ggplot2::ggsave(paste0(out$folders$output_path_plots,de_plot_path,"\\Volcano_plot_", k, "_",
    #                        eb_robust, "_", eb_trend, "_", thresh, "_", tag_name, ".pdf"),
    #                 plot = volcano_plot, device = "pdf", width = 7, height = 7)
    # 
    # ggplot2::ggsave(paste0(out$folders$output_path_plots,de_plot_path,"\\Volcano_plot_2_", k, "_",
    #                        eb_robust, "_", eb_trend, "_", thresh, "_", tag_name, ".pdf"),
    #                 plot = volcano_plot_2, device = "pdf", width = 7, height = 7)
    # 
    # if(interactive == TRUE){
    #   ma <- plotly::ggplotly(ma)
    #   htmlwidgets::saveWidget(plotly::as_widget(ma),
    #                           paste0(out$folders$output_path_plots,de_plot_path,"\\MA_plot_", k, "_",
    #                                  eb_robust, "_", eb_trend, "_", thresh, "_", tag_name, ".html"))
    #   
    #   ma_2 <- ma_plot_2 + ggplot2::theme(legend.title = element_blank())
    #   ma_2 <- plotly::ggplotly(ma_2) %>%
    #     plotly::layout(legend = list(title = list(text = "DE, Class")))
    #   htmlwidgets::saveWidget(plotly::as_widget(ma_2),
    #                           paste0(out$folders$output_path_plots,de_plot_path,"\\MA_plot_2_", k, "_",
    #                                  eb_robust, "_", eb_trend, "_", thresh, "_", tag_name, ".html"))
    #   
    #   vo <- plotly::ggplotly(volcano_plot)
    #   htmlwidgets::saveWidget(plotly::as_widget(vo),
    #                           paste0(out$folders$output_path_plots,de_plot_path,"\\Volcano_plot_", k, "_",
    #                                  eb_robust, "_", eb_trend, "_", thresh, "_", tag_name, ".html"))
    #   
    #   vo_2 <- volcano_plot_2 + ggplot2::theme(legend.title = element_blank())
    #   vo_2 <- plotly::ggplotly(vo_2) %>%
    #     plotly::layout(legend = list(title = list(text = "DE, Class")))
    #   htmlwidgets::saveWidget(plotly::as_widget(vo_2),
    #                           paste0(out$folders$output_path_plots,de_plot_path,"\\Volcano_plot_2_", k, "_",
    #                                  eb_robust, "_", eb_trend, "_", thresh, "_", tag_name, ".html"))
    # }
    # 
  }
  
  # Differential expressed lipids barplot
  output_2 <- limma::decideTests(fit_c, method = decide_met, adjust.method = "BH", p.value = thresh)
  # ndiff <- colSums(abs(output_2))
  # message("Number of DE lipids per comparison:")
  # print(knitr::kable(ndiff, format = "rst"))
  # 
  # ndiff_long <- data.frame(contrast = colnames(out$contrast_matrix), nDE = ndiff)
  # message("Number of UP/DW regulated lipids per comparison:")
  # print(head(summary(output_2)))
  # 
  # numde_plot <- ggplot2::ggplot(ndiff_long, aes(x = contrast, y = nDE)) +
  #   geom_bar(stat = "identity", width = 0.1) +
  #   geom_text(aes(label = nDE), vjust = 2.5, color = "SteelBlue", size = 3.5) +
  #   ggtitle(paste0("Number of DE lipids ", tag_name)) +
  #   coord_flip()
  # 
  # aux_df <- as.data.frame(summary(output_2))
  # aux_df <- aux_df[which(aux_df$Var1 != "NotSig"),]
  # colnames(aux_df) <- c("DE", "Comparison", "NumberOfLipids")
  # 
  # numde2_plot <- ggplot2::ggplot(aux_df, aes(x = Comparison, y = NumberOfLipids, fill = DE)) +
  #   geom_bar(stat = "identity", position = "dodge", width = 0.1) +
  #   scale_color_manual(name = "Gene Status",
  #                      labels = c("DW regulated", "UP regulated"),
  #                      values = c("green","red")) +
  #   geom_text(aes(label = NumberOfLipids), vjust = 2.5, color = "SteelBlue", size = 3.5) +
  #   ggtitle(paste0("Number of UP/DW regulated lipids ", tag_name)) +
  #   coord_flip()
  # 
  # ggplot2::ggsave(paste0(out$folders$output_path_plots,de_plot_path,"\\DE_barplot_",
  #                        tag_name, ".pdf"),
  #                 plot = numde_plot, device = "pdf", width = 7, height = 7)
  # 
  # ggplot2::ggsave(paste0(out$folders$output_path_plots,de_plot_path,"\\DE_barplot_2_",
  #                        tag_name, ".pdf"),
  #                 plot = numde2_plot, device = "pdf", width = 7, height = 7)
  
  # Updating
  aux_out <- list(test_result = output_2, limma_result = output_limma)
  
  message("---> DIFFERENTIAL EXPRESSION ADVISE-LIPIDOMICS PIPELINE END <---")
  
  res <- purrr::flatten(list(out,aux_out))
  
  ####AGGIUNTA 
  showNotification(tagList(icon("check"), HTML("&nbsp;Differential expression completed!")), type = "message")
  return(res)
  ##### FINE AGGIUNTA
  
}