#' comparison_advise_lipidomics
#'
#' @description \code{comparison_advise_lipidomics}
#'
#' @param out List. It is the result from the \code{diffexp_advise_lipidomics} function.
#' @param gradient Logical value. If TRUE, the Venn Diagram is colored taking into account
#' the number of common lipids by a gradient. Default = FALSE.
#'
#' @return Venn diagram and Upset plot of common lipids among different comparisons
#' (comparison of all the selected contrasts).
#'
#' @import dplyr
#' @import shiny
#' @importFrom SummarizedExperiment colData assay
#' @import limma
#' @importFrom purrr flatten
#' @import statmod
#' @importFrom tibble rownames_to_column
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette
#' @importFrom ComplexHeatmap make_comb_mat UpSet
#' @import ggVennDiagram
#' @rawNamespace import(ggplot2, except = last_plot)
#'
#' @export
#'
#' @details The Venn Diagramm can be in classical form or in gradient form.
#'
#' @note Last change 17/12/2021
#' 



comparison_advise_lipidomics <- function(out,
                                         gradient = FALSE){
  
  message("---> COMPARISON PLOTS ADVISE-LIPIDOMICS PIPELINE START <---")
  
  # Differential expressed lipids Venn diagram
  if(ncol(out$test_result) > 2 & ncol(out$test_result) < 6){
    
    # out_list <- list()
    # for(k in 1:ncol(out$test_result)){
    #   aux_col <- as.data.frame(out$test_result[,k])
    #   idx <- which(as.data.frame(aux_col) != 0)
    #   aux_col[idx,] <- rownames(aux_col)[idx]
    #   out_list[[k]] <- aux_col[aux_col!="0"]
    # }
    
    aux_con <- tibble::rownames_to_column(as.data.frame(out$test_result)) %>%
      dplyr::mutate_all(list(~dplyr::case_when(. != 0 ~ rowname,  TRUE ~ "0"))) %>%
      dplyr::select(-rowname)
    aux_con <- lapply(as.list(aux_con), function(x) x[x != "0"])
    venn <- ggVennDiagram::Venn(aux_con)
    data_venn <- ggVennDiagram::process_data(venn)
    
    # Table of Venn Diagram
    table_venn <- stringi::stri_list2matrix(data_venn@region$item, byrow = FALSE)
    colnames(table_venn) <-  gsub("\\.+","\U2229",data_venn@region$name)
    
    if(gradient == TRUE){
      
      data_venn@setEdge$id <- rep("3",length(data_venn@setEdge$id))
      vd_plot <- ggVennDiagram::ggVennDiagram(aux_con,
                                              category.names = colnames(out$test_result),
                                              label_percent_digit = 1,
                                              label_alpha = 0,
                                              label_color = "black",
                                              edge_size = 0) +
        ggplot2::geom_sf(aes(color = id), data = ggVennDiagram::venn_setedge(data_venn), show.legend = FALSE) +
        ggplot2::scale_fill_gradient(low = "white", high = "blue") +
        ggplot2::scale_x_continuous(expand = expansion(mult = .2))
      
      
      # ggplot2::ggsave(paste0(out$folders$output_path_plots,de_plot_path,"\\venn_diagram_gradient_",
      #                        tag_name, ".pdf"),
      #                 plot = vd_plot, device = "pdf", width = 7, height = 7)
      
    } else {
      
      region_label <- data_venn@region %>%
        dplyr::filter(.data$component == "region") %>%
        dplyr::mutate(percent = paste(round(.data$count*100/sum(.data$count),
                                            digits = 1),"%", sep="")) %>%
        dplyr::mutate(both = paste(.data$count,paste0("(",.data$percent,")"),sep = "\n"))
      
      vd_plot <- ggplot() +
        geom_sf(aes(fill = id), data = venn_region(data_venn), show.legend = FALSE) +
        geom_sf(color = "black", size = 1, data = venn_setedge(data_venn), show.legend = FALSE) +
        geom_sf_text(aes(label = name), data = venn_setlabel(data_venn)) +
        geom_sf_text(aes(label = both), data = region_label) +
        scale_x_continuous(expand = expansion(mult = .2)) +
        theme_void()
      
      # ggplot2::ggsave(paste0(out$folders$output_path_plots,de_plot_path,"\\venn_diagram_color_",
      #                        tag_name, ".pdf"),
      #                 plot = vd_plot, device = "pdf", width = 7, height = 7)
      
    }
    
  } else {
    
    message("Venn diagram for DE lipids is plotted from two to five contrasts!")
    
  }
  
  
  # Differential expressed lipids upset plot
  aux_color <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(12,"Paired"))
  comb_mat <- ComplexHeatmap::make_comb_mat(out$test_result)
  
  #pdf(file = paste0(out$folders$output_path_plots,de_plot_path,"\\upset_plot_", tag_name, ".pdf"))
  upset_plot <- ComplexHeatmap::UpSet(comb_mat, comb_col = aux_color(length(comb_degree(comb_mat))))
  plot(upset_plot)
  #dev.off()
  
  
  message("---> COMPARISON PLOTS ADVISE-LIPIDOMICS PIPELINE END <---")
  
}