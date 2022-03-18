#' expdesign_advise_lipidomics
#' @description \code{enrichment_advise_lipidomics} is the function for enrichment analysis.
#'
#' @param out The result limma_result.
#' @param rank_c Character string. Value to plot. Possible options are: "logFC" , "P.Value", "adj.P.Val", "B"
#' @param min_elem Numeric. Internal value.
#' @param thresh p-value threshold. Possible options are 0.001, 0.01, 0.05.
#' @param k The name of the contrast.
#'
#' @return write something
#'
#' @import dplyr
#' @import shiny
#' @rawNamespace import(ggplot2, except = last_plot)
#' @importFrom SummarizedExperiment rowData
#' @importFrom tidyr gather unite
#' @importFrom ggpubr ggarrange rremove text_grob annotate_figure
#' @importFrom fgsea fgsea
#'
#' @export
#'
#' @details WRITE SOMETHING
#'
#' @note Last change 17/11/2021

enrichment_advise_lipidomics <- function(out,
                                         rank_c = "logFC",
                                         min_elem = 2,
                                         thresh = 0.05,
                                         k = ""
                                         ){
  # Creation lipid sets
  lipid <- SummarizedExperiment::rowData(out$sumexp_data)
  lipid <- as.data.frame(lipid) %>%
    tidyr::gather("variable", "value", Class, total_sn, total_un) %>%
    dplyr::filter(!is.na(value)) %>%
    tidyr::unite("set", variable, value, sep = "_") %>%
    dplyr::group_by(set) %>%
    dplyr::filter(dplyr::n() >= min_elem)
  
  set <- split(as.character(lipid$Lipids), lipid$set)
  if (length(unlist(set)) == 0) {
    stop('Lacking of lipid sets!')
  }
  
  # Application of fast gene set enrichment analysis

    result <- out$limma_result[[k]] %>% dplyr::arrange(-dplyr::across(rank_c))
    stats <- as.data.frame(result)[,rank_c]
    names(stats) <- as.character(as.data.frame(result)[, "Lipids"])
    

    fun_gsea <- fgsea::fgsea(pathways = set, stats = stats, minSize = min_elem, nperm=5000000)
    fun_gsea <- fun_gsea %>%
      dplyr::arrange(padj) %>%
      dplyr::rename(set = pathway)
    
    aux_fgsea <- lipid %>%
      dplyr::left_join(fun_gsea, by = "set") %>%
      dplyr::left_join(result, by = "Lipids") %>%
      dplyr::select(Lipids,set,pval,padj,logFC,DE,Class) %>%
      dplyr::mutate(signif = case_when(padj > thresh ~ "No Sign", TRUE ~ "Sign" ))
    
    aux_class <- aux_fgsea[grep("Class",lipid$set),]
    aux_tsn <- aux_fgsea[grep("total_sn",lipid$set),]
    aux_tun <- aux_fgsea[grep("total_un",lipid$set),]
    
    sig_plot_1 <- ggplot2::ggplot(aux_class,aes(Class,logFC, fill = signif)) +
      geom_boxplot() + geom_jitter() +
      geom_hline(yintercept = 0, lty = 2) +
      coord_flip() + xlab("Class")
    sig_plot_2 <- ggplot2::ggplot(aux_tsn,aes(set,logFC, fill = signif)) +
      geom_boxplot() + geom_jitter() +
      geom_hline(yintercept = 0, lty = 2) +
      coord_flip() + xlab("Total Chain Length")
    sig_plot_3 <- ggplot2::ggplot(aux_tun,aes(set,logFC, fill = signif)) +
      geom_boxplot() + geom_jitter() +
      geom_hline(yintercept = 0, lty = 2) +
      coord_flip() + xlab("Total Unsaturation")
    sig_plot <- ggpubr::ggarrange(sig_plot_1 + ggpubr::rremove("xlab"),
                                  sig_plot_2 + ggpubr::rremove("xlab"),
                                  sig_plot_3,
                                  ncol = 1, common.legend = TRUE, legend = "right") %>%
      ggpubr::annotate_figure(top = ggpubr::text_grob(paste0("Lipid Set Enrichment Analysis: ", k)))


  return(sig_plot)
}