#' sumexp_advise_lipidomics
#'
#' @description \code{sumexp_advise_lipidomics} is the function for the creation of the
#' SummarizedExperiment object from concentration matrix, target file information, and
#' metadata (already stored in 'out' list)
#'
#' @param out List. It is the result from the \code{recovery_advise_lipidomics} function.
#'
#' @return res: a list with results from the recovery step, updated with two SummerizedExperiment
#' objects with concentration matrix, target file information, lipid features and metadata.
#'
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom purrr flatten
#' @importFrom stringr str_extract
#' @importFrom tibble column_to_rownames
#' @import shiny
#' @import dplyr
#'
#' @export
#'
#' @details The conversion step creates the SummarizeExperiment object from the concentration
#' matrix (i.e. the feature matrix), the target file data, and laboratory metadata. Moreover,
#' lipid features are parsed from the name of the lipid species. Two SummarizedExperiment objects
#' are stored: one with each single replicate, one with replicates averaged per sample.
#'
#' @note Last change 17/12/2021

sumexp_advise_lipidomics <- function(out){

  message("---> SUMMERIZED EXPERIMENT DATA ADVISE-LIPIDOMICS PIPELINE START <---")

  message("Feature matrix initialization...")
  
  if(shiny::isRunning()){
    showNotification(tagList(icon("cogs"), HTML("&nbsp;Feature matrix initialization...")), type = "default")
  }
  
  # Replicates
  aux_assay <- out$concentration_matrix_filt
  rownames(aux_assay) = aux_assay[,1]
  aux_assay <- aux_assay[,-1]

  if(out$replicates == TRUE){
    
    # Samples (averaged replicates)
    aux_assay_mean <- as.data.frame(t(aux_assay))
    aux_assay_mean$Sample <- rownames(aux_assay_mean)
    aux_assay_mean$Sample <- sapply(strsplit(aux_assay_mean$Sample,"_"), `[`, 1)
    aux_assay_mean$Sample <- factor(aux_assay_mean$Sample, levels = unique(aux_assay_mean$Sample))
    aux_assay_mean <- aux_assay_mean %>%
      dplyr::group_by(Sample) %>%
      dplyr::summarise(dplyr::across(where(is.double), mean, na.rm = T)) %>%
      tibble::column_to_rownames("Sample") %>% t() %>% as.data.frame()
  }

  message("Column annotations initialization...")
  
  if(shiny::isRunning()){
    showNotification(tagList(icon("cogs"), HTML("&nbsp;Column annotations initialization...")), type = "default")
  }
  
  
  # Replicates
  aux_coldata <- as.data.frame(out$targets$targetfile_lipidomics)
  rownames(aux_coldata) <-  aux_coldata$SampleID
  aux_coldata <- aux_coldata[!(rownames(aux_coldata) %in% out$sample_filtered), ]
  
  # removing na columns
  aux_coldata = aux_coldata[,colSums(is.na(aux_coldata))<nrow(aux_coldata)]

  if(out$replicates == TRUE){
    # Samples (averaged replicates)
    aux_coldata_mean <- aux_coldata
    aux_coldata_mean$SampleID <- sapply(strsplit(aux_coldata_mean$SampleID,"_"), `[`, 1)
    aux_coldata_mean$SampleID <- factor(aux_coldata_mean$SampleID, levels = unique(aux_coldata_mean$SampleID))
    aux_coldata_mean <- aux_coldata_mean %>%
      dplyr::distinct(SampleID, .keep_all = TRUE)
    rownames(aux_coldata_mean) <- aux_coldata_mean$SampleID
    #removing all na columns
    aux_coldata_mean = aux_coldata_mean[,colSums(is.na(aux_coldata_mean))<nrow(aux_coldata_mean)]
  }
  
  message("Row annotations initialization...")
  
  if(shiny::isRunning()){
    showNotification(tagList(icon("cogs"), HTML("&nbsp;Row annotations initialization...")), type = "default")
  }
  
  aux_rowdata <- strsplit(out$concentration_matrix_filt[,1], split = ")")
  aux_rowdata_ion <- unlist(lapply(aux_rowdata,"[",2))
  aux_rowdata_species <- unlist(lapply(aux_rowdata,"[",1)) %>% strsplit(split = "[(_/:]") ####MODIFICA aggiunto "/" in split
  aux_rowdata_species <- t(sapply(aux_rowdata_species, function(x) `length<-`(unlist(x), max(sapply(aux_rowdata_species, length)))))
  aux_rowdata_class <- aux_rowdata_species[,1]

  #---Subfunction---#
  fun_h <- function(h){
    hh_1 <- stringr::str_extract(h, "[[:alpha:]]+")
    hh_2 <- stringr::str_extract(h, "[[:digit:]]+")
    hh <- cbind(hh_1,hh_2)
    return(hh)
  }

  aux_list <- list()
  for(k in 2:ncol(aux_rowdata_species)){
    aux_list[[k]] <- t(sapply(aux_rowdata_species[,k],fun_h))
  }
  aux_rowdata <- as.data.frame(do.call(cbind, aux_list))
  col_names <- c("n_acyl","sn","ether","unsat")
  if(ncol(aux_rowdata) %% 4 == 0){
    aux_t <- ncol(aux_rowdata)/4
    colnames(aux_rowdata) <- paste(rep(col_names, times = aux_t),
                                   rep(1:aux_t, each = 4), sep = "_")
  } else {
    message("There is an issue on number of variable after parsing.")
    
    if(shiny::isRunning()){
      showNotification(tagList(icon("times-circle"), HTML("&nbsp;There is an issue on number of variable after parsing")), type = "error")
    }
  }

  aux_sn <- c(grep("sn",colnames(aux_rowdata), value = TRUE))
  aux_un <- c(grep("unsat",colnames(aux_rowdata), value = TRUE))
  aux_rowdata <- as.data.frame(aux_rowdata) %>%
    dplyr::mutate_at(c(aux_sn,aux_un), as.numeric) %>%
    dplyr::mutate(total_sn = rowSums(dplyr::across(all_of(aux_sn)),na.rm = TRUE),
                  total_un = rowSums(dplyr::across(all_of(aux_un)),na.rm = TRUE))
  aux_rowdata <- cbind(Lipids = out$concentration_matrix_filt[,1], Class = aux_rowdata_class, 
                       aux_rowdata, Ion = aux_rowdata_ion)
  rownames(aux_rowdata) <-  aux_rowdata$Lipids

  message("Metadata initialization...")
  
  if(shiny::isRunning()){
    showNotification(tagList(icon("cogs"), HTML("&nbsp;Metadata initialization...")), type = "default")
  }
  
  aux_metadata <- list(AnalysisLab = out$analysis$lab_analyst,
                       DataAnalyst = out$analysis$data_analyst,
                       AnalysisDate = out$analysis$analysis_date,
                       PipelineVersion = out$analysis$pipe_ver)

  # Replicates
  sumexp = SummarizedExperiment::SummarizedExperiment(assays = aux_assay,
                                                      rowData = aux_rowdata,
                                                      colData = aux_coldata,
                                                      metadata = aux_metadata)
  
  if(out$replicates == TRUE){
    # Samples (averaged replicates)
    sumexp_mean <- SummarizedExperiment::SummarizedExperiment(assays = aux_assay_mean,
                                                              rowData = aux_rowdata,
                                                              colData = aux_coldata_mean,
                                                              metadata = aux_metadata)
  } else {
    sumexp_mean = NA
  }
  
  # Updating
  aux_out <- list(sumexp_data = sumexp, sumexp_data_mean = sumexp_mean)

  message("---> SUMMERIZED EXPERIMENT DATA ADVISE-LIPIDOMICS PIPELINE END <---")

  res <- purrr::flatten(list(out,aux_out))
  
  if(shiny::isRunning()){
    showNotification(tagList(icon("check"), HTML("&nbsp;Summarized experiment data saved!")), type = "message")
  }
  
  return(res)

}
