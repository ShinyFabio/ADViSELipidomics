#' filter_sumexp
#'
#' @description \code{filter_sumexp} is the function for the filtering of the
#' SummarizedExperiment object.
#'
#' @param data A list containing both the SumExp objects and other information like replicates and data_type.
#' @param filt_targ Character. A column from the colData where to apply the filter (e.g. "Product").
#' @param value_targ Character. The values to be filtered on the colData (e.g. c("CTRL", "BC")).
#' @param filt_class Character. The lipid classes to be filtered (e.g. c("Cer", "PC", "PE")).
#'
#' @importFrom SummarizedExperiment SummarizedExperiment colData assay rowData
#' @importFrom S4Vectors metadata
#' @importFrom tibble column_to_rownames rownames_to_column
#' @import shiny
#' @import dplyr
#'
#'
#'
#' @note Last change 23/06/2022



filter_sumexp = function(data, 
                        filt_targ = NULL, 
                        value_targ = NULL, 
                        filt_class = NULL){

  
  #### DATA
  coldata = data$sumexp_data %>% colData() %>% as.data.frame()
  assay = data$sumexp_data %>% assay()
  rowdata = data$sumexp_data %>% rowData() %>% as.data.frame()
  metadata = metadata(data$sumexp_data)
  
  #se non ci sono replicati il sumexp_all$sumexp_data_mean = NA
  if(data$replicates == TRUE){
    coldata_mean = data$sumexp_data_mean %>% colData() %>% as.data.frame()
    assay_mean = data$sumexp_data_mean %>% assay()
    rowdata_mean = data$sumexp_data_mean %>% rowData() %>% as.data.frame()
  }
  

  #### Filtering based on target file
  if(!is.null(value_targ) && value_targ != ""){
    message("Filtering samples...")
    if(shiny::isRunning()){
      showNotification(tagList(icon("cogs"), HTML("&nbsp;Filtering samples...")), type = "default")
    }
    sampleid_before = coldata$SampleID #just for the message
    
    coldata = coldata %>% dplyr::filter(!!sym(filt_targ) %in% value_targ)
    sampleid = coldata$SampleID
    assay = assay %>% dplyr::select(sampleid)
    
    if(data$replicates == TRUE){
      coldata_mean = coldata_mean %>% dplyr::filter(!!sym(filt_targ) %in% value_targ)
      sampleid_mean = coldata_mean$SampleID
      assay_mean = assay_mean %>% dplyr::select(sampleid_mean)
    }
    
    message("You correctly filtered ", length(sampleid), " out of ", length(sampleid_before), " samples (considering replicates).")
    if(shiny::isRunning()){
      showNotification(tagList(icon("check"), type = "message", 
                               HTML("&nbsp;You correctly filtered ", length(sampleid), " out of ", 
                                    length(sampleid_before), " samples (considering replicates).")))
    }
  }
  
  
  #### Filtering based on lipid class
  if(!is.null(filt_class) && filt_class != ""){
    
    message("Filtering lipids in ",paste0(filt_class, collapse = ", "), " classes...")
    if(shiny::isRunning()){
      showNotification(tagList(icon("cogs"), HTML("&nbsp;Filtering lipids in ",paste0(filt_class, collapse = ", "), " classes...")), type = "default")
    }
    lipids_before = rowdata$Lipids #just for the message
    
    
    rowdata = rowdata %>% dplyr::filter(Class %in% filt_class)
    lipids = rowdata$Lipids
    assay = assay %>% tibble::rownames_to_column("Lipids") %>% dplyr::filter(Lipids %in% lipids) %>% 
      tibble::column_to_rownames("Lipids")
    
    if(data$replicates == TRUE){
      rowdata_mean = rowdata_mean %>% dplyr::filter(Class %in% filt_class)
      lipids_mean = rowdata_mean$Lipids
      assay_mean = assay_mean %>% tibble::rownames_to_column("Lipids") %>% dplyr::filter(Lipids %in% lipids_mean) %>% 
        tibble::column_to_rownames("Lipids")
    }
    
    message("You correctly filtered ", length(lipids), " out of ", length(lipids_before), " lipids.")
    if(shiny::isRunning()){
      showNotification(tagList(icon("check"), type = "message", 
                               HTML("&nbsp;You correctly filtered ", length(lipids), " out of ", 
                                    length(lipids_before), " lipids.")))
    }
  }
  
   
  
  #### Recreate sumexp
  # Replicates
  sumexp = SummarizedExperiment::SummarizedExperiment(assays = assay,
                                                      rowData = rowdata,
                                                      colData = coldata,
                                                      metadata = metadata)
  
  if(data$replicates == TRUE){
    # Samples (averaged replicates)
    sumexp_mean <- SummarizedExperiment::SummarizedExperiment(assays = assay_mean,
                                                              rowData = rowdata_mean,
                                                              colData = coldata_mean,
                                                              metadata = metadata)
  } else {
    sumexp_mean = NA
  }
  
  out = list(sumexp_data = sumexp, sumexp_data_mean = sumexp_mean, replicates = data$replicates, data_type = data$data_type)
  
  return(out)
}
  