#' workbench_advise_lipidomics
#'
#' @description \code{workbench_advise_lipidomics} is the function for the creation of the
#' SummarizedExperiment object from a Metabolomics Workbench output.
#'
#' @param filtered List. It is the result from the \code{na_advise_lipidomics} function.
#' @param coldata the colData object from the Metabolomics Workbench output.
#' @param metadata metadata to be saved in the final summarized experiment.
#' @param data_type information about the data type (Area, Concentration, Peak Heigth...).
#' @param id the id of the selected metabolomics workbench experiment .
#'
#' @return res: a list with two SummerizedExperiment
#'
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom stringi stri_list2matrix
#' @importFrom stringr str_extract str_split_fixed
#' @import shiny
#' @import dplyr
#' @importFrom tibble column_to_rownames
#'
#' @export
#'
#' @details This function take the filtered matrix (assay), the target file (coldata) and 
#' the metadata to build a SummarizedExperiment. The rowData of the SumExp is also builded.
#'
#' @note Last change 28/12/2021


workbench_advise_lipidomics <-  function(filtered, coldata, metadata, data_type, id){
  
  #### assay ####
  
  message("Feature matrix initialization...")
  showNotification(tagList(icon("cogs"), HTML("&nbsp;Feature matrix initialization...")), type = "default")

  aux_assay <- filtered$concentration_matrix_filt[,-1]
  
  
  ####### for replicates
  if(id == "ST001115" || id == "ST000608"){

    # Samples (averaged replicates)
    aux_assay_mean <- as.data.frame(t(aux_assay))
    aux_assay_mean$Sample <- rownames(aux_assay_mean)
    aux_assay_mean$Sample <- sapply(strsplit(aux_assay_mean$Sample,"_"), `[`, 1)
    aux_assay_mean$Sample <- factor(aux_assay_mean$Sample, levels = unique(aux_assay_mean$Sample))
    aux_assay_mean <- aux_assay_mean %>%
      dplyr::group_by(Sample) %>%
      dplyr::summarise(dplyr::across(where(is.double), mean, na.rm = T)) %>%
      tibble::column_to_rownames("Sample") %>% t() %>% as.data.frame() #####aggiuta di as.data.frame()
  }
  ######
  
  

  #### rowData ####
  
  message("Row annotations initialization...")
  
  showNotification(tagList(icon("cogs"), HTML("&nbsp;Row annotations initialization...")), type = "default")

  if(id == "ST000608"){
    test = strsplit(filtered$concentration_matrix_filt[,1], split = "\\(")
    aux_rowdata_class = unlist(lapply(test,"[",1))
    
    aux_rowdata_species = unlist(lapply(test,"[",2))
    aux_rowdata_species = gsub("\\)", "", aux_rowdata_species)
    aux_rowdata_species = strsplit(aux_rowdata_species, split = "[/:]")
    aux_rowdata_species <- stringi::stri_list2matrix(aux_rowdata_species, byrow = TRUE)
    aux_rowdata_species = cbind(aux_rowdata_class, aux_rowdata_species)
    aux_rowdata_ion = unlist(lapply(test,"[",44)) #44 perchè sicuro non c'è nulla quindi tutti NA
  }else{
    aux_rowdata <- strsplit(filtered$concentration_matrix_filt[,1], split = ")")
    aux_rowdata_ion <- unlist(lapply(aux_rowdata,"[",2))
    aux_rowdata_species <- unlist(lapply(aux_rowdata,"[",1)) %>% strsplit(split = "[(/:]")
    aux_rowdata_species <- stringi::stri_list2matrix(aux_rowdata_species, byrow = TRUE)
    aux_rowdata_class <- aux_rowdata_species[,1]
  }
  

  
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
    showNotification(tagList(icon("times-circle"), HTML("&nbsp;There is an issue on number of variable after parsing")), type = "error")
  }
  
  aux_sn <- c(grep("sn",colnames(aux_rowdata), value = TRUE))
  aux_un <- c(grep("unsat",colnames(aux_rowdata), value = TRUE))
  aux_rowdata <- as.data.frame(aux_rowdata) %>%
    dplyr::mutate_at(c(aux_sn,aux_un), as.numeric) %>%
    dplyr::mutate(total_sn = rowSums(dplyr::across(all_of(aux_sn)),na.rm = TRUE),
                  total_un = rowSums(dplyr::across(all_of(aux_un)),na.rm = TRUE))
  aux_rowdata <- cbind(Lipids = filtered$concentration_matrix_filt[,1], Class = aux_rowdata_class, 
                       aux_rowdata, Ion = aux_rowdata_ion)
  rownames(aux_rowdata) <-  aux_rowdata$Lipids
  
  
  
  #### colData ####
  
  message("Column annotations initialization...")
  showNotification(tagList(icon("cogs"), HTML("&nbsp;Column annotations initialization...")), type = "default")


  aux_coldata = coldata %>% as.data.frame() %>% dplyr::rename(SampleID = local_sample_id)
  
  
  #######for replicates #####
  if(id == "ST001115" || id == "ST000608"){
    sample_corrected = stringr::str_split_fixed(aux_coldata$SampleID, "_",n = 3)
    sample_corrected = paste(paste(sample_corrected[,1], sample_corrected[,2], sep = "-"), sample_corrected[,3], sep = "_")
    aux_coldata$SampleID = sample_corrected
  }
  ####
  
  
  rownames(aux_coldata) <-  aux_coldata$SampleID
  aux_coldata <- aux_coldata[!(rownames(aux_coldata) %in% filtered$sample_filtered), ]
  aux_coldata = aux_coldata[,colSums(is.na(aux_coldata))<nrow(aux_coldata)]
  
  
  if(id == "ST001073"){
    aux_coldata = aux_coldata %>% dplyr::mutate(Treatment = dplyr::case_when(
      Saline_treatment == "Yes " ~ "Saline",
      Wnt3a_treatment == "Yes " ~ "Wnt3a",
      Saline_treatment == "No " & Wnt3a_treatment == "No " ~ "None"), .keep = "unused")
    aux_coldata$Day_post_crush = as.numeric(levels(aux_coldata$Day_post_crush))[aux_coldata$Day_post_crush]
  }

  
  if(id == "ST001115" || id == "ST000608"){
    # Samples (averaged replicates)
    aux_coldata_mean <- aux_coldata
    aux_coldata_mean$SampleID <- sapply(strsplit(aux_coldata_mean$SampleID,"_"), `[`, 1)
    aux_coldata_mean$SampleID <- factor(aux_coldata_mean$SampleID, levels = unique(aux_coldata_mean$SampleID))
    aux_coldata_mean <- aux_coldata_mean %>%
      dplyr::distinct(SampleID, .keep_all = TRUE)
    rownames(aux_coldata_mean) <- aux_coldata_mean$SampleID
    aux_coldata_mean = aux_coldata_mean[,colSums(is.na(aux_coldata_mean))<nrow(aux_coldata_mean)]
  }


  message("Metadata initialization...")
  showNotification(tagList(icon("cogs"), HTML("&nbsp;Metadata initialization...")), type = "default")

  

  if(id == "ST001115" || id == "ST000608"){
    # Samples (averaged replicates)
    tec_rep = TRUE
    sumexp_mean <- SummarizedExperiment::SummarizedExperiment(assays = aux_assay_mean,
                                                              rowData = aux_rowdata,
                                                              colData = aux_coldata_mean,
                                                              metadata = metadata)
  } else {
    sumexp_mean = NA
    tec_rep = FALSE
  }
  
  
  # Replicates
  sumexp = SummarizedExperiment::SummarizedExperiment(assays = aux_assay,
                                                      rowData = aux_rowdata,
                                                      colData = aux_coldata,
                                                      metadata = metadata)
  
  

  # Updating
  res <- list(sumexp_data = sumexp, sumexp_data_mean = sumexp_mean, replicates = tec_rep, data_type = data_type)
  
  message("---> SUMMERIZED EXPERIMENT DATA ADVISE-LIPIDOMICS PIPELINE END <---")
  
  showNotification(tagList(icon("check"), HTML("&nbsp;Summarized experiment data saved!")), type = "message")
  return(res)

}
