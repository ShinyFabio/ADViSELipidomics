#' filter_advise_lipidomics
#'
#' @description \code{filter_advise_lipidomics} is the function for filtering lipid data on
#' several features.
#'
#' @param out List. It is the result from the \code{read_advise_lipidomics} function.
#' @param ca_bound Numerical vector. It is the range of the carbon atoms for each lipid
#' species: first number is the lower bound, second number is the upper bound. Default = c(14,24).
#' @param db_bound Numerical vector. It is the range of the double bonds for each lipid
#' species: first number is the lower bound, second number is the upper bound. Default = c(0,6).
#' @param data_type Character. Where come from the data. Can be "LipidSearch" or "LIQUID".
#'
#' @return res: a list with results from the reading step, updated with the filtered data, divided
#' into non-labeled data (lipid_filtered) and deuterated data (lipid_deuterated).
#'
#' @import dplyr
#' @import shiny
#' @importFrom purrr flatten
#'
#' @export
#'
#' @details The filter step is composed of different filters applied on non-labeled files:
#' a) retention time, the retention times of the lipid species should be in the range reported
#' in the internal standard file for each class;
#' b) carbon atom, the carbon atoms of the lipid species should be in the range reported
#' in the internal standard file for each class, and should be an even number;
#' c) double bond, the double bonds of the lipid species should be in the range reported
#' in the internal standard file for each class;
#' d) duplicated area, the duplicated areas for a lipid species should be ingnored, conserving
#' only the minimum area for one lipid species (duplication is checkd on m/z values).
#' After this step, filtered data are stored separately from raw data.
#'
#' @note Last change 17/12/2021

filter_advise_lipidomics <- function(out,
                                     ca_bound = c(14,24),
                                     db_bound= c(0,6),
                                     data_type = "LipidSearch"){

  message("---> FILTERING DATA ADVISE-LIPIDOMICS PIPELINE START <---")

  # I) Check on retention time range
  message("Check on retention time range...")
  
  if(shiny::isRunning()){
    showNotification(tagList(icon("cogs"), HTML("&nbsp;Check on retention time range...")), type = "default")
  }

  if(data_type == "LipidSearch"){
    rt_range <- out$targets$internal_standard %>% dplyr::select(Class,Ion,MinRt,MaxRt,InternalStandardLipidIon)
    lipid_filtered <- lapply(out$lipid_data[grepl("nonlabeled",names(out$lipid_data), fixed = TRUE)],
                             function(x) merge(x, rt_range, all.x=TRUE) %>%
                               dplyr::filter(TopRT >= MinRt & TopRT <= MaxRt))
    
    lipid_deuterated <- lapply(out$lipid_data[grepl("deuterated",names(out$lipid_data), fixed = TRUE)],
                               function(x) merge(x, rt_range, by.x = "LipidIon",
                                                 by.y = "InternalStandardLipidIon") %>%
                                 dplyr::select(LipidIon,Area))
    
  }else{
    ##### LIQUID
    rt_range <- out$targets$internal_standard %>% dplyr::select(Class,Adduct,MinRt,MaxRt) %>% dplyr::rename(Adduct_std = Adduct)
    # Filtering lipids with 0 intensity
    lipid_filtered = lapply(out$lipid_data, function(x)
      x %>% dplyr::filter(Intensity != 0) %>% 
        dplyr::mutate(Class =  strsplit(Common_Name, split = "\\(") %>% lapply("[",1) %>% unlist())
    )

    # Filtering lipids based on Adduct and Apex_RT
    lipid_filtered = lapply(lipid_filtered, function(x) merge(x, rt_range, all.x=TRUE) %>%
                              dplyr::filter(Adduct == Adduct_std) %>%
                              dplyr::filter(Apex_RT >= MinRt & Apex_RT <= MaxRt))
  }


  # II) Check on number of carbon atoms: range and even number
  message("Check on number of carbon atoms...")
  
  if(shiny::isRunning()){
    showNotification(tagList(icon("cogs"), HTML("&nbsp;Check on number of carbon atoms...")), type = "default")
  }
  
  if(data_type == "LipidSearch"){
      #--- Subfunction---#
    fun_a <- function(a){
      aa <- a$FattyAcid
      aa <- gsub("[a-zA-Z()]", "", aa)
      aa <- strsplit(aa, split = "[:_]")
      aa <- list(ca = lapply(aa, function(x) c(x[seq(length(x)) %% 2 == 1])),
                 db = lapply(aa, function(x) c(x[seq(length(x)) %% 2 == 0])) )
      return(aa)
    }
  }else{
    #--- Subfunction---#
    fun_a <- function(a){
      aa <- a$Common_Name
      aa <- gsub("[a-zA-Z()-]", "", aa)
      aa <- strsplit(aa, split = "[:/]")
      aa <- list(ca = lapply(aa, function(x) c(x[seq(length(x)) %% 2 == 1])),
                 db = lapply(aa, function(x) c(x[seq(length(x)) %% 2 == 0])) )
      return(aa)
    }
  }


  ca_list <- lapply(lipid_filtered, function (x) fun_a(x)$ca)
  db_list <- lapply(lipid_filtered, function (x) fun_a(x)$db)

  ## II.1) Range

  #---Subfunction---#
  fun_b <- function(b){
    bb <- lapply(b, function(x) (as.numeric(x) >= ca_bound[1] & as.numeric(x) <= ca_bound[2]))
    bb <- which(sapply(bb, function(x) any(x == FALSE)))
    return(bb)
  }

  ca_list_idx <- lapply(ca_list, fun_b)
  for(k in 1:length(lipid_filtered)){
    if (length(ca_list_idx[[k]]) != 0){
      aux = lipid_filtered[[k]]
      lipid_filtered[[k]] <- aux[-c(ca_list_idx[[k]]),]
      rownames(lipid_filtered[[k]]) <- 1:nrow(lipid_filtered[[k]])
    }
  }

  ca_list_2 <- lapply(lipid_filtered, function (x) fun_a(x)$ca)
  db_list_2 <- lapply(lipid_filtered, function (x) fun_a(x)$db)

  ## II.2) Even Number

  #---Subfunction---#
  fun_c <- function(c){
    cc <- lapply(c, function(x) sum(as.numeric(x)))
    cc <- which(sapply(cc, function(x) (x %% 2 == 1)))
    return(cc)
  }

  ca_list_idx_2 <- lapply(ca_list_2, fun_c)
  for(k in 1:length(lipid_filtered)){
    if (length(ca_list_idx_2[[k]]) != 0){
      aux = lipid_filtered[[k]]
      lipid_filtered[[k]] <- aux[-c(ca_list_idx_2[[k]]),]
      rownames(lipid_filtered[[k]]) <- 1:nrow(lipid_filtered[[k]])
    }
  }

  ca_list_3 <- lapply(lipid_filtered, function (x) fun_a(x)$ca)
  db_list_3 <- lapply(lipid_filtered, function (x) fun_a(x)$db)

  # III) Check on insaturation (double bonds)
  message("Check on number of double bonds...")
  
  if(shiny::isRunning()){
    showNotification(tagList(icon("cogs"), HTML("&nbsp;Check on number of double bonds...")), type = "default")
  }
  
  #---Subfunction---#
  fun_d <- function(d){
    dd <- lapply(d, function(x) (as.numeric(x) >= db_bound[1] & as.numeric(x) <= db_bound[2]))
    dd <- which(sapply(dd, function(x) any(x == FALSE)))
    return(dd)
  }

  db_list_idx <- lapply(db_list_3,fun_d)
  for(k in 1:length(lipid_filtered)){
    if (length(db_list_idx[[k]]) != 0){
      aux = lipid_filtered[[k]]
      lipid_filtered[[k]] <- aux[-c(db_list_idx[[k]]),]
      rownames(lipid_filtered[[k]]) <- 1:nrow(lipid_filtered[[k]])
    }
  }

  
  if(data_type == "LipidSearch"){
    #IV) Check on same species per class (observed m/z values)
    message("Check on observed m/z values...")
    
    if(shiny::isRunning()){
      showNotification(tagList(icon("cogs"), HTML("&nbsp;Check on observed m/z values...")), type = "default")
    }
    
    #---Subfunction---#
    fun_e <- function(e){
      ee <- e %>% dplyr::select(Class,ObsMz,Area)
      ee$ObsMz <- as.integer(ee$ObsMz)
      ee$idx = rownames(ee)
      return(ee)
    }
    
    mz_value <- lapply(lipid_filtered,fun_e)
    mz_value_list <- lapply(mz_value, function(x) split(x,x$Class))
    
    #---Subfunction---#
    fun_f <- function(f){
      dup <- f$ObsMz[duplicated(f$ObsMz)]
      if(length(dup) != 0){
        all_dup = list()
        max_dup = list()
        for(k in 1:length(dup)){
          all_dup[[k]] <- f[f$ObsMz == dup[k],]
          max_dup[[k]] <- all_dup[[k]] %>% slice_max(Area,n = 1,with_ties = FALSE)
        }
        all_dup_df <- as.data.frame(do.call(rbind, all_dup))
        max_dup_df <- as.data.frame(do.call(rbind, max_dup))
        min_dup_idx <- setdiff(all_dup_df$idx, max_dup_df$idx)
        return(as.numeric(min_dup_idx))
      }
    }
    
    min_idx <- list()
    for(kk in 1:length(mz_value_list)){
      if(!is.null(unlist(lapply(mz_value_list[[kk]],fun_f)))){
        min_idx[[kk]] <- unlist(lapply(mz_value_list[[kk]],fun_f))
        aux <- lipid_filtered[[kk]]
        lipid_filtered[[kk]] <- aux[-c(min_idx[[kk]]),]
      } else {
        min_idx[[kk]] = 0
      }
    }
    
    # Updating
    aux_out <- list(lipid_filtered = lipid_filtered, lipid_deuterated = lipid_deuterated)
  }else{
    #### LIQUID 
    #IV) Check on replicates lipids
    message("Check on replicates lipids...")
    if(shiny::isRunning()){
      showNotification(tagList(icon("cogs"), HTML("&nbsp;Check on replicates lipids...")), type = "default")
    }
    
    #---Subfunction---#
    fun_f <- function(f){
      dup <- f$Common_Name[duplicated(f$Common_Name)] %>% unique()
      not_dup = setdiff(f$Common_Name, dup)
      if(length(dup) != 0){
        all_dup = list()
        max_dup = list()
        for(k in 1:length(dup)){
          all_dup[[k]] <- f[f$Common_Name == dup[k],]
          max_dup[[k]] <- all_dup[[k]] %>% dplyr::slice_max(Intensity,n = 1,with_ties = FALSE)
        }
        max_dup_df <- as.data.frame(do.call(rbind, max_dup))
        max_dup_df <- rbind(max_dup_df, dplyr::filter(f, Common_Name %in% not_dup))
      }
    }
    
    #Filter duplicated lipids. Only lipids with maximum intensity will be picked.
    lipid_filtered = lapply(lipid_filtered, fun_f)
    
    # Updating
    aux_out <- list(lipid_filtered = lipid_filtered)
  }
 

  message("---> FILTERING DATA ADVISE-LIPIDOMICS PIPELINE END <---")
  
  res <- purrr::flatten(list(out,aux_out))

  if(shiny::isRunning()){
    showNotification(tagList(icon("check"), HTML("&nbsp;Filtering data completed!")), type = "message")
  }
  
  return(res)

}
