#' Create contrast list
#' 
#' @description This function make a contrast list with maximum 2 variables.
#' 
#' 
#' @param out The SummarizedExperiment result.
#' @param design_vars Character string. A design variable.
#' @param batch_method Method for the batch variable. "given" or....
#' @param batch_type Character. "remove" or "fit"
#' @param batch_vars Character. Batch effect variable.
#' 
#' @import dplyr
#' @importFrom tidyr unite
#' 

contrast_advise_lipidomics <- function(out,
                                       design_vars = "",
                                       batch_vars = "",
                                       batch_type = "fit",
                                       batch_method = "given",
                                       ref = ""){
  
  message("---> CONTRAST GENERATION ADVISE-LIPIDOMICS PIPELINE START <---")
  
  #---Subfunction---#
  fun_i <- function(i,k,av,ref){
    if(length(av) == 1){
      ii <- dplyr::distinct(i, dplyr::across(names(av)), .keep_all = TRUE)
      #ordering by ref
      levels = c(ii[,1][!ii[,1] %in% ref], ref)
      ii <-  dplyr::slice(ii, match(levels, ii[,1]))
    } else if (length(av) == 2){
      ii <- i %>% 
        dplyr::filter(dplyr::across(names(av)[1]) == k) %>%  ###MODIFICA CON ACROSS  
        dplyr::distinct(dplyr::across(names(av)[2]), .keep_all = TRUE)  #MODIFICA CON ACROSS
      #ordering by ref
      levels = c(ii[,2][!ii[,2] %in% ref], ref)
      ii <-  dplyr::slice(ii, match(levels, ii[,2]))
    }
    ii <- ii %>% tidyr::unite(col = "unif_1", sep = "_") %>% t() 
    ii <- tryCatch({combn(ii, m = 2)}, 
             error = function(e){
               message("A variable has too few levels. Try with a different variable.")
               showNotification(tagList(icon("times-circle"), HTML("&nbsp;A variable has too few levels. Try with a different variable.")), type = "error")
             }
    )
      ii = ii %>% t() %>% as.data.frame() %>%
        tidyr::unite(col = "unif_2", sep = "vs") %>%
        dplyr::mutate(unif_3 = stringr::str_replace(unif_2,"vs","-")) %>%
        tidyr::unite(col = "contrast", sep = "=")
      return(ii)
  }
  
  if(batch_type == "fit" & batch_method == "given"){
    tot_vars = c(design_vars,batch_vars)
    tot_vars = tot_vars[tot_vars != ""]
  }else{
    tot_vars = design_vars
    tot_vars = tot_vars[tot_vars != ""]
  }

  if(length(tot_vars) > 2){
    stop("Too many variables (design and batch)!")
  } else {
    aux_tot <- as.data.frame(colData(out)) %>% dplyr::select(all_of(tot_vars)) #MODIFICA era out$sumexp_data
    av <- lapply(aux_tot,unique)
    con_list <- list()
    
    if(length(av) == 1){
      con_list <- fun_i(aux_tot,0,av,ref)
      #check
      if(length(con_list[[1]]) == 1){
        con_list <- NULL
      }
    } else if(length(av) == 2){
      for(k in av[[1]]){
        con_list[[k]] <- fun_i(aux_tot,k,av,ref)
        #check
        if(nrow(con_list[[1]]) == 1){
          con_list <- NULL
        }
      }
    }
    
  }


  message("---> CONTRAST GENERATION ADVISE-LIPIDOMICS PIPELINE END <---")
  
  return(con_list)

}