#' calibplot_advise_lipidomics
#'
#' @description \code{calibplot_advise_lipidomics} is the function for creating the calibration
#' curves by the calibration matrix, and plotting the curves for the all the internal standard lipid species.
#'
#' @param out List. It is the result from the \code{filter_advise_lipidomics} function.
#' @param cal_plot_path Character string. It is the name of folder in which the plots of the
#' calibration curves will be stored. Default = "CalibrationPlot".
#' @param cal_mat Dataframe. It is the calibration matrix, generated after the calibration step.
#' @param intercept_flag Logical value. If TRUE, linear model will be not forced to have the
#' intercept at zero. Default = FALSE.
#' @param lm_robust_flag Logical value. If TRUE, linear model will be performed in robust method
#' (by using the function by MASS package). Default = FALSE.
#' @param lin_calibration Logical value. If TRUE, a filter on concentration range of linearity is
#' applied, considering the min and max values in the internal standard file. Default  = TRUE.
#' @param plot_calibration Logical value. If TRUE calibration plots are drawn.
#'
#' @return res: a list with results from filtering step, updated with the calibration matrix
#' and the linear model coefficients.
#'
#' @import dplyr
#' @rawNamespace import(ggplot2, except = last_plot)
#' @import shiny
#' @importFrom MASS rlm
#' @importFrom purrr flatten
#' @importFrom reshape2 melt dcast
#' @importFrom stats lm coef
#'
#' @export
#'
#' @details The second part of the calibration step calculates the coefficient values
#' (intercept and slope) from the application of a linear regression on the values of
#' the calibration matrix. The linear regression can be normal or robust, furthermore
#' the intercept can be force to be zero. All the curves (lines) can be plotted and
#' stored in the suitable folder.
#'
#' @note Last change 02/11/2021

calibplot_advise_lipidomics <- function(out,
                                        cal_plot_path = "CalibrationPlot",
                                        cal_mat,
                                        intercept_flag = TRUE,
                                        lm_robust_flag = FALSE,
                                        lin_calibration = TRUE,
                                        plot_calibration = TRUE
                                        ){

  if (!dir.exists(paste0(out$folders$output_path_plots,cal_plot_path))){
    dir.create(file.path(paste0(out$folders$output_path_plots,cal_plot_path)))
  }

  message("Checking calibration linearity...")
  
  if(lin_calibration == TRUE){
    
    if(shiny::isRunning()){
      showNotification(tagList(icon("cogs"), HTML("&nbsp;Checking calibration linearity...")), type = "default")
    }
    
    aux_cal <- cbind(rownames(cal_mat),cal_mat) %>% dplyr::select(-LipidIon)
    aux_lin <- out$targets$internal_standard[,c("InternalStandardLipidIon",
                                                "MinLinearity","MaxLinearity")]
    colnames(aux_cal)[1] <- "InternalStandardLipidIon"
    aux_cal <- reshape2::melt(aux_cal, id.vars="InternalStandardLipidIon")
    aux_cal <- base::merge(aux_cal,aux_lin)
    aux_cal[as.numeric(aux_cal$variable) < aux_cal$MinLinearity |
              as.numeric(as.character(aux_cal$variable)) > aux_cal$MaxLinearity, "value"] <- NA
    aux_cal <- reshape2::dcast(aux_cal,InternalStandardLipidIon ~ variable)
    cal_mat <- aux_cal[,-1]
    rownames(cal_mat) <- aux_cal[,1]
  }
  
  message("---> CALIBRATION COEFFICIENTS AND PLOTS ADVISE-LIPIDOMICS PIPELINE START <---")

  if(intercept_flag == TRUE){

    message("Calibration without intercept at zero...")
    
    if(shiny::isRunning()){
      showNotification(tagList(icon("info"), HTML("&nbsp;Calibration without intercept at zero...")), type = "default")
    }
    
    coef_df <- data.frame(matrix(nrow = nrow(cal_mat),ncol = 2))
    colnames(coef_df) <- c("q","m")

    for(k in 1:nrow(cal_mat)){

      aux = t(rbind(as.numeric(cal_mat[k,]),as.numeric(colnames(cal_mat))))
      colnames(aux) = c("Area",rownames(cal_mat[k,]))

      # Coefficients
      if(lm_robust_flag == FALSE){
        model <- stats::lm(aux[,1]~aux[,2])
        coef_df[k,] <- stats::coef(model)
        met <- "lm"
        metgg = "lm"
      } else {
        model <- MASS::rlm(aux[,1]~aux[,2])
        coef_df[k,] <- stats::coef(model)
        met <-  "rlm"
        metgg = MASS::rlm
      }

      # Plots
      if(plot_calibration == TRUE){
        g <- ggplot(as.data.frame(aux),aes(x = aux[,2], y = aux[,1])) +
          stat_smooth(method = metgg, formula = y~x, se = FALSE, na.rm = TRUE) +
          geom_point(na.rm = TRUE) +
          labs(title = colnames(aux)[2]) +
          xlab("Concentration (ng/ml)") + ylab("Mean Area")
        colnames(aux)[2] <- gsub(":","-",colnames(aux)[2])
        colnames(aux)[2] <- gsub("/","o",colnames(aux)[2])
        ggsave(paste0(out$folders$output_path_plots, cal_plot_path, "\\Calibration_", met,
                      "_", colnames(aux)[2],".png"), plot = g, device = "png")
      }

    }

    coef_df <- cbind(InternalStandardLipidIon = rownames(cal_mat),coef_df)

  } else {

    if(intercept_flag == FALSE){

      message("Calibration with intercept at zero...")
      
      if(shiny::isRunning()){
        showNotification(tagList(icon("info"), HTML("&nbsp;Calibration with intercept at zero...")), type = "default")
      }
      
      
      coef_df <- data.frame(matrix(nrow = nrow(cal_mat),ncol = 1))
      colnames(coef_df) <- "m"

      for(k in 1:nrow(cal_mat)){

        aux = t(rbind(as.numeric(cal_mat[k,]),as.numeric(colnames(cal_mat))))
        colnames(aux) = c("Area",rownames(cal_mat[k,]))

        # Coefficients
        if(lm_robust_flag == FALSE){
          model <- stats::lm(aux[,1] ~ 0 + aux[,2])
          coef_df[k,] <- stats::coef(model)
          met = "lm"
          metgg = "lm"
        } else {
          model <- MASS::rlm(aux[,1] ~ 0 + aux[,2])
          coef_df[k,] <- stats::coef(model)
          met = "rlm"
          metgg = MASS::rlm
        }

        # Plots
        if(plot_calibration == TRUE){
          g <- ggplot(as.data.frame(aux),aes(x = aux[,2], y = aux[,1])) +
            stat_smooth(method = metgg, formula = y~0+x, se = FALSE, na.rm = TRUE) +
            geom_point(na.rm = TRUE) +
            labs(title = colnames(aux)[2]) +
            xlab("Concentration (ng/ml)") + ylab("Mean Area")
          colnames(aux)[2] <- gsub(":","-",colnames(aux)[2])
          colnames(aux)[2] <- gsub("/","o",colnames(aux)[2])
          ggsave(paste0(out$folders$output_path_plots, cal_plot_path,"\\Calibration_", met,
                        "_", colnames(aux)[2],".png"), plot = g, device = "png")
        }

      }

      coef_df <- cbind(InternalStandardLipidIon = rownames(cal_mat),coef_df)

    }

  }



  # Updating
  aux_out <- list(calibration_matrix = cal_mat, coefficients = coef_df)

  message("---> CALIBRATION COEFFICIENTS AND PLOTS ADVISE-LIPIDOMICS PIPELINE END <---")

  res <- purrr::flatten(list(out,aux_out))
  
  if(shiny::isRunning()){
    if(plot_calibration == TRUE){
      showNotification(tagList(icon("check"), HTML("&nbsp;Calibration coefficients calculated and plots done!")), type = "message")
      showNotification(tagList(icon("info"), HTML("&nbsp;Plots saved in ", cal_plot_path)), type = "default")
    }else{
      showNotification(tagList(icon("check"), HTML("&nbsp;Calibration coefficients calculated!")), type = "message")
    }
  }

  return(res)

}
