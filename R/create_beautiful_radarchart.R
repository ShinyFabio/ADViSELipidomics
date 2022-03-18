#' Create a customized spider plot
#' 
#' @description This function make a customized spider plot using the {fmsb} package.
#' 
#' 
#' @param data A dataframe with 3 or 4 rows and some column (columns will be the variables in the spider plot). 
#' The rows must contain a row with maximum values and a row with mininum values, plus one or two rows with a "Codice_azienda".
#' @param caxislabels Range of values inside the spider plot (eg. caxislabels = c(0, 2, 5, 7, 10)).
#' @param color Function that provides colors to the plot. For example color = grDevices::hcl.colors(2, palette = "Dynamic").
#' 
#' 
#' @importFrom fmsb radarchart
#' 
#' @examples \dontrun{
#' 
#' create_beautiful_radarchart(
#'  radardatamm[c("Max", "Min", "SA_10_massa_ultima","AV_02_massa_ultima"),], 
#'  caxislabels = c(0, 2, 4, 6, 8, 10), 
#'  color = grDevices::hcl.colors(2, palette = "Dynamic")
#'  )
#' 
#' }
#' 


create_beautiful_radarchart <- function(data, color = "#00AFBB", 
                                        vlabels = colnames(data),
                                        caxislabels = NULL, title = NULL, ...){
  radarchart(
    data, axistype = 1,
    # Customize the polygon
    pcol = color, pfcol = scales::alpha(color, 0.5), plwd = 2, plty = 1,
    # Customize the grid
    cglcol = "grey", cglty = 1, cglwd = 0.8,
    # Customize the axis
    axislabcol = "grey", 
    # Variable labels
    vlcex = 1.3, vlabels = vlabels, seg = 5,
    caxislabels = caxislabels, title = title, ...
  )
}