#' filter_advise_lipidomics
#'
#' @description \code{linebreaks} draw multiple break line (br())
#'
#' @param n number of breaks
#'
#' @import shiny
#'
#' @examples \dontrun{
#' linebreaks(10)   #add 10 breaks
#' }


linebreaks <- function(n){
  HTML(strrep(br(), n))
  }
