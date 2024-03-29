# Theme: Poor Man's Flatly modified (ADViSE theme)------------------------------------------------------------------------------
#' @title theme_ADViSE
#' @description Poor Man's Flatly theme for a shinydashboard application edited
#' @importFrom dashboardthemes shinyDashboardThemeDIY
#'
#' @return Object produced by shinyDashboardThemeDIY
#' @seealso \code{\link{shinyDashboardThemeDIY}}
#' @export
#' 
theme_ADViSE <- dashboardthemes::shinyDashboardThemeDIY(
  
  ### general
  appFontFamily = "Arial"
  ,appFontColor = "rgb(33,37,41)"
  ,primaryFontColor = "rgb(245,245,245)"
  ,infoFontColor = "rgb(245,245,245)"
  ,successFontColor = "rgb(33,37,41)"
  ,warningFontColor = "rgb(33,37,41)"
  ,dangerFontColor = "rgb(33,37,41)"
  ,bodyBackColor = "rgb(255,255,255)"
  
  ### header
  ,logoBackColor = "rgb(42,170,226)"            #header
  
  ,headerButtonBackColor = "rgb(42,170,226)"    #header
  ,headerButtonIconColor = "rgb(44,47,118)"       #sidebar
  ,headerButtonBackColorHover = "rgb(20,154,128)"
  ,headerButtonIconColorHover = "rgb(0,0,0)"
  
  ,headerBackColor = "rgb(42,170,226)"             #header
  ,headerBoxShadowColor = ""
  ,headerBoxShadowSize = "0px 0px 0px"
  
  ### sidebar
  ,sidebarBackColor = "rgb(44,47,118)"            #sidebar
  ,sidebarPadding = 0
  
  ,sidebarMenuBackColor = "inherit"
  ,sidebarMenuPadding = 0
  ,sidebarMenuBorderRadius = 0
  
  ,sidebarShadowRadius = ""
  ,sidebarShadowColor = "0px 0px 0px"
  
  ,sidebarUserTextColor = "rgb(255,255,255)"
  
  ,sidebarSearchBackColor = "rgb(255,255,255)"
  ,sidebarSearchIconColor = "rgb(44,47,118)"         #sidebar
  ,sidebarSearchBorderColor = "rgb(255,255,255)"
  
  ,sidebarTabTextColor = "rgb(255,255,255)"
  ,sidebarTabTextSize = 14
  ,sidebarTabBorderStyle = "none"
  ,sidebarTabBorderColor = "none"
  ,sidebarTabBorderWidth = 0
  
  ,sidebarTabBackColorSelected = "rgb(30,43,55)"
  ,sidebarTabTextColorSelected = "rgb(42,170,226)"     #header
  ,sidebarTabRadiusSelected = "0px"
  
  ,sidebarTabBackColorHover = "rgb(44,47,118)"          #sidebar
  ,sidebarTabTextColorHover = "rgb(42,170,226)"        #header
  ,sidebarTabBorderStyleHover = "none"
  ,sidebarTabBorderColorHover = "none"
  ,sidebarTabBorderWidthHover = 0
  ,sidebarTabRadiusHover = "0px"
  
  ### boxes
  ,boxBackColor = "rgb(245,245,245)"
  ,boxBorderRadius = 0
  ,boxShadowSize = "0px 0px 0px"
  ,boxShadowColor = ""
  ,boxTitleSize = 19
  ,boxDefaultColor = "rgb(52,152,219)"
  ,boxPrimaryColor = "rgb(44,47,118)"        #sidebar
  ,boxInfoColor = "rgb(52,152,219)"
  ,boxSuccessColor = "rgb(24, 188, 156)"
  ,boxWarningColor = "rgb(243,156,18)"
  ,boxDangerColor = "rgb(231,76,60)"
  
  ,tabBoxTabColor = "rgb(44,47,118)"          #sidebar
  ,tabBoxTabTextSize = 14
  ,tabBoxTabTextColor = "rgb(24, 188, 156)"
  ,tabBoxTabTextColorSelected = "rgb(255, 255, 255)"
  ,tabBoxBackColor = "rgb(255,255,255)"
  ,tabBoxHighlightColor = "rgb(255,255,255)"
  ,tabBoxBorderRadius = 10
  
  ### inputs
  ,buttonBackColor = "rgb(44,47,118)"         #sidebar
  ,buttonTextColor = "rgb(255,255,255)"
  ,buttonBorderColor = "rgb(44,47,118)"        #sidebar
  ,buttonBorderRadius = 5
  
  ,buttonBackColorHover = "rgb(30,43,55)"
  ,buttonTextColorHover = "rgb(255,255,255)"
  ,buttonBorderColorHover = "rgb(30,43,55)"
  
  ,textboxBackColor = "rgb(255,255,255)"
  ,textboxBorderColor = "rgb(206,212,218)"
  ,textboxBorderRadius = 5
  ,textboxBackColorSelect = "rgb(255,255,255)"
  ,textboxBorderColorSelect = "rgb(89,126,162)"
  
  ### tables
  ,tableBackColor = "rgb(255,255,255)"
  ,tableBorderColor = "rgb(236,240,241)"
  ,tableBorderTopSize = 1
  ,tableBorderRowSize = 1
  
)