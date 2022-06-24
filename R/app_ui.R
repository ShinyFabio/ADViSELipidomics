#' The application User-Interface
#' 
#' @param request Internal parameter for `{shiny}`. 
#'     DO NOT REMOVE.
#' @import shiny
#' @import shinydashboard
#' @import shinyWidgets
#' @import dashboardthemes
#' @importFrom shinyFiles shinyDirChoose 
#' @importFrom DT DTOutput
#' @import shinyBS
#' @importFrom bsplus use_bs_tooltip bs_embed_tooltip
#' @importFrom DataEditR dataEditUI dataOutputUI
#' @importFrom shinycssloaders withSpinner
#' @noRd


app_ui <- function(request) {
  tagList(
    #shinythemes::themeSelector(),
    # Leave this function for adding external resources
    golem_add_external_resources(),
    
    bsplus::use_bs_tooltip(), 
    
    tags$head(tags$link(rel = "stylesheet", type = "text/css", href = "www/custom_notifications.css")),
    
    #custom infobox with height fixed (for the citation)
    tags$head(tags$style(HTML('.info-box {height: 90px;}'))),
    
    tags$head(tags$link(rel = "stylesheet", type = "text/css", href = "www/custom_dashboardheader_title.css")),

    # Your application UI logic 

    dashboardPage(
      dashboardHeader(title = span(
                                   tagList(tags$img(src = "www/NewLogoAL.png", width = "32px",style="margin-right: 4px;"), "ADViSELipidomics")),
                      tags$li(class = "dropdown", actionBttn("jumptohome", icon = icon("home"), style = "stretch", size = "lg"))),
      
      
    dashboardSidebar(
      sidebarUserPanel(name = textOutput("nome"),
                       subtitle = actionButton('change','Change', style='padding:0px; height:18px; font-size:90%'),
                       image = "www/userimage.png"
      ),
      
      sidebarMenu(id = "sidebarmenu",
        menuItem("Home", tabName = "home", icon = icon("home")),
        menuItem("Data Import & Preprocessing", tabName = "rawsub",icon = icon("file-import")),
        menuItem("SumExp Visualization", tabName = "seedatatab", icon = icon("table")),
        menuItem("Exploratory Analysis", tabName = "explormenu", icon = icon("chart-bar"), #startExpanded = T,
          menuSubItem("Plots", tabName = "graphtab"),
          menuSubItem("Clustering", tabName = "clusttab"),
          menuSubItem("Dimensionality Reduction", tabName = "dimredtab")
        ),
        menuItem("Statistical Analysis", tabName = "statanalmenu", icon = icon("balance-scale-right"), #startExpanded = T,
          menuSubItem("Differential Analysis", tabName = "diffexptab"),
          menuSubItem("Enrichment Analysis", tabName = "enrmenu")
        )
      )
    ),
    
    dashboardBody(
      
      #custom theme made with dashboardthemes (check R files)
      theme_ADViSE,
      
      #### Home ####
      tabItems(
        tabItem(
          tabName = "home",
          fluidRow(
          shinydashboard::box(width = 12, status = "primary",
            column(3, br(), tags$img(src = "www/advise_logo.png", width = "400px")),
            column(7,br(), br(), br(), 
              HTML("<h1 style = 'text-align: right;font-size: 53px;color: #0e3d51;'>
                     <strong>&nbsp;&nbsp;Welcome to ADViSELipidomics!</strong>
                     </h1>")),
            column(2, tags$img(src = "www/NewLogoAL.png", width = "140px"), style = "text-align: right"))),
          br(),fluidRow(column(12,wellPanel(
             h3(strong("ADViSELipidomics")," is a Shiny app with a complete workflow preprocessing and analyzing 
                lipidomics data from different sources. The user can upload the data files together with lipid 
                and sample details, select filters and statistical methods to apply to the dataset, and obtain the 
                results as tables and interactive plots.",style = "color: #0e3d51;")
          ))),
          linebreaks(3),
          fluidRow(column(width = 12, style = "text-align:center",
                          actionButton("gotoimport", "Start!", icon = icon("rocket"),class = "btn btn-primary btn-lg", 
                                       width = "290px", style='padding:10px; font-size:310%; font-weight: bold;'))),
          linebreaks(5),
          fluidRow(
            shinydashboard::infoBox(width = 3, icon = icon("github"), title = "Bug reports", 
            subtitle = "If you encounter a bug or a problem click here and go to the GitHub page.", 
            color = "aqua", href = "https://github.com/ShinyFabio/ADViSELipidomics/issues", fill = T),
          shinydashboard::infoBox(
            width = 6, icon = icon("newspaper"), title = "Citation", color = "aqua", fill = T,
            subtitle = "Please cite our article (E. Del Prete, A. M. Campos, 
            F. Della Rocca, C. Gallo, A. Fontana, G. Nuzzo, C. Angelini, ADViSELipidomics: a workflow for analyzing lipidomics data, 
            Submitted to Bioinformatics, 2022) when you publish using this tool."),
          shinydashboard::infoBox(width = 3, icon = icon("book"), title = "Guide", 
                                  subtitle = "If you need help in the use of ADViSELipidomics, click here.", 
                                  color = "aqua", href = "https://shinyfabio.github.io/ADViSELipidomics_book/", fill = T)
          ),
          linebreaks(2),
          column(8, offset = 2, wellPanel(
            h4("This work was supported by the project 'Antitumor Drugs and Vaccines from the Sea 
               (ADViSE)' project (CUP B43D18000240007-SURF 17061BP000000011) funded by POR Campania
               FESR 2014-2020 'Technology Platform for Therapeutic Strategies against Cancer'-Action 1.2.1 and 1.2.2.")))

        ),
        
        tabItem(tabName = "rawsub",
          fluidPage(
            fluidRow(
              
                box(width = 3, title = "What is your data source?", status = "primary", solidHeader = TRUE,
                    column(
                      10,offset = 1,
                      
                  radioGroupButtons("sel_inputdata_type", "Choose one :",
                                    choiceValues = list("lipidsearch", "liquid", "excel", "sumexp", "mw"),
                                    choiceNames = list(
                                      paste(tags$img(src = "www/lipidsearch_icon.png", height = "27px", width = "27px"), HTML("<b style=font-size:16px>&nbsp;Lipid Search</b>")),
                                      paste(tags$img(src = "www/liquid_icon.png", height = "29px", width = "29px"), HTML("<b style=font-size:16px>&nbsp;LIQUID</b>")),
                                      paste(tags$img(src = "www/excel_icon.png", height = "23px", width = "23px"), HTML("<b style=font-size:16px>&nbsp;Excel files</b>")),
                                      paste(icon("database",  style='font-size:16px;'), HTML("<b style=font-size:16px>&nbsp;Summarized Experiment</b>")),
                                      paste(tags$img(src = "www/mw_icon.png", height = "27px", width = "27px"), HTML("<b style=font-size:16px>&nbsp;Metabolomics Workbench</b>"))
                                    ),
                                    direction = "vertical",justified = TRUE, status = "primary",size = "lg"
                  ))),
                  
                  conditionalPanel(condition = "input.sel_inputdata_type == 'lipidsearch'",
                    column(2,
                      box(width = NULL, title = "Internal standard", status = "primary", solidHeader = TRUE,
                          helpText("Here you can decide if use internal standard or not. If not, the calibration step will be
                               skipped and your final output will be an area matrix, otherwise lipids area are 
                               calibrated and the final output will be a concentration matrix."),
                          awesomeRadio("type_lipsearch",
                                      label = "Do you have internal standards?", 
                                      choices = c("Yes","No"), selected = "Yes")
                      )
                    ))
                
            ),
            
            
            hr(),
            
            

              ##### Sumexp Output #####
              conditionalPanel(
                condition = "input.sel_inputdata_type == 'sumexp'",
                box(width = 5, title = "Summarized Experiment file", status = "primary", solidHeader = TRUE,
                    fileInput("finalsumexpinput", "Import a summarizeExperiment file (.rds)", accept = ".rds")
                )
                
              ),
              
              
              ##### LIQUID Output ####
              conditionalPanel(
                condition = "input.sel_inputdata_type == 'liquid'",
                
                ###### Box Filtering #####
                column(
                  5,
                  box(width = NULL, title = "1. Importing & Filtering", status = "primary", solidHeader = TRUE,
                      
                      # Step 1.
                      div(icon("circle"), HTML("&nbsp;Step. 1: Import targetfile and internal reference files (.xlsx)"), style = "text-align: left; font-size: 18px; font-weight: bold;"),
                      br(),
                      fluidRow(
                        column(10, fileInput("targetfile_liquid", "Select the Targetfile Lipidomics (.xlsx)", accept = ".xlsx")),
                        column(2, br(), 
                               tags$head(tags$style("#edit_target_liquid-upinternalmodal .modal-dialog{ width:1300px}")),
                               tags$head(tags$style("#edit_target_liquid-upinternalmodal .modal-body{ min-height:1000px}")),
                               mod_edit_data_ui("edit_target_liquid"))
                        
                      ),
                      fluidRow(
                        column(10, fileInput("internalstdpath_liquid", "Select the Internal Reference file (.xlsx)", accept = ".xlsx")),
                        column(2, 
                               br(),
                               tags$head(tags$style("#edit_internal_liquid-upinternalmodal .modal-dialog{ width:1300px}")),
                               tags$head(tags$style("#edit_internal_liquid-upinternalmodal .modal-body{ min-height:1000px}")),
                               mod_edit_data_ui("edit_internal_liquid")) 
                      ),
                      hr(),
                      
                      #Step 2.
                      conditionalPanel(
                        condition = "output.checktargets_liquid == false",
                        div(icon("circle"), HTML("&nbsp;Step. 2: Choose the data folder and read all the data."), style = "text-align: left; font-size: 18px; font-weight: bold;"),
                        br(),
                        fluidRow(
                          div(
                            column(4, shinyFiles::shinyDirButton("datafolder_liquid", label = "Browse...", title = "Please select the data folder", multiple = FALSE, icon = icon("folder-open"))),
                            column(4, actionButton("readdatabttn_liquid", "Read Data", icon("cogs"))),
                            column(4, 
                                   conditionalPanel(
                                     condition = "output.checkdatafiles_liquid == false",
                                     actionButton("showqualcheck_liquid", "Quality check", icon("eye")))
                            ), style = "text-align: center;")
                        ),
                        tags$head(tags$style("#viewqualcheckplot_liquid .modal-dialog{ width:1400px}")),
                        tags$head(tags$style("#viewqualcheckplot_liquid .modal-body{ min-height:800px}")),
                        shinyBS::bsModal("viewqualcheckplot_liquid", trigger = "showqualcheck_liquid", size = "large", title = "Quality check on intensity",
                                         plotlyOutput("qualcheckplot_liquid",height = "600px")),
                        hr()
                      ),
                      
                      
                      #Step 3.
                      conditionalPanel(
                        condition = "output.checkdatafiles_liquid == false",
                        div(icon("circle"), HTML("&nbsp;Step. 3: Select a range for the carbon number and the double bound number"), style = "text-align: left; font-size: 18px; font-weight: bold;"),
                        br(),
                        sliderInput("ca_bound_liquid", div("Range carbon number", style = "font-size: 15px; font-weight: bold;"), min = 0, max = 50, value = c(14, 24)),
                        sliderInput("db_bound_liquid", div("Range double bound number", style = "font-size: 15px; font-weight: bold;"), min = 0, max = 10, value = c(0, 6)),
                        fluidRow(
                          column(6,
                                 div(actionButton("gofilterbttn_liquid", "Filter Data", class = "btn-primary btn-lg", icon("cogs"), width = "160px", style='height:40px; padding:5px; font-size:140%; font-weight: bold;'), 
                                     style = "text-align:center;")),
                          column(6,
                                 conditionalPanel(
                                   condition = "output.checkendfiltering_liquid == false",
                                   div(actionButton("viewfiltered_liquid", "Check filtered data", class = "btn-primary btn-lg", icon("eye"), style='width:230px; height:40px; padding:5px; font-size:140%; font-weight: bold;'), 
                                       style = "text-align:center;"),
                                   tags$head(tags$style("#viewfilteredmodal_liquid .modal-dialog{ width:1100px}")),
                                   shinyBS::bsModal("viewfilteredmodal_liquid", trigger = "viewfiltered_liquid", size = "large",
                                                    selectInput("selcol_lipidfilt_liquid", "Select a sample", choices = ""),
                                                    div(DT::DTOutput("dtlipidfilterd_liquid"), style = "overflow-x: scroll;"))
                                 )
                          )
                        ) #end of fluidrow
                        
                      )
                      
                  ) #end of box filtering
                ), #end of column for the first part "filtering"
                
                
                ###### Imputation step 
                conditionalPanel(
                  condition = "output.checkendfiltering_liquid == false",
                  mod_imputation_step_ui("imputation_step_liquid", offset = 1)
                )
                
              ),
              
              ##### MW output ####
              conditionalPanel(
                condition = "input.sel_inputdata_type == 'mw'",
                box(width = 5, title = "Metabolomics Workbench", status = "primary", solidHeader = TRUE,
                fluidRow(
                  column(8,selectInput("selmwid", "Select a MW study ID", choices = c("ST001073","ST001115","ST000608"))),
                  column(4,br(),actionButton("loadmw", "Load!",icon("cloud-download-alt")))
                ),
                conditionalPanel(
                  condition = "output.checkmw == false",
                  selectInput("typedata_assay_MW", choices = c("Concentration", "Area", "Peak intensity"), selected = "Area",
                              label = tags$span("Select the assay measure", 
                                                tags$i(class = "glyphicon glyphicon-info-sign", style = "color:#0072B2;",
                                                       title = "This information will be used for the plot axis labels."))
                  ),
                  hr(),
                  ###Filtering and imputation
                  conditionalPanel(
                    condition = "output.check_naassay_MW == true",
                    fluidRow(
                      column(5, h5(strong("Filter and impute NA? ")), style="padding-right: 0px; width: 175px;"),
                      column(4, prettyToggle("togglefilterna_ass_MW", label_on = "Yes!", label_off = "No..", outline = TRUE, plain = TRUE, bigger = TRUE, icon_on = icon("thumbs-up"), icon_off = icon("thumbs-down")),
                             style="padding-top: 8px;padding-left: 0px;")
                    ),
                    conditionalPanel(
                      condition = "input.togglefilterna_ass_MW == true",
                      fluidRow(
                        column(
                          8,
                          sliderInput("na_filt_lip_ass_MW", "Max missing data percentage allowable on lipids", min = 0, max = 1, step = 0.1, value = 0.3),
                          sliderInput("na_filt_sam_ass_MW", "Max missing data percentage allowable on samples", min = 0, max = 1, step = 0.1, value = 0.6),
                          selectInput("imput_method_ass_MW", "Imputation method", choices = c("mean", "median", "none", "knn", "irmi"), selected = "knn")
                        ),
                        column(
                          4,
                          br(), br(), br(), br(),
                          div(actionButton("gofilterna_MW", "Check filtered NA", class = "btn-primary btn-lg", icon("eye"), style="white-space: normal; height:60px; width:140px;"),
                              style = "text-align:center;")
                        )
                      )
                    ),
                    
                    tags$head(tags$style("#modalfilteredna_MW .modal-dialog{ width:1300px}")),
                    tags$head(tags$style("#modalfilteredna_MW .modal-body{ min-height:1000px}")),
                    shinyBS::bsModal("modalfilteredna_MW", trigger = "gofilterna_MW",
                                     radioGroupButtons("vim_or_dt_filtermodal_MW", "Choose what to show :",
                                                       choiceValues = list("table", "vim"),
                                                       choiceNames = list(
                                                         paste(shiny::icon("table",  style='font-size:16px;'), HTML("<b style=font-size:16px>&nbsp;Table</b>")),
                                                         paste(shiny::icon("buromobelexperte",  style='font-size:16px;'), HTML("<b style=font-size:16px>&nbsp;NA plots</b>"))),
                                                       justified = TRUE, status = "primary"),
                                     conditionalPanel(condition = "input.vim_or_dt_filtermodal_MW == 'table'",
                                                      div(DT::DTOutput("dt_filteredna_MW"), style = "overflow-x: scroll;"),
                                                      br(),br(),
                                                      fluidPage(
                                                        fluidRow(
                                                          column(3, div(shinydashboard::valueBoxOutput("nadim1_MW",width = NULL)), style="padding-right: 0px;"),
                                                          column(1, tags$img(src = "www/right_arrow.png", width = "90px")),
                                                          column(3, div(shinydashboard::valueBoxOutput("nadim2_MW",width = NULL)), style="padding-right: 0px;padding-left: 0px;")
                                                        )
                                                      )
                                                      
                                     ),
                                     conditionalPanel(condition = "input.vim_or_dt_filtermodal_MW == 'vim'",
                                                      switchInput("vimopz_filt_MW", label = "Combined", value = TRUE),
                                                      plotOutput("vimplot_filteredna_MW", height = "1000px")
                                     )
                    ) #end of bsModal
                    
                  ),
                  hr(),
                  fluidRow(
                    column(6,
                      div(actionButton("make_sumexp_MW", "Create Summarized Experiment", icon = icon("cogs")), style = "text-align: center")
                    ),
                    column(
                      6,
                      conditionalPanel(
                        condition = "output.checkstepmw == false",
                        div(actionButton("gotosumexp_from_MW", "See the results", icon("eye")), style = "text-align: center")
                        )
                      )
                  )
                  
                )
              )
              ),
              
              ##### File excel output #####
              conditionalPanel(condition = "input.sel_inputdata_type == 'excel'",
                box(width = 5, title = "Excel files", status = "primary", solidHeader = TRUE,
                ####ColData
                fileInput("coldatainput", "Target file (.xlsx)",accept = ".xlsx"),
                
                ### ASSAY
                conditionalPanel(condition = "output.checkcold == false",
                  hr(),
                  fileInput("assayinput", "Data matrix (.xlsx)", accept = ".xlsx"),
                  conditionalPanel(condition = "output.checkassay == false",
                      selectInput("typedata_assay", choices = c("Concentration", "Area", "Peak intensity"),
                                  label = tags$span("Select the assay measure", 
                                                    tags$i(class = "glyphicon glyphicon-info-sign", style = "color:#0072B2;",
                                                           title = "This information will be used for the plot axis labels."))
                                  ),
                    ###Filtering and imputation
                    conditionalPanel(condition = "output.check_naassay == true",
                      hr(),
                      fluidRow(
                        column(5, h5(strong("Filter and impute NA? ")), style="padding-right: 0px; width: 175px;"),
                        column(4, prettyToggle("togglefilterna_ass", label_on = "Yes!", label_off = "No..", outline = TRUE, plain = TRUE, bigger = TRUE, icon_on = icon("thumbs-up"), icon_off = icon("thumbs-down")),
                               style="padding-top: 8px;padding-left: 0px;")
                      ),
                      conditionalPanel(
                        condition = "input.togglefilterna_ass == true",
                        fluidRow(
                          column(
                            8,
                            sliderInput("na_filt_lip_ass", "Max missing data percentage allowable on lipids", min = 0, max = 1, step = 0.1, value = 0.3),
                            sliderInput("na_filt_sam_ass", "Max missing data percentage allowable on samples", min = 0, max = 1, step = 0.1, value = 0.6),
                            selectInput("imput_method_ass", "Imputation method", choices = c("mean", "median", "none", "knn", "irmi"), selected = "knn")
                          ),
                          column(
                            4,
                            br(), br(), br(), br(),
                            div(actionButton("gofilterna_csv", "Check filtered NA", class = "btn-primary btn-lg", icon("eye"), style="white-space: normal; height:60px; width:140px;"),
                                style = "text-align:center;")
                          )
                        )
                      ),
                      
                      tags$head(tags$style("#modalfilteredna_csv .modal-dialog{ width:1300px}")),
                      tags$head(tags$style("#modalfilteredna_csv .modal-body{ min-height:1000px}")),
                      shinyBS::bsModal("modalfilteredna_csv", trigger = "gofilterna_csv",
                                       radioGroupButtons("vim_or_dt_filtermodal_csv", "Choose what to show :",
                                                         choiceValues = list("table", "vim"),
                                                         choiceNames = list(
                                                           paste(shiny::icon("table",  style='font-size:16px;'), HTML("<b style=font-size:16px>&nbsp;Table</b>")),
                                                           paste(shiny::icon("buromobelexperte",  style='font-size:16px;'), HTML("<b style=font-size:16px>&nbsp;NA plots</b>"))),
                                                         justified = TRUE, status = "primary"),
                                       conditionalPanel(condition = "input.vim_or_dt_filtermodal_csv == 'table'",
                                                        div(DT::DTOutput("dt_filteredna_csv"), style = "overflow-x: scroll;"),
                                                        br(),br(),
                                                        fluidPage(
                                                          fluidRow(
                                                            column(3, div(shinydashboard::valueBoxOutput("nadim1_csv",width = NULL)), style="padding-right: 0px;"),
                                                            column(1, tags$img(src = "www/right_arrow.png", width = "90px")),
                                                            column(3, div(shinydashboard::valueBoxOutput("nadim2_csv",width = NULL)), style="padding-right: 0px;padding-left: 0px;")
                                                          )
                                                        )
                                                        
                                       ),
                                       conditionalPanel(condition = "input.vim_or_dt_filtermodal_csv == 'vim'",
                                                        switchInput("vimopz_filt_csv", label = "Combined", value = TRUE),
                                                        plotOutput("vimplot_filteredna_csv", height = "1000px")
                                       )
                      ) #end of bsModal
                    ),
                    hr(),
                    fluidRow(
                      column(6, br(),
                             div(actionButton("make_sumexpcsv", "Create Summarized Experiment", icon = icon("cogs")), style = "text-align: center")
                      ),
                      column(
                        6,
                        conditionalPanel(
                          condition = "output.checkstepcsv == false",
                          div(
                            actionButton("gotosumexp_from_csv", "See the results", icon("eye")),
                            br(),br(),
                            downloadButton("downloadsumexp_csv", "Download"),
                            style = "text-align:center;")
                          
                        )
                      )
                    )
                    
                  )
                )
              ) #end of box
            ),
            
              
              ###### Lipidsearch output ####
              conditionalPanel(
                condition = "input.sel_inputdata_type == 'lipidsearch'",
                
                fluidRow(
                  mod_reading_step_ui("reading_step_lipsearch"),
                  
                  ##### Box Calibration #####
                  conditionalPanel(
                    condition = "output.checkendfiltering == false && input.type_lipsearch == 'Yes'", 
                    
                    mod_calibration_step_ui("calibration_step_lipsearch_with")
                    
                  ),
                  
                  conditionalPanel(
                    condition = "output.checkendfiltering == false && input.type_lipsearch == 'No'",
                    mod_imputation_step_ui("imputation_step_lipsearch_without", offset = 1)
                  )
                ),
                
                
                #### Box NA Filtering ####
                fluidRow(
                  conditionalPanel(
                    condition = "output.checkstepg == false && input.type_lipsearch == 'Yes'",
                    
                    mod_imputation_step_ui("imputation_step_lipsearch_with", offset = 0)        
                    
                  ) #end of conditional panel stepg
                )
                

              )
              


          ) #end of fluidpage
      ), #end of tabItem
      

      ##### Summ EXp tab ######
      tabItem(tabName = "seedatatab",
        
        box(width = 12, status = "primary", #style = "overflow-y: scroll; overflow-x: scroll;",
            fluidRow(column(
              width = 1,
              column(
                width = 7, style="padding-left: 0px;",
                dropdownButton(
                  circle = TRUE, status = "danger", icon = icon("cog"), width = "300px",
                  #tooltip only works if the text is in one line
                  tooltip = tooltipOptions(title = "Here you can decide which part of the Summarized Experiment to show and also generate a filtered version of it (both filtering on lipids and samples)."),
                  selectInput("sumexpselectobj", "Select what to show", choices = c("colData", "rowData", "assays",  "metadata"), multiple = FALSE),
                  conditionalPanel(condition = "output.check_replicates == true",
                                   h5(strong("Summarized data")),
                                   materialSwitch("summ_viewtable", value = FALSE, status = "primary")
                  ),
                  div(actionButton("filt_yourdata", "Filter data", icon("edit"), style ="font-size:120%"), style = "text-align:center")
                )
                )
              )),
            
            shinyBS::bsModal(
              "filt_yourdata_modal", "Filter data", trigger = "filt_yourdata", size = "large",
              fluidRow(
                column(
                  6,
                  box(
                    width = NULL, status = "primary", title = "Filter samples", solidHeader = TRUE,
                    fluidRow(
                      column(6, selectInput("filt_by_target", "Select a variable", choices = "")),
                      column(6, selectInput("filt_by_target_values", "Select the value(s)", choices = "", multiple = TRUE))
                    )
                  )
                ),
                column(
                  6,
                  box(width = NULL, status = "primary", title = "Filter lipids", solidHeader = TRUE,
                    selectInput("filt_by_lipids", "Filter by lipid class", choices = "", multiple = TRUE)
                  )
                )
              ),
              
              fluidRow(
                column(4,offset = 2, actionButton("gofilt_sumexp", "Apply filter", icon("cogs"), style='font-size:140%; font-weight: bold;')),
                column(4, offset = 2, conditionalPanel(
                  condition = "output.check_filt == false", downloadButton("down_filtsumexp", "Download", style='font-size:140%; font-weight: bold;')
                ))
              )
                             
            ),
            DT::DTOutput("dtsumexp")
            
            ),
        
      ), #end of tabitem  seedatatab


      ##### Plots #####

      tabItem(tabName = "graphtab",
        tabsetPanel(
          tabPanel("Lipids",
          fluidRow(
            column(
              2,
              box(width = NULL, status = "primary",
                  radioButtons("selplotlip1", choices = c("Lipid class distribution", "Lipid class proportion", "Lipid species distribution", "Lipid boxplot"),
                    label = tags$span("Data information",
                      tags$i(style = "color:#0072B2;",class = "glyphicon glyphicon-info-sign",
                             title = "Lipid class distribution is the number of lipids that belong to each class. 
                             Lipid class proportion is the relative proportion among lipids classes for each sample.
                             Lipid species distribution is the abudance of each lipid species for each condition in a given
                             lipidic class. Lipid boxplot shows the amount of a selected lipid compared to a variable."))
                      ),

                conditionalPanel(
                  condition = "input.selplotlip1 == 'Lipid class distribution'",
                  radioButtons("selplotlip", "Plot type", choices = list("Pie chart" = 1, "Barplot" = 2, "Spiderplot" = 3), selected = 2)
                ),
                conditionalPanel(
                  condition = "input.selplotlip1 == 'Lipid class proportion'",
                  awesomeCheckbox("summ_taxabar", "Summarize data"),
                  selectInput("annot_taxa", "Annotation", choices = "")
                ),
                
                conditionalPanel(
                  condition = "input.selplotlip1 == 'Lipid species distribution'",
                  selectInput("class_spec_dist", "Select a lipid class", choices = ""),
                  selectInput("fill_spec_dist", "Select a fill variable", choices = ""),
                  awesomeCheckbox("summ_spec_dist", "Summarise by fill variable", value = TRUE),
                  awesomeCheckbox("facet_spec_dist", "Faceting with a variable", value = FALSE),
                  conditionalPanel(
                    condition = "input.facet_spec_dist == true",
                    selectInput("facet_spec_var", "Faceting variable", choices = "")),
                  awesomeCheckbox("flip_spec", "Reverse axis", value = FALSE)
                  
                ),
                
                #boxplot
                conditionalPanel(
                  condition = "input.selplotlip1 == 'Lipid boxplot'",
                  selectInput("select_lipidplot", "Select a lipid", choices = ""),
                  selectInput("select_fillboxplot", "Select fill column", choices = ""),
                  fluidRow(
                    column(6,awesomeCheckbox("log_lipidboxplot", "Log scale", value = TRUE)),
                    column(6, conditionalPanel(condition = "output.check_replicates == true",
                                               awesomeCheckbox("summ_lipidboxplot", "Summarize data", value = FALSE)))
                  ),
                  awesomeCheckbox("addpoints_lipidboxplot", "Add points", value = FALSE)
                )
              )
            ),
            
            conditionalPanel(condition = "input.selplotlip1 == 'Lipid class distribution'",
              column(width=9, offset = 1, shinycssloaders::withSpinner(uiOutput("lipidsplot")))
            ),
            conditionalPanel(condition = "input.selplotlip1 == 'Lipid class proportion'",
              column(10, shinycssloaders::withSpinner(uiOutput("taxabarplot_ui"))),
            ),
            conditionalPanel(
              condition = "input.selplotlip1 == 'Lipid species distribution'",
              column(10, shinycssloaders::withSpinner(plotly::plotlyOutput("lipspec_barplot",height = "650px")))
            ),
            conditionalPanel(
              condition = "input.selplotlip1 == 'Lipid boxplot'",
              column(10, shinycssloaders::withSpinner(plotly::plotlyOutput("lipcondplot")))
              )
          )
        ),


        tabPanel("Scatterplot",
          sidebarLayout(
            sidebarPanel(
              width = 3,
              selectInput("sel_sample1_scatt", "Select a sample", choices = ""),
              selectInput("sel_sample2_scatt", "Select another sample", choices = ""),
              fluidRow(
                column(6, awesomeCheckbox("log_scatt", "Log scale", value = TRUE)),
                column(6, 
                  conditionalPanel(condition = "output.check_replicates == true",
                    awesomeCheckbox("summ_scatt", "Summarize data", value = FALSE))
                  )
              )
            ),
            mainPanel(
              width = 9,
              shinycssloaders::withSpinner(plotly::plotlyOutput("scattsampleplot"))
            )
          )
        ), #end of tabpanel scatterplot
        
        tabPanel("Heatmap",
          sidebarLayout(
            sidebarPanel(width = 3,
              div(actionButton("makeheatmap", label = "Make Heatmap", icon = icon("cogs"),class = "btn btn-primary btn-lg", width = "190px", style='padding:5px; font-size:140%; font-weight: bold;'), align= "center"),
              br(),
              h4(strong("Preprocessing")),
              awesomeCheckbox("filter_heatmap", "Filter data", value = FALSE),
              
              conditionalPanel(
                condition = "input.filter_heatmap == true",
                fluidRow(
                  column(6, selectInput("filtheatmapcol", "Column filtering", choices = "")),
                  column(6, selectInput("filtheatmapval", "Value filtering", choices = "", multiple = TRUE))
                ),
                
                fluidRow(
                  conditionalPanel(
                    condition = "output.checkadd2filt_heat == 'twovar'",
                    column(6, selectInput("filtheatmapcol2", "Column filtering", choices = "")),
                    column(6, selectInput("filtheatmapval2", "Value filtering", choices = "", multiple = TRUE)))
                ),
                fluidRow(
                  column(5,offset = 1, bsButton("add2filter_heat", label = HTML("&nbsp;Add"), style="success", icon("plus"))),
                  column(6, actionButton("gofilter_heat", "Apply filtering"))
                )
              ),

              hr(),
              awesomeCheckbox("logheat", "Log2 scale", value = FALSE),
              fluidRow(
                column(6,
                  selectInput("selscaleheat", "Standardize data:", 
                              choices = c("None" = "none", "By row" = "row", "By column" = "column"), 
                              selected = "column")
                ),
                conditionalPanel(condition = "input.selscaleheat == 'none'",
                  column(6,textInput("unitlegend_ht", "Unit measure", value = "ug/ml"))
                )
              ),
              hr(),
              h4(strong("Annotations")),
              fluidRow(
                column(6, selectInput("selectannot_row", "Row annotation:", choices = "")),
                column(6, selectInput("selectannot_col", "Column annotation:", choices = ""))
              ),
              
              hr(),
              h4(strong("Dendrogram options")),
              ###dendrogramm on column or row?
              h5(strong("Where to show dendrograms")),
              fluidRow(
                column(6, materialSwitch(inputId = "rowdend", label = "Row",  value = TRUE, status = "primary", width = "90%")),
                column(6, materialSwitch(inputId = "columndend", label = "Column",  value = TRUE, status = "primary", width = "90%"))
              ),

              conditionalPanel(condition = "input.rowdend == 1 || input.columndend == 1",
                fluidRow(
                  column(6,
                    selectInput("seldistheat", "Distance function:", choices = c("euclidean", "maximum", "canberra"), selected = "euclidean") #, "minkowski","manhattan",
                  ),
                  column(6,
                    selectInput("selhclustheat", "Clustering method:", choices = c("ward.D2", "complete", "average" , "median"), selected = "complete") #, "centroid","mcquitty","ward.D2", "single", 
                  )
                )
              ),
              
              conditionalPanel(condition = "input.rowdend == 0",
                               h5(strong("Order data by annotation?")),
                               awesomeCheckbox("heatsort", label = "Order", value = TRUE)
              ),

              conditionalPanel(condition = "input.rowdend == 1",
                hr(),
                uiOutput("sliderheatrow")
              ),
              
              conditionalPanel(condition = "input.columndend == 1",
                               hr(),
                               uiOutput("sliderheatcol")
            )
          ), #end of sidebarpanel
          
          mainPanel(width = 9, 
            InteractiveComplexHeatmap::InteractiveComplexHeatmapOutput("heatmap_output", layout = "1|(2-3)", width1 = 1000, height1 = 600)
            )
          )
        ), #end of tabpanel heatmap
        
        
        ##### Quality plots #####
        tabPanel(
          "Quality plots",
            sidebarLayout(
              sidebarPanel(
                width = 2,
                awesomeRadio("type_qualplots", "Quality plots on:", choices = c("Samples", "Lipids")),
                
                ###quality on samples
                conditionalPanel(
                  condition = "input.type_qualplots == 'Samples'",
                  awesomeRadio("selplot_samp", "Plot type:", choices = c("Barplot", "Boxplot", "Density plot")),
                  conditionalPanel(condition = "input.selplot_samp == 'Barplot' || input.selplot_samp == 'Boxplot'",
                    selectInput("densplot_sampfill", "Select fill column", choices = "")
                  ),
                  conditionalPanel(condition = "input.selplot_samp == 'Barplot'",
                    awesomeCheckbox("plotsamp_log", "Log scale", value = TRUE)
                  ),
                  
                  conditionalPanel(condition = "input.selplot_samp == 'Boxplot'",
                                   awesomeCheckbox("addpoints_boxdens", "Add points", value = FALSE)
                  ),
                  conditionalPanel(condition = "output.check_replicates == true",
                                   awesomeCheckbox("summ_qualplot", "Summarize data", value = FALSE)
                  )
                ),
                
                ###quality on lipids
                conditionalPanel(
                  condition = "input.type_qualplots == 'Lipids'",
                  awesomeRadio("selplot_lip", "Plot type:", choices = c("Barplot", "Boxplot")),
                  
                  #barplot
                  conditionalPanel(
                    condition = "input.selplot_lip == 'Barplot'",
                    awesomeCheckbox("concbar_log", "Log scale", value = TRUE),
                    selectInput("filtclass_concbar", "Filter by lipid class", choices = "", multiple = TRUE),
                    awesomeCheckbox("flip_concbar", "Reverse axis", value = FALSE),
                  ),
                  
                  #boxplots
                  conditionalPanel(condition = "input.selplot_lip == 'Boxplot'",
                                   selectInput("concbox_samp", "Select a sample", choices = ""),
                                   selectInput("concbox_sampfill", "Select fill column", choices = ""),
                                   awesomeCheckbox("concbox_log", "Log scale", value = TRUE),
                                   awesomeCheckbox("concbox_points", "Add points", value = FALSE),
                                   hr(),
                                   
                                   div(actionButton("second_concbox", "Add another boxplot"), align = "center"),
                                   conditionalPanel(
                                     condition = "input.second_concbox != 0",
                                     br(),
                                     selectInput("concbox_samp2", "Select a sample", choices = ""),
                                     selectInput("concbox_sampfill2", "Select fill column", choices = ""),
                                     awesomeCheckbox("concbox_log2", "Log scale", value = TRUE),
                                     awesomeCheckbox("concbox_points2", "Add points", value = FALSE)
                                   )
                  )
                )

              ), #end of sidebarpanel
              
              mainPanel(
                width = 10,
                conditionalPanel(
                  condition = "input.type_qualplots == 'Samples'",
                  conditionalPanel(condition = "input.selplot_samp == 'Barplot'",
                                   shinycssloaders::withSpinner(plotly::plotlyOutput("concbarplot_samp", height = "650px"))
                  ),
                  conditionalPanel(condition = "input.selplot_samp == 'Boxplot'",
                                   shinycssloaders::withSpinner(plotly::plotlyOutput("boxplot_density",height = "800px"))
                  ),
                  conditionalPanel(condition = "input.selplot_samp == 'Density plot'",
                                   shinycssloaders::withSpinner(plotly::plotlyOutput("densplot_samp"))
                  )
                ),
                
                conditionalPanel(
                  condition = "input.type_qualplots == 'Lipids'",
                  #barplot
                  conditionalPanel(condition = "input.selplot_lip == 'Barplot'",
                                   shinycssloaders::withSpinner(uiOutput("concbarplot_lip_UI"))
                  ),
                  
                  #boxplot
                  conditionalPanel(
                    condition = "input.selplot_lip == 'Boxplot'",
                    shinycssloaders::withSpinner(plotly::plotlyOutput("concboxplot_samp")),
                    conditionalPanel(condition = "input.second_concbox != 0",
                                     hr(),
                                     shinycssloaders::withSpinner(plotly::plotlyOutput("concboxplot_samp2"))
                    )
                  )
                ),

              ) #end of mainpanel
            ) #end of sidebarlayout

        ) #end of tabpanel quality plot
        )
      ), #end of tabitem graphtab
        

      #### PCA ####
      tabItem(tabName = "dimredtab",
        tabsetPanel(
          tabPanel("PCA",
            sidebarLayout(
              sidebarPanel(
                width = 2,
                fluidRow(
                  column(6,style="padding-right: 0px;",
                         bsplus::bs_embed_tooltip(awesomeCheckbox("scalepcalc", "Scale data", value = TRUE), 
                                            title = "If set to TRUE, PCA will be performed on correlation matrix, 
                                            otherwhise will be performed on covariance matrix.")),
                  column(6,awesomeCheckbox("logpcalc", "Log2 data", value = TRUE))),
                conditionalPanel(condition = "output.check_replicates == true",
                                 awesomeCheckbox("summarize_pca", "Summarize data", value = FALSE),
                ),
                div(actionButton("gopca", "Perform PCA", icon("cogs")), style = "text-align:center;"),
                hr(),
                selectInput("firstPC", "First PC", choices = ""),
                selectInput("secondPC", "Second PC", choices = ""),
                conditionalPanel(
                  condition = "input.selbiplotlc == 'Biplot'",
                  awesomeCheckbox("filt_biplot", "Filter first N lipids (biplot)", value = TRUE),
                  conditionalPanel(
                    condition = "input.filt_biplot == true",
                    bsplus::bs_embed_tooltip(
                      sliderInput("N_filt_biplot", "N lipids", value = 15, min = 1, max = 50),
                      "Biplot will show the first N lipids (arrows) for each PC in descending order of loadings."
                    )
                  )
                )
              ),
              
              mainPanel(width = 10,
                        tabsetPanel(
                          tabPanel("Plot",
                                   br(),
                                   fluidPage(
                                     fluidRow(
                                       column(3,box(width = NULL, status = "primary",
                                                    awesomeRadio("selbiplotlc", "Select plot type", choices = c("Plot", "Biplot"))
                                       )),
                                       column(4, box(width = NULL, status = "primary",
                                                     selectInput("colbiplotlc", "Select fill column", choices = ""))),
                                       column(4, box(width = NULL, status = "primary",
                                                     selectInput("shpbiplotlc", "Select shapes column", choices = "")
                                       ))
                                     ),
                                     fluidRow(shinycssloaders::withSpinner(plotly::plotlyOutput("biplotlc", height = "500px")))
                                   )
                          ),
                          
                          tabPanel("Scree plot",
                                   br(),
                                   shinycssloaders::withSpinner(plotly::plotlyOutput("screeplotlc"))
                          ),
                          
                          
                          tabPanel("Loadings",
                                   br(),
                                   fluidPage(
                                     fluidRow(
                                       column(4, box(width = NULL, status = "primary", uiOutput("sliderpclc"))),
                                       column(4, box(width = NULL, status = "primary", uiOutput("pcantopui")))),
                                     fluidRow(shinycssloaders::withSpinner(plotly::plotlyOutput("loadingslc", height = "550px"))))
                          ),
                          
                          
                          
                          tabPanel("Plot 3D",
                                   br(),
                                   fluidPage(
                                     fluidRow(
                                       column(4, box(width = NULL, status = "primary",
                                                     selectInput("col3dlc", "Select fill column", choices = "")))
                                     ),
                                     fluidRow(shinycssloaders::withSpinner(plotly::plotlyOutput("pca3dlc", height = "500px"))))
                          )
                          
                        )
              ) #end of mainpanel pca
            )
          ), #end of tabpanel pca
          
          ##### PLS-DA #####
          tabPanel("PLS-DA",
                   sidebarLayout(
                     sidebarPanel(
                       width = 2,
                       conditionalPanel(condition = "output.check_replicates == true",
                       awesomeCheckbox("summ_plsda", "Summarize data", value = FALSE)),
                       selectInput("selcol_pls", "Select a variable", choices = ""),
                       uiOutput("plsda_ncompui"),
                       div(actionButton("gopls", "Perform PLS-DA", icon("cogs")), style = "text-align:center;"),
                       hr(),
                       div(actionButton("tunepls", "Tune parameters", icon("wrench")), style = "text-align:center;"),
                       shinyBS::bsModal("view_tunepls", trigger = "tunepls", size = "large", title = "Tune parameters",
                                        sidebarLayout(
                                          sidebarPanel(
                                            numericInput("selfolds", "Folds", value = 5),
                                            numericInput("selnrepeat", "Number of times of Cross-Validation", value = 10),
                                            div(actionButton("gotunepls", "Run tuning", icon("cogs")), style = "text-align:center;")),
                                          mainPanel(
                                            shinycssloaders::withSpinner(plotOutput("perfplsda")),
                                            hr(),
                                            h4(strong("Choice of number of components:")),
                                            verbatimTextOutput("bestncomp_pls")
                                          )
                                        )
                       )
                     ),
                     
                     mainPanel(width = 10,
                               fluidRow(
                                 column(6,
                                   box(width = NULL, status = "primary", 
                                       fluidRow(
                                         column(
                                           width = 1,
                                           dropdownButton(
                                             awesomeRadio("pls_back_ellipse", "Plot with:", choices = c("Background", "Ellipse", "Nothing")),
                                             awesomeCheckbox("indname_pls", "Use rownames instead of shapes", value = T),
                                             circle = TRUE, status = "danger", icon = icon("cog"), width = "300px",
                                             tooltip = tooltipOptions(title = "Click to see options!")
                                           )            
                                         )),
                                     shinycssloaders::withSpinner(plotOutput("plotindiv_pls"))),
                                     ),
                                   
                                 column(6, 
                                   box(width = NULL, status = "primary", br(),br(),
                                       shinycssloaders::withSpinner(plotOutput("plotvar_pls")))
                                   )
                               )
                     )
                   )
          ),
          
          ##### sPLS-DA #####
          tabPanel(
            "sPLS-DA",
            sidebarLayout(
              sidebarPanel(
                width = 2,
                conditionalPanel(condition = "output.check_replicates == true",
                awesomeCheckbox("summ_splsda", "Summarize data", value = FALSE)),
                selectInput("selcol_spls", "Select a variable", choices = ""),
                uiOutput("splsda_ncompui"),
                bsplus::bs_embed_tooltip(
                  textInput("keepx_splsda", "KeepX", value = "50,50", placeholder = "Enter values separated by a comma..."),
                  "The number of variables to select on each component for sparse PLS-DA. You can click on 'Tune parameters' 
                  to get the best value for this option."),
                div(actionButton("gospls", "Perform sPLS-DA", icon("cogs")), style = "text-align:center;"),
                hr(),
                div(actionButton("tunespls", "Tune parameters", icon("wrench")), style = "text-align:center;"),
                shinyBS::bsModal("view_tunespls", trigger = "tunespls", size = "large", title = "Tune parameters",
                                 sidebarLayout(
                                   sidebarPanel(
                                     numericInput("selfolds_spls", "Folds", value = 5),
                                     numericInput("selnrepeat_spls", "Number of times of Cross-Validation", value = 10),
                                     div(actionButton("gotunespls", "Run tuning", icon("cogs")), style = "text-align:center;")
                                     ),
                                   mainPanel(
                                     shinycssloaders::withSpinner(plotOutput("plot_bestsparse")),
                                     hr(),
                                     h4(strong("Best number of components:")),
                                     verbatimTextOutput("bestncomp_spls"),
                                     h4(strong("Best KeepX:")),
                                     verbatimTextOutput("bestkeepx_spls")
                                   )
                                 )
                )
              ),
              
              mainPanel(width = 10,
                fluidRow(
                  column(
                    6,
                    box(width = NULL, status = "primary", 
                        fluidRow(
                          column(
                            width = 1,
                            dropdownButton(
                              awesomeRadio("spls_back_ellipse", "Plot with:", choices = c("Background", "Ellipse", "Nothing")),
                              awesomeCheckbox("indname_spls", "Use rownames instead of shapes", value = T),
                              circle = TRUE, status = "danger", icon = icon("cog"), width = "300px",
                              tooltip = tooltipOptions(title = "Click to see options!")
                            )
                          )),
                        shinycssloaders::withSpinner(plotOutput("plotindiv_spls"))),
                  ),
                  column(6, 
                         box(width = NULL, status = "primary", br(),br(),
                             shinycssloaders::withSpinner(plotOutput("plotvar_spls")))
                  )
                )
              )
            )
          )

        ) #end of tabsetpanel


      ), #end of tabitem pca




      ####### Clustering ######
      tabItem(tabName = "clusttab",
              
        sidebarLayout(
          sidebarPanel(width = 2,
            conditionalPanel(condition = "output.check_replicates == true",
            awesomeCheckbox("cluster_summar", "Summarize data", value = FALSE)),
            awesomeCheckbox("clust_scale", "Scale data", value = TRUE),
            #awesomeCheckbox("clust_pca", "Perform PCA", value = FALSE),
            awesomeRadio("clust_transp", "Clustering on:", choices = c("Lipids", "Samples")),
            radioGroupButtons("selclustmethod", "Clustering type", choices = c("Hierarchical", "Partitioning"), justified = TRUE, status = "primary"),
            conditionalPanel(condition = "input.selclustmethod == 'Partitioning'",
              selectInput("selclusthmorfo", "Clustering algorithm", choices = c("K-means", "PAM", "Clara"), selected = "K-means", multiple = FALSE),
            ),
            div(actionButton("tuneclust", "Tune parameters", icon("wrench")), style = "text-align: center;"),
            shinyBS::bsModal("view_tuneclust", trigger = "tuneclust", size = "large", title = "Tune parameters",
              shinycssloaders::withSpinner(plotOutput("numclustergraph", height = "800px"))
            )
          ),
          
          mainPanel(width = 10,
            fluidPage(
              fluidRow(
                column(4, 
                  box(width = NULL, status = "primary",
                      sliderInput("selnumclust", "Cluster number:", min = 1, max = 10, value = 2))),
                
                column(4, 
                  conditionalPanel(
                    condition = "input.selclustmethod == 'Hierarchical'",
                    box(width=NULL, status = "primary",
                        selectInput("selhclustmeth", "Agglomeration method", choices = c("single", "complete", "ward.D", "ward.D2"), selected = "ward.D2"))
                  )
                )
              ),
              conditionalPanel(
                condition = "output.check_clustering == true",
                h3(strong("Error in clustering! Probably there are too many missing values.", style = "color:red"))
              ),
              conditionalPanel(
                condition = "output.check_clustering == false",
                fluidRow(shinycssloaders::withSpinner(plotOutput("plotcluster", height = "600px")))
              )
            )
          )
        )
        
      ),
      
      

      ####### EXP design ######
      tabItem(tabName = "diffexptab",
        tabsetPanel(
          
          #### build DE
          tabPanel("Build DA",
            sidebarLayout(
              sidebarPanel(
                width = 3,
                #contrast list
                h4(strong("DA options")),
                conditionalPanel(condition = "output.check_replicates == true",
                  column(6, style="padding-left: 0px;", awesomeCheckbox("expdes_summar", "Summarize data", value = FALSE)),
                  conditionalPanel(
                    condition = "input.expdes_summar == false",
                    column(
                      6,
                      bsplus::bs_embed_tooltip(
                        title = "If you have technical replicates, you can also incorporate the replicate 
                        effect in the model by checking this box.",
                        awesomeCheckbox("expdes_repeffect", "Replicates effect", value = FALSE))
                      )
                  )
                ),
                

                awesomeRadio("expdes_bs_norm", "Normalization between replicates or samples", choices = "", inline = TRUE),
                #awesomeCheckbox("expdes_bs_norm", "Normalization between replicates or samples", value = FALSE),
                hr(),
                h4(strong("Variables")),
                
                fluidRow(
                  column(8, selectInput("expdes_design_vars", "Select primary variable", choices = "")),
                  column(4, br(), bsButton("add2var", label = HTML("&nbsp;Add"), style="success", icon("plus")))
                ),
                fluidRow(conditionalPanel(condition = "output.checkadd2var == 'twovar'",
                                          column(8, selectInput("expdes_design_vars2", "Select secondary variable", choices = "")))
                ),
                hr(),
                
                h4(strong("Batch effect")),
                awesomeCheckbox("expdes_batch_effect", "Batch effect", value = TRUE),
                conditionalPanel(
                  condition = "input.expdes_batch_effect == true",
                  fluidRow(
                    column(
                      6, 
                      bsplus::bs_embed_tooltip(
                        selectInput("expdes_batch_type", "Batch type", choices = c("remove", "fit")),
                        title = "ADViSELipidomics copes with the batch effects by either fitting the model with the batch variables 
                        or removing the batch effect before fitting the model. If you select 'fit', the software requires your batch variable
                        also in the contrasts list. Remember that you can generate a contrasts list up to two variables.")),
                    column(
                      6, 
                      bsplus::bs_embed_tooltip(
                        title = "Algorithm used for the batch effect. If 'limma',it will be used the 'removeBatchEffect' function, otherwhise
                        the 'ComBat' function from SVA package (parametric or non-parametric)",
                        selectInput("batch_meth", "Batch method", choices = "")))
                  ),
                  conditionalPanel(
                    condition = "input.expdes_batch_type == 'remove' || input.batch_meth == 'given'",
                    fluidRow(
                      column(8, selectInput("expdes_batch_var1", "Primary batch variable", choices = "")),
                      column(4, br(), bsButton("add2batch", label = HTML("&nbsp;Add"), style="success", icon("plus")))
                    ),
                    fluidRow(
                      conditionalPanel(condition = "output.checkadd2batch == 'twobatch'",
                                       column(8, selectInput("expdes_batch_var2", "Second batch variable", choices = ""))
                      )
                    )
                  )
                ),
                hr(),
                
                div(actionButton("open_writecontrast", "Write contrasts", icon("pen")), style = "text-align: center;"),
                br(),
                
                tags$head(tags$style("#view_writecontrast .modal-dialog{ width:1100px}")),
                shinyBS::bsModal("view_writecontrast", trigger = "open_writecontrast", size = "large", title = "Write the contrasts",
                                   fluidPage(
                                     column(
                                       3,
                                       box(width = NULL, status = "primary", title = "First term", solidHeader = TRUE,
                                           selectInput("firstelement", "Level 1", choices = ""),
                                           conditionalPanel(condition = "output.checkntotvars == 1",
                                                            selectInput("secondelement1", "Level 2", choices = "")
                                           ),
                                           conditionalPanel(condition = "output.checkntotvars == 2",
                                                            selectInput("secondelement2", "Level 2", choices = "")
                                           )
                                       ),
                                       conditionalPanel(condition = "output.checkntotvars == 2",
                                                        box(width = NULL, status = "primary", title = "Second term", solidHeader = TRUE,
                                                            radioButtons("fixedlevel", "Fixed Level", choices = c("Level 1", "Level 2"), inline = TRUE), #, checkbox = TRUE
                                                            conditionalPanel(condition = "input.fixedlevel == 'Level 1'",
                                                                             br(),
                                                                             h5(strong("Level 1: fixed")), 
                                                                             br(),
                                                                             selectInput("thirdelement2", "Level 2", choices = "")
                                                            ),
                                                            conditionalPanel(condition = "input.fixedlevel == 'Level 2'",
                                                                             selectInput("thirdelement1", "Level 1", choices = ""),
                                                                             h5(strong("Level 2: fixed"))
                                                            )
                                                        )
                                                        
                                       ),
                                       hr(),
                                       # Row selection
                                       numericInput("row.selection", "Select row to be deleted", min = 1, max = 100, value = ""),
                                       # Add button
                                       actionButton("add.button", "Add", icon("plus")), 
                                       # Delete button 
                                       actionButton("delete.button", "Delete", icon("minus"))
                                     ),
                                     column(8,
                                            DTOutput("tablewritten"),
                                            br(), 
                                            hr(), 
                                            br(),
                                            fluidRow(
                                              column(9,verbatimTextOutput("printconlist")),
                                              column(3, actionButton("submit_written", "Submit", icon("file-import"), style = "width: 150px;height: 50px;font-size: 20px;"))
                                            )
                                     )
                                   )

                ),
                
                prettyRadioButtons("expdes_thresh", "p-value threshold", choices = c(0.001, 0.01, 0.05), inline = TRUE, selected = 0.05),
                bsplus::bs_embed_tooltip(title = "A method for the decideTests function (limma package) used to identify 'differentially expressed' lipids.",
                selectInput("expdes_decide_met", "Select method", choices = c("separate", "global", "hierarchical", "nestedF"), selected = "separate")),
                fluidRow(
                  column(
                    6,
                    div(actionButton("runDE", "Run DA", icon("cogs"), style = "font-size:140%;"), style = "text-align:center;"),
                  ),
                  column(
                    6,
                    conditionalPanel(
                      condition = "output.checktoptable == false",
                      div(actionButton("opentoptable", "Check TopTable", icon("eye")), style = "text-align:center;")
                    )
                  )
                ),
                tags$head(tags$style("#view_toptable .modal-dialog{ width:1300px}")),
                tags$head(tags$style("#view_toptable .modal-body{ min-height:800px}")),
                shinyBS::bsModal("view_toptable", trigger = "opentoptable", size = "large", title = "Top Table",
                                 fluidRow(
                                   column(3,
                                          box(width = NULL, status = "primary", 
                                              selectInput("sel_toptable", "Select contrast", choices = ""),
                                              downloadButton("downloadtoptable", "Download"))),
                                   column(9, box(width = NULL, status = "primary",div(DT::DTOutput("toptable"), style = "overflow-x: scroll;")))
                                 )
                                 ),
                
                
                hr(),
                selectInput("expdes_colmaplot", "Select contrast", choices = ""),
                sliderInput("expdes_lfc", "Select logFC", min = 1, max = 10, value = 2, step = 0.5),
                hr(),
                h4(strong("Plot options")),
                div(actionButton("addmorevolcanos", "Add another plot"), align = "center")
                
              ),
              
              mainPanel(width = 9,
                        mod_ma_volcano_plot_ui("ma_volcano_plot1"),
                        
                        conditionalPanel(condition = "input.addmorevolcanos != 0",
                                         mod_ma_volcano_plot_ui("ma_volcano_plot2")
                        )
              ) #end of mainpanel
              
              
            ) #end of sidebarLayout
          ),
          
          ##### Comparisons of DE (venn etc)
          tabPanel("Comparisons",
            sidebarLayout(
              sidebarPanel(
                width = 2,
                radioGroupButtons("sel_venn_upset", "Choose a graph :",
                                  choiceValues = list("Venn Diagram", "Upset plot"),
                                  choiceNames = list(
                                    paste(tags$img(src = "www/venn_icon.png", height = "21px", width = "21px"), HTML("<b style=font-size:16px>&nbsp;Venn Diagram</b>")),
                                    paste(tags$img(src = "www/upset_icon.png", height = "20px", width = "20px"), HTML("<b style=font-size:16px>&nbsp;Upset Plot</b>"))
                                  ),
                                  direction = "vertical",justified = TRUE, status = "primary"
                ),

                conditionalPanel(condition = "input.sel_venn_upset == 'Venn Diagram'",
                  selectInput("venncontrast", "Select contrasts (min 2 max 5)", choices = "", multiple = TRUE),
                  awesomeCheckbox("venngrad", "Gradient", value = FALSE)
                ) 
              ),
              
              mainPanel(width = 10,
                conditionalPanel(condition = "input.sel_venn_upset == 'Venn Diagram'",
                  shinycssloaders::withSpinner(uiOutput("venndiag"))) ,
                
                conditionalPanel(condition = "input.sel_venn_upset == 'Upset plot'",
                  fluidPage(
                    column(width = 10, offset = 1,br(), shinycssloaders::withSpinner(plotOutput("upsetplot")))
                  )
                )
              )
            )
            
            
          ) #end of tabpanel comparisons
        )
      ),
      

      
      ##### Enrichment ####
      tabItem(tabName = "enrmenu",
        sidebarLayout(
          sidebarPanel(width = 2,
            selectInput("rank_c", "rank_c", choices = c("logFC" , "P.Value", "adj.P.Val", "B")),
            prettyRadioButtons("enrich_thresh", "p-value threshold", choices = c(0.001, 0.01, 0.05), selected = 0.05),
            selectInput("enrich_selcont", "Select contrast", choices = "")
          ),
          
          mainPanel(width = 10,
            shinycssloaders::withSpinner(plotOutput("enrich_plot", height = "900px"))
          )
        )
        )
      
      
      
      
      )
    )
    ) #end of dashboardpage
    

  )
}

#' Add external Resources to the Application
#' 
#' This function is internally used to add external 
#' resources inside the Shiny application. 
#' 
#' @import shiny
#' @importFrom golem add_resource_path activate_js favicon bundle_resources
#' @noRd
golem_add_external_resources <- function(){
  
  add_resource_path(
    'www', app_sys('app/www')
  )
 
  tags$head(
    favicon(),
    bundle_resources(
      path = app_sys('app/www'),
      app_title = 'ADViSELipidomics'
    )
    # Add here other external resources
    # for example, you can add shinyalert::useShinyalert() 
  )
}

