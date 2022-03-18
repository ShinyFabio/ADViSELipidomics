#' imputation_step UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList
#' @importFrom stats prcomp
mod_imputation_step_ui <- function(id, offset){
  ns <- NS(id)
  tagList(
    
    #fluidRow(
    #step filtering and imputation
    column(5, offset = offset,
           box(width = NULL, title = "3. MISSING DATA & SUMMARIZED EXPERIMENT", status = "primary", solidHeader = TRUE,
               # Filtering
               div(icon("circle"), HTML("&nbsp;Step. 1: Missing data filtering"), style = "text-align: left; font-size: 18px; font-weight: bold;"),
               br(),
               fluidRow(
                 column(5, h4("Do you want to filter data/NA? "), style="padding-right: 0px; width: 275px;"),
                 column(4, prettyToggle(ns("togglefilterna"), label_on = "Yes!", 
                                        label_off = "No..", outline = TRUE, plain = TRUE, bigger = TRUE, 
                                        icon_on = icon("thumbs-up"), icon_off = icon("thumbs-down")),
                        style="padding-top: 8px;padding-left: 0px;")
               ),
               
               fluidRow(
                 conditionalPanel(condition = "input.togglefilterna == true", ns = ns,
                                  column(
                                    7 ,offset = 1,
                                    sliderInput(ns("na_filt_lip"), "Max missing data percentage allowable on lipids", min = 0, max = 1, step = 0.1, value = 0.3),
                                    sliderInput(ns("na_filt_sam"), "Max missing data percentage allowable on samples", min = 0, max = 1, step = 0.1, value = 0.6)
                                  ),
                                  
                                  column(
                                    4,
                                    br(), br(), br(),
                                    div(actionButton(ns("gofilterna"), "Check filtered NA", class = "btn-primary btn-lg", icon("eye"), style="white-space: normal; height:60px; width:140px;"),
                                        style = "text-align:center;"),
                                  )
                 ),
                 tags$head(tags$style(paste0("#", ns("modalfilteredna")," .modal-dialog{ width:1300px}"))),
                 tags$head(tags$style(paste0("#", ns("modalfilteredna")," .modal-body{ min-height:1000px}"))),
                 shinyBS::bsModal(ns("modalfilteredna"), trigger = ns("gofilterna"),
                                  radioGroupButtons(ns("vim_or_dt_filtermodal"), "Choose what to show :",
                                                    choiceValues = list("table", "vim"),
                                                    choiceNames = list(
                                                      paste(shiny::icon("table",  style='font-size:16px;'), HTML("<b style=font-size:16px>&nbsp;Table</b>")),
                                                      paste(shiny::icon("buromobelexperte",  style='font-size:16px;'), HTML("<b style=font-size:16px>&nbsp;NA plots</b>"))),
                                                    justified = TRUE, status = "primary"),
                                  conditionalPanel(condition = "input.vim_or_dt_filtermodal == 'table'", ns = ns,
                                                   div(DT::DTOutput(ns("dt_filteredna")), style = "overflow-x: scroll;"),
                                                   br(),br(),
                                                   fluidPage(
                                                     fluidRow(
                                                       column(3, div(shinydashboard::valueBoxOutput(ns("nadim1"),width = NULL)), style="padding-right: 0px;"),
                                                       column(1, tags$img(src = "www/right_arrow.png", width = "90px")),
                                                       column(3, div(shinydashboard::valueBoxOutput(ns("nadim2"),width = NULL)), style="padding-right: 0px;padding-left: 0px;")
                                                     )
                                                   )
                                                   
                                  ),
                                  conditionalPanel(condition = "input.vim_or_dt_filtermodal == 'vim'", ns = ns,
                                                   switchInput(ns("vimopz_filt"), label = "Combined", value = TRUE),
                                                   plotOutput(ns("vimplot_filteredna"), height = "1000px")
                                  )
                 ) #end of bsModal
                 
               ),
               
               hr(),
               
               # Imputation
               div(icon("circle"), HTML("&nbsp;Step. 2: Imputation"), style = "text-align: left; font-size: 18px; font-weight: bold;"),
               br(),
               
               fluidRow(
                 column(5, h4("Do you want to impute missing data?"), style="padding-right: 0px; width: 330px;"),
                 column(4, prettyToggle(ns("toggleimputena"), label_on = "Yes!", label_off = "No..", outline = TRUE, plain = TRUE, bigger = TRUE, icon_on = icon("thumbs-up"), icon_off = icon("thumbs-down")),
                        style="padding-top: 8px;padding-left: 0px;")
               ),
               
               fluidRow(
                 conditionalPanel(condition = "input.toggleimputena == true", ns = ns,
                                  column(6, offset = 1,
                                         selectInput(ns("imput_method"), "Imputation method", choices = c("mean", "median", "none", "knn","irmi"), selected = "knn"),
                                  )
                                  ,
                                  column(4, offset = 1,
                                         div(actionButton(ns("goimputena"), "Check imputed NA", class = "btn-primary btn-lg", icon("eye"), style="white-space: normal; height:60px; width:140px;"),
                                             style = "text-align:center;")
                                  )
                 ),
                 
                 tags$head(tags$style(paste0("#", ns("modalimputena")," .modal-dialog{ width:1300px}"))),
                 tags$head(tags$style(paste0("#"), ns("modalimputena")," .modal-body{ min-height:1000px}"))),
                 shinyBS::bsModal(ns("modalimputena"), trigger = ns("goimputena"),
                                  radioGroupButtons(ns("vim_or_dt_imputemodal"), "Choose what to show :",
                                                    choiceValues = list("table", "vim"),
                                                    choiceNames = list(
                                                      paste(shiny::icon("table",  style='font-size:16px;'), HTML("<b style=font-size:16px>&nbsp;Table</b>")),
                                                      paste(shiny::icon("buromobelexperte",  style='font-size:16px;'), HTML("<b style=font-size:16px>&nbsp;NA plots</b>"))),
                                                    justified = TRUE, status = "primary"),
                                  conditionalPanel(condition = "input.vim_or_dt_imputemodal == 'table'", ns = ns,
                                                   div(DT::DTOutput(ns("dt_imputena")), style = "overflow-x: scroll;")
                                  ),
                                  conditionalPanel(condition = "input.vim_or_dt_imputemodal == 'vim'", ns = ns,
                                                   switchInput(ns("vimopz_imput"), label = "Combined", value = TRUE),
                                                   plotOutput(ns("vimplot_imputena"), height = "1000px")
                                  )
                 ), #end of bsModal
               
               hr(),
               
               
               # finish and sumexp
               div(icon("circle"), HTML("&nbsp;Step. 3: Creation of Summarized Experiment"), style = "text-align: left; font-size: 18px; font-weight: bold;"),
               fluidRow(
                 column(7, br(),
                        div(
                          actionButton(ns("gofinish"), "Finish!", class = "btn-primary btn-lg", icon("rocket"), width = "310px", style='height:45px; padding:5px; font-size:140%; font-weight: bold;'),
                          style = "text-align:center;"), br()
                 ),
                 column(5,
                        conditionalPanel(
                          condition = "output.checksteph == false", ns = ns,
                          div(
                            actionButton(ns("gotosumexp"), "See the results", icon("eye")),
                            br(),br(),
                            downloadButton(ns("downloadsumexp"), "Download"),
                            style = "text-align:center;")
                        )
                 )
                 
               )
           ) #end of box
    ) #end of column
  #)
  )
}
    
#' imputation_step Server Functions
#'
#' @noRd 
mod_imputation_step_server <- function(id,parent, stepg, data_type){
  moduleServer( id, function(input, output, session){
    ns <- session$ns
 
    ##### FILTERING
    
    filterstep = eventReactive(input$gofilterna,{
      na_advise_lipidomics(
        out = stepg(),
        na_filter_lip = as.numeric(input$na_filt_lip),
        na_filter_sam = as.numeric(input$na_filt_sam),
        imputation_met = "none"
      )
    })
    
    

    output$dt_filteredna = DT::renderDT({
      req(filterstep())
      filterstep()$concentration_matrix_filt
    })
    
    output$vimplot_filteredna = renderPlot({
      req(filterstep())
      filterstep()$concentration_matrix_filt %>% dplyr::select(-1) %>% VIM::aggr(combined = input$vimopz_filt, cex.axis = 0.7, cex.lab = 0.7)
    })
    
    #box table dimensions after na filtering
    output$nadim1 = shinydashboard::renderInfoBox({
      dim = dim(filterstep()$concentration_matrix[,-1])
      shinydashboard::infoBox(
        title = div(HTML(paste0("Table dimension", br(), "before filtering")), style = "color:white; font-size:100%;"),
        value = div(paste0(dim[1], " x ", dim[2]), style = "font-size:140%"),
        icon = icon("table"), color = "yellow", fill = TRUE)
    })
    
    output$nadim2 = shinydashboard::renderInfoBox({
      dim = dim(filterstep()$concentration_matrix_filt[,-1])
      shinydashboard::infoBox(
        title = div(HTML(paste("Table dimension", br(), "after filtering")), style = "color:white; font-size:100%;"),
        value = div(paste0(dim[1], " x ", dim[2]), style = "font-size:140%"),
        icon = icon("table"), color = "green", fill = TRUE)
    })
    
    
    ##### IMPUTATION

    imputestep = eventReactive(input$goimputena,{
      if(input$togglefilterna == TRUE){
        na_lip = as.numeric(input$na_filt_lip)
        na_sam = as.numeric(input$na_filt_sam)
      }else{
        na_lip = 1
        na_sam = 1
      }

      na_advise_lipidomics(
        out = stepg(),
        na_filter_lip = na_lip,
        na_filter_sam = na_sam,
        imputation_met = input$imput_method
      )
    })


    output$dt_imputena = renderDT({
      req(imputestep())
      imputestep()$concentration_matrix_filt
    })

    output$vimplot_imputena = renderPlot({
      req(imputestep())
      imputestep()$concentration_matrix_filt %>% dplyr::select(-1) %>% VIM::aggr(combined = input$vimopz_imput, cex.axis = 0.7,cex.lab = 0.7)
    })

    
    
    
    #steph()$concentration_matrix_filt
    #steph()$sample_filtered
    
    
    
    steph = eventReactive(input$gofinish,{
      req(stepg())
      
      if(input$togglefilterna == FALSE && input$toggleimputena == FALSE){
        data = na_advise_lipidomics(
          out = stepg(),
          na_filter_lip = 1,
          na_filter_sam = 1,
          imputation_met = "none"
        )
      }else if(input$togglefilterna == TRUE && input$toggleimputena == FALSE){
        data = na_advise_lipidomics(
          out = stepg(),
          na_filter_lip = as.numeric(input$na_filt_lip),
          na_filter_sam = as.numeric(input$na_filt_sam),
          imputation_met = "none"
        )
      }else if(input$togglefilterna == FALSE && input$toggleimputena == TRUE){
        data = na_advise_lipidomics(
          out = stepg(),
          na_filter_lip = 1,
          na_filter_sam = 1,
          imputation_met = input$imput_method
        )
      }else{
        data = na_advise_lipidomics(
          out = stepg(),
          na_filter_lip = as.numeric(input$na_filt_lip),
          na_filter_sam = as.numeric(input$na_filt_sam),
          imputation_met = input$imput_method
        )
      }

      sum = sumexp_advise_lipidomics(out = data)
      
      g1 = sum$sumexp_data %>% SummarizedExperiment::assay() %>% t()
      g = try(stats::prcomp(g1))
      if(class(g) == "try-error"){
        shinyWidgets::sendSweetAlert(session, title = "Too many NAs", type = "warning", width = "600px",
                                     text = "If you will proceed, PCA won't work.")
      }
      return(sum)
    })
    
    #check all data correctly loaded
    output$checksteph = reactive(
      tryCatch({
        is.null(steph())
      },
      shiny.silent.error = function(e) {
        return(TRUE)
      })
    )
    outputOptions(output, "checksteph", suspendWhenHidden = FALSE)
    
    
    
    #### Download handler for the download button
    output$downloadsumexp <- downloadHandler(
      #put the file name with also the file extension
      filename = function() {
        paste0("summ_EXP_", Sys.Date(), ".rds")
      },
      
      # This function should write data to a file given to it by the argument 'file'.
      content = function(file) {
        summ = list(sumexp_data = steph()$sumexp_data, sumexp_data_mean = steph()$sumexp_data_mean, replicates = steph()$replicates, data_type = data_type)
        saveRDS(summ, file)
      }
    )
    
    observeEvent(input$gotosumexp, {
      shinydashboard::updateTabItems(session = parent, "sidebarmenu", "seedatatab")
    })
    
    
    return(reactive({
      req(steph())
      return(steph())
    }))
    
    
  })
}
    
## To be copied in the UI
# mod_imputation_step_ui("imputation_step_ui_1")
    
## To be copied in the server
# mod_imputation_step_server("imputation_step_ui_1")
