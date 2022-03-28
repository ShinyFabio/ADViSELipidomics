#' reading_step UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList
#' @importFrom shinyBS bsModal
#' @importFrom DT DTOutput

mod_reading_step_ui <- function(id){
  ns <- NS(id)
  tagList(
    
    ###### Box Filtering #####
    column(
      5,
      box(width = NULL, title = "1. IMPORTING & FILTERING", status = "primary", solidHeader = TRUE,
          
          # Step 1.
          div(icon("circle"), HTML("&nbsp;Step. 1: Import targetfile and internal reference files (.xlsx)"), style = "text-align: left; font-size: 18px; font-weight: bold;"),
          br(),
          fluidRow(
            column(10, fileInput(ns("targetfilepath"), "Select the Targetfile Lipidomics (.xlsx)", accept = ".xlsx")),
            column(2, br(), 
                   mod_edit_data_ui(ns("edit_target")))
            
          ),
          fluidRow(
            column(10, fileInput(ns("internalstdpath"), "Select the Internal Reference file (.xlsx)", accept = ".xlsx")),
            column(2, 
                   br(),
                   mod_edit_data_ui(ns("edit_internal"))) 
          ),
          hr(),
          
          #Step 2.
          conditionalPanel(
            condition = "output.checktargets == false", ns = ns,
            div(icon("circle"), HTML("&nbsp;Step. 2: Choose the data folder and read all the data."), style = "text-align: left; font-size: 18px; font-weight: bold;"),
            br(),
            fluidRow(
              div(
                column(4, shinyFiles::shinyDirButton(ns("datafolder"), label = "Browse...", title = "Please select the data folder", multiple = FALSE, icon = icon("folder-open"))),
                column(4, actionButton(ns("readdatabttn"), "Read Data", icon("cogs"))),
                column(4, 
                       conditionalPanel(
                         condition = "output.checkdatafiles == false", ns = ns,
                         actionButton(ns("showqualcheck"), "Quality check", icon("eye")))
                ), style = "text-align: center;")
            ),
            tags$head(tags$style(paste0("#", ns("viewqualcheckplot"), " .modal-dialog{ width:1400px}"))),
            tags$head(tags$style(paste0("#"), ns("viewqualcheckplot")," .modal-body{ min-height:800px}"))),
            shinyBS::bsModal(ns("viewqualcheckplot"), trigger = ns("showqualcheck"), size = "large", title = "Quality check on area",
                             conditionalPanel(
                               condition = "output.checkint_std == 'Yes'", ns = ns,
                               radioGroupButtons(ns("selqualplot"), "Data", choices = c("nonlabeled", "deuterated"), individual = TRUE, 
                                                 checkIcon = list(yes = tags$i(class = "fa fa-circle",  style = "color: steelblue"),
                                                                  no = tags$i(class = "fa fa-circle-o", style = "color: steelblue")))
                             ), 
                             plotlyOutput(ns("qualcheckplot"),height = "600px")),
            hr(),
          
          
          #Step 3.
          conditionalPanel(
            condition = "output.checkdatafiles == false", ns = ns,
            div(icon("circle"), HTML("&nbsp;Step. 3: Select a range for the carbon number and the double bound number"), style = "text-align: left; font-size: 18px; font-weight: bold;"),
            br(),
            sliderInput(ns("ca_bound"), div("Range carbon number", style = "font-size: 15px; font-weight: bold;"), min = 0, max = 50, value = c(14, 24)),
            sliderInput(ns("db_bound"), div("Range double bound number", style = "font-size: 15px; font-weight: bold;"), min = 0, max = 10, value = c(0, 6)),
            fluidRow(
              column(6,
                     div(actionButton(ns("gofilterbttn"), "Filter Data", class = "btn-primary btn-lg", icon("cogs"), width = "160px", style='height:40px; padding:5px; font-size:140%; font-weight: bold;'), 
                         style = "text-align:center;")),
              column(6,
                     conditionalPanel(
                       condition = "output.checkendfiltering == false", ns = ns,
                       div(actionButton(ns("viewfiltered"), "Check filtered data", class = "btn-primary btn-lg", icon("eye"), style='width:230px; height:40px; padding:5px; font-size:140%; font-weight: bold;'), 
                           style = "text-align:center;"),
                       tags$head(tags$style(paste0("#", ns("viewfilteredmodal")," .modal-dialog{ width:1100px}"))),
                       shinyBS::bsModal(ns("viewfilteredmodal"), trigger = ns("viewfiltered"), size = "large",
                                        selectInput(ns("selcol_lipidfilt"), "Select a sample", choices = ""),
                                        div(DT::DTOutput(ns("dtlipidfilterd")), style = "overflow-x: scroll;"))
                                    )
                             )
                           ) #end of fluidrow
                           
          )
          
      ) #end of box filtering
    ) #end of column for the first part "filtering"
 
  )
}
    
#' reading_step Server Functions
#'
#' @importFrom shinyFiles shinyDirButton getVolumes shinyDirChoose parseDirPath
#' @importFrom shinyWidgets show_alert
#' @importFrom readxl read_xlsx
#' @importFrom lubridate as_date
#' @importFrom tools file_ext
#' @importFrom DT renderDT
#' @importFrom purrr pluck
#' @noRd 
mod_reading_step_server <- function(id, analysis, int_std){
  moduleServer( id, function(input, output, session){
    ns <- session$ns
    ##### STEP 1.
    ####import target and internal files xlsx
    
    #targetfile
    targetfile_to_edit = reactive({
      req(input$targetfilepath)
      ext <- tools::file_ext(input$targetfilepath$name)
      if(ext != "xlsx"){
        shinyWidgets::show_alert("Invalid file!", "Please upload a .xlsx file", type = "error")
      }
      validate(need(ext == "xlsx", "Invalid file! Please upload a .xlsx file"))
      x = readxl::read_xlsx(input$targetfilepath$datapath, na = c("", "NA"))
      #x$Exp_date = lubridate::as_date(x$Exp_date)
      return(x)
    })
    
    
    targetfile_edit = mod_edit_data_server("edit_target", data_input = targetfile_to_edit)
    
    
    #internal standard
    internalstd_to_edit = reactive({
      req(input$internalstdpath)
      ext <- tools::file_ext(input$internalstdpath$name)
      if(ext != "xlsx"){
        shinyWidgets::show_alert("Invalid file!", "Please upload a .xlsx file", type = "error")
      }
      validate(need(ext == "xlsx", "Invalid file! Please upload a .xlsx file"))
      readxl::read_xlsx(input$internalstdpath$datapath, na = c("", "NA"))
    })
    
    internalstd_edit = mod_edit_data_server("edit_internal", data_input = internalstd_to_edit)
    
    
    
    
    #unisco
    targets = reactive({
      req(targetfile_edit(), internalstd_edit())
      list(
        targetfile_lipidomics = targetfile_edit(),
        internal_standard = internalstd_edit()
      )
    })
    
    
    
    
    ###### STEP 2.
    #se tutti i file sono stati caricati
    output$checktargets = reactive(
      return(is.null(targets()))
    )
    outputOptions(output, "checktargets", suspendWhenHidden = FALSE)
    
    
    ####select the data folder where to read data
    
    # this makes the directory at the base of your computer.
    volumes = c(Home = fs::path_home(), shinyFiles::getVolumes()())
    
    shinyFiles::shinyDirChoose(input, 'datafolder', roots = volumes, session = session)
    
    data_path = reactive({
      if(length(input$datafolder) != 1 ) {
        shinyFiles::parseDirPath(volumes,input$datafolder)
      }else{
        NULL
      }
    })
    
    
    
    stepb = eventReactive(input$readdatabttn,{
      req(data_path(), targets(), analysis())
      
      #creo il primo stepa (primo out)
      stepa = list(targets = targets(), analysis = analysis())
      withProgress(message = "Reading data...", value=0, {
        
        read_advise_lipidomics(out = stepa, datapath = data_path(), target_file = stepa$targets$targetfile_lipidomics)
      })
      
    })
    
    #stepb prende il percorso dei lipid_data (es. AF-1C-M_deuterated_1.txt) e il targetfile e carica tutti i file.
    #usare così:
    #stepb()$lipid_data
    #stepb()$replicates  dove dice TRUE se ci sono replicati tecnici o FALSE se non ci sono
    
    
    #se tutti i data sono stati caricati
    output$checkint_std = reactive(int_std())
    outputOptions(output, "checkint_std", suspendWhenHidden = FALSE)
    
    
    ####quality check delle aree
    output$qualcheckplot = renderPlotly({
      req(stepb())
      
      if(int_std() == "No"){
        lipid_nonfilt = stepb()$lipid_data
      }else{
        lipid_nonfilt = stepb()$lipid_data[grepl(input$selqualplot, names(stepb()$lipid_data), fixed = TRUE)]
      }

      area = lapply(lipid_nonfilt, function(x) log2(sum(x$Area))) %>% tibble::as_tibble() %>% t() %>% as.data.frame() %>% tibble::rownames_to_column()
      colnames(area) = c("Samples", "Area")
      
      hh = ggplot(area) + geom_col(aes(x= Samples, y = Area)) + ylab("Log2(Area)") + 
        theme(axis.text.x = element_text(angle = 270, hjust = 0, size = 8))
      ggplotly(hh)
    })
    
    
    #######filtering options
    
    #se tutti i data sono stati caricati
    output$checkdatafiles = reactive(
      return(is.null(stepb()))
    )
    outputOptions(output, "checkdatafiles", suspendWhenHidden = FALSE)
    
    stepc = eventReactive(input$gofilterbttn,{
      req(stepb())
      filter_advise_lipidomics(
        out = stepb(),
        ca_bound = input$ca_bound,
        db_bound = input$db_bound
      )
    })
    
    #stepc() prende gli internal standard in targets() e i lipid_data in readed_data() e li filtra.
    #usare così:
    #stepc()$lipid_filtered
    #stepc()$lipid_deuterated
    
    
    ###check end of filtering box
    output$checkendfiltering = reactive(
      return(is.null(stepc()))
    )
    outputOptions(output, "checkendfiltering", suspendWhenHidden = FALSE)
    
    
    observeEvent(stepc(),{
      updateSelectInput(session, "selcol_lipidfilt", choices = names(stepc()$lipid_filtered))
    })
    
    output$dtlipidfilterd = DT::renderDT({
      req(stepc())
      stepc()$lipid_filtered %>% purrr::pluck(input$selcol_lipidfilt)
    })
    
    
    
    return(reactive({
      req(stepc())
      return(stepc())
    }))
    
    
  })
}
    
## To be copied in the UI
# mod_reading_step_ui("reading_step_ui_1")
    
## To be copied in the server
# mod_reading_step_server("reading_step_ui_1")
