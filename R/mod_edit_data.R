#' edit_data UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#' @param data_input put the internal/default data
#'
#' @noRd 
#'
#' @import shiny
#' @import DataEditR
#' @import shinyBS
#' @importFrom shinyjs useShinyjs
#' @importFrom miniUI gadgetTitleBar miniTitleBarCancelButton miniTitleBarButton
#' @importFrom readr write_excel_csv
#' @importFrom bsplus use_bs_tooltip bs_embed_tooltip
#' @importFrom cicerone Cicerone use_cicerone
#' 
#' 



mod_edit_data_ui <- function(id){
  ns <- NS(id)
  tagList(
    fluidPage(
    useShinyjs(),
    cicerone::use_cicerone(),

    bsplus::use_bs_tooltip(),    #mi serve per il bs_embed_tooltip()
    actionButton(ns("editinternal"), icon("edit")),
    tags$head(tags$style(paste0("#", ns("upinternalmodal")," .modal-dialog{ width:1300px}"))),
    tags$head(tags$style(paste0("#", ns("upinternalmodal"), " .modal-body{ min-height:1000px}"))),
    shinyBS::bsModal(ns("upinternalmodal"),
      title = miniUI::gadgetTitleBar(title = "Data Editor", 
                                     right = bsplus::bs_embed_tooltip(miniUI::miniTitleBarButton(ns("done"), "Done", primary = TRUE), title = "Apply editing"),
                                     left = miniUI::miniTitleBarCancelButton(ns("cancel"), label = "Cancel", primary = FALSE)
                                     ), 
      trigger = ns("editinternal"), size = "large", 

      fluidRow(
        column(5, style = "padding-left: 5px;",
          div(bsplus::bs_embed_tooltip(dataSelectUI(ns("select1")), title = "Select columns"), 
              id = ns("idselect1"),style = "display:inline-block;"), 
          div(bsplus::bs_embed_tooltip(dataFilterUI(ns("filter1")), title = "Filter rows by condition"),
              id = ns("idfilter1"),style = "display:inline-block;"), 
          div(bsplus::bs_embed_tooltip(dataOutputUI(ns("output-active")), title = "Download edited data"),
              id = ns("idoutput-active"),style = "display:inline-block;"),
           #il bsTooltip con lo switch non funziona e quindi uso questo
          div(
            bsplus::bs_embed_tooltip(switchInput(ns("cut"), value = FALSE), title = "Enable/disable editing"), 
            style = "display:inline-block;", id = ns("idcut")),
            
        ),
         column(1,offset = 6, 
           div(bsplus::bs_embed_tooltip(shinyWidgets::actionBttn(inputId = ns("help_button"), size = "sm",style = "material-circle",icon("question"), color = "primary"),
                                        title = "Help"), style = "text-align: right;")
         )
      ), 
      
      fluidRow(column(12, style = "padding-left: 5px;", dataEditUI(ns("edit1"))))
                     
    ) #end of shinybsmodal
    )
  )
}
    
#' edit_data Server Functions
#'
#' @noRd 

mod_edit_data_server <- function(id, data_input){
  moduleServer(id, function(input, output, session){
    ns <- session$ns
 
    guided_tour = cicerone::Cicerone$
      new()$
      step(
        ns("idselect1"),
        "Select columns",
        "Here you can select which columns (variables) you need (geen = taken, red = not taken)."
      )$
      step(
        ns("idfilter1"),
        "Filter by rows",
        "Here you can filter rows by conditions. Use just one condition for each variable."
      )$
      step(
        ns("idoutput-active"),
        "Download edited data",
        "Here you can download the edited data as csv file."
      )$
      step(
        ns("idcut"),
        "Enable editing",
        "Button to enable/disable editing. ADViSELipidomics will use the edited file only if this button is on 'ON'."
      )$
      step(
        ns("done"),
        "Done",
        "Close and apply changes if there are."
      )$
      step(
        ns("cancel"),
        "Cancel",
        "Close without saving changes."
      )
    
    

     observeEvent(input$help_button, {
       guided_tour$init()$start()
     })
     
    
    # DATA STORAGE
    values <- reactiveValues(
      data = NULL, # original data
      data_active = NULL, # displayed data
      rows = NULL,
      columns = NULL
    )
    
    
    # RESET FILTERS
    observeEvent(data_input(), {
      # RESET FILTERS
      values$rows <- NULL
      values$columns <- NULL
      # BIND ROWS/COLUMNS
      values$data <- data_input()
      
    })
    
    # FILTERS ALWAYS RESET ON DATA SYNC
    
    # DATA SELECT
    data_select <- dataSelectServer("select1",
                                    data = reactive(values$data),
                                    hover_text = "select columns")
    
    # DATA FILTER
    data_filter <- dataFilterServer("filter1",
                                    data = reactive(values$data),
                                    hover_text = "filter rows")
    
    # UPDATE FILTERS
    observe({
      values$rows <- data_filter$rows()
      values$columns <- data_select$columns()
    })
    
    # DATA FILTERING
    observe({
      # ENTIRE DATA
      if (length(values$rows) == 0 & length(values$columns) == 0) {
        values$data_active <- values$data
        # DATA SUBSET
      } else {
        # ROWS
        if (length(values$rows) != 0 & length(values$columns) == 0) {
          values$data_active <- values$data[values$rows,
                                            ,
                                            drop = FALSE]
          # COLUMNS
        } else if (length(values$rows) == 0 & length(values$columns) != 0) {
          values$data_active <- values$data[,
                                            values$columns,
                                            drop = FALSE]
          # ROWS & COLUMNS
        } else if (length(values$rows) != 0 & length(values$columns) != 0) {
          values$data_active <- values$data[values$rows,
                                            values$columns,
                                            drop = FALSE]
        }
      }
    })
    
    # DATAEDIT - ENTIRE DATASET
    data_update <- dataEditServer("edit1",
                                  data = reactive({values$data_active}),
                                  col_names = FALSE)
    
    # UPDATE ACTIVE DATA
    observe({
      values$data_active <- data_update()
    })
    
    
    
    # DATA OUTPUT - DATA ACTIVE
    dataOutputServer("output-active",
                     data = reactive({values$data_active}),
                     hover_text = "save selection \n to file",
                     write_fun = openxlsx::write.xlsx,
                     save_as = paste0(Sys.Date(),"_edited.xlsx"))
    

    # DONE
    finaldata2 = eventReactive(input$done, {
      # DATA ACTIVE
      if (input$cut == TRUE){
        values$data_active
      }
    }, ignoreNULL = FALSE)
    
        # CANCEL
    observeEvent(input$cancel, {
      toggleModal(session = session, "upinternalmodal", toggle = "close")
    })
    
    # DONE
    observeEvent(input$done, {
      toggleModal(session = session, "upinternalmodal", toggle = "close")
    })
    
    addTooltip(session = session, id = "cut", title = "enable/disable editing")  
    
    return(reactive({
      if (is.null(finaldata2())) {
        finaldata = data_input()
      }else{
        finaldata = finaldata2()
      }
      return(finaldata)
    }))
    
    
  })
}
    
## To be copied in the UI
# mod_edit_data_ui("edit_data_ui_1")
    
## To be copied in the server
# mod_edit_data_server("edit_data_ui_1")
