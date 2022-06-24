#' calibration_step UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList
#' @importFrom bsplus bs_embed_tooltip
mod_calibration_step_ui <- function(id){
  ns <- NS(id)
  tagList(
    
    column(5, offset = 1,
           box(width = NULL, title = "2. CALIBRATION", status = "primary", solidHeader = TRUE,
               
               #Step 1.
               div(icon("circle"), HTML("&nbsp;Step. 1: Import calibration files (.xlsx)"), style = "text-align: left; font-size: 18px; font-weight: bold;"),
               br(),
               fluidRow(
                 column(10, fileInput(ns("calibdeuteratedpath"), "Select the Calibration Deuterated (.xlsx)", accept = ".xlsx")),
                 column(2, br(), actionButton(ns("editcalibdeu"), icon("edit"))),
                 shinyBS::bsModal(ns("upcalibdeumodal"), "Edit data", trigger = ns("editcalibdeu"), size = "large",
                                  DataEditR::dataEditUI(ns("editedcalibeu")),
                                  DataEditR::dataOutputUI(ns("output_calibeu_edit"))
                 )
               ),
               fluidRow(
                 column(10, fileInput(ns("calibnonlabelpath"), "Select the Calibration NonLabeled (.xlsx)", accept = ".xlsx")),
                 column(2, br(), actionButton(ns("editcaliblab"), icon("edit"))),
                 shinyBS::bsModal(ns("upcaliblabmodal"), "Edit data", trigger = ns("editcaliblab"), size = "large",
                                  DataEditR::dataEditUI(ns("editedcalilab")),
                                  DataEditR::dataOutputUI(ns("output_calilab_edit"))
                 )
               ),
               hr(),
               
               #Step 2.
               conditionalPanel(condition = "output.checkfilterdata == false", ns = ns,
                                div(icon("circle"), HTML("&nbsp;Step. 2: Select the folder with the concentration files"), style = "text-align: left; font-size: 18px; font-weight: bold;"),
                                br(),
                                fluidRow(
                                  div(
                                    column(5, shinyFiles::shinyDirButton(ns("calibrationfolder"), label = "Browse...", title = "Please select the concentration folder", multiple = FALSE, icon = icon("folder-open"))),
                                    column(7, actionButton(ns("gocalibrationfiles"), "Read the concentration files", icon("cogs"))),
                                    style = "text-align: center;")),
                                hr()
               ),
               
               #Step 3.
               conditionalPanel(condition = "output.checkcalibmatrix == false", ns = ns,
                                div(icon("circle"), HTML("&nbsp;Step. 3: Choose a folder where to save the outputs"), style = "text-align: left; font-size: 18px; font-weight: bold;"),
                                br(),
                                fluidRow(
                                  div(
                                    column(5, shinyFiles::shinyDirButton(ns("mainfolder"), label = "Browse...", title = "Please select the main folder", multiple = FALSE, icon = icon("folder-open"))),
                                    column(7, actionButton(ns("createfolders"), "Create the folder structure", icon("cogs"))),
                                    style = "text-align: center;")),
                                hr()
               ),
               
               #Step 4.
               conditionalPanel(condition = "input.createfolders > 0 && output.checkfolders == 1", ns = ns,
                                div(icon("circle"), HTML("&nbsp;Step. 4: Choose the calibration options"), style = "text-align: left; font-size: 18px; font-weight: bold;"),
                                br(),
                                fluidRow(
                                  column(6, offset = 1,
                                         fluidRow(checkboxInput(ns("intercept_flag"), "Intercept at zero", value = TRUE)),
                                         fluidRow(checkboxInput(ns("lm_robust_flag"), "Use robust method linear model", value = FALSE)),
                                         fluidRow(checkboxInput(ns("lin_calibration"), "Filter on concentration range of linearity", value = TRUE))
                                  ),
                                  column(4,
                                         div(
                                           actionButton(ns("gocalibplot"), "Run calibration", icon("cogs")),
                                           conditionalPanel(condition = "output.checkstepf == false", ns = ns,
                                                            br(), br(),
                                                            actionButton(ns("viewcalibplot"), "See the calibration plots", icon("eye")),
                                           ),
                                           style = "text-align: center;"),
                                         
                                         shinyBS::bsModal(ns("modalcalibnplots"), "Calibration plot", trigger = ns("viewcalibplot"), size = "large",
                                                          fluidRow(
                                                            column(4, DT::DTOutput(ns("dtlistcalibplot"))),
                                                            column(8, imageOutput(ns("phcalibplot")))
                                                          ))
                                  )
                                ),
                                hr()
               ),
               
               
               #Step 5.
               conditionalPanel(
                 condition = "output.checkstepf == false", ns = ns,
                 div(icon("circle"), HTML("&nbsp;Step. 5: Apply recovery"), style = "text-align: left; font-size: 18px; font-weight: bold;"),
                 fluidRow(
                   column(4, offset = 1,
                          br(),
                          div(
                            actionButton(ns("gorecovery"), "Apply recovery", class = "btn-primary btn-lg", icon("cogs")),
                            style = "text-align:center;")
                   ),
                   column(5, offset = 1,
                          conditionalPanel(condition = "output.checkstepg == false", ns = ns,
                                           actionButton(ns("checkrecovery"), "Check concentration matrix", class = "btn-primary btn-lg", icon("eye")),
                                           br(), br(),
                                           div(
                                             bsplus::bs_embed_tooltip(
                                               downloadButton(ns("table_not_lin"), "Download LOL"), 
                                               title = "Download the table of Lipids out Of Linearity"),
                                             style = "text-align:center;"),
                                           tags$head(tags$style("#viewrecoverymodal .modal-dialog{ width:1300px}")),
                                           tags$head(tags$style("#viewrecoverymodal .modal-body{ min-height:1000px}")),
                                           shinyBS::bsModal("viewrecoverymodal", trigger = ns("checkrecovery"),
                                                            radioGroupButtons(ns("vimordtrecovery"), "Choose what to show :",
                                                                              choiceValues = list("table", "vim"),
                                                                              choiceNames = list(
                                                                                paste(shiny::icon("table",  style='font-size:16px;'), HTML("<b style=font-size:16px>&nbsp;Table</b>")),
                                                                                paste(shiny::icon("buromobelexperte",  style='font-size:16px;'), HTML("<b style=font-size:16px>&nbsp;NA plots</b>"))),
                                                                              justified = TRUE, status = "primary"),
                                                            conditionalPanel(condition = "input.vimordtrecovery == 'table'", ns = ns,
                                                                             div(DT::DTOutput(ns("dtrecovery")), style = "overflow-x: scroll;")
                                                            ),
                                                            conditionalPanel(condition = "input.vimordtrecovery == 'vim'", ns = ns,
                                                                             plotOutput(ns("vimplotrecovery"), height = "1000px")
                                                            )
                                           ) #end of bsModal recovery
                          )
                   )
                 )
               )
           ) #end of box
    ) #end of column box calibration
    
 
  )
}
    
#' calibration_step Server Functions
#'
#' @noRd 
mod_calibration_step_server <- function(id, stepc){
  moduleServer( id, function(input, output, session){
    ns <- session$ns
 
    calibdeu_to_edit = reactive({
      req(input$calibdeuteratedpath)
      ext <- tools::file_ext(input$calibdeuteratedpath$name)
      if(ext != "xlsx"){
        shinyWidgets::show_alert("Invalid file!", "Please upload a .xlsx file", type = "error")
      }
      validate(need(ext == "xlsx", "Invalid file! Please upload a .xlsx file"))
      readxl::read_xlsx(input$calibdeuteratedpath$datapath, na = c("", "NA"))
    })
    
    calibdeu_edit <- DataEditR::dataEditServer("editedcalibeu", data = calibdeu_to_edit())
    
    DataEditR::dataOutputServer("output_calibeu_edit", data = calibdeu_edit)
    
    
    #internal standard
    caliblab_to_edit = reactive({
      req(input$calibnonlabelpath)
      ext <- tools::file_ext(input$calibnonlabelpath$name)
      if(ext != "xlsx"){
        shinyWidgets::show_alert("Invalid file!", "Please upload a .xlsx file", type = "error")
      }
      validate(need(ext == "xlsx", "Invalid file! Please upload a .xlsx file"))
      readxl::read_xlsx(input$calibnonlabelpath$datapath, na = c("", "NA"))
    })
    
    caliblab_edit <- DataEditR::dataEditServer("editedcalilab", data = caliblab_to_edit())
    
    DataEditR::dataOutputServer("output_calilab_edit", data = caliblab_edit)
    
    
    
    #carico e aggiorno la lista dei targets aggiungendo gli ultimi due file
    stepd = reactive({
      req(calibdeu_edit(), caliblab_edit())
      new = list(
        calibration_deuterated = calibdeu_edit(),
        calibration_nonlabeled = caliblab_edit()
      )
      stepcx = stepc()
      stepcx$targets = list(
        stepcx$targets,
        new
      ) %>% purrr::flatten()
      return(stepcx)
      
    })
    
    
    
    ####PART 2. Select the calibration curve folder
    
    output$checkfilterdata = reactive(
      return(is.null(stepd()))
    )
    outputOptions(output, "checkfilterdata", suspendWhenHidden = FALSE)
    
    
    volumes = c(Home = fs::path_home(), shinyFiles::getVolumes()())
    
    shinyFiles::shinyDirChoose(input, 'calibrationfolder', roots = volumes, session = session)
    
    calibration_path = reactive({
      if(length(input$calibrationfolder) != 1 ) {
        shinyFiles::parseDirPath(volumes,input$calibrationfolder)   
      }else{
        NULL
      }
    })
    
    
    
    calibration_matrix = eventReactive(input$gocalibrationfiles,{
      #calibration matrix deuterated
      cal_mat_de = caliblist_advise_lipidomics(
        info = "deuterated",
        out = stepd(),
        calibration_targetfile = stepd()$targets$calibration_deuterated,
        calibration_path = paste0(calibration_path(), "/")
      )
      #cal_mat_de <- base::Reduce(f = dplyr::full_join, x = Filter(Negate(is.null), cal_mat_de))
      
      #calibration matrix nonlabeled
      cal_mat_nl <- caliblist_advise_lipidomics(
        info = "non labeled",
        out = stepd(),
        calibration_targetfile = stepd()$targets$calibration_nonlabeled,
        calibration_path = paste0(calibration_path(), "/")
      )
      #unisco le due matrix
      #cal_mat_nl <- base::Reduce(f = dplyr::full_join, x = Filter(Negate(is.null), cal_mat_nl))
      cal_mat <- dplyr::bind_rows(cal_mat_nl, cal_mat_de)
      cal_mat = data.frame(cal_mat[,order(as.numeric(colnames(cal_mat)))],
                           row.names = cal_mat$LipidIon, check.names = FALSE)
      return(cal_mat)
    })
    
    
    
    
    
    
    
    ##### STEP 3. Select a folder where to save the outputs
    
    #se tutti i data sono stati caricati
    output$checkcalibmatrix = reactive(
      return(is.null(calibration_matrix()))
    )
    outputOptions(output, "checkcalibmatrix", suspendWhenHidden = FALSE)
    
    
    shinyFiles::shinyDirChoose(input, 'mainfolder', roots = volumes, session = session)
    
    
    input_path = reactive({
      if(length(input$mainfolder) != 1 ) {
        shinyFiles::parseDirPath(volumes,input$mainfolder)   #"input_path" same of config_advise_lipidomics()
      }else{
        NULL
      }
    })
    
    
    #here we have list with all the paths (both input and output fo)
    folders = reactive({
      if(!is.null(input_path())){
        list(
          input_path(),
          output_path = paste0(input_path(), "/Output_Lipidomics/"),
          target_path = paste0(input_path(), "/TargetFiles_Lipidomics"),
          calibration_path = paste0(input_path(), "/CalibrationCurveRuns/"),
          output_path_plots = paste0(input_path(), "/Output_Lipidomics/Plots/"),
          output_path_tables = paste0(input_path(), "/Output_Lipidomics/Tables/"),
          output_path_analysis = paste0(input_path(), "/Output_Lipidomics/Analysis/")
        )
      }
    })
    checkfolders1 = reactiveVal(0)
    observeEvent(input$createfolders, {
      if(!is.null(folders())){
        showNotification(tagList(icon("cogs"), HTML("&nbsp;Creating Output Data Structure...")), type = "default")
        
        if (!dir.exists(folders()$output_path)){
          dir.create(file.path(folders()$output_path))
        } else {
          showNotification(tagList(icon("info"), HTML("&nbsp;Note: the output folder already exists")), type = "default")
        }
        showNotification(tagList(icon("info"), HTML("&nbsp;Results will be saved in:", folders()$output_path)), type = "default")
        
        
        if (!dir.exists(folders()$output_path_plots)){
          dir.create(file.path(folders()$output_path_plots))
        }
        
        if (!dir.exists(folders()$output_path_tables)){
          dir.create(file.path(folders()$output_path_tables))
        }
        
        if (!dir.exists(folders()$output_path_analysis)){
          dir.create(file.path(folders()$output_path_analysis))
        }
        
        showNotification(tagList(icon("check"), HTML("&nbsp;Data Structure successfully created!")), type = "message")
        checkfolders1(1)
        output$checkfolders = reactive(checkfolders1())
        outputOptions(output, "checkfolders", suspendWhenHidden = FALSE)
        
        
      }else{
        showNotification(tagList(icon("exclamation"), HTML("&nbsp;Please, first select a folder from the Browse button")), type = "warning")
        checkfolders1(0)
      }
      
    })
    
    
    ###### STEP 4. Choose the calibration options
    
    #se tutti i data sono stati caricati
    output$checkfolder = reactive(
      return(is.null(folders()))
    )
    outputOptions(output, "checkfolder", suspendWhenHidden = FALSE)
    
    
    stepf = eventReactive(input$gocalibplot,{
      req(calibration_matrix())
      
      ###unisco prima folders nella lista
      stepe = list(stepd(), folders = list(folders())) %>% purrr::flatten()
      
      calibplot_advise_lipidomics(
        out = stepe,
        cal_plot_path = "/CalibrationPlot",
        cal_mat = calibration_matrix(),
        intercept_flag = !input$intercept_flag,
        lm_robust_flag = input$lm_robust_flag,
        lin_calibration = input$lin_calibration,
        plot_calibration = TRUE
      )
    })
    
    #stepf() prende il path stepd()$folders$output_path_plots, crea la cartella CalibrationPlot e crea i plot. Inoltre
    #dalla calibration_matrix() restituisce i seguenti valori:
    #stepf()$calibration_matrix
    #stepf()$coefficients
    #ex calibplot_data()
    
    
    
    list_calib_plot = reactive({
      req(stepf())
      x =list.files(paste0(stepf()$folders$output_path_plots, "/CalibrationPlot"))
      x2 = stringr::str_sub(x, end=-5)
      x2 = data.frame("Plot" = x2) #, Method = "NULL"
      #x3 = x2 %>% dplyr::mutate(Method = ifelse(str_detect(x2$Plot, "_lm_") == TRUE, "Linear model", "Robust linear model"))
      x2
    })
    
    output$dtlistcalibplot = renderDT({
      req(list_calib_plot())
      list_calib_plot()
    }, selection = "single", server = FALSE, rownames = FALSE, option = list(pageLength = 15))
    
    
    
    output$phcalibplot = renderImage({
      req(input$dtlistcalibplot_rows_selected)
      nroww=input$dtlistcalibplot_rows_selected
      #x = dplyr::select(list_calib_plot(), "Plot")
      x = list_calib_plot()
      path = paste0(stepf()$folders$output_path_plots, "CalibrationPlot/", x[nroww,], ".png")
      
      list(src = path, width = 500)
    }, deleteFile = FALSE)
    
    
    
    ####### STEP 5. Make the recovery percentage
    
    
    #se tutti i data sono stati caricati
    output$checkstepf = reactive(
      return(is.null(stepf()))
    )
    outputOptions(output, "checkstepf", suspendWhenHidden = FALSE)
    
    
    stepg = eventReactive(input$gorecovery,{
      
      #qui mi serve out$calibration_matrix,out$coefficients    calibplot_data()$
      #poi out$targets$internal_standard                       targets()$internal_standard
      #out$lipid_data                                          readed_data()$lipid_data
      #out$lipid_filtered                                      filtered_data()$lipid_filtered

      recovery_advise_lipidomics(
        out = stepf(),
        intercept_flag = !input$intercept_flag
      )
    })
    
    #stepg ha in più
    #stepg()$concentration_matrix
    #stepg()$recovery_percentage
    
    
    #### Download handler for the download button
    output$table_not_lin <- downloadHandler(
      #put the file name with also the file extension
      filename = function() {
        paste0("Table_LOL_", Sys.Date(), ".xlsx")
      },
      
      # This function should write data to a file given to it by the argument 'file'.
      content = function(file) {
        openxlsx::write.xlsx(stepg()$table_not_lin, file)
      }
    )
    
    
    #se tutti i data sono stati caricati
    output$checkstepg = reactive(
      return(is.null(stepg()))
    )
    outputOptions(output, "checkstepg", suspendWhenHidden = FALSE)

    
    
    output$dtrecovery = renderDT({
      req(stepg())
      stepg()$concentration_matrix
    })
    
    output$vimplotrecovery = renderPlot({
      req(stepg())
      stepg()$concentration_matrix %>% dplyr::select(-1) %>% VIM::aggr(combined = T, cex.axis = 0.7)
    })
    
    
    return(reactive({
      req(stepg())
      return(stepg())
    }))
    
  })
}
    
## To be copied in the UI
# mod_calibration_step_ui("calibration_step_ui_1")
    
## To be copied in the server
# mod_calibration_step_server("calibration_step_ui_1")
