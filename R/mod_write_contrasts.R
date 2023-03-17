#' write_contrasts UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_write_contrasts_ui <- function(id){
  ns <- NS(id)
  tagList(


    conditionalPanel(ns = ns,
      "output.check_if_batch_rem_out == false",
      h4(strong("Design variables")),
      selectInput(ns("expdes_design_vars_all"), "Select primary variable", choices = "", multiple = TRUE),
      hr(),

      h4(strong("Batch effect")),
      awesomeCheckbox(ns("expdes_batch_effect"), "Batch effect", value = TRUE),
      conditionalPanel(
        condition = "input.expdes_batch_effect == true", ns= ns,
        fluidRow(
          column(
            6,
            selectInput(
              ns("expdes_batch_type"), choices = c("Explicity removal" = "remove","Fit batch vars with design" = "fit"),
              label = tags$span("Batch type",
                tags$i(style = "color:#0072B2;",class = "glyphicon glyphicon-info-sign",
                       title = "ADViSELipidomics copes with the batch effects by either fitting the model with the batch variables
                    or removing the batch effect before fitting the model. If you select 'fit', the software requires your batch variable
                    also in the contrasts list.")))
            ),

          column(
            6,
            selectInput(ns("batch_meth"),choices = "",
                        label = tags$span("Batch method",
                                          tags$i(style = "color:#0072B2;",class = "glyphicon glyphicon-info-sign",
                                                 title = "Algorithm used for the batch effect. If 'limma',it will be used the 'removeBatchEffect' function, 
                                                 otherwhise the 'ComBat' function from SVA package (parametric or non-parametric).")))
            )
        ),
        conditionalPanel( ns= ns,
          condition = "input.expdes_batch_type == 'remove' || input.batch_meth == 'given'",
          fluidRow(
            column(8, selectInput(ns("expdes_batch_var1"), "Primary batch variable", choices = "")),
            column(4, br(), bsButton(ns("add2batch"), label = HTML("&nbsp;Add"), style="success", icon("plus")))
          ),
          fluidRow(
            conditionalPanel(condition = "output.checkadd2batch == 'twobatch'", ns= ns,
                             column(8, selectInput(ns("expdes_batch_var2"), "Second batch variable", choices = ""))
            )
          )
        )
      )
    ),

    conditionalPanel(
      "output.check_if_batch_rem_out == true", ns= ns,
      helpText("Batch already removed. Design variables and batch variables taken from the 'remove batch' step.")
    ),
    hr(),

    h4(strong("Contrasts")),
    awesomeRadio(ns("input_of_contrast"), "How do you want to load the contrasts?", choices = c("Write contrasts", "Upload a txt file"),inline =TRUE),
    conditionalPanel(ns= ns,
      "input.input_of_contrast == 'Write contrasts'",
      column(12, actionButton(ns("open_writecontrast"), "Write contrasts", icon("pen")), style = "text-align: center;"), br(), br(),
    ),
    conditionalPanel(ns= ns,
      "input.input_of_contrast == 'Upload a txt file'",
      column(12, fileInput(ns("contrast_txt_file"),"Load a contrast file", accept = ".txt"), style = "text-align: center;")

    ),

    hr(),

    tags$head(tags$style(paste0("#",ns("view_writecontrast")," .modal-dialog{ width:1100px}"))),
    shinyBS::bsModal(ns("view_writecontrast"), trigger = ns("open_writecontrast"), size = "large", title = "Write the contrasts",
                     fluidRow(
                       column(
                         3,
                         box(width = NULL, status = "primary", title = "First term", solidHeader = TRUE,
                             selectInput(ns("firstelement"), "Level 1", choices = ""),
                             conditionalPanel(condition = "output.checkntotvars == 1", ns= ns,
                                              selectInput(ns("secondelement1"), "Level 2", choices = "")
                             ),
                             conditionalPanel(condition = "output.checkntotvars == 2", ns= ns,
                                              selectInput(ns("secondelement2"), "Level 2", choices = "")
                             )
                         ),
                         conditionalPanel(
                           condition = "output.checkntotvars == 2", ns= ns,
                           box(width = NULL, status = "primary", title = "Second term", solidHeader = TRUE,
                               radioButtons(ns("fixedlevel"), "Fixed Level", choices = c("Level 1", "Level 2"), inline = TRUE), #, checkbox = TRUE
                               conditionalPanel(condition = "input.fixedlevel == 'Level 1'", ns= ns,
                                                br(),
                                                h5(strong("Level 1: fixed")),
                                                br(),
                                                selectInput(ns("thirdelement2"), "Level 2", choices = "")
                               ),
                               conditionalPanel(condition = "input.fixedlevel == 'Level 2'", ns= ns,
                                                selectInput(ns("thirdelement1"), "Level 1", choices = ""),
                                                h5(strong("Level 2: fixed"))
                               )
                           )

                         ),
                         hr(),
                         # Row selection
                         numericInput(ns("row.selection"), "Select row to be deleted", min = 1, max = 100, value = ""),
                         # Add button
                         actionButton(ns("add.button"), "Add", icon("plus")),
                         # Delete button
                         actionButton(ns("delete.button"), "Delete", icon("minus"))
                       ),
                       column(9,
                              DTOutput(ns("tablewritten")),
                              br(),
                              hr(),
                              br(),
                              fluidRow(
                                column(9,verbatimTextOutput(ns("printconlist"))),
                                column(3, actionButton(ns("submit_written"), "Submit", icon("file-import"), style = "width: 150px;height: 50px;font-size: 20px;"))
                              )
                       )
                     )

    ),
  )
}

#' write_contrasts Server Functions
#'
#' @noRd
mod_write_contrasts_server <- function(id, sumexp, is_batch_removed, Batch_Options){
  moduleServer( id, function(input, output, session){
    ns <- session$ns

    #TRUE if batch already removed
    output$check_if_batch_rem_out = reactive({
      is_batch_removed()
    })
    outputOptions(output,"check_if_batch_rem_out", suspendWhenHidden = FALSE)



    #update batch method
    observeEvent(input$expdes_batch_type, {
      if(input$expdes_batch_type == "fit"){
        updateSelectInput(session, "batch_meth", choices =  c("Known batch" = "given", "Unknown batch" = "estimate"))
      }else if(input$expdes_batch_type == "remove"){
        updateSelectInput(session, "batch_meth", choices = c("limma", "combat_param", "combat_nonparam"))
      }
    })


    observeEvent(sumexp(), {

      #var replicate effect
      updateSelectInput(session, "repeff_vars", choices = colnames(SummarizedExperiment::colData(sumexp())))

      #var designs
      datavar1 = sumexp() %>% SummarizedExperiment::colData() %>%
        as.data.frame() %>% dplyr::select(!where(is.double)) %>% dplyr::select(where(function(x) length(unique(x))>1))

      datavar2 = sumexp() %>% SummarizedExperiment::colData() %>%
        as.data.frame() %>% dplyr::select(where(function(x) length(unique(x))>1))

      updateSelectInput(session, "expdes_design_vars_all", choices = colnames(datavar2))

      #batch
      updateSelectInput(session, "expdes_batch_var1", choices = colnames(SummarizedExperiment::colData(sumexp())))
      updateSelectInput(session, "expdes_batch_var2", choices = colnames(SummarizedExperiment::colData(sumexp())))
    })


    #update design and batch vars if batch already removed
    observeEvent(is_batch_removed(),{
      if(is_batch_removed()){

        #update design variables
        updateSelectInput(session, "expdes_design_vars_all", choices = Batch_Options()$design_vars, selected = Batch_Options()$design_vars)


        #batch
        if(length(Batch_Options()$batch_vars) == 2){
          updateSelectInput(session, "expdes_batch_var1", choices = Batch_Options()$design_vars[1])
          updateSelectInput(session, "expdes_batch_var2", choices = Batch_Options()$design_vars[2])
        }else{
          updateSelectInput(session, "expdes_batch_var1", choices = Batch_Options()$design_vars[1])
          updateSelectInput(session, "expdes_batch_var2", choices = "")
        }
      }
    })



    #Batch variable
    output$checkadd2batch = reactive({
      if(input$add2batch %%2 == 0){"onebatch"} else{"twobatch"}
    })
    outputOptions(output, "checkadd2batch", suspendWhenHidden = FALSE)

    observeEvent(input$add2batch,{
      if(input$add2batch %%2 == 1){
        updateButton(session, "add2batch",label = HTML("&nbsp;Remove"), style = "danger", icon("minus"))
      }else{
        updateButton(session, "add2batch", label = HTML("&nbsp;Add"), style="success", icon("plus"))
      }
    })



    varsde = reactive({
      if(is.null(is_batch_removed())){return(NULL)}
      if(is_batch_removed()) return(Batch_Options()$design_vars)

      validate(need(input$expdes_design_vars_all, "Select at least one design variable."))
      return(input$expdes_design_vars_all)
    })



    #all batches here
    batches = reactive({
      if(is.null(is_batch_removed())){return(NULL)}

      if(is_batch_removed()){
        Batch_Options()$batch_vars
      }else{
        if(input$expdes_batch_effect == TRUE){
          if(input$add2batch %%2 == 1){
            batchs = c(input$expdes_batch_var, input$expdes_batch_var2)
          }else{ batchs = input$expdes_batch_var1 }
        }else{ batchs = "" }
        batchs
      }
    })


    ntotvars = reactive({
      req(varsde())
      if(input$batch_meth == "given"){
        totvar = c(varsde(),batches())
      }else{
        totvar = c(varsde())
      }
      length(totvar[totvar != ""])
    })


    observeEvent(ntotvars(),{
      if(ntotvars() > 2){
        updateAwesomeRadio(session, "input_of_contrast", choices = c("Upload a txt file"),inline = TRUE)
      }else{
        updateAwesomeRadio(session, "input_of_contrast", choices = c("Write contrasts", "Upload a txt file"),inline = TRUE)
      }
    })


    output$checkntotvars = reactive({
      ntotvars()
    })
    outputOptions(output, "checkntotvars", suspendWhenHidden = FALSE)


    #### write table ##


    firstvar = eventReactive(input$expdes_design_vars_all,{
      req(sumexp()) #firstvariab
      validate(need(input$expdes_design_vars_all, "Select at least one design variable"))
      sumexp() %>% SummarizedExperiment::colData() %>% as.data.frame() %>%
        dplyr::select(input$expdes_design_vars_all[1]) %>% unique()
    })



    observeEvent(firstvar(),{
      updateSelectInput(session, "firstelement", label = paste0("Level 1 for ", input$expdes_design_vars_all[1]), choices = firstvar())
      updateSelectInput(session, "secondelement1", label = paste0("Level 2 for ", input$expdes_design_vars_all[1]), choices = firstvar())
    })

    #second var
    secondvar = reactive({
      req(sumexp(), ntotvars()) #input$secondvariab
      if(ntotvars() == 2){
        if(input$batch_meth == "given"){
          sumexp() %>% SummarizedExperiment::colData() %>% as.data.frame() %>%
            dplyr::select(input$expdes_batch_var1) %>% unique()
        }else{
          sumexp() %>% SummarizedExperiment::colData() %>% as.data.frame() %>%
            dplyr::select(input$expdes_design_vars_all[2]) %>% unique()
        }
      }else{
        NULL
      }
    })

    observeEvent(secondvar(),{
      updateSelectInput(session, "secondelement2", label = paste("Level 2 for", colnames(secondvar())), choices = secondvar())
    })


    #third var
    thirdvar = reactive({
      req(sumexp())
      if(input$fixedlevel == "Level 2"){firstvar()}else{secondvar()}
    })


    observeEvent(thirdvar(),{
      updateSelectInput(session, "thirdelement1", label = paste("Level 1 for", colnames(thirdvar())), choices = thirdvar())
      updateSelectInput(session, "thirdelement2", label = paste("Level 2 for",colnames(thirdvar())), choices = thirdvar())
    })

    observeEvent(secondvar(),{
      req(firstvar())
      updateRadioButtons(session, "fixedlevel", choiceNames  = c(input$expdes_design_vars_all[1], colnames(thirdvar())), choiceValues = c("Level 1", "Level 2"))
    })

    #thirdvariab
    values <- reactiveValues()
    values$df <- na.omit(data.frame(First_level = NA,
                                    Second_level = NA, Third_level = NA, Fourth_level = NA))


    observeEvent(ntotvars(),{
      if(ntotvars() == 1){
        values$df <- na.omit(data.frame(First_level = NA, Second_level = NA))
      }else{
        values$df <- na.omit(data.frame(First_level = NA, Second_level = NA,
                                        Third_level = NA, Fourth_level = NA))
      }
    })

    observeEvent(input$add.button,{
      cat("addEntry\n")
      print(input$firstelement)
      print(input$secondelement)
      if(ntotvars() == 2){
        print(input$thirdelement2)
        if(input$fixedlevel == "Level 2"){
          newRow <- data.frame(input$firstelement, input$secondelement2, input$thirdelement1, input$secondelement2)
        }else{
          newRow <- data.frame(input$firstelement, input$secondelement2, input$firstelement, input$thirdelement2)
        }
      }else{
        newRow <- data.frame(input$firstelement, input$secondelement1)
      }
      colnames(newRow)<-colnames(values$df)
      values$df <- rbind(values$df,newRow)
    })

    observeEvent(input$delete.button,{
      cat("deleteEntry\n")
      if(is.na(input$row.selection)){
        values$df <- values$df[-nrow(values$df), ]
      } else {
        values$df <- values$df[-input$row.selection, ]
      }
    })

    output$tablewritten1 = renderDT({
      values$df[,c(1,2)]
    },options = list(dom = 't'))
    output$tablewritten2 = renderDT({
      values$df[,c(3,4)]
    },options = list(dom = 't'))

    output$tablewritten = renderDT({
      values$df
    },options = list(dom = 't'))


    writtenlist = reactive({
      validate(need(nrow(values$df) > 0, "Dataframe empty. Waiting for new data..."))
      if(ntotvars() == 2){
        d1= values$df[,c(1,2)] %>% tidyr::unite(col = "unif_1", sep = "_")
        d2= values$df[,c(3,4)] %>% tidyr::unite(col = "unif_2", sep = "_")
        cbind(d1,d2) %>% tidyr::unite(col = "unif_3", sep = "vs") %>% dplyr::mutate(unif_4 = unif_3 ) %>%
          dplyr::mutate(unif_4 = stringr::str_replace(unif_4,"vs","-")) %>%
          tidyr::unite(col = "contrast", sep = "=")
      }else{
        values$df %>% tidyr::unite(col = "unif_1", sep = "vs") %>% dplyr::mutate(unif_2 = unif_1 ) %>%
          dplyr::mutate(unif_2 = stringr::str_replace(unif_2,"vs","-")) %>%
          tidyr::unite(col = "contrast", sep = "=")
      }
    })

    output$printconlist = renderPrint({
      req(writtenlist())
      writtenlist()
    })



    written_contrast = eventReactive(input$submit_written,{
      req(writtenlist())

      if(ntotvars() == 2){
        tocheck = stringr::str_split_fixed(writtenlist()$contrast, pattern = "=", n = 2)[,1] %>%
          stringr::str_split_fixed(pattern = "vs", n = 2) %>% as.vector()

        if(input$batch_meth == "given"){
          totvar = c(varsde(),batches())
        }else{
          totvar = c(varsde())}
        totvar = totvar[totvar != ""]

        united = sumexp() %>% SummarizedExperiment::colData() %>%
          as.data.frame() %>% tidyr::unite("merged", totvar, sep = "_")

        check = tocheck %in% united$merged  #check if are all present
        if(FALSE %in% check){
          no_comb = tocheck[!check]
          shinyWidgets::show_alert(title = "One or more combinations don't exist!", type = "error",
                                   text = paste0("The combination(s) '", no_comb,
                                                 "' doesn't exist in your dataset. Please change the levels and submit again."))

        }
        validate(need(!(FALSE %in% check), "A combination is not present"))
      }
      x = writtenlist()
      colnames(x) <- "X1"
      x%>% tibble::as_tibble()
      return(x)
    })

    observeEvent(input$submit_written,{
      showNotification(tagList(icon("check"), HTML("&nbsp;Written contrasts submitted!")), type = "message")
    })

    batch_type = reactive({
      if(is.null(is_batch_removed())){return(NULL)}
      if(is_batch_removed()) return("already")
      return(ifelse(input$expdes_batch_effect == FALSE, "none", input$expdes_batch_type))
    })



    ### load txt file contrast
    file_contrast = reactive({
      req(input$contrast_txt_file)
      ext <- tools::file_ext(input$contrast_txt_file$name)
      if(ext != "txt"){
        shinyWidgets::show_alert("Invalid file!", "Please upload a .txt file", type = "error")
        return(NULL)
      }

      readr::read_csv(input$contrast_txt_file$datapath, col_names = FALSE,show_col_types = FALSE) %>% as.data.frame()
    })

    observeEvent(file_contrast(),{
      if(!is.null(file_contrast())){
        sendmessages("Contrasts file loaded!", type = "success")
      }
    })


    final_contrast = reactive({
      if(input$input_of_contrast == "Write contrasts"){
        req(written_contrast())
        data = file_cont = written_contrast()
      }else{
        req(file_contrast())
        data = file_cont = file_contrast()
      }

      list(file_cont = data,
           varsde = varsde(),
           batches = batches(),
           batch_meth = input$batch_meth,
           batch_type = batch_type()
           )
    })

    return(reactive({
      req(final_contrast())
      return(final_contrast())
    }))


  })
}

## To be copied in the UI
# mod_write_contrasts_ui("write_contrasts_1")

## To be copied in the server
# mod_write_contrasts_server("write_contrasts_1")
