#' The application server-side
#' 
#' @param input,output,session Internal parameters for {shiny}. 
#'     DO NOT REMOVE.
#' @import shiny
#' @rawNamespace import(ggplot2, except = last_plot)
#' @rawNamespace import(plotly, except = rename)
#' @importFrom stats prcomp kmeans dist hclust
#' @importFrom shinyFiles shinyDirChoose getVolumes parseDirPath
#' @importFrom fs path_home
#' @importFrom DT renderDT datatable JS
#' @importFrom stringr str_sub str_split_fixed str_replace
#' @importFrom purrr flatten pluck
#' @import shinyBS
#' @importFrom DataEditR dataEditServer dataOutputServer
#' @importFrom lubridate as_date
#' @importFrom VIM aggr
#' @importFrom shinyWidgets sendSweetAlert show_alert
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @import ggfortify
#' @importFrom tidyr gather separate pivot_longer unite
#' @importFrom forcats fct_rev
#' @import utils
#' @importFrom tibble rownames_to_column as_tibble column_to_rownames
#' @import ggVennDiagram
#' @importFrom stringi stri_list2matrix
#' @importFrom ComplexHeatmap make_comb_mat comb_size comb_degree UpSet upset_right_annotation draw column_order decorate_annotation anno_barplot
#' @importFrom grid gpar unit grid.text
#' @importFrom mixOmics plsda background.predict plotIndiv plotVar perf splsda tune.splsda color.jet color.mixo
#' @importFrom gridExtra grid.arrange
#' @importFrom cluster pam clara
#' @importFrom factoextra fviz_nbclust fviz_cluster eclust fviz_dend fviz_silhouette
#' @importFrom metabolomicsWorkbenchR do_query
#' @importFrom openxlsx write.xlsx
#' @importFrom graphics par
#' @import scales
#' @noRd




app_server <- function( input, output, session ) {
  # Your application server logic 
  
  options(scipen = 999)
  
  observeEvent(input$jumptohome, {
    shinydashboard::updateTabItems(session, "sidebarmenu", "home")
  })
  
  observeEvent(input$gotoimport, {
    shinydashboard::updateTabItems(session, "sidebarmenu", selected = "rawsub")
  })
  
  
  # the modal dialog where the user can enter the query details.
  query_modal <- modalDialog(
    h3(strong("Welcome to ADViSELipidomics"), style = "text-align: center"),
    br(),
    h4("Before start, please enter your name and your company."),
    textInput("indata_analyst", list(HTML("&nbsp;"), icon("user"), HTML("&nbsp;Enter your name")), value = "Name"),
    textInput("inlab_analyst", list(HTML("&nbsp;"), icon("building"), HTML("&nbsp;Enter your company")), value = "Company"),
    easyClose = F,
    footer = tagList(
      actionButton("run", "Run")
    )
  )

  
  
  # Show the model on start up ...
  observeEvent(input$sidebarmenu,{
    req(input$sidebarmenu == "rawsub")
    showModal(query_modal)
  },ignoreInit = TRUE, once = TRUE)
  

  # ... or when user wants to change query
  observeEvent(input$change,{
    showModal(query_modal)
  })

  observeEvent(input$run, {
    removeModal()
  })
  
  output$nome = renderText({
    if(is.null(input$indata_analyst)){
      "Name"
    }else{
      input$indata_analyst
    }
  })


  analysis = reactive({
    list(
      lab_analyst = input$inlab_analyst,
      data_analyst = input$indata_analyst,
      pipe_ver = utils::packageVersion("ADViSELipidomics"),
      analysis_date = date()
    )
  })

  
  
  # this makes the directory at the base of your computer.
  volumes = c(Home = fs::path_home(), shinyFiles::getVolumes()())
  
  
  ######## STEP FILTERING ###########
  
  stepc = mod_reading_step_server("reading_step_lipsearch", analysis = analysis, int_std = reactive(input$type_lipsearch))
  

  
  ######### STEP CALIBRATION ##########
  
  ###check end of filtering box
  output$checkendfiltering = reactive(
    return(is.null(stepc()))
  )
  outputOptions(output, "checkendfiltering", suspendWhenHidden = FALSE)
  
  
  ####PART 1. import calibration files
  
  stepgwith = mod_calibration_step_server("calibration_step_lipsearch_with", stepc = stepc)
  
  stepg = reactive({
    req(stepc())
    if(input$type_lipsearch == 'No'){
      data_list <- lapply(stepc()$lipid_filtered, function(x) x %>%
                            dplyr::select(LipidIon,Area))
      
      data_updated = list()
      for(k in 1:length(data_list)){
        name = gsub("_nonlabeled","",names(data_list[k]))
        data_updated[k] = list(data_list[[k]] %>% dplyr::rename(!!name := "Area"))
      }
      
      merged = data_updated %>% purrr::reduce(dplyr::full_join, by = "LipidIon") %>% 
        dplyr::rename("Lipids" = "LipidIon")
      
      #return stepg$concentration_matrix with everything else
      list(stepc(), concentration_matrix = list(concentration_matrix = merged)) %>% purrr::flatten()
    }else{
      stepgwith()
    }

  })
  

  
  #check data correctly loaded
  output$checkstepg = reactive(
    return(is.null(stepg()))
  )
  outputOptions(output, "checkstepg", suspendWhenHidden = FALSE)
  
  observeEvent(stepg(),{
    num_na = sum(is.na(stepg()$concentration_matrix))
    not_na = sum(!is.na(stepg()$concentration_matrix))
    
    if (num_na/(num_na+not_na) > 0.5){
      shinyWidgets::sendSweetAlert(session, title = "Too many NAs", type = "warning", width = "600px",
                                   text = tags$div(
                                     p(tags$div("More than 50% of the calibration matrix data are NAs.", style = "font-weight: bold;")),
                                     p(div("You can:",style = "text-align: left;")),
                                     p(tags$div(icon("caret-right"), " Press the ", strong("'Check concentration matrix'"), " button to visualize them.", style = "text-align: left;")),
                                     p(tags$div(icon("caret-right"), " Filter and imputate them in the ", strong("next step."), style = "text-align: left;")),
                                   )
      )
    }
    
  })


  ######### STEP NA filtering and imputation ####


  stephwith = mod_imputation_step_server("imputation_step_lipsearch_with", parent = session, stepg = stepg, data_type = "Concentration")
  
  stephwithout = mod_imputation_step_server("imputation_step_lipsearch_without", parent = session, stepg = stepg, data_type = "Area")
  
  steph = reactive({
    if(input$type_lipsearch == 'No'){
      stephwithout()
    }else{
      stephwith()
    }
  })



  ######## Create sumexp object from LIQUID ########
  
  ##### STEP 1.
  ####import target and internal files xlsx
  
  #targetfile
  targetfile_to_edit_liquid = reactive({
    req(input$targetfile_liquid)
    ext <- tools::file_ext(input$targetfile_liquid$name)
    if(ext != "xlsx"){
      shinyWidgets::show_alert("Invalid file!", "Please upload a .xlsx file", type = "error")
    }
    validate(need(ext == "xlsx", "Invalid file! Please upload a .xlsx file"))
    x = readxl::read_xlsx(input$targetfile_liquid$datapath)
    #x$Exp_date = lubridate::as_date(x$Exp_date)
    return(x)
  })
  
  
  targetfile_edit_liquid = mod_edit_data_server("edit_target_liquid", data_input = targetfile_to_edit_liquid)
  
  
  #internal standard
  internalstd_to_edit_liquid = reactive({
    req(input$internalstdpath_liquid)
    ext <- tools::file_ext(input$internalstdpath_liquid$name)
    if(ext != "xlsx"){
      shinyWidgets::show_alert("Invalid file!", "Please upload a .xlsx file", type = "error")
    }
    validate(need(ext == "xlsx", "Invalid file! Please upload a .xlsx file"))
    readxl::read_xlsx(input$internalstdpath_liquid$datapath)
  })
  
  internalstd_edit_liquid = mod_edit_data_server("edit_internal_liquid", data_input = internalstd_to_edit_liquid)
  
  
  
  targets_liquid = reactive({
    req(targetfile_edit_liquid(), internalstd_edit_liquid())
    list(
      targetfile_lipidomics = targetfile_edit_liquid(),
      internal_standard = internalstd_edit_liquid()
    )
  })
  
  
  ###### STEP 2.
  #check data correctly loaded
  output$checktargets_liquid = reactive(
    return(is.null(targets_liquid()))
  )
  outputOptions(output, "checktargets_liquid", suspendWhenHidden = FALSE)
  
  
  ####select the data folder where to read data
  
  # this makes the directory at the base of your computer.
  volumes = c(Home = fs::path_home(), shinyFiles::getVolumes()())
  
  shinyFiles::shinyDirChoose(input, 'datafolder_liquid', roots = volumes, session = session)
  
  data_path_liquid = reactive({
    if(length(input$datafolder_liquid) != 1 ) {
      shinyFiles::parseDirPath(volumes,input$datafolder_liquid)
    }else{
      NULL
    }
  })
  
  
  stepb_liquid = eventReactive(input$readdatabttn_liquid,{
    req(data_path_liquid(), targets_liquid(), analysis())
    
    stepa = list(targets = targets_liquid(), analysis = analysis())
    withProgress(message = "Reading data...", value=0, {
      
      read_advise_lipidomics(
        out = stepa,
        datapath = data_path_liquid(),
        target_file = stepa$targets$targetfile_lipidomics,
        data_type = "Liquid")
    })
    
  })
  
  #stepb take the lipid_data path (ex. AF-1C-M_deuterated_1.txt) and targetfile and load all files.
  #stepb_liquid()$lipid_data
  #stepb_liquid()$replicates  TRUE if there are replicates otherwise FALSE
  
  
  #check data correctly loaded
  output$checkdatafiles_liquid = reactive(
    return(is.null(stepb_liquid()))
  )
  outputOptions(output, "checkdatafiles_liquid", suspendWhenHidden = FALSE)
  
  
  ####quality check
  output$qualcheckplot_liquid = renderPlotly({
    req(stepb_liquid())
    
    lipid_nonfilt = stepb_liquid()$lipid_data
    area = lapply(lipid_nonfilt, function(x) log2(sum(x$Intensity))) %>% tibble::as_tibble() %>% t() %>% as.data.frame() %>% tibble::rownames_to_column()
    colnames(area) = c("Samples", "Intensity")
    
    hh = ggplot(area) + geom_col(aes(x= Samples, y = Intensity)) + ylab("Log2(Intensity)") + 
      theme(axis.text.x = element_text(angle = 270, hjust = 0, size = 8))
    ggplotly(hh)
  })
  
  
  #######filtering options
  

  
  stepc_liquid = eventReactive(input$gofilterbttn_liquid,{
    req(stepb_liquid())
    filter_advise_lipidomics(
      out = stepb_liquid(),
      ca_bound = input$ca_bound_liquid,
      db_bound = input$db_bound_liquid,
      data_type = "Liquid"
    )
  })
  
  #stepc_liquid() take internal standard in targets() and lipid_data in readed_data() and filters them.
  #stepc_liquid()$lipid_filtered
  #stepc_liquid()$lipid_deuterated

  
  ###check end of filtering box
  output$checkendfiltering_liquid = reactive(
    return(is.null(stepc_liquid()))
  )
  outputOptions(output, "checkendfiltering_liquid", suspendWhenHidden = FALSE)
  
  
  observeEvent(stepc_liquid(),{
    updateSelectInput(session, "selcol_lipidfilt_liquid", choices = names(stepc_liquid()$lipid_filtered))
  })
  
  output$dtlipidfilterd_liquid = DT::renderDT({
    req(stepc_liquid())
    stepc_liquid()$lipid_filtered %>% purrr::pluck(input$selcol_lipidfilt_liquid)
  })
  
  
  
  ####PART 1.
  
  stepg_liquid = reactive({
    req(stepc_liquid())
    
    data_list_liquid <- lapply(stepc_liquid()$lipid_filtered, function(x) x %>%
                                 dplyr::select(Common_Name, Intensity))
    
    data_updated = list()
    for(k in 1:length(data_list_liquid)){
      data_updated[k] = list(data_list_liquid[[k]] %>% dplyr::rename(!!names(data_list_liquid)[k] := "Intensity"))
    }
    
    merged = data_updated %>% purrr::reduce(dplyr::full_join, by = "Common_Name") %>% 
      dplyr::rename("Lipids" = "Common_Name")
    
    #return stepg_liquid$concentration_matrix with everything else
    list(stepc_liquid(), concentration_matrix = list(concentration_matrix = merged)) %>% purrr::flatten()
    
    
  })
  
  
  
  #check data correctly loaded
  output$checkstepg_liquid = reactive(
    return(is.null(stepg_liquid()))
  )
  outputOptions(output, "checkstepg_liquid", suspendWhenHidden = FALSE)
  
  observeEvent(stepg_liquid(),{
    num_na = sum(is.na(stepg_liquid()$concentration_matrix))
    not_na = sum(!is.na(stepg_liquid()$concentration_matrix))
    
    if (num_na/(num_na+not_na) > 0.5){
      shinyWidgets::sendSweetAlert(session, title = "Too many NAs", type = "warning", width = "600px",
                                   text = tags$div(
                                     p(tags$div("More than 50% of the calibration matrix data are NAs.", style = "font-weight: bold;")),
                                     p(div("You can:",style = "text-align: left;")),
                                     p(tags$div(icon("caret-right"), " Press the ", strong("'Check concentration matrix'"), " button to visualize them.", style = "text-align: left;")),
                                     p(tags$div(icon("caret-right"), " Filter and imputate them in the ", strong("next step."), style = "text-align: left;")),
                                   )
      )
    }
    
  })
  
  
  stephliquid = mod_imputation_step_server("imputation_step_liquid", parent = session, stepg = stepg_liquid, data_type = "Peak Intensity")
  
  
  
  ######### Create sumexp object from Metabolomics Workbench #####
  
  mwdata = eventReactive(input$loadmw,{
    req(input$selmwid)
    showNotification(tagList(icon("cogs"), HTML("&nbsp;Downloading data from Metabolomics Workbench...")), type = "default")
    query = metabolomicsWorkbenchR::do_query('study','study_id', input$selmwid,'SummarizedExperiment')
    showNotification(tagList(icon("check"), HTML(paste0("&nbsp;Download of ",input$selmwid, " completed."))), type = "message")
    return(query)
  })
  
  #se tutti i data sono stati caricati
  output$checkmw = reactive(
    return(is.null(mwdata()))
  )
  outputOptions(output, "checkmw", suspendWhenHidden = FALSE)
  
  
  #assay
  mwassay = reactive({
    req(mwdata())
    if(input$selmwid == "ST000608"){
      asss = list()
      for(i in 1:2){
        roww = mwdata()[[i]] %>% SummarizedExperiment::rowData() %>% as.data.frame() %>% dplyr::select(metabolite_id, metabolite_name)
        aux_assay = mwdata()[[i]] %>% assay() %>% tibble::rownames_to_column("metabolite_id")
        aux_assay = dplyr::left_join(roww, aux_assay, by = "metabolite_id") %>% 
          dplyr::select(-metabolite_id) %>% dplyr::rename(Lipids = metabolite_name)
        rownames(aux_assay) = aux_assay$Lipids
        asss[[i]] = aux_assay
      }
      tot_ass = rbind(asss[[1]],asss[[2]])
      lip1 = stringr::str_replace(tot_ass$Lipids,";.*","")
      lip2 = stringr::str_replace(lip1,"_.*","")
      lip3 = stringr::str_replace(lip2, "PE\\(P-", "PEP\\(")
      index_L = grep("(?=^PE|PC|PS|PI|PG.*$)(?=.*\\b0:0\\b)",lip3, perl = T)
      lip3[index_L] = gsub("0:0/|/0:0", "", lip3[index_L])
      lip3[index_L] = paste0("L",lip3[index_L])
      tot_ass$Lipids = lip3
      aux_assay = tot_ass %>% dplyr::group_by(Lipids) %>% dplyr::summarise(dplyr::across(where(is.double), mean, na.rm = T)) %>% 
        as.data.frame()
    }else{
      roww = mwdata() %>% SummarizedExperiment::rowData() %>% as.data.frame() %>% dplyr::select(metabolite_id, metabolite_name)
      aux_assay = mwdata() %>% SummarizedExperiment::assay() %>% tibble::rownames_to_column("metabolite_id")
      aux_assay = dplyr::left_join(roww, aux_assay, by = "metabolite_id") %>% 
        dplyr::select(-metabolite_id) %>% dplyr::rename(Lipids = metabolite_name)
      #removing QC
      if(input$selmwid == "ST001073"){
        aux_assay = aux_assay %>% dplyr::select(!dplyr::starts_with("QC"))
      }
    }
    
    rownames(aux_assay) = aux_assay$Lipids
    if(input$selmwid == "ST001115" || input$selmwid == "ST000608"){
      sample_corrected = stringr::str_split_fixed(colnames(aux_assay[,-1]), "_",n = 3)
      sample_corrected = paste(paste(sample_corrected[,1], sample_corrected[,2], sep = "-"), sample_corrected[,3], sep = "_")
      colnames(aux_assay) = c("Lipids", sample_corrected)
    }
    aux_assay
  })
  
  #check na
  output$check_naassay_MW = reactive({
    req(mwassay())
    x = mwassay() %>% dplyr::select(where(is.double))
    num_na = sum(is.na(x))
    #not_na = sum(!is.na(x))
    if(num_na != 0){
      shinyWidgets::sendSweetAlert(session, title = "Missing values!", type = "warning", width = "600px",
                                   text = div("There are some missing values. You should filter and imputate 
                     them in order to prevent errors.", style = "font-weight: bold;"))
      TRUE
    }
  })
  outputOptions(output, "check_naassay_MW", suspendWhenHidden = FALSE)
  
  
  #### Check filtering
  filterstep_MW = eventReactive(input$gofilterna_MW,{
    req(mwassay())
    #filtered assay
    filt2 = list(concentration_matrix = as.data.frame(mwassay()))
    na_advise_lipidomics(
      out = filt2,
      na_filter_lip = as.numeric(input$na_filt_lip_ass_MW),
      na_filter_sam = as.numeric(input$na_filt_sam_ass_MW),
      imputation_met = "none"
    )
  })
  
  output$dt_filteredna_MW = renderDT({
    req(filterstep_MW())
    filterstep_MW()$concentration_matrix_filt
  })
  
  output$vimplot_filteredna_MW = renderPlot({
    req(filterstep_MW())
    filterstep_MW()$concentration_matrix_filt %>% dplyr::select(-1) %>% 
      VIM::aggr(combined = input$vimopz_filt_MW, cex.axis = 0.7, cex.lab = 0.7)
  })
  
  #box table dimensions after na filtering
  output$nadim1_MW = shinydashboard::renderInfoBox({
    dim = dim(filterstep_MW()$concentration_matrix[,-1])
    shinydashboard::infoBox(
      title = div(HTML(paste0("Table dimension", br(), "before filtering")), style = "color:white; font-size:100%;"),
      value = div(paste0(dim[1], " x ", dim[2]), style = "font-size:140%"),
      icon = icon("table"), color = "yellow", fill = TRUE)
  })
  
  output$nadim2_MW = shinydashboard::renderInfoBox({
    dim = dim(filterstep_MW()$concentration_matrix_filt[,-1])
    shinydashboard::infoBox(
      title = div(HTML(paste("Table dimension", br(), "after filtering")), style = "color:white; font-size:100%;"),
      value = div(paste0(dim[1], " x ", dim[2]), style = "font-size:140%"),
      icon = icon("table"), color = "green", fill = TRUE)
  })
  
  
  
  #### Create final output SumExp 
  sumexp_mw = eventReactive(input$make_sumexp_MW,{
    req(mwassay())
    #filtered assay
    filt2 = list(concentration_matrix = as.data.frame(mwassay()))
    if(input$togglefilterna_ass_MW == FALSE){
      filt = na_advise_lipidomics(
        out = filt2,
        na_filter_lip = 1,
        na_filter_sam = 1,
        imputation_met = "none"
      )
    }else{
      filt = na_advise_lipidomics(
        out = filt2,
        na_filter_lip = as.numeric(input$na_filt_lip_ass_MW),
        na_filter_sam = as.numeric(input$na_filt_sam_ass_MW),
        imputation_met = input$imput_method_ass_MW
      )
    }


    anal = list(
      lab_analyst = input$inlab_analyst,
      data_analyst = input$indata_analyst,
      pipe_ver = utils::packageVersion("ADViSELipidomics"), #"0.3.0", #golem::get_golem_version(),
      analysis_date = date()
    )
    
    if(input$selmwid == "ST000608"){
      coldata = mwdata()[[1]] %>% SummarizedExperiment::colData() %>% as.data.frame() %>% 
        dplyr::mutate(Sample_Type = gsub(" ","",Sample_Type))
    }
    
    if(input$selmwid == "ST001073"){
      #removing QC
      coldata = mwdata() %>% SummarizedExperiment::colData() %>% 
        as.data.frame() %>% dplyr::filter(!grepl("QC",local_sample_id)) %>% 
        dplyr::mutate(Crush = gsub(" ","",Crush))
    }
    
    if(input$selmwid == "ST001115"){
      coldata = mwdata() %>% SummarizedExperiment::colData() %>% 
        as.data.frame() %>% dplyr::mutate(Fraction = gsub(" ","",Fraction))
    }

    
    workbench_advise_lipidomics(filtered = filt, 
                                coldata = coldata,
                                metadata = anal,
                                data_type = input$typedata_assay_MW,
                                id = input$selmwid)
    
  })
  
  
  #check data correctly loaded
  output$checkstepmw = reactive(
    return(is.null(sumexp_mw()))
  )
  outputOptions(output, "checkstepmw", suspendWhenHidden = FALSE)
  
  
  observeEvent(input$gotosumexp_from_MW, {
    shinydashboard::updateTabItems(session, "sidebarmenu", "seedatatab")
  })
  
######### Create sumexp object from csv #####
  
  ###### COLDATA
  cold = reactive({
    req(input$coldatainput)
    ext <- tools::file_ext(input$coldatainput$name)
    if(ext != "xlsx"){
      shinyWidgets::show_alert("Invalid file!", "Please upload a .xlsx file", type = "error")
    }
    validate(need(ext == "xlsx", "Invalid file! Please upload a .xlsx file"))
    readxl::read_xlsx(input$coldatainput$datapath)
  })
  
  
  observeEvent(cold(),{
    if(!"SampleID" %in% colnames(cold())){
      showNotification(tagList(icon("times-circle"), HTML("&nbsp;Required 'SampleID' column. Check if exists.")), type = "error")
    }
    if(length(cold()) == 1){
      showNotification(tagList(icon("times-circle"), HTML("&nbsp;Something wrong in reading. Check the delimiter.")), type = "error")
    }
  })
  
  
  #check data correctly loaded
  output$checkcold = reactive(
    if("SampleID" %in% colnames(cold()) & length(cold()) > 1){
      showNotification(tagList(icon("check"), HTML("&nbsp;Target file loaded correctly!")), type = "message")
      FALSE #instead of is.null(), so if it's null returns TRUE
    }
  )
  outputOptions(output, "checkcold", suspendWhenHidden = FALSE)
  
  cold2 = reactive({
    req(cold())
    validate(need("SampleID" %in% colnames(cold()), "Required 'SampleID' column. Check if exists."))
    validate(need(length(cold()) > 1, "Something wrong in reading. Check the delimiter."))
    cold()
  })
  
  ###### ASSAY ##
  
  ass = reactive({
    req(input$assayinput)
    ext <- tools::file_ext(input$assayinput$name)
    if(ext != "xlsx"){
      shinyWidgets::show_alert("Invalid file!", "Please upload a .xlsx file", type = "error")
    }
    validate(need(ext == "xlsx", "Invalid file! Please upload a .xlsx file"))
    readxl::read_xlsx(input$assayinput$datapath)
  })
  
  
  observeEvent(ass(),{
    if(!"Lipids" %in% colnames(ass())){
      showNotification(tagList(icon("times-circle"), HTML("&nbsp;Required 'Lipids' column. Check if exists.")), type = "error")
    }
    if(length(ass()) == 1){
      showNotification(tagList(icon("times-circle"), HTML("&nbsp;Something wrong in reading. Check the delimiter.")), type = "error")
    }
  })
  
  
  #check data correctly loaded
  output$checkassay = reactive(
    if("Lipids" %in% colnames(ass()) & length(ass()) > 1){
      showNotification(tagList(icon("check"), HTML("&nbsp;Assay loaded correctly!")), type = "message")
      FALSE #al posto di is.null() che se è null restituisce TRUE
    }
  )
  outputOptions(output, "checkassay", suspendWhenHidden = FALSE)
  
  #check
  output$check_naassay = reactive({
    x = ass() %>% as.data.frame() %>% dplyr::select(where(is.double))
    num_na = sum(is.na(x))
    if(num_na != 0){
      shinyWidgets::sendSweetAlert(session, title = "Missing values!", type = "warning", width = "600px",
                                   text = div("There are some missing values. You should filter and imputate 
                     them in order to prevent errors.", style = "font-weight: bold;"))
      TRUE
    }
  })
  outputOptions(output, "check_naassay", suspendWhenHidden = FALSE)
  
  
  ass2 = reactive({
    req(ass())
    validate(need("Lipids" %in% colnames(ass()), "Required 'Lipids' column. Check if exists."))
    validate(need(length(ass()) > 1, "Something wrong in reading. Check the delimiter."))
    ass()
  })
  
  
  #### Check filtering NA
  filterstep_csv = eventReactive(input$gofilterna_csv,{
    req(ass2())
    #filtered assay
    filt2 = list(concentration_matrix = as.data.frame(ass2()))
    na_advise_lipidomics(
      out = filt2,
      na_filter_lip = as.numeric(input$na_filt_lip_ass),
      na_filter_sam = as.numeric(input$na_filt_sam_ass),
      imputation_met = "none"
    )
  })
  
  
  output$dt_filteredna_csv = renderDT({
    req(filterstep_csv())
    filterstep_csv()$concentration_matrix_filt
  })
  
  output$vimplot_filteredna_csv = renderPlot({
    req(filterstep_csv())
    filterstep_csv()$concentration_matrix_filt %>% dplyr::select(-1) %>% 
      VIM::aggr(combined = input$vimopz_filt_csv, cex.axis = 0.7, cex.lab = 0.7)
  })
  
  #box table dimensions after na filtering
  output$nadim1_csv = shinydashboard::renderInfoBox({
    dim = dim(filterstep_csv()$concentration_matrix[,-1])
    shinydashboard::infoBox(
      title = div(HTML(paste0("Table dimension", br(), "before filtering")), style = "color:white; font-size:100%;"),
      value = div(paste0(dim[1], " x ", dim[2]), style = "font-size:140%"),
      icon = icon("table"), color = "yellow", fill = TRUE)
  })
  
  output$nadim2_csv = shinydashboard::renderInfoBox({
    dim = dim(filterstep_csv()$concentration_matrix_filt[,-1])
    shinydashboard::infoBox(
      title = div(HTML(paste("Table dimension", br(), "after filtering")), style = "color:white; font-size:100%;"),
      value = div(paste0(dim[1], " x ", dim[2]), style = "font-size:140%"),
      icon = icon("table"), color = "green", fill = TRUE)
  })
  
  
  #### Final output SumExp
  sumexp_csv = eventReactive(input$make_sumexpcsv,{
    req(ass2(), cold2())
    #filtered assay
    filt2 = list(concentration_matrix = as.data.frame(ass2()))
    if(input$togglefilterna_ass == FALSE){
      filt = na_advise_lipidomics(
        out = filt2,
        na_filter_lip = 1,
        na_filter_sam = 1,
        imputation_met = "none"
      )
    }else{
      filt = na_advise_lipidomics(
        out = filt2,
        na_filter_lip = as.numeric(input$na_filt_lip_ass),
        na_filter_sam = as.numeric(input$na_filt_sam_ass),
        imputation_met = input$imput_method_ass
      )
    }
    
    #replicates
    tec_rep2 = stringr::str_count(cold2()$SampleID, "_")
    if(sum(tec_rep2) != 0){
      tec_rep = TRUE
    }else{
      tec_rep = FALSE
    }
    
    anal = list(
      lab_analyst = input$inlab_analyst,
      data_analyst = input$indata_analyst,
      pipe_ver = utils::packageVersion("ADViSELipidomics"), #"0.5.0", #golem::get_golem_version(),
      analysis_date = date()
    )
    
    out = purrr::flatten(
      list(filt,
           replicates = tec_rep,
           analysis = list(analysis = anal),
           targ = list(targets = list(targetfile_lipidomics = cold2()))
      )
    )
    
    sumexp_advise_lipidomics(out)
  })
  
  
  #check data correctly loaded
  output$checkstepcsv = reactive(
    return(is.null(sumexp_csv()))
  )
  outputOptions(output, "checkstepcsv", suspendWhenHidden = FALSE)
  
  observeEvent(input$gotosumexp_from_csv, {
    shinydashboard::updateTabItems(session, "sidebarmenu", "seedatatab")
  })
  
  
  #### Download handler for the download button
  output$downloadsumexp_csv <- downloadHandler(
    #put the file name with also the file extension
    filename = function() {
      paste0("summ_EXP_", Sys.Date(), ".rds")
    },
    
    # This function should write data to a file given to it by the argument 'file'.
    content = function(file) {
      summ = list(sumexp_data = sumexp_csv()$sumexp_data, 
                  sumexp_data_mean = sumexp_csv()$sumexp_data_mean, 
                  replicates = sumexp_csv()$replicates, 
                  data_type = input$typedata_assay)
      saveRDS(summ, file)
    }
  )
  
  
  
#____________________________________________________________________________
  
  output$checksteph = reactive(
    tryCatch({
      is.null(steph())
    },
    shiny.silent.error = function(e) {
      TRUE
    })
  )


  ####selection for the sumexp object (from file or step or csv) ####
  
  #import sumexp
  sumexp_import = reactive({
    req(input$finalsumexpinput)
    ext <- tools::file_ext(input$finalsumexpinput$name)
    if(ext != "rds"){
      shinyWidgets::show_alert("Invalid file!", "Please upload a .rds file", type = "error")
    }
    validate(need(ext == "rds", "Invalid file! Please upload a .rds file"))
    file = readRDS(file = input$finalsumexpinput$datapath)
    if(!is.null(file)){
      showNotification(tagList(icon("check"), HTML("&nbsp;Summarize Experiment successfully loaded!")), type = "message")
      return(file)
    }
  })
  

  
  sumexp_all = reactive({
    checksumimport = tryCatch({sumexp_import()
      FALSE
    },shiny.silent.error = function(e) {TRUE})
    
    checkcsv = tryCatch({sumexp_csv()
      FALSE
    },shiny.silent.error = function(e) {TRUE})
    
    checksteph = tryCatch({steph()
      FALSE
    },shiny.silent.error = function(e) {TRUE})
    
    checkmw = tryCatch({sumexp_mw()
      FALSE
    },shiny.silent.error = function(e) {TRUE})
    
    checkliquid = tryCatch({stephliquid()
      FALSE
    },shiny.silent.error = function(e) {TRUE})
    
    
    
    if(input$sel_inputdata_type == "lipidsearch" && checksteph == FALSE){
      if(input$type_lipsearch == "Yes"){
        type = "Concentration"
      }else{
        type = "Area"
      }
      list(sumexp_data = steph()$sumexp_data, sumexp_data_mean = steph()$sumexp_data_mean, replicates = steph()$replicates, data_type = type)
    }else if(input$sel_inputdata_type == "excel" && checkcsv == FALSE){
      list(sumexp_data = sumexp_csv()$sumexp_data, sumexp_data_mean = sumexp_csv()$sumexp_data_mean, replicates = sumexp_csv()$replicates, data_type = input$typedata_assay)
    }else if(input$sel_inputdata_type == "sumexp" && checksumimport == FALSE){
      sumexp_import()
    }else if(input$sel_inputdata_type == "mw" && checkmw == FALSE){
      sumexp_mw()
    }else if(input$sel_inputdata_type == "liquid" && checkliquid == FALSE){
      list(sumexp_data = stephliquid()$sumexp_data, 
           sumexp_data_mean = stephliquid()$sumexp_data_mean, 
           replicates = stephliquid()$replicates, 
           data_type = "Peak Intensity")
    }
  })
  
  


  sumexpdata = reactive({
    req(sumexp_all())
    sumexp_all()$sumexp_data
  })
  
  
  sumexpdatamean = reactive({
    req(sumexp_all())
    sumexp_all()$sumexp_data_mean
  })
  
  #check data correctly loaded
  output$check_replicates = reactive({
    req(sumexp_all())
    sumexp_all()$replicates  #TRUE if there are replicates
  })
  outputOptions(output, "check_replicates", suspendWhenHidden = FALSE)


  output$dtsumexp = renderDT({
    req(sumexpdata())
    if(input$summ_viewtable == TRUE && sumexp_all()$replicates == TRUE){
      x = sumexpdatamean()
    }else{
      x = sumexpdata()
    }
    if(input$sumexpselectobj == "colData"){
      SummarizedExperiment::colData(x) %>% as.data.frame()
    }else if(input$sumexpselectobj == "assays"){
      SummarizedExperiment::assay(x) %>% as.data.frame()
    }else if(input$sumexpselectobj == "rowData"){
      row = SummarizedExperiment::rowData(x) %>% as.data.frame()
      temp = gsub(").*", ")", row$Lipids)
      temp = gsub("_","/", temp)
      row$Lipids = temp
      render <- c(
        "function(data, type, row){",
        "  if(type === 'display' && data){",
        "    var a = '<a href=\"http://www.swisslipids.org/#/search/' + data + '\">' + data + '</a>';",
        "    return a;",
        "  } else {",
        "    return data;",
        "  }",
        "}"
      )
      DT::datatable(row, rownames = T, 
                options = list(
                  columnDefs = list(
                    list(targets = 1, render = DT::JS(render)),
                    list(targets = "_all", className = "dt-center")
                  )
                )
      )
    }else{
      x@metadata %>% as.data.frame()
    }
  },options = list(scrollX = TRUE, scrollY = "700px"))



  


  #######PCA ########


  pcadata = eventReactive(input$gopca,{
    req(sumexpdata())
    if(input$summarize_pca == TRUE && sumexp_all()$replicates == TRUE){
      data = SummarizedExperiment::assay(sumexpdatamean()) %>% t()
    }else{
      data = SummarizedExperiment::assay(sumexpdata()) %>% t()
    }
    
    if(input$logpcalc == TRUE){
      data = log2(data)
    }
    g = try(stats::prcomp(data, scale = input$scalepcalc))
    if(class(g) == "try-error"){
      showNotification("Error in performing PCA. Too many NAs.", type = "error")
      return(NULL)
    }else{
      showNotification("PCA performed.")
      return(g)
    }
  })

  observeEvent(pcadata(),{
    pcs = colnames(pcadata()$rotation)[1:10] #take first 10 pcs
    pcs  = pcs[!is.na(pcs)] #if there are less than 10pcs, remove the NA generated
    updateSelectInput(session, "firstPC", choices = pcs)
    #updateSelectInput(session, "secondPC", choices = pcs)
  })
  
  observeEvent(input$firstPC,{
    pcs = colnames(pcadata()$rotation)[1:10] #take first 10 pcs
    pcs  = pcs[!is.na(pcs)] #if there are less than 10pcs, remove the NA generated
    updateSelectInput(session, "secondPC", choices = pcs[!pcs %in% input$firstPC])
  })

  
  #slider for PCs choice
  output$sliderpclc <- renderUI({
    req(pcadata())
    pca = pcadata()
    sliderInput("selpcslc", "Number of Principal Components (PC)", min=1, max=length(pca$sdev), value=2, step = 1)
  })


  output$pcantopui = renderUI({
    req(pcadata())
    max = pcadata()$rotation[,1] %>% length()
    numericInput("ntop_loadings", paste0("ntop ", "(from 1 to ", max, " )"), value = 25, min = 1, max = max)
  })

  ###plot loadings
  output$loadingslc = plotly::renderPlotly({
    req(pcadata(), input$selpcslc)
    pca = pcadata()
    loadpca = as.data.frame(pca$rotation[, input$selpcslc]) #invece di loadings ci sono i rotation
    loadpca = tibble::rownames_to_column(loadpca)

    pcasdev = as.data.frame(round(pca$sdev^2/sum(pca$sdev^2)*100, 2))
    colnames(loadpca) = c("Compounds", paste0("PC", input$selpcslc))
    ###ordering
    loadpca = loadpca %>% dplyr::arrange(desc(abs(loadpca[,2])))
    loadpca$Compounds = factor(loadpca$Compounds, levels = loadpca$Compounds)

    #ntop
    loadpca = loadpca[1:input$ntop_loadings,]

    loadplot = ggplot(loadpca, aes(x = Compounds, y = loadpca[,2], fill = Compounds)) + geom_col() +
      labs(y = paste0("PC", input$selpcslc, " ", "(", pcasdev[as.numeric(input$selpcslc), ], "%", ")"), title = "Loadings") +
      theme(axis.text.x = element_text(angle = 315, hjust = 0))
    plotly::ggplotly(loadplot)
  })

  ###screeplot
  output$screeplotlc <- plotly::renderPlotly({
    req(pcadata())
    pca = pcadata()
    var = cumsum(100*pca$sdev^2/sum(pca$sdev^2))
    var = as.data.frame(cbind(var)) %>% tibble::rownames_to_column()
    colnames(var) = c("Principal_components", "Explained_variation")
    var$Principal_components = factor(var$Principal_components, levels = c(1:length(var$Principal_components)))


    screegg = ggplot(var, aes(Principal_components,Explained_variation)) +
      geom_line(colour = "red", group = 1, linetype = "dashed", size = 1) + geom_point(size = 4, colour = "red") +
      labs(x = "Principal components", y = "Explained variation (%)", title = "Screeplot") +
      scale_y_continuous(limits = c(0, 100), breaks = c(seq(0, 100, by = 10)))
    plotly::ggplotly(screegg)

  })

  
  pcadata_info = eventReactive(input$gopca,{
    req(pcadata())
    if(input$summarize_pca == TRUE && sumexp_all()$replicates == TRUE){
      coll = SummarizedExperiment::colData(sumexpdatamean()) %>% as.data.frame() %>% tibble::rownames_to_column() 
    }else{
      coll = SummarizedExperiment::colData(sumexpdata()) %>% as.data.frame() %>% tibble::rownames_to_column()
    }
    

    tjoin = pcadata()$x %>% as.data.frame() %>% tibble::rownames_to_column()
    
    dplyr::inner_join(tjoin, coll, by = "rowname") %>% tibble::column_to_rownames("rowname")
    
  })
  

  observeEvent(pcadata_info(), {
    data = colnames(dplyr::select(pcadata_info(), -dplyr::starts_with("PC")))
    if("Product_Batch" %in% data){
      sel = "Product_Batch"
    }else{sel = data[1]}
    updateSelectInput(session, "shpbiplotlc", choices = data, selected = sel)
    updateSelectInput(session, "colbiplotlc", choices = data, selected = sel)
    updateSelectInput(session, "col3dlc", choices = data)
  })

  ###biplot
  output$biplotlc = plotly::renderPlotly({
    req(pcadata(), input$colbiplotlc)
    data_info = pcadata_info()
    data_info[,input$shpbiplotlc] = as.factor(data_info[,input$shpbiplotlc])
    data_info[,input$colbiplotlc] = as.factor(data_info[,input$colbiplotlc])
    x_axs = as.numeric(gsub("PC", "", input$firstPC))
    y_axs = as.numeric(gsub("PC", "", input$secondPC))
    
    if(input$selbiplotlc == "Biplot"){
      temp = autoplot(pcadata(), x = x_axs, y = y_axs, data = data_info, shape = input$shpbiplotlc, colour = input$colbiplotlc, loadings = TRUE, loadings.colour = 'blue',
                      loadings.label = TRUE, loadings.label.size = 4, title = "Biplot") + theme(legend.title = element_blank())
    }else{
      temp = autoplot(pcadata(), x = x_axs, y = y_axs, title = "Plot", data = data_info, shape = input$shpbiplotlc, colour = input$colbiplotlc) +
        theme(legend.title = element_blank())
    }

    plotly::ggplotly(temp) %>% plotly::layout(legend = list(title = list(text = paste(input$colbiplotlc, input$shpbiplotlc, sep = ", "))))
  })


  ### plot 3D
  output$pca3dlc = plotly::renderPlotly({
    req(pcadata_info())
    data_info = pcadata_info()
    data_info[,input$col3dlc] = as.factor(data_info[,input$col3dlc])
    data_info %>% plotly::plot_ly(x = ~PC1, y = ~PC2, z= ~PC3, color = ~base::get(input$col3dlc), type = "scatter3d", mode = "markers") %>% 
      plotly::layout(legend = list(title = list(text = input$col3dlc)))
  })


####### Plots #####


  #####lipid class distribution plots
  output$lipidsplot = renderUI({
    req(sumexpdata())

    lipids = SummarizedExperiment::rowData(sumexpdata())$Lipids %>% stringr::str_split_fixed(pattern = "\\(", 2)  #\\( perchè è un carattere speciale
    lipids = lipids[,1] %>% as.data.frame()
    colnames(lipids) = "Lipid"

    if(input$selplotlip == 1){
      output$pie1 = plotly::renderPlotly({
        lipids %>% plotly::plot_ly(labels = ~Lipid, type= "pie", textposition = 'inside', textinfo = 'label+value',
                                   marker = list(line = list(color = '#FFFFFF', width = 1)),
                                   showlegend = FALSE, width = "600px") %>% plotly::layout(title = "Lipid class distribution")
      })
      plotly::plotlyOutput("pie1")
    } else if(input$selplotlip == 2){
      output$bar1 = plotly::renderPlotly({
        getPalette = grDevices::colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))

        dd = ggplot(data= lipids) + geom_bar(mapping = aes(x = Lipid, fill = Lipid)) +
          scale_y_continuous(breaks = scales::pretty_breaks()) + ggtitle("Lipid class distribution") +
          xlab("Lipids") + labs(fill="Lipid Classes") + theme(axis.title = element_text(face="bold", size = 13), axis.text.x = element_text(angle = 315, hjust = 0, size = 11)) +
          scale_fill_manual(values = getPalette(length(unique(lipids$Lipid))))

      })
      plotly::plotlyOutput("bar1", width = "80%")
    } else{
      # Spider plot
      output$assaggispider = renderPlot({
        tab = table(lipids) %>% as.data.frame()
        gggg = data.frame(rbind(tab$Freq))
        colnames(gggg) = paste(tab$lipids)
        maxx = max(tab$Freq)
        max2 = maxx
        while(max2%%5 != 0){
          max2 = max2+1
        }

        g2 = rbind("Max" = max2, "Min" = 0, gggg)
        rownames(g2) <- c("Max", "Min", "freq")

        lab = c(0, max2*2/10, max2*4/10, max2*6/10, max2*8/10, max2)

        create_beautiful_radarchart(g2,  caxislabels = lab, color = grDevices::hcl.colors(2, palette = "Dynamic"), title = "Lipid class distribution")

      }, width = 800, height = 600)
      plotOutput("assaggispider")
    }
  })
  
  
  
  ##### taxa barplot lipid classes
  
  data_for_taxa = reactive({
    req(sumexpdata())
    if(input$summ_taxabar == TRUE && sumexp_all()$replicates == TRUE){
      sumexpdatamean()
    }else{
      sumexpdata()
    }
  })
  
  observeEvent(data_for_taxa(),{
    updateSelectInput(session, "annot_taxa", choices = c("No", colnames(SummarizedExperiment::colData(data_for_taxa()))))
  })
  
  output$taxabarplot = plotly::renderPlotly({
    req(data_for_taxa())
    size = ifelse(input$summ_taxabar == TRUE, 10, 8)

    data1 = data_for_taxa()

    taxabar = data1 %>% assay() %>% tibble::rownames_to_column("Lipids") %>% 
      dplyr::mutate(Lipids = gsub("\\(.*", "", Lipids)) %>% dplyr::rename(Class = Lipids)
    
    taxa_long = taxabar %>% tidyr::pivot_longer(cols = 2:length(taxabar), names_to = "SampleID", values_to = "Concentration")
    taxa_long2 = taxa_long %>% dplyr::group_by(Class, SampleID) %>% dplyr::summarise(dplyr::across(Concentration, sum))
    
    if(input$annot_taxa != "No"){
      
      coldata = data_for_taxa() %>% SummarizedExperiment::colData() %>% as.data.frame() %>% 
        dplyr::select(SampleID, input$annot_taxa)
      taxa_long2 = dplyr::left_join(taxa_long2, coldata, by = "SampleID")

      taxa = ggplot(taxa_long2, aes(x = SampleID, y = Concentration, fill = Class)) + geom_col(position = "fill") + 
        scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) + ggtitle("Lipid class proportion") + 
        facet_grid(~get(input$annot_taxa), scales = "free", switch = "x") + ylab(paste(sumexp_all()$data_type, "proportion")) +
        theme(axis.text.x = element_text(angle = 315, hjust = 0, size = size, margin=margin(t=30)),legend.title = element_blank())
    }else{
      taxa = ggplot(taxa_long2, aes(x = SampleID, y = Concentration, fill = Class)) + geom_col(position = "fill") + 
        scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) + ylab(paste(sumexp_all()$data_type, "proportion")) +
        ggtitle("Lipid class proportion") + 
        theme(axis.text.x = element_text(angle = 315, hjust = 0, size = size),legend.title = element_blank())
    }

    plotly::ggplotly(taxa)
  })


  ##### lipid species distribution barplot #####
  
  
  observeEvent(sumexpdatamean(),{
    #lipids class
    updateSelectInput(session, "class_spec_dist", choices = unique(SummarizedExperiment::rowData(sumexpdatamean())$Class))
    
    #fill var
    coln = sumexpdatamean() %>% SummarizedExperiment::colData() %>% as.data.frame() %>% dplyr::select(!where(is.numeric)) %>% colnames()
    if("Product_Batch" %in% coln){
      sel = "Product_Batch"
    }else{sel = coln[1]}
    
    updateSelectInput(session, "fill_spec_dist", choices = coln, selected = sel)
  })
  
  
  output$lipspec_barplot = renderPlotly({
    req(sumexpdatamean())
    
    ass_mean = SummarizedExperiment::assay(sumexpdatamean())
    
    rowd = SummarizedExperiment::rowData(sumexpdatamean()) %>% as.data.frame() %>% dplyr::select(Lipids, Class)
    ass_mean2 = ass_mean %>% as.data.frame() %>% tibble::rownames_to_column("Lipids") %>% dplyr::left_join(rowd, by = "Lipids")
    
    
    #filter by class
    cold_mean = SummarizedExperiment::colData(sumexpdatamean()) %>% as.data.frame()
    
    ass_mean_filt = ass_mean2 %>% dplyr::filter(Class == input$class_spec_dist) %>% dplyr::select(-Class) %>% 
      tibble::column_to_rownames("Lipids") %>% t() %>% as.data.frame() %>% tibble::rownames_to_column("SampleID") %>% 
      dplyr::left_join(dplyr::select(cold_mean, SampleID, input$fill_spec_dist), by = "SampleID") %>% 
      tidyr::pivot_longer(cols = !c(input$fill_spec_dist, SampleID), names_to = "Lipids", values_to = "Value")
    

    if(input$summ_spec_dist == TRUE){
      #summarized
      temp2 = ggplot(ass_mean_filt, aes_string(x="Lipids", y = "Value", fill = input$fill_spec_dist)) +
        geom_bar(position = position_dodge(), stat = "summary",fun = "mean") +
        stat_summary(geom = "errorbar", fun.data = mean_se, position = position_dodge2(padding = 1.6, preserve = "single")) +
        theme(axis.text.x = element_text(angle = 315, hjust = 0), legend.title = element_blank()) + ylab(sumexp_all()$data_type)
    }else{
      #not summarized
      temp2 = ggplot(ass_mean_filt, aes_string(x = "Lipids", y = "Value", fill = input$fill_spec_dist, label = "SampleID")) +
        geom_col(position = position_dodge2(preserve = "single")) + ylab(sumexp_all()$data_type) +
        theme(axis.text.x = element_text(angle = 315, hjust = 0), legend.title = element_blank())
    }



    plotly::ggplotly(temp2)
  })

  ######boxplot lipide (es. cer) contro condizioni (sulf, bc)



  observeEvent(sumexpdata(), {
    updateSelectInput(session, "select_lipidplot", choices = SummarizedExperiment::rowData(sumexpdata())$Lipids)
    
    if("Product_Batch" %in% colnames(SummarizedExperiment::colData(sumexpdata()))){
      sel = "Product_Batch"
    }else{sel = colnames(SummarizedExperiment::colData(sumexpdata()))[1]}
    updateSelectInput(session, "select_fillboxplot", choices = colnames(SummarizedExperiment::colData(sumexpdata())), selected = sel)
  })


  output$lipcondplot = renderPlotly({
    req(sumexpdata())
    if(input$summ_lipidboxplot == TRUE && sumexp_all()$replicates == TRUE){
      data1 = sumexpdatamean()
    }else{
      data1 = sumexpdata()
    }
    
    if(input$log_lipidboxplot == TRUE){
      data = SummarizedExperiment::assay(data1) %>% log2()
    }else{
      data = SummarizedExperiment::assay(data1)
    }
    
    ass = data %>% t() %>% as.data.frame() %>% tibble::rownames_to_column("sample")
    cold = SummarizedExperiment::colData(data1) %>% as.data.frame() %>% tibble::rownames_to_column("sample")

    joinedt1 = dplyr::inner_join(ass, cold, by = "sample")

    Concentration = joinedt1[, colnames(joinedt1) %in% input$select_lipidplot]

    jok = data.frame(dplyr::select(joinedt1, input$select_fillboxplot), Concentration)
    jok[,input$select_fillboxplot] = as.factor(jok[,input$select_fillboxplot])
    
    if(input$addpoints_lipidboxplot == FALSE){
      gg = ggplot(jok, aes_string(x = input$select_fillboxplot, y = "Concentration", fill = input$select_fillboxplot)) +
        geom_boxplot() + ylab(input$select_lipidplot)
    }else{
      gg = ggplot(jok, aes_string(x = input$select_fillboxplot, y = "Concentration", fill = input$select_fillboxplot)) +
        geom_boxplot() + geom_point() + geom_jitter(width = 0.2) + ylab(input$select_lipidplot) 
    }
    plotly::ggplotly(gg)
  })

  
  #####scatterplot sample vs sample
  samples_for_update = reactive({
    req(sumexpdata())
    if(input$summ_scatt == TRUE && sumexp_all()$replicates == TRUE){
      SummarizedExperiment::assay(sumexpdatamean()) %>% colnames()
    }else{
      SummarizedExperiment::assay(sumexpdata()) %>% colnames()
    }
  })

  observeEvent(samples_for_update(), {
    req(samples_for_update())
    updateSelectInput(session, "sel_sample1_scatt", choices = samples_for_update())
    updateSelectInput(session, "sel_sample2_scatt", choices = samples_for_update())
  })


  output$scattsampleplot = renderPlotly({
    req(samples_for_update())
    xlab = ifelse(input$log_scatt == FALSE, input$sel_sample1_scatt, paste0("Log2(", input$sel_sample1_scatt, ")"))
    ylab = ifelse(input$log_scatt == FALSE, input$sel_sample1_scatt, paste0("Log2(", input$sel_sample2_scatt, ")"))

    if(input$summ_scatt == TRUE && sumexp_all()$replicates == TRUE){
      assay = SummarizedExperiment::assay(sumexpdatamean())
    }else{
      assay = SummarizedExperiment::assay(sumexpdata()) 
    }
    assay = assay %>% tibble::rownames_to_column("Lipid")
    

    if(input$log_scatt == TRUE){
      assay = assay %>% dplyr::mutate(dplyr::across(where(is.numeric), log2))
      }
    
    gg = ggplot(assay, aes_string(x = paste0("`",input$sel_sample1_scatt,"`"), y = paste0("`",input$sel_sample2_scatt,"`"), color = "Lipid")) +
      geom_point() + theme(legend.title = element_blank()) + labs(x = xlab, y = ylab)
    plotly::ggplotly(gg) %>% plotly::layout(legend = list(title = list(text = "Lipids")))
  })

  
  
  ############ HEATMAP
  
  dataforheatmap = reactive({
    req(sumexp_all())
    if(sumexp_all()$replicates == TRUE){
      sumexpdatamean()
    }else{
      sumexpdata()
    }
  })
  #slider for columns
  output$sliderheatcol <- renderUI({
    req(dataforheatmap())
    len = SummarizedExperiment::rowData(dataforheatmap())$Class %>% unique() %>% length() #numero di lipidi (poi si ruota)
    sliderInput("slidercolheat", "Column cluster number:", min=2, max = len, value=2, step = 1)
  })
  
  #slider for rows
  output$sliderheatrow <- renderUI({
    req(dataforheatmap())
    len = SummarizedExperiment::colData(dataforheatmap())$SampleID %>% unique() %>% length() #numero di sample (poi si ruota)
    sliderInput("sliderrowheat", "Row cluster number:", min = 2, max = len, value = 2, step = 1)
  })
  
  observeEvent(dataforheatmap(),{
    if("Product_Batch" %in% colnames(SummarizedExperiment::colData(dataforheatmap()))){
      sel = "Product_Batch"
    }else{sel = colnames(SummarizedExperiment::colData(dataforheatmap()))[1]}
    updateSelectInput(session, "selectannot_row", choices = colnames(SummarizedExperiment::colData(dataforheatmap())), selected = sel)
    
    #rowdata without NA column
    rowdata = SummarizedExperiment::rowData(dataforheatmap())[,colSums(is.na(SummarizedExperiment::rowData(dataforheatmap())))<nrow(SummarizedExperiment::rowData(dataforheatmap()))] %>%
      colnames()
    if("Class" %in% rowdata){
      sel2 = "Class"
    }else{sel2 = rowdata[1]}
    updateSelectInput(session, "selectannot_col", choices = rowdata, selected = sel2)
    
    
    #filtering
    filtcols = dataforheatmap() %>% SummarizedExperiment::colData() %>% as.data.frame() %>% 
      dplyr::select(where(~length(unique(.))>1)) %>% colnames()
    updateSelectInput(session, "filtheatmapcol",  choices = c("none", filtcols))
  })
  

  observeEvent(input$filtheatmapcol,{
    if(input$filtheatmapcol != "none"){
      vals = dataforheatmap() %>% SummarizedExperiment::colData() %>% as.data.frame() %>% 
        dplyr::pull(input$filtheatmapcol) %>% unique()%>% as.character()
      updateSelectInput(session, "filtheatmapval",  choices = vals)
    }
  })

  dataheat = reactive({
    req(dataforheatmap())
    make_heatmap(
      data = dataforheatmap(),
      filter = c(input$filtheatmapcol, input$filtheatmapval),
      add_rowannot = input$selectannot_row,
      add_colannot = input$selectannot_col,
      log_data = input$logheat,
      scale_data = input$selscaleheat,
      order_data = input$heatsort,
      dist_method = input$seldistheat,
      clust_method = input$selhclustheat,
      row_dend = input$rowdend,
      row_nclust = input$sliderrowheat,
      col_dend = input$columndend,
      col_nclust = input$slidercolheat,
      padding = c(4,2,2,15),
      unit_legend = input$unitlegend_ht
    )
  })

  
  observeEvent(input$makeheatmap,{
    InteractiveComplexHeatmap::makeInteractiveComplexHeatmap(input, output, session, dataheat(), heatmap_id  = "heatmap_output")
  })
  
  
  
  
  
  
  
#######Quality plots ######
  
  observeEvent(sumexpdata(), {
    #### SAMPLE
    #barplot density plot e boxplot
    if("Product_Batch" %in% colnames(SummarizedExperiment::colData(sumexpdata()))){
      sel = "Product_Batch"
    }else{sel = colnames(SummarizedExperiment::colData(sumexpdata()))[1]}
    updateSelectInput(session, "densplot_sampfill", choices = colnames(SummarizedExperiment::colData(sumexpdata())), selected = sel)

    ### LIPIDS
    #boxplot lipids
    updateSelectInput(session, "concbox_samp", choices = substr(SummarizedExperiment::colData(sumexpdata())$SampleID, 1, nchar(SummarizedExperiment::colData(sumexpdata())$SampleID)-2))
    updateSelectInput(session, "concbox_sampfill", choices = c("Class", colnames(SummarizedExperiment::colData(sumexpdata()))))

    updateSelectInput(session, "concbox_samp2", choices = substr(SummarizedExperiment::colData(sumexpdata())$SampleID, 1, nchar(SummarizedExperiment::colData(sumexpdata())$SampleID)-2))
    updateSelectInput(session, "concbox_sampfill2", choices = c("Class", colnames(SummarizedExperiment::colData(sumexpdata()))))
  })
  
  
  #se ci sono replicati metto di default la spunta a summarize data
  observeEvent(sumexp_all(),{
    if(sumexp_all()$replicates == TRUE){
      updateAwesomeCheckbox(session, "summ_qualplot", value = TRUE)
    }
  })
  
  


  
  ####### Samples _____________________________________________________
  
  
  #Concentration Barplot (Samples)
  output$concbarplot_samp = renderPlotly({
    req(sumexpdata(), input$densplot_sampfill)
    
    if(input$summ_qualplot == TRUE && sumexp_all()$replicates == TRUE){
      se_data1 = sumexpdatamean()
    }else{
      se_data1 <- sumexpdata()
    }
    se_data = SummarizedExperiment::assay(se_data1)
    
    sum_data <- data.frame(Sum_conc = apply(se_data, 2, function(x) sum(x)))
    if(input$plotsamp_log == TRUE){
      sum_data = log2(sum_data)
      y_sm = paste0("Log2(", sumexp_all()$data_type, " Sum)")
    } else {
      y_sm = paste0(sumexp_all()$data_type, " Sum)")
    }
    sum_data$SampleID <- factor(rownames(sum_data), levels = rownames(sum_data))
    sum_data <- merge(x = sum_data,
                      y = data.frame(SummarizedExperiment::colData(se_data1)[, c("SampleID", input$densplot_sampfill)]),
                      by = "SampleID", all.x = TRUE)
    sum_data[,input$densplot_sampfill] = as.factor(sum_data[,input$densplot_sampfill])

    sample_barplot <- ggplot2::ggplot(sum_data, aes_string(x = "SampleID", y = "Sum_conc", fill = input$densplot_sampfill)) +
      geom_bar(stat = "identity") +
      labs(x = "Samples", y = y_sm, title = paste0(sumexp_all()$data_type, " Barplot (Samples)")) +
      theme(axis.text.x = element_text(angle = -90, vjust = 0.5))

    plotly::ggplotly(sample_barplot)

  })

  
  #### sample boxplot sample with log conc #####
  output$boxplot_density = renderPlotly({
    req(sumexpdata())
    #se_data <- SummarizedExperiment::assay(sumexpdata())
    
    if(input$summ_qualplot == TRUE && sumexp_all()$replicates == TRUE){
      se_data1 = sumexpdatamean() 
    }else{
      se_data1 = sumexpdata()
    }
    se_data = SummarizedExperiment::assay(se_data1)

    box_data = log2(se_data)
    
    box_data = as.data.frame(t(box_data)) %>% tibble::rownames_to_column("Sample")
    

    box_data$temp <- SummarizedExperiment::colData(se_data1)[, input$densplot_sampfill]
    colnames(box_data)[colnames(box_data) == 'temp'] <- input$densplot_sampfill
    

    aux_name = box_data %>% dplyr::ungroup() %>% dplyr::select(-c("Sample", input$densplot_sampfill)) %>% colnames()
    box_data <- box_data %>% tidyr::gather(Lipid,Conc, all_of(aux_name)) #diventa long
    box_data <- cbind(box_data, Class = sapply(strsplit(box_data$Lipid,"\\("), `[`, 1))
    box_data[,input$densplot_sampfill] = as.factor(box_data[,input$densplot_sampfill])
    
    if(input$addpoints_boxdens == TRUE){
      au2 = ggplot(box_data, aes_string(y = "Conc",  x ="Sample", fill = input$densplot_sampfill)) + geom_boxplot() + 
        geom_point() + labs( y = "log2(Concentration)", title = "Boxplot") + coord_flip()
    }else{
      au2 = ggplot(box_data, aes_string(y = "Conc",  x ="Sample", fill = input$densplot_sampfill)) + geom_boxplot() + 
      labs( y = paste0("log2(",sumexp_all()$data_type,")"), title = "Boxplot") + coord_flip()
    }

    plotly::ggplotly(au2)
    
  })
  

  ##### density plot ##
  
  output$densplot_samp = plotly::renderPlotly({
    req(sumexpdata())

    if(input$summ_qualplot == TRUE && sumexp_all()$replicates == TRUE){
      se_data = SummarizedExperiment::assay(sumexpdatamean())
    }else{
      se_data = SummarizedExperiment::assay(sumexpdata()) 
    }
    
    box_data <- log2(se_data)
    box_data = as.data.frame(t(box_data)) %>% tibble::rownames_to_column("Sample")

    box_data = box_data %>% tidyr::gather(Lipid, Conc, all_of(colnames(box_data[,-1])))
    
    au2 = ggplot(box_data) + geom_density(aes_string(x = "Conc",  fill = "Sample"), alpha = 0.2) + 
      labs( x = paste0("log2(",sumexp_all()$data_type,")"), y = "density", title = "Density boxplot")
    
    plotly::ggplotly(au2)
  })
  
  
  ####### LIPIDS #_______________________________________________________________
  
    
  # Concentration Barplot (Lipids)
  output$concbarplot_lip = renderPlotly({
    req(sumexpdata())

    se_data <- SummarizedExperiment::assay(sumexpdata())
    if(input$concbar_log == TRUE){
      cv_data = log2(se_data)
      y_cv = paste0("% CV Log2(",sumexp_all()$data_type,")")
    } else {
      cv_data = se_data
      y_cv = paste0("% CV ",sumexp_all()$data_type)
    }
    cv_data <- data.frame(cv_conc = apply(cv_data, 1, function(x) sd(as.numeric(x))/mean(as.numeric(x)) * 100))
    cv_data$Lipid <- factor(rownames(cv_data), levels = rownames(cv_data))
    cv_data$Class <- SummarizedExperiment::rowData(sumexpdata())$Class

    plot = ggplot2::ggplot(cv_data, aes(x = forcats::fct_rev(Lipid), y = cv_conc, fill = Class)) +
      geom_bar(stat = "identity") + labs(x = "Lipids", y = y_cv, title = paste0(sumexp_all()$data_type, " Barplot (Lipids)")) +
      theme(axis.text.x = element_text(angle = -90, size = 7))
    plotly::ggplotly(plot)

  })

  
  #Plots: Concentration Boxplot (Lipids)
  output$concboxplot_samp = renderPlotly({
    req(sumexpdata())

    se_data <- SummarizedExperiment::assay(sumexpdata())

    if(input$concbox_log == TRUE){
      box_data = log2(se_data)
      y_box = paste0("Log2(",sumexp_all()$data_type," Sum)")
    } else {
      box_data = se_data
      y_box = paste0(sumexp_all()$data_type," Sum")
    }

    box_data <- as.data.frame(t(box_data)) %>% tibble::rownames_to_column("Sample") %>%
      tidyr::separate(Sample, into = c("Sample", "rep"), sep = "_", fill = "right") %>% dplyr::select(-rep) %>%
      dplyr::group_by(Sample)
    if(input$concbox_sampfill != "Class"){
      box_data$temp <- as.factor(SummarizedExperiment::colData(sumexpdata())[, input$concbox_sampfill])
      colnames(box_data)[colnames(box_data) == 'temp'] <- input$concbox_sampfill
    }

    aux_name <- colnames(box_data)[-c(1,dim(box_data)[2])]
    box_data <- box_data %>% tidyr::gather(Lipid,Conc,all_of(aux_name))
    box_data <- cbind(box_data, Class = sapply(strsplit(box_data$Lipid,"\\("), `[`, 1))

    bx_data2 = dplyr::filter(box_data, Sample == input$concbox_samp)
    
    if(input$concbox_points == FALSE){
     plot = ggplot2::ggplot(bx_data2, aes_string(x = "Class", y = "Conc", fill = input$concbox_sampfill)) +geom_boxplot() +
      labs(y = y_box, title = paste0(sumexp_all()$data_type, " Boxplot (Lipids)")) 
    }else{
      plot = ggplot2::ggplot(bx_data2, aes_string(x = "Class", y = "Conc", fill = input$concbox_sampfill)) + geom_boxplot() +
        geom_point() + geom_jitter(width = 0.2) + labs( y = y_box, title = paste0(sumexp_all()$data_type, " Boxplot (Lipids)")) 
    }

    
    plotly::ggplotly(plot)

  })

  #second boxplot
  output$concboxplot_samp2 = renderPlotly({
    req(sumexpdata())

    se_data <- SummarizedExperiment::assay(sumexpdata())

    if(input$concbox_log2 == TRUE){
      box_data = log2(se_data)
      y_box = paste0("Log2(",sumexp_all()$data_type," Sum)")
    } else {
      box_data = se_data
      y_box = paste0(sumexp_all()$data_type," Sum")
    }

    box_data <- as.data.frame(t(box_data)) %>% tibble::rownames_to_column("Sample") %>%
      tidyr::separate(Sample, into = c("Sample", "rep"), sep = "_", fill = "right") %>% dplyr::select(-rep) %>%
      dplyr::group_by(Sample)
    if(input$concbox_sampfill2 != "Class"){
      box_data$temp <- SummarizedExperiment::colData(sumexpdata())[, input$concbox_sampfill2]
      colnames(box_data)[colnames(box_data) == 'temp'] <- input$concbox_sampfill2
    }

    aux_name <- colnames(box_data)[-c(1,dim(box_data)[2])]
    box_data <- box_data %>% tidyr::gather(Lipid,Conc,all_of(aux_name))
    box_data <- cbind(box_data, Class = sapply(strsplit(box_data$Lipid,"\\("), `[`, 1))

    bx_data2 = dplyr::filter(box_data, Sample == input$concbox_samp2)
    
    if(input$concbox_points2 == FALSE){
      plot = ggplot2::ggplot(bx_data2, aes_string(x = "Class", y = "Conc", fill = input$concbox_sampfill2)) +geom_boxplot() +
      labs( y = y_box, title = paste0(sumexp_all()$data_type, " Boxplot (Lipids)"))
    }else{
      plot = ggplot2::ggplot(bx_data2, aes_string(x = "Class", y = "Conc", fill = input$concbox_sampfill2)) +geom_boxplot() +
        geom_point() + geom_jitter(width = 0.2) + labs( y = y_box, title = paste0(sumexp_all()$data_type, " Boxplot (Lipids)"))
    }
    
    plotly::ggplotly(plot)

  })

  
  #variables
  output$checkadd2var = reactive({
    if(input$add2var %%2 == 0){
      "onevar"
    }else{"twovar"}
  })
  outputOptions(output, "checkadd2var", suspendWhenHidden = FALSE)
  
  observeEvent(input$add2var,{
    if(input$add2var %%2 == 1){
      updateButton(session, "add2var",label = HTML("&nbsp;Remove"), style = "danger", icon("minus")) 
    }else{
      updateButton(session, "add2var", label = HTML("&nbsp;Add"), style="success", icon("plus"))
    }
  })
  
  
  #Batch variable
  output$checkadd2batch = reactive({
    if(input$add2batch %%2 == 0){
      "onebatch"
    }else{"twobatch"}
  })
  outputOptions(output, "checkadd2batch", suspendWhenHidden = FALSE)
  
  observeEvent(input$add2batch,{
    if(input$add2batch %%2 == 1){
      updateButton(session, "add2batch",label = HTML("&nbsp;Remove"), style = "danger", icon("minus")) 
    }else{
      updateButton(session, "add2batch", label = HTML("&nbsp;Add"), style="success", icon("plus"))
    }
  })

#### Clustering #####
  
  sumexpde_forclust = reactive({
    req(sumexp_all())
    if(input$cluster_summar == FALSE){
      data = sumexp_all()$sumexp_data
    }else{
      data = sumexp_all()$sumexp_data_mean
    }
    if(input$clust_transp == "Samples"){
      data = data %>% SummarizedExperiment::assay() %>% t()
    }else{
      data = data %>% SummarizedExperiment::assay()
    }

    if(input$clust_scale == TRUE){
      scale(data)
    }else{data}
  })
  

  
  # plots for clusters number choice
  output$numclustergraph = renderPlot({
    req(sumexpde_forclust())
    
    if(input$selclusthmorfo == "K-means"){
      p1 = factoextra::fviz_nbclust(sumexpde_forclust(), stats::kmeans, method = "gap_stat")
      p2 = factoextra::fviz_nbclust(sumexpde_forclust(), stats::kmeans, method = "wss")
      p3 = factoextra::fviz_nbclust(sumexpde_forclust(), stats::kmeans, method = "silhouette")
    } else if(input$selclusthmorfo == "PAM"){
      p1 = factoextra::fviz_nbclust(sumexpde_forclust(), cluster::pam, method = "gap_stat")
      p2 = factoextra::fviz_nbclust(sumexpde_forclust(), cluster::pam, method = "wss")
      p3 = factoextra::fviz_nbclust(sumexpde_forclust(), cluster::pam, method = "silhouette")
    } else{
      p1 = factoextra::fviz_nbclust(sumexpde_forclust(), cluster::clara, method = "gap_stat")
      p2 = factoextra::fviz_nbclust(sumexpde_forclust(), cluster::clara, method = "wss")
      p3 = factoextra::fviz_nbclust(sumexpde_forclust(), cluster::clara, method = "silhouette")
    }
    
    if(input$selclustmethod == "Partitioning"){
      gridExtra::grid.arrange(p1, p2, p3, ncol = 2)
    } else{
      meth = c("single","complete","ward.D","ward.D2")
      d = stats::dist(sumexpde_forclust())
      par(mfrow=c(2,2))
      for(i in seq(1,4)){
        hs = stats::hclust(d, method = meth[i])
        plot(hs$height, pch=16, main = meth[i], ylab = "Height")
      }
    }
  })
  
  #data cluster
  output$plotcluster = renderPlot({
    req(sumexpde_forclust())
    if(input$selclustmethod == "Partitioning"){
      if(input$selclusthmorfo == "K-means"){
        clust = stats::kmeans(sumexpde_forclust(), centers = input$selnumclust, nstart = 25)
      } else if (input$selclusthmorfo == "PAM"){
        clust = cluster::pam(sumexpde_forclust(), k = input$selnumclust)
      } else {
        clust = cluster::clara(sumexpde_forclust(), k = input$selnumclust)
      }
      factoextra::fviz_cluster(clust, data = sumexpde_forclust(), ellipse.type = "t", palette = "jco", ggtheme = theme_minimal())
    } else {
      hcluster = factoextra::eclust(sumexpde_forclust(), "hclust", hc_method = input$selhclustmeth, k = input$selnumclust)
      p1 = factoextra::fviz_dend(hcluster, palette = "jco", rect = TRUE, show_labels = T, cex = 0.7, ggtheme = theme_minimal())
      p2 = factoextra::fviz_silhouette(hcluster)
      gridExtra::grid.arrange(p1, p2, ncol = 2)
    }
  })
  
  
  
    ##### EXP Design #####

  #update batch method
  observeEvent(input$expdes_batch_type, {
    if(input$expdes_batch_type == "fit"){
      updateSelectInput(session, "batch_meth", choices = c("given", "estimate"))
    }else if(input$expdes_batch_type == "remove"){
      updateSelectInput(session, "batch_meth", choices = c("limma", "combat_param", "combat_nonparam"))
    }
  })


  sumexpde_forcoldata = reactive({
    req(sumexp_all())
    if(input$expdes_summar == TRUE && sumexp_all()$replicates == TRUE){
      sumexpdatamean()
    }else{
      sumexpdata()
    }
  })


  observeEvent(sumexpde_forcoldata(), {
    datavar1 = sumexpde_forcoldata() %>% SummarizedExperiment::colData() %>%
      as.data.frame() %>% dplyr::select(!where(is.double)) %>% dplyr::select(where(function(x) length(unique(x))>1))
    
    datavar2 = sumexpde_forcoldata() %>% SummarizedExperiment::colData() %>%
      as.data.frame() %>% dplyr::select(where(function(x) length(unique(x))>1))
    
    updateSelectInput(session, "expdes_design_vars", choices = colnames(datavar1))
    updateSelectInput(session, "expdes_design_vars2", choices = colnames(datavar2))
    #batch
    updateSelectInput(session, "expdes_batch_var1", choices = colnames(SummarizedExperiment::colData(sumexpde_forcoldata())))
    updateSelectInput(session, "expdes_batch_var2", choices = colnames(SummarizedExperiment::colData(sumexpde_forcoldata())))
  })

  varsde = reactive({
    req(input$expdes_design_vars)
    if(input$add2var %%2 == 1){
      vars = c(input$expdes_design_vars, input$expdes_design_vars2)
    }else{
      vars = input$expdes_design_vars
    }
    vars
  })
  
  #update batch method
  observeEvent(varsde(), {
    if(length(varsde()) == 2){
      updateSelectInput(session, "expdes_batch_type", choices = c("remove"))
    }else{
      updateSelectInput(session, "expdes_batch_type", choices = c("remove", "fit"))
    }
  })

  batch_type1 = reactive({
    if(input$expdes_batch_effect == FALSE){
      "none"
    }else{
      input$expdes_batch_type
    }
  })

  #all batches here
  batches = reactive({
    if(input$expdes_batch_effect == TRUE){
      if(input$add2batch %%2 == 1){
        batchs = c(input$expdes_batch_var, input$expdes_batch_var2)
      }else{ batchs = input$expdes_batch_var1 }
    }else{ batchs = "" }
    batchs
  })


  ntotvars = reactive({
    req(varsde())
    if(input$batch_meth == "given"){
      totvar = c(varsde(),batches())
    }else{
      totvar = c(varsde())
    }
    totvar = totvar[totvar != ""]
    length(totvar)
  })

  output$checkntotvars = reactive(
    ntotvars()
  )
  outputOptions(output, "checkntotvars", suspendWhenHidden = FALSE)

  #### write table ##
  #first var
  firstvar = eventReactive(input$expdes_design_vars,{
    req(sumexpde_forcoldata()) #firstvariab
    sumexpde_forcoldata() %>% SummarizedExperiment::colData() %>% as.data.frame() %>%
      dplyr::select(input$expdes_design_vars) %>% unique()
  })
  
  observeEvent(firstvar(),{
    updateSelectInput(session, "firstelement", label = paste0("Level 1 for ", input$expdes_design_vars), choices = firstvar())
    updateSelectInput(session, "secondelement1", label = paste0("Level 2 for ", input$expdes_design_vars), choices = firstvar())
  })

  #second var
  secondvar = reactive({
    req(sumexpde_forcoldata(), ntotvars()) #input$secondvariab
    if(ntotvars() == 2){
      if(input$batch_meth == "given"){
        sumexpde_forcoldata() %>% SummarizedExperiment::colData() %>% as.data.frame() %>% 
          dplyr::select(input$expdes_batch_var1) %>% unique()
      }else{
        sumexpde_forcoldata() %>% SummarizedExperiment::colData() %>% as.data.frame() %>% 
          dplyr::select(input$expdes_design_vars2) %>% unique()
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
    req(sumexpde_forcoldata())
    if(input$fixedlevel == "Level 2"){
      firstvar()
    }else{
      secondvar()
    }
  })


  observeEvent(thirdvar(),{
    updateSelectInput(session, "thirdelement1", label = paste("Level 1 for", colnames(thirdvar())), choices = thirdvar())
    updateSelectInput(session, "thirdelement2", label = paste("Level 2 for",colnames(thirdvar())), choices = thirdvar())
  })

  observeEvent(secondvar(),{
    req(firstvar())
    updateRadioButtons(session, "fixedlevel", choiceNames  = c(input$expdes_design_vars, colnames(thirdvar())), choiceValues = c("Level 1", "Level 2"))
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


  

  file_contrast = eventReactive(input$submit_written,{
    req(writtenlist())
    
    if(ntotvars() == 2){
    tocheck = stringr::str_split_fixed(writtenlist()$contrast, pattern = "=", n = 2)[,1] %>% 
      stringr::str_split_fixed(pattern = "vs", n = 2) %>% as.vector()
    
    if(input$batch_meth == "given"){
      totvar = c(varsde(),batches())
      }else{
      totvar = c(varsde())}
    totvar = totvar[totvar != ""]
    
    united = sumexpde_forcoldata() %>% SummarizedExperiment::colData() %>% 
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

  
  observeEvent(sumexp_all(),{
    if(sumexp_all()$data_type != "Concentration"){
      updateAwesomeRadio(session, "expdes_bs_norm", choices = c("none", "scale", "quantile"), inline = TRUE)
    }else{
      updateAwesomeRadio(session, "expdes_bs_norm", choices = c("none", "scale"), inline = TRUE)
    }
  })
  
  
  #expdes_summar
  expdesign = eventReactive(input$runDE,{
    req(sumexp_all(), file_contrast(), varsde())

    #se non ci sono repliche faccio finta che sia mediato
    if(sumexp_all()$replicates == FALSE){
      exp_summ = TRUE
    }else{
      exp_summ = input$expdes_summar
    }

    exp1 = expdesign_advise_lipidomics(
      out = sumexp_all(),
      design_vars = varsde(),
      rep_mean = exp_summ,
      file_contrast = file_contrast(),
      batch_type = batch_type1(),
      batch_vars = batches(),
      batch_method = input$batch_meth
      )
    tryCatch({
      diffexp_advise_lipidomics(
        out = exp1,
        rep_mean = exp_summ,
        rep_effect = input$expdes_repeffect,
        bs_norm = input$expdes_bs_norm,
        batch_vars = batches(),
        batch_type = batch_type1(),
        batch_method = input$batch_meth,
        thresh = input$expdes_thresh,
        decide_met = input$expdes_decide_met
      )
    },error = function(err){
      if(err$message == "No residual degrees of freedom in linear model fits"){
        shinyWidgets::show_alert("No residual degrees of freedom in linear model fits.
                               We suggest to use the model with replicates.", type = "error")
      }else{
        shinyWidgets::show_alert(paste0(err), type = "error")
      }

    })

  })


  
  output$checktoptable = reactive(
    is.null(expdesign()$limma_result)
  )
  outputOptions(output, "checktoptable", suspendWhenHidden = FALSE)
  


  observeEvent(expdesign(), {
    updateSelectInput(session, "expdes_colmaplot", choices = colnames(expdesign()$contrast_matrix))
    updateSelectInput(session, "sel_toptable", choices = colnames(expdesign()$contrast_matrix))
  })
  
  output$toptable = DT::renderDT({
    req(expdesign(), input$sel_toptable)
    expdesign()$limma_result[[input$sel_toptable]]
  })
  
  #### Download handler for the download button
  output$downloadtoptable <- downloadHandler(
    #put the file name with also the file extension
    filename = function() {
      paste0("TopTable", Sys.Date(), ".xlsx")
    },
    
    # This function should write data to a file given to it by the argument 'file'.
    content = function(file) {
      openxlsx::write.xlsx(expdesign()$limma_result[[input$sel_toptable]], file)
    }
  )
  


  expdes_data_forplot = reactive({
    req(expdesign(), input$expdes_colmaplot)
    expdesign()$limma_result[[input$expdes_colmaplot]]
  })


  mod_ma_volcano_plot_server("ma_volcano_plot1",
                             data_input = expdes_data_forplot,
                             contrast = reactive(input$expdes_colmaplot),
                             logfc_val = reactive(input$expdes_lfc),
                             thresh = reactive(input$expdes_thresh))


  mod_ma_volcano_plot_server("ma_volcano_plot2",
                             data_input = expdes_data_forplot,
                             contrast = reactive(input$expdes_colmaplot),
                             logfc_val = reactive(input$expdes_lfc),
                             thresh = reactive(input$expdes_thresh))


  observeEvent(expdesign(),{
    updateSelectInput(session, "venncontrast", choices = colnames(expdesign()$test_result), selected = colnames(expdesign()$test_result))
  })

# Differential expressed lipids Venn diagram
output$venndiag = renderUI({
  req(expdesign())
  test_result = expdesign()$test_result[,input$venncontrast]
  validate(
    need(ncol(test_result) >= 2 && ncol(test_result) <6, "Venn diagram for DC lipids is plotted from two to five contrasts!")
  )
    aux_con <- tibble::rownames_to_column(as.data.frame(test_result)) %>%
      dplyr::mutate_all(list(~dplyr::case_when(. != 0 ~ rowname,  TRUE ~ "0"))) %>%
      dplyr::select(-rowname)
    aux_con <- lapply(as.list(aux_con), function(x) x[x != "0"])
    venn <- ggVennDiagram::Venn(aux_con)
    data_venn <- ggVennDiagram::process_data(venn)

    output$vennplot = renderPlot({
      if(input$venngrad == TRUE){

      data_venn@setEdge$id <- rep("3",length(data_venn@setEdge$id))
      vd_plot <- ggVennDiagram::ggVennDiagram(aux_con,
                                              category.names = colnames(test_result),
                                              label_percent_digit = 1,
                                              label_alpha = 0,
                                              label_color = "black",
                                              edge_size = 0) +
        ggplot2::geom_sf(aes(color = id), data = ggVennDiagram::venn_setedge(data_venn), show.legend = FALSE) +
        ggplot2::scale_fill_gradient(low = "white", high = "blue") +
        ggplot2::scale_x_continuous(expand = expansion(mult = .2))

    } else {
      region_label <- data_venn@region %>%
        dplyr::filter(.data$component == "region") %>%
        dplyr::mutate(percent = paste(round(.data$count*100/sum(.data$count),
                                            digits = 1),"%", sep="")) %>%
        dplyr::mutate(both = paste(.data$count,paste0("(",.data$percent,")"),sep = "\n"))

      vd_plot <- ggplot() +
        geom_sf(aes(fill = id), data = ggVennDiagram::venn_region(data_venn), show.legend = FALSE) +
        geom_sf(color = "black", size = 1, data = ggVennDiagram::venn_setedge(data_venn), show.legend = FALSE) +
        geom_sf_text(aes(label = name), data = ggVennDiagram::venn_setlabel(data_venn)) +
        geom_sf_text(aes(label = both), data = region_label) +
        scale_x_continuous(expand = expansion(mult = .2)) +
        theme_void()

    }
    vd_plot
    })

    output$dtvenn = DT::renderDT({
      table_venn <- stringi::stri_list2matrix(data_venn@region$item, byrow = FALSE)
      colnames(table_venn) <-  gsub("\\.+","\U2229",data_venn@region$name)
      #table_venn = gsub(").*", ")", table_venn)
      #table_venn = gsub("_","/", table_venn)
      render <- c(
        "function(data, type, row){",
        "  if(type === 'display' && data){",
        "    var newstr = data.replace(/_/g, '/').replace(/\\).*/g, ')');",
        "    var a = '<a href=\"http://www.swisslipids.org/#/search/' + newstr + '\">' + data + '</a>';",
        "    return a;",
        "  } else {",
        "    return data;",
        "  }",
        "}"
      )
      DT::datatable(table_venn, rownames = F,
                options = list(
                  columnDefs = list(
                    list(targets = "_all", render = DT::JS(render)),
                    list(targets = "_all", className = "dt-center")
                  )
                )
      )
    },options = list("pageLength" = 15))

    fluidRow(
      column(5, br(), br(), br(), plotOutput("vennplot")),
      column(7, div(DT::DTOutput("dtvenn"), style = "overflow-y: scroll; overflow-x: scroll;"))
    )
})



output$upsetplot = renderPlot({
  req(expdesign())
  validate(
    need(length(colnames(expdesign()$test_result)) > 1, "Upset Plot for DC lipids is plotted from more than one contrast!")
  )
  validate(
    need(all(expdesign()$test_result == 0) == FALSE, "No DC Lipids into the contrasts.")
  )

  comb_mat <- ComplexHeatmap::make_comb_mat(expdesign()$test_result)
  cs = ComplexHeatmap::comb_size(comb_mat)
  cm_degree = ComplexHeatmap::comb_degree(comb_mat)
  aux_color <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(12,"Paired"))

  ht = ComplexHeatmap::UpSet(comb_mat, pt_size = grid::unit(5,"mm"), lwd = 3,
             comb_col = aux_color(length(cm_degree)),
             comb_order = order(cm_degree, -cs),
             top_annotation = HeatmapAnnotation(
               "Intersection size" = ComplexHeatmap::anno_barplot(cs,
                                                  ylim = c(0, max(cs)*1.1),
                                                  border = FALSE,
                                                  gp = grid::gpar(fill = aux_color(length(cm_degree))),
                                                  height = grid::unit(10, "cm")
               ),
               annotation_name_side = "left",
               annotation_name_rot = 0),
             right_annotation = ComplexHeatmap::upset_right_annotation(comb_mat, add_numbers = TRUE),
  )
  ht = ComplexHeatmap::draw(ht)
  od = ComplexHeatmap::column_order(ht)
  ComplexHeatmap::decorate_annotation("Intersection size", {
    grid::grid.text(cs[od], x = seq_along(cs), y = grid::unit(cs[od], "native") + grid::unit(4, "pt"),
              default.units = "native", just = c("center", "bottom"),
              gp = grid::gpar(fontsize = 10, col = "#404040"), rot = 0)
  })
})



##### Enrichment ####
observeEvent(expdesign(),{
  updateSelectInput(session, "enrich_selcont", choices = names(expdesign()$limma_result))
})

output$enrich_plot = renderPlot({
  req(expdesign())
  enrichment_advise_lipidomics(
    out = expdesign(),
    rank_c = input$rank_c,
    thresh = input$enrich_thresh,
    k = input$enrich_selcont
  )

})


####### PLS-DA ######


sumexpde_forplsda = reactive({
  req(sumexp_all())
  if(input$summ_plsda == TRUE && sumexp_all()$replicates == TRUE){
    sumexpdatamean()
  }else{
    sumexpdata()
  }
})

observeEvent(sumexpde_forplsda(),{
  data = sumexpde_forplsda() %>% SummarizedExperiment::colData() %>% as.data.frame()
  #selectin columns with more than one level
  data = data[, sapply(data, function(col) length(unique(col))) > 1]
  if("Product_Batch" %in% colnames(data)){
    sel = "Product_Batch"
  }else{
    sel= colnames(data)[1]
  }
  updateSelectInput(session, "selcol_pls", choices = colnames(data), selected = sel)
})


#number of allowed components (usually K-1)
output$plsda_ncompui = renderUI({
  req(sumexpde_forplsda())
  Y = as.factor(colData(sumexpde_forplsda())[,input$selcol_pls])
  max = length(unique(Y))
  numericInput("plsda_ncomp", paste0("Number of components ", "(from 2 to ", max, " )"), value = 2, min = 1, max = max, step = 1)
})


#pls-da
plsda = eventReactive(input$gopls,{
  req(sumexpde_forplsda(), input$plsda_ncomp)
  X = sumexpde_forplsda() %>% SummarizedExperiment::assay() %>% as.matrix() %>% t()
  Y = as.factor(SummarizedExperiment::colData(sumexpde_forplsda())[,input$selcol_pls])
  if(nrow(X)<30){
    showNotification(tagList(icon("exclamation-circle"), HTML("&nbsp;There are less than 30 samples. PLS-DA could not be statistically significant.")), type = "warning")
  }
  mixOmics::plsda(X, Y, ncomp = input$plsda_ncomp)
})


#plotIndiv
output$plotindiv_pls = renderPlot({
  req(plsda())
  validate(need(input$plsda_ncomp > 1, "PLS-DA plots require at least 2 components."))
  if(input$pls_back_ellipse == "Background"){
    background <- mixOmics::background.predict(plsda(), comp.predicted = 2, dist = "max.dist")
    len = length(unique(plsda()$Y))
    if(len > 8){
      palette = grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, name  = "Dark2"))
      colors = palette(len)
    }else{
      colors = RColorBrewer::brewer.pal(len, name  = "Dark2")
      colors = colors[1:len]
    }
    mixOmics::plotIndiv(plsda(), comp = 1:2, ind.names = input$indname_pls, legend = TRUE, background = background,
                        col = colors)
  }else if(input$pls_back_ellipse == "Ellipse"){
    mixOmics::plotIndiv(plsda(), ellipse = T, legend = T, ind.names = input$indname_pls)
  }else{
    mixOmics::plotIndiv(plsda(), legend = T, ind.names = input$indname_pls)
  }
})

output$plotvar_pls = renderPlot({
  req(plsda())
  validate(need(input$plsda_ncomp > 1, "PLS-DA plots require at least 2 components."))
  mixOmics::plotVar(plsda())
})

perf.plsda = eventReactive(input$gotunepls,{
  req(sumexpde_forplsda(), input$selcol_pls, input$selfolds)
  X = sumexpde_forplsda() %>% SummarizedExperiment::assay() %>% as.matrix() %>% t()
  Y = as.factor(SummarizedExperiment::colData(sumexpde_forplsda())[,input$selcol_pls])

  plsda.res = mixOmics::plsda(X, Y, ncomp = length(unique(Y)))
  mixOmics::perf(plsda.res, validation = "Mfold", folds = input$selfolds,
       progressBar = TRUE, auc = TRUE, nrepeat = input$selnrepeat)#, cpus = parallel::detectCores(logical = F)-1)
})


output$perfplsda = renderPlot({
  req(perf.plsda())
  plot(perf.plsda(), col = mixOmics::color.mixo(5:7), sd = TRUE, legend.position = "horizontal")
})

output$bestncomp_pls = renderPrint({
  req(perf.plsda())
  perf.plsda()$choice.ncomp
  
})









####sPLS-DA #####


sumexpde_forsplsda = reactive({
  req(sumexp_all())
  if(input$summ_splsda == TRUE && sumexp_all()$replicates == TRUE){
    sumexpdatamean()
  }else{
    sumexpdata()
  }
})

observeEvent(sumexpde_forsplsda(),{
  data = sumexpde_forsplsda() %>% SummarizedExperiment::colData() %>% as.data.frame()
  #selectin columns with more than one level
  data = data[, sapply(data, function(col) length(unique(col))) > 1]
  if("Product_Batch" %in% colnames(data)){
    sel = "Product_Batch"
  }else{
    sel= colnames(data)[1]
  }
  updateSelectInput(session, "selcol_spls", choices = colnames(data), selected = sel)
})


#number of allowed components (usually K-1)
output$splsda_ncompui = renderUI({
  req(sumexpde_forsplsda())
  Y = as.factor(colData(sumexpde_forsplsda())[,input$selcol_spls])
  max = length(unique(Y))
  numericInput("splsda_ncomp", paste0("Number of components ", "(from 2 to ", max, " )"), value = 2, min = 1, max = max, step = 1)
})


#spls-da
splsda = eventReactive(input$gospls,{
  req(sumexpde_forsplsda(), input$splsda_ncomp, input$keepx_splsda)
  X = sumexpde_forsplsda() %>% SummarizedExperiment::assay() %>% as.matrix() %>% t()
  Y = as.factor(SummarizedExperiment::colData(sumexpde_forsplsda())[,input$selcol_spls])
  keepx = gsub(" ", "", input$keepx_splsda)
  keepx = strsplit(keepx,",", fixed = FALSE)[[1]] 
  keepx = as.numeric(keepx)
  validate(need(anyNA(keepx) == FALSE, "Error in KeepX. Please verify KeepX contains only numbers separated by a comma."))
  if(nrow(X)<30){
    showNotification(tagList(icon("exclamation-circle"), HTML("&nbsp;There are less than 30 samples. sPLS-DA could not be statistically significant.")), type = "warning")
  }
  mixOmics::splsda(X, Y, ncomp = input$splsda_ncomp, keepX = keepx)
})


#plotIndiv
output$plotindiv_spls = renderPlot({
  req(splsda())
  if(input$spls_back_ellipse == "Background"){
    background <- mixOmics::background.predict(splsda(), comp.predicted = 2, dist = "max.dist")
    len = length(unique(splsda()$Y))
    if(len > 8){
      palette = grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, name  = "Dark2"))
      colors = palette(len)
    }else{
      colors = RColorBrewer::brewer.pal(len, name  = "Dark2")
      colors = colors[1:len]
    }
    
    mixOmics::plotIndiv(splsda(), comp = 1:2, ind.names = input$indname_spls, legend = TRUE, background = background,
                        col = colors)
  }else if(input$spls_back_ellipse == "Ellipse"){
    mixOmics::plotIndiv(splsda(), ellipse = T, legend = T, ind.names = input$indname_spls)
  }else{
    mixOmics::plotIndiv(splsda(), legend = T, ind.names = input$indname_spls)
  }
})



output$plotvar_spls = renderPlot({
  req(splsda())
  mixOmics::plotVar(splsda())
})




tune.splsda = eventReactive(input$gotunespls,{
  req(sumexpde_forsplsda(), input$selnrepeat_spls,input$selcol_spls)
  X = sumexpde_forsplsda() %>% SummarizedExperiment::assay() %>% as.matrix() %>% t()
  Y = as.factor(SummarizedExperiment::colData(sumexpde_forsplsda())[,input$selcol_spls])
  
  list.keepX <- c(5:10,  seq(20, 100, 10))
  mixOmics::tune.splsda(X, Y, ncomp = length(unique(Y)),
                        validation = 'Mfold', 
                        folds = input$selfolds_spls, 
                        dist = 'max.dist', progressBar = TRUE,
                        measure = "BER", test.keepX = list.keepX,
                        nrepeat = input$selnrepeat_spls)
})


output$plot_bestsparse = renderPlot({
  req(tune.splsda())
  n = dim(tune.splsda()$choice.ncomp$values)[2]
  plot(tune.splsda(), col = mixOmics::color.jet(n))
})


output$bestncomp_spls = renderPrint({
  req(tune.splsda())
  tune.splsda()$choice.ncomp$ncomp
})

output$bestkeepx_spls = renderPrint({
  req(tune.splsda())
  tune.splsda()$choice.keepX
})





}
