#' ma_volcano_plot UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#' @param data_input data to be plotted. Use expdes_data_forplot().
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
#' @importFrom shinyWidgets awesomeRadio
#' @rawNamespace import(plotly, except = rename)
#' @importFrom shinycssloaders withSpinner
#' 
mod_ma_volcano_plot_ui <- function(id){
  ns <- NS(id)
  tagList(
    fluidPage(
      fluidRow(
      box(width = NULL, status = "primary", solidHeader = TRUE, 
          column(2,awesomeRadio(ns("expdes_selplot"), "Plot type:", choices = c("MA Plot", "Volcano Plot"))),
          column(2,selectInput(ns("expdes_fillplot"), "Select fill column", choices = c("DC", "Class"))),
          column(2,selectInput(ns("expdes_lipidplot"), "Label on lipids", choices = "", multiple = TRUE))
      )
    ),
    fluidRow(
      conditionalPanel(condition = "input.expdes_lipidplot == 'None' || input.expdes_lipidplot == ''", ns = ns,
                       shinycssloaders::withSpinner(plotly::plotlyOutput(ns("plot_expdes")))
      ),
      conditionalPanel(condition = "input.expdes_lipidplot != 'None' && input.expdes_lipidplot != ''", ns = ns,
                       shinycssloaders::withSpinner(plotOutput(ns("plot_expdes_labeled")))
      )
    )
    )
    
 
  )
}
    
#' ma_volcano_plot Server Functions
#'
#' @noRd 
#' 
#' @rawNamespace import(ggplot2, except = last_plot)
#' @rawNamespace import(plotly, except = rename)
#' @importFrom dplyr filter
#' @importFrom ggrepel geom_label_repel
#' 
mod_ma_volcano_plot_server <- function(id, data_input, contrast, logfc_val, thresh){
  moduleServer( id, function(input, output, session){
    ns <- session$ns
    
    observeEvent(data_input(),{
      updateSelectInput(session, "expdes_lipidplot", choices = data_input()$Lipids) #, selected = "None"
    })
    

    output$plot_expdes = plotly::renderPlotly({
      req(data_input())
      if(input$expdes_selplot == "MA Plot"){
        if(input$expdes_fillplot == "Class"){
          #MA plot con class
          ma_plot2 <- ggplot(data_input(), aes(x = AveExpr, y = logFC, col = as.factor(Class), label = Lipids)) +
            geom_point(aes(shape = as.factor(DE))) +
            geom_hline(yintercept = 0) +
            geom_hline(yintercept = logfc_val(), linetype = "dotted", size = 1) +  #qui ci stava lfc
            geom_hline(yintercept = -logfc_val(), linetype = "dotted", size = 1) +  #anche qui lfc
            ggtitle(paste("MA Plot:", gsub("vs"," vs ",contrast()))) + theme(legend.title = element_blank())
          plotly::ggplotly(ma_plot2) %>% plotly::layout(legend = list(title = list(text = "DC, Class")))
          
        }else{
          #MA plot con DE
          ma_plot <- ggplot(data_input(), aes(x = AveExpr, y = logFC, col = as.factor(DE), label = Lipids)) +
            geom_point() + scale_color_manual(name = "DC Lipids", labels = c("DW-reg", "Unchanged", "UP-reg"),
                                              values = c("-1" = "#00a65a","0" = "black","1" = "red")) +  #qui ci stava col_de
            geom_hline(yintercept = 0) +
            geom_hline(yintercept = logfc_val(), linetype = "dotted", size = 1) +  #qui ci stava lfc
            geom_hline(yintercept = -logfc_val(), linetype = "dotted", size = 1) +  #anche qui lfc
            ggtitle(paste("MA Plot:", gsub("vs"," vs ",contrast()))) + theme(legend.title = element_blank())
          plotly::ggplotly(ma_plot) %>% plotly::layout(legend = list(title = list(text = "DC")))
        }
      }else{
        if(input$expdes_fillplot == "Class"){
          #Volcano plot con class
          volcano_plot2 <- ggplot(data_input(), aes(x = logFC, y = -log10(adj.P.Val), col = as.factor(Class), label = Lipids)) +
            geom_point(aes(shape = as.factor(DE))) +
            geom_vline(xintercept = 0) +
            geom_hline(yintercept = -log10(as.numeric(thresh())), linetype = "dotted", size = 1) +
            geom_vline(xintercept = logfc_val(), linetype = "dotted", size = 1) +
            geom_vline(xintercept = -logfc_val(), linetype = "dotted", size = 1) +
            ggtitle(paste("Volcano Plot:", gsub("vs"," vs ",contrast()))) +
            labs(shape = "DC", colour = "Class") + theme(legend.title = element_blank())
          plotly::ggplotly(volcano_plot2) %>% plotly::layout(legend = list(title = list(text = "DC, Class")))
        }else{
          #Volcano plot con DE
          volcano_plot1 <- ggplot(data_input(), aes(x = logFC, y = -log10(adj.P.Val), col = as.factor(DE), label = Lipids)) +
            geom_point() +
            scale_color_manual(name = "DC Lipids", labels = c("DW-reg", "Unchanged", "UP-reg"), values = c("-1" = "#00a65a","0" = "black","1" = "red")) +
            geom_vline(xintercept = 0) + 
            geom_hline(yintercept = -log10(as.numeric(thresh())), linetype = "dotted", size = 1) +
            geom_vline(xintercept = logfc_val(), linetype = "dotted", size = 1) +
            geom_vline(xintercept = -logfc_val(), linetype = "dotted", size = 1) +
            ggtitle(paste("Volcano Plot:", gsub("vs"," vs ",contrast()))) + theme(legend.title = element_blank())
          plotly::ggplotly(volcano_plot1) %>% plotly::layout(legend = list(title = list(text = "DC")))
        }
      }
    })
    
    
    #plot con label
    output$plot_expdes_labeled = renderPlot({
      req(data_input())
      filtered = data_input() %>% dplyr::filter(Lipids %in% input$expdes_lipidplot)
      if(input$expdes_selplot == "MA Plot"){
        if(input$expdes_fillplot == "Class"){
          #MA plot con class
          ggplot(data_input(), aes(x = AveExpr, y = logFC, col = as.factor(Class), label = Lipids)) +
            geom_point(aes(shape = as.factor(DE))) +
            geom_hline(yintercept = 0) +
            geom_hline(yintercept = logfc_val(), linetype = "dotted", size = 1) +  #qui ci stava lfc
            geom_hline(yintercept = -logfc_val(), linetype = "dotted", size = 1) +  #anche qui lfc
            ggtitle(paste("MA Plot:", gsub("vs"," vs ",contrast()))) + ggrepel::geom_label_repel(aes(label = Lipids), data = filtered)
          
        }else{
          #MA plot con DE
          ggplot(data_input(), aes(x = AveExpr, y = logFC, col = as.factor(DE), label = Lipids)) +
            geom_point() + scale_color_manual(name = "Lipids", labels = c("DW-reg", "Unchanged", "UP-reg"),
                                              values = c("-1" = "#00a65a","0" = "black","1" = "red")) +  #qui ci stava col_de
            geom_hline(yintercept = 0) +
            geom_hline(yintercept = logfc_val(), linetype = "dotted", size = 1) +  #qui ci stava lfc
            geom_hline(yintercept = -logfc_val(), linetype = "dotted", size = 1) +  #anche qui lfc
            ggtitle(paste("MA Plot:", gsub("vs"," vs ",contrast()))) + ggrepel::geom_label_repel(aes(label = Lipids), data = filtered)
        }
      }else{
        if(input$expdes_fillplot == "Class"){
          #Volcano plot con class
          ggplot(data_input(), aes(x = logFC, y = -log10(adj.P.Val), col = as.factor(Class), label = Lipids)) +
            geom_point(aes(shape = as.factor(DE))) +
            geom_vline(xintercept = 0) +
            geom_hline(yintercept = -log10(as.numeric(thresh())), linetype = "dotted", size = 1) +
            geom_vline(xintercept = logfc_val(), linetype = "dotted", size = 1) +
            geom_vline(xintercept = -logfc_val(), linetype = "dotted", size = 1) +
            ggtitle(paste("Volcano Plot:", gsub("vs"," vs ",contrast()))) +
            labs(shape = "DC", colour = "Class") + ggrepel::geom_label_repel(aes(label = Lipids), data = filtered)
        }else{
          #Volcano plot con DE
          ggplot(data_input(), aes(x = logFC, y = -log10(adj.P.Val), col = as.factor(DE), label = Lipids)) +
            geom_point() +
            scale_color_manual(name = "Lipids", labels = c("DW-reg", "Unchanged", "UP-reg"), values = c("-1" = "#00a65a","0" = "black","1" = "red")) +
            geom_vline(xintercept = 0) + 
            geom_hline(yintercept = -log10(as.numeric(thresh())), linetype = "dotted", size = 1) +
            geom_vline(xintercept = logfc_val(), linetype = "dotted", size = 1) +
            geom_vline(xintercept = -logfc_val(), linetype = "dotted", size = 1) +
            ggtitle(paste("Volcano Plot:", gsub("vs"," vs ",contrast()))) + ggrepel::geom_label_repel(aes(label = Lipids), data = filtered)
        }
      }
    })
  })
}
    
## To be copied in the UI
# mod_ma_volcano_plot_ui("ma_volcano_plot_ui_1")
    
## To be copied in the server
# mod_ma_volcano_plot_server("ma_volcano_plot_ui_1")
