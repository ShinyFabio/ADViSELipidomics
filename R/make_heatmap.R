#' Make a Heatmap using ComplexHeatmap package.
#' 
#' @description This function render a Heatmap using ComplexHeatmap package.
#' 
#' @param data A summarizedexperiment object. Tested only on sumexp averaged
#' @param add_rowannot A column used for the left (row) annotation (e.g. "Product_Batch")
#' @param add_colannot A column used for the bottom (column) annotation (e.g. "Class")
#' @param scale_data Scaling data option. Data can be scaled (using scale()) by row ("row"), by column ("column") or not scaled ("none"). The result is the z-score.
#' @param row_dend Boolean. Render dendogram on rows (TRUE or FALSE).
#' @param col_dend Boolean. Render dendogram on columns (TRUE or FALSE).
#' @param row_nclust Numeric. Number of cluster in the row dendrogram.
#' @param col_nclust Numeric. Number of cluster in the column dendrogram.
#' @param dist_method Choose a distance method (e.g "euclidean"). See ?stats::dist().
#' @param clust_method Choose a clustering method (e.g "ward.D"). See ?stats::hclust().
#' @param col_lab A label for the column (i.e. "Lipids").
#' @param unit_legend An unit label for the legend. If data are scaled will be "Z-score" otherwise will be the input.
#' @param col_label_size Size of the column labels. By default col_label_size = 10.
#' @param padding padding for the heatmap in mm. By default padding = c(2,2,2,15).
#' @param filter Character vector. The filtering on column by value. For example: c(input$filtheatmapcol, input$filtheatmapval).
#' @param log_data Logical value. Set to TRUE to perform a log2 scaling.
#' @param order_data Logical value. Set to TRUE if you want to performin the ordering based on the add_rowannot value.
#'
#' @importFrom dplyr select arrange across left_join
#' @importFrom stats dist hclust as.dendrogram setNames
#' @importFrom dendextend color_branches
#' @importFrom ComplexHeatmap HeatmapAnnotation Heatmap draw ht_opt
#' @importFrom SummarizedExperiment colData assay rowData
#' @importFrom grid unit
#' @importFrom grDevices colorRampPalette rainbow
#' @importFrom RColorBrewer brewer.pal
#' @importFrom circlize colorRamp2
#' @importFrom tibble column_to_rownames
#'


make_heatmap = function(data,
                        #filter = "none",
                        add_rowannot = "Product_Batch",
                        add_colannot = "Class",
                        log_data = FALSE,
                        scale_data = "column",
                        order_data = TRUE,
                        row_dend = TRUE, 
                        row_nclust = 2, 
                        col_dend = TRUE, 
                        col_nclust = 3, 
                        dist_method = "euclidean", 
                        clust_method = "ward.D2", 
                        col_lab = "Lipids", 
                        unit_legend = "ug/ml",
                        col_label_size = 8, 
                        padding = c(2,2,2,15)){

coldata = data$coldata
assay = data$assay
# if(!"none" %in% filter){
#   coldata = coldata %>% dplyr::filter(!!sym(filter[1]) == filter[2])
#   assay = assay %>% dplyr::select(rownames(coldata))
#   # if(dim(assay)[2] <= 1){
#   #   if(shiny::isRunning()){
#   #     shinyWidgets::show_alert(title ="Warning!",type = "warning",
#   #                              text = "Your filter returns a matrix with only one column. Dendrograms and scaling on column are disabled.")
#   #   }
#   #   message("Your filter returns a matrix with only one column. Dendrograms and scaling on column are disabled.")
#   # }
# }
  
temp = assay %>% t()



if(order_data == TRUE){
  order = assay %>% t() %>% as.data.frame() %>% tibble::rownames_to_column("SampleID")
  cold  = coldata %>% dplyr::select(SampleID, add_rowannot)
  joined = dplyr::left_join(order, cold, by = "SampleID")
  if(add_rowannot == "SampleID"){
    temp = joined %>% dplyr::arrange( dplyr::across(add_rowannot)) %>%
      tibble::column_to_rownames("SampleID") %>% as.matrix()
  }else{
    temp = joined %>% dplyr::arrange( dplyr::across(add_rowannot)) %>% dplyr::select(-add_rowannot) %>%
      tibble::column_to_rownames("SampleID") %>% as.matrix()
  }

}



if(log_data == TRUE){
  temp = log2(temp)
}


#scale none, row, column
if(scale_data == "column"){
  temp = scale(temp) # scale and center columns
  #ora se ho una colonna con tutti 0, scale restituisce una colonna con tutti NaN. Qui sostituisco le colonne
  #con tutti NaN con tutti 0. QUesta parte non dovrebbe servire
  for(i in seq(1:length(temp[1,]))){
    if(mean(is.na(temp[,i])) == 1){
      temp[,i] = 0
    }
  }
  unit_legend = "Z-score"
} else if(scale_data == "row"){
  temp = t(scale(t(temp))) # scale and center rows
  unit_legend = "Z-score"
}

if(scale_data == "none"){
  legend_col = c("#ffffff", "#ff8080", "#ff0000")
}else{
  max_legend = max(abs(temp), na.rm = TRUE)
  legend_col = circlize::colorRamp2(c(-max_legend, 0, max_legend), c("blue", "white", "red"))
}


#dendrogram = none', 'row', 'column' or 'both' 
if(row_dend == TRUE){
  row_dend2 = tryCatch({
    temp %>% stats::dist(method = dist_method) %>% stats::hclust(method = clust_method) %>% stats::as.dendrogram()
  },error = function(err){
    str(err)
    if(err$message == "NA/NaN/Inf in foreign function call (arg 10)"){
      if(shiny::isRunning()){
        shinyWidgets::show_alert(title ="Error in row dendrogramm!",
                                 text = "There are too many NAs. Distance matrix cannot be calculated. 
                                 Try to remove the dendrogram on rows or impute NAs.", type = "error")
      }
      message("There are too many NAs. Distance matrix cannot be calculated. Try to remove the dendrogram on rows or impute NAs.")
    }else{
      if(shiny::isRunning()){
        shinyWidgets::show_alert(title ="Error in row dendrogramm!",
                                 text = "Something wrong in the clustering. Maybe too many NAs? Try to remove the dendrogram on rows or impute NAs.", type = "error")
      }
      message("Something wrong in the clustering. Maybe too many NAs? Try to remove the dendrogram on rows or impute NAs.")
      cat(crayon::bgYellow(crayon::red(err$message)))
    }
  })
  if(is.null(row_dend2)){return(NULL)}
  row_dend2 = temp %>% stats::dist(method = dist_method) %>% stats::hclust(method = clust_method) %>% stats::as.dendrogram()
  row_dend2 = dendextend::color_branches(row_dend2, k = row_nclust)
  row_split = row_nclust
} else {
  row_dend2 = FALSE
  row_split = NULL
  
}

if(col_dend == TRUE){
  col_dend2 = tryCatch({
    temp %>% t() %>% stats::dist(method = dist_method) %>% stats::hclust(method = clust_method) %>% stats::as.dendrogram()
  },error = function(err){
    if(err$message == "NA/NaN/Inf in foreign function call (arg 10)"){
      if(shiny::isRunning()){
        shinyWidgets::show_alert(title ="Error in column dendrogramm!",
                                 text = "There are too many NAs. Distance matrix cannot be calculated. 
                                 Try to remove the dendrogram on columns or impute NAs.", type = "error")
      }
      message("There are too many NAs. Distance matrix cannot be calculated. Try to remove the dendrogram on columns or impute NAs.")
    }else{
      if(shiny::isRunning()){
        shinyWidgets::show_alert(title ="Error in column dendrogramm!",
                                 text = "Something wrong in the clustering. ? Try to remove the dendrogram on columns or impute NAs.", type = "error")
      }
      message("Something wrong in the clustering. Maybe too many NAs? Try to remove the dendrogram on columns or impute NAs.")
      cat(crayon::bgYellow(crayon::red(err$message)))
    }
  })
  if(is.null(col_dend2)){return(NULL)}
  col_dend2 = dendextend::color_branches(col_dend2, k = col_nclust)
  col_split = col_nclust
} else {
  col_dend2 = FALSE
  col_split = NULL
}


getPalette = grDevices::colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))

#row annotation
annotdata_row = coldata %>% dplyr::select(add_rowannot) 
if(order_data == TRUE){
  annotdata_row = annotdata_row %>% dplyr::arrange(dplyr::across(add_rowannot))
}
leng_row = annotdata_row %>% table() %>% length()
colorannot_row = stats::setNames(grDevices::rainbow(n = leng_row), c(row.names(table(annotdata_row)))) #oppure getPalette
colorannot_row = stats::setNames(list(colorannot_row), paste(add_rowannot))
row_ha = ComplexHeatmap::HeatmapAnnotation(df = annotdata_row, which = "row", col = colorannot_row, border = TRUE)

#col annotation
annotdata_col = data$rowdata %>% as.data.frame() %>% dplyr::select(add_colannot)
leng_col = annotdata_col %>% table() %>% length()
colorannot_col = stats::setNames(getPalette(leng_col), c(row.names(table(annotdata_col))))
colorannot_col = stats::setNames(list(colorannot_col), paste(add_colannot))
col_ha = ComplexHeatmap::HeatmapAnnotation(df = annotdata_col, which = "column", col = colorannot_col, border = TRUE)

#add space between annotation and heatmap
ComplexHeatmap::ht_opt(ROW_ANNO_PADDING = grid::unit(4, "mm"), COLUMN_ANNO_PADDING = grid::unit(4, "mm"))


ht = ComplexHeatmap::Heatmap(temp, name = unit_legend, rect_gp = grid::gpar(col = "white", lwd = 1), 
                             row_title = "Sample", column_title = "Lipids", 
                             row_names_gp = grid::gpar(fontsize = 10), column_names_gp = grid::gpar(fontsize = col_label_size), #size testo
                             cluster_rows = row_dend2, cluster_columns = col_dend2, 
                             left_annotation = row_ha, bottom_annotation = col_ha,
                             column_split = col_split, row_split = row_split,
                             row_gap = grid::unit(2, "mm"), column_gap = grid::unit(2, "mm"), #space between divisions
                             col = legend_col
)
ht = ComplexHeatmap::draw(ht, padding = grid::unit(padding, "mm"))
}
