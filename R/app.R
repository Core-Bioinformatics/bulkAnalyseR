#' @import shiny
#' @import edgeR
#' @import ggplot2
bulkApp <- function(...){
  
  ui <- navbarPage(
    "Interaction of BCL11A and CHD8 in triple negative breast cancer", 
    theme = shinythemes::shinytheme("flatly"),
    tabPanel("RNAseq",
             tabsetPanel(
               DEpanelUI("DE"),
               enrichmentPanelUI("Enrichment")
             )),
    peaksPanelUI("chip", "ChIPseq", c("control BCL11A IP" = "control_11AIP",
                                      "control CHD8 IP" = "control_CHD8IP",
                                      "BCL11A KD BCL11A IP" = "11AKD_11AIP",
                                      "BCL11A KD CHD8 IP" = "11AKD_CHD8IP",
                                      "CHD8 KD BCL11A IP" = "CHD8KD_11AIP",
                                      "CHD8 KD CHD8 IP" = "CHD8KD_CHD8AIP")),
    peaksPanelUI("atac", "ATACseq", c("control" = "control",
                                      "BCL11A" = "BCL11A",
                                      "CHD8" = "CHD8"))
  )
  
  server <- function(input, output, session){
    getPlotData.DE <- DEpanelServer("DE")
    enrichmentPanelServer("Enrichment", getPlotData.DE)
    peaksPanelServer("chip", ChIPseqdata)
    peaksPanelServer("atac", ATACseqdata)
  }
  
  shinyApp(ui, server, ...)
  
}