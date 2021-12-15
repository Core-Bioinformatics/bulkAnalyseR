#' @import shiny
#' @import ggplot2
#' @importFrom rlang .data
bulkApp <- function(...){
  
  anno <- AnnotationDbi::select(
    getExportedValue('org.Mm.eg.db', 'org.Mm.eg.db'),
    keys = rownames(expression.matrix),
    keytype = 'ENSEMBL',
    columns = 'SYMBOL'
  ) %>%
    dplyr::distinct(ENSEMBL, .keep_all = TRUE) %>%
    dplyr::mutate(NAME = ifelse(is.na(SYMBOL), ENSEMBL, SYMBOL))
  
  ui <- navbarPage(
    "Interaction of BCL11A and CHD8 in triple negative breast cancer", 
    theme = shinythemes::shinytheme("flatly"),
    tabPanel("RNAseq",
             tabsetPanel(
               sampleSelectPanelUI("SampleSelect"),
               QCpanelUI("QC", metadata),
               DEpanelUI("DE", metadata),
               DEplotPanelUI("DEplot"),
               DEsummaryPanelUI("DEsummary", metadata),
               enrichmentPanelUI("Enrichment"),
               patternPanelUI("Patterns", metadata),
               crossPanelUI("Cross", metadata),
               GRNpanelUI("GRN"),
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
    # expression.matrix <- reactiveVal(expression.matrix)
    # metadata <- reactiveVal(metadata)
    filteredInputs <- sampleSelectPanelServer("SampleSelect", expression.matrix, metadata)
    expression.matrix <- reactive(filteredInputs()[["expression.matrix"]])
    metadata <- reactive(filteredInputs()[["metadata"]])
    QCpanelServer("QC", expression.matrix, metadata, anno)
    DEresults <- DEpanelServer("DE", expression.matrix, metadata, anno)
    DEplotPanelServer("DEplot", DEresults, anno)
    DEsummaryPanelServer("DEsummary", expression.matrix, metadata, DEresults, anno)
    enrichmentPanelServer("Enrichment", DEresults, organism = "mmusculus")
    patternPanelServer("Patterns", expression.matrix, metadata, anno)
    crossPanelServer("Cross", expression.matrix, metadata, anno)
    GRNpanelServer("GRN", expression.matrix, anno)
    peaksPanelServer("chip", ChIPseqdata)
    peaksPanelServer("atac", ATACseqdata)
  }
  
  shinyApp(ui, server, ...)
  
}
