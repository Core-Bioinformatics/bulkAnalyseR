enrichmentPanelUI <- function(id){
  ns <- NS(id)
  
  tabPanel(
    'Enrichment',
    sidebarLayout(
      
      # Sidebar panel for inputs ----
      sidebarPanel(
        #file upload
        fileInput(ns('fileIn'), 'Choose CSV File', accept = '.csv'),
        
        checkboxGroupInput(ns('gprofilerSources'), 'Select data sources', 
                           choices = c('GO:BP', 'GO:MF', 'GO:CC', 'KEGG', 'REAC', 
                                       'TF', 'MIRNA', 'CORUM', 'HP', 'HPA', 'WP'), 
                           selected = c('GO:BP', 'GO:MF', 'GO:CC', 'KEGG', 'REAC', 'TF', 'MIRNA')),
        
        #download file name and button
        textInput(ns('fileName'), 'File name for download', value ='EnrichmentSet.csv'),
        downloadButton(ns('download'), 'Download')
      ),
      
      #main panel for displaying outputs
      mainPanel(
        plotOutput(ns('plot'), click = ns('plot_click')),
        tableOutput(ns('data')),
      )
    )
  )
}

enrichmentPanelServer <- function(id, getPlotData.DE){
  # check whether inputs (other than id) are reactive or not
  stopifnot({
    is.reactive(getPlotData.DE)
  })
  
  moduleServer(id, function(input, output, session){
    getenrichmentData <- reactive({
      req(input[['fileIn']])
      tryCatch(
        {
          data <- utils::read.csv(input[['fileIn']]$datapath)
        },
        error = function(e) {
          # return a safeError if a parsing error occurs
          stop(safeError(e))
        }
      )
      enrichment <- gprofiler2::gost(query = data$gene_id,
                                     organism = 'mmusculus',
                                     correction_method = 'fdr',
                                     custom_bg = getPlotData.DE()$gene_id,
                                     sources = input[['gprofilerSources']])
      enrichment$result[, c('query', 'significant', 'p_value', 'term_size',
                            'query_size', 'intersection_size', 'precision', 'recall',
                            'term_id', 'source', 'term_name', 'effective_domain_size')]
    })
    
    #plot enrichment data
    source <- NULL; p_value <- NULL
    output[['plot']] <- renderPlot({
      ggplot(getenrichmentData(), aes(x = source, y = p_value, colour = source)) + 
        geom_point() + 
        theme_bw()
    })
    
    #define clicking on enrichment data table
    output[['data']] <- renderTable({
      req(input[['plot_click']])
      nearPoints(df = getenrichmentData(), coordinfo = input[['plot_click']], maxpoints = 5)
    })
    
    #download enrichment
    output[['download']] <- downloadHandler(
      filename = function(){
        paste(input[['fileName']])
      },
      content = function(file){
        utils::write.csv(x = getenrichmentData(), file, row.names = FALSE)
      }
    )
  })
}

enrichmentPanelApp <- function(){
  shinyApp(
    ui = fluidPage(enrichmentPanelUI('RNA')),
    server = function(input, output, session){
      enrichmentPanelServer('RNA', reactive(gene.df[,'gene_id', drop = FALSE]))
    }
  )
}