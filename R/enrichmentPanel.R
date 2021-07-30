enrichmentPanelUI <- function(id){
  ns <- NS(id)
  
  tabPanel(
    'Enrichment',
    sidebarLayout(
      
      # Sidebar panel for inputs ----
      sidebarPanel(
        
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
    
    #Run enrichment
    getenrichmentData <- reactive({
      inputdata = getPlotData.DE()
      data = inputdata$de
      listofgenes = inputdata$genelist
      enrichment <- gprofiler2::gost(query = data$gene_id,
                                     organism = 'mmusculus',
                                     correction_method = 'fdr',
                                     custom_bg = listofgenes,
                                     sources = input[['gprofilerSources']])
      enrichment$result[, c('query', 'significant', 'p_value', 'term_size',
                            'query_size', 'intersection_size', 'precision', 'recall',
                            'term_id', 'source', 'term_name', 'effective_domain_size')]
    })
    
    #Jitter plot and save coordinates
    getenrichmentPlot <- reactive({
      jitter.plot = ggplot(getenrichmentData(), aes(x = source, y = p_value, colour = source)) + 
        geom_jitter()
      jitter.build <- ggplot_build(jitter.plot)
      x=jitter.build$data[[1]]$x
      df = getenrichmentData()
      df$jitter = x
      df$`-log10(pVal)`= -log10(df$p_value)
      return(df)
    })
    
    #Plot enrichment data
    source <- NULL; p_value <- NULL
    output[['plot']] <- renderPlot({
      ggplot(getenrichmentPlot(), aes(x = jitter, y = `-log10(pVal)`, colour = source)) + 
        geom_point() + 
        theme_bw()+ scale_x_continuous(breaks=seq(1,length(unique(getenrichmentPlot()$source)),1),labels = unique(getenrichmentPlot()$source))+xlab('')
    })
    
    #Define clicking on enrichment data table
    output[['data']] <- renderTable({
      req(input[['plot_click']])
      nearPoints(df = getenrichmentPlot(), coordinfo = input[['plot_click']], maxpoints = 5)
    })
    
    #Download enrichment
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
