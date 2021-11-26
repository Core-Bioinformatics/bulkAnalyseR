enrichmentPanelUI <- function(id){
  ns <- NS(id)
  
  tabPanel(
    'Enrichment',
    shinyWidgets::dropdownButton(
      checkboxGroupInput(ns('gprofilerSources'), 'Select data sources', 
                         choices = c('GO:BP', 'GO:MF', 'GO:CC', 'KEGG', 'REAC', 
                                     'TF', 'MIRNA', 'CORUM', 'HP', 'HPA', 'WP'), 
                         selected = c('GO:BP', 'GO:MF', 'GO:CC', 'KEGG', 'REAC', 'TF', 'MIRNA')),
      actionButton(ns('goEnrichment'), label = 'Start enrichment analysis'),
      textInput(ns('fileName'), 'File name for data download', value ='EnrichmentSet.csv'),
      downloadButton(ns('download'), 'Download Data'),
      textInput(ns('plotFileName'), 'File name for plot download', value ='EnrichmentPlot.png'),
      downloadButton(ns('downloadPlot'), 'Download Plot'),
      status = "info",
      icon = icon("gear", verify_fa = FALSE), 
      tooltip = shinyWidgets::tooltipOptions(title = "Click to see inputs!")
    ),
    plotOutput(ns('plot'), click = ns('plot_click')),
    tableOutput(ns('data'))
  )
}

enrichmentPanelServer <- function(id, DEresults, organism){
  # check whether inputs (other than id) are reactive or not
  stopifnot({
    is.reactive(DEresults)
    !is.reactive(organism)
  })
  
  moduleServer(id, function(input, output, session){
    
    #Run enrichment
    getenrichmentData <- eventReactive({
      DEresults()
      input[["goEnrichment"]]
    }, 
    {
      inputdata = DEresults()
      enrichment <- gprofiler2::gost(query = inputdata$DEtableSubset$gene_id,
                                     organism = organism,
                                     correction_method = 'fdr',
                                     custom_bg = inputdata$DEtable$gene_id,
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
      x = jitter.build$data[[1]]$x
      df = getenrichmentData()
      df$jitter = x
      df$`-log10(pVal)`= -log10(df$p_value)
      return(df)
    })
    
    #Plot enrichment data
    source <- NULL; p_value <- NULL
    plotenrichmentPlot <- reactive({
      plotdata <- getenrichmentPlot()
      return(      ggplot(plotdata, aes(x = jitter, y = `-log10(pVal)`, colour = source)) + 
                     geom_point() + 
                     theme_bw()+ 
                     scale_x_continuous(breaks = seq(1, length(unique(plotdata$source)), 1), 
                                        labels = unique(plotdata$source)) + 
                     xlab("")
      )
    })
    output[['plot']]<-renderPlot({plotenrichmentPlot()})

    
    #Define clicking on enrichment data table
    output[['data']] <- renderTable({
      req(input[['plot_click']])
      nearPoints(df = getenrichmentPlot()[,c('term_name','source','term_id','-log10(pVal)','intersection_size','jitter')], coordinfo = input[['plot_click']], maxpoints = 5)
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
    
    output[['downloadPlot']] <- downloadHandler(
      filename = function() { input[['plotFileName']] },
      content = function(file) {
        device <- function(..., width, height) grDevices::png(..., width = width, height = height, res = 300, units = "in")
        ggsave(file, plot = plotenrichmentPlot(), device = device)
      }
    )
  })
}

# enrichmentPanelApp <- function(){
#   shinyApp(
#     ui = fluidPage(enrichmentPanelUI('RNA')),
#     server = function(input, output, session){
#       enrichmentPanelServer('RNA', reactive(gene.df[,'gene_id', drop = FALSE]))
#     }
#   )
# }
