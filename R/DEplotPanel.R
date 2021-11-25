DEplotPanelUI <- function(id){
  ns <- NS(id)
  
  tabPanel(
    'Volcano and MA plots',
    shinyWidgets::dropdownButton(
      selectInput(ns('plotType'), 'Type of plot:', c('Volcano', 'MA')),
      shinyWidgets::switchInput(
        inputId = ns('autoLabel'),
        label = "Auto labels", 
        labelWidth = "80px",
        onLabel = 'On',
        offLabel = 'Off',
        value = FALSE,
        onStatus = FALSE
      ),
      shinyWidgets::switchInput(
        inputId = ns('allGenes'),
        label = "Showing on click:", 
        labelWidth = "80px",
        onLabel = 'All genes',
        offLabel = 'Only DE genes',
        value = FALSE,
        onStatus = FALSE
      ),
      selectInput(ns("geneName"), "Genes to highlight:", multiple = TRUE, choices = character(0)),
      status = "info",
      icon = icon("gear", verify_fa = FALSE), 
      tooltip = shinyWidgets::tooltipOptions(title = "Click to see inputs!")
    ),
    plotOutput(ns('plot'), click = ns('plot_click')),
    tableOutput(ns('data')) 
  )
}

DEplotPanelServer <- function(id, DEresults){
  
  # check whether inputs (other than id) are reactive or not
  stopifnot({
    is.reactive(DEresults)
  })
  
  moduleServer(id, function(input, output, session){
    
    #Set up server-side search for gene names
    # observe({
    #   input[["geneName"]]
    #   updateSelectizeInput(session, "geneName", choices = DEresults()$DEtable$gene_name, server = TRUE)
    # })
    
    #Define plot (MA or volcano)
    DEplot <- reactive({
      results = DEresults()
      
      if(input[['plotType']] == 'Volcano'){
        myplot <- volcano_plot(
          genes.de.results = results$DEtable,
          pval.threshold = results$pvalThreshold, 
          lfc.threshold = results$lfcThreshold,
          add.labels.auto = input[["autoLabel"]],
          n.labels.auto = c(10, 10),
          add.labels.custom = length(input[["geneName"]]) > 0,
          genes.to.label = input[["geneName"]]
        )
      }
      
      data = results$DEtable
      lfcThreshold = results$lfcThreshold
      pvalThreshold = results$pvalThreshold
      plotdata.DE = results$DEtableSubset
      
      plotdata.logFC = data[((abs(data$log2FC) > lfcThreshold)), ]
      plotdata.pval = data[(data$pvalAdj < pvalThreshold), ]
      plotdata.mygene = data[data$gene_name %in% input[['geneName']], ]
      max.M = max(abs(data$log2FC))
      
      if (input[['plotType']] == 'MA'){
        myplot <- ggplot(data, aes(x = avgExp, y = log2FC)) +
          geom_point(color = 'black', alpha = 0.1) +
          ylim(-max.M, max.M) +
          geom_point(data = plotdata.DE, color = 'red', alpha=0.5) +
          geom_point(data = plotdata.mygene, color = 'green', alpha = 1, size = 2) +
          ggrepel::geom_text_repel(data = plotdata.mygene, label = plotdata.mygene$gene_name) +
          theme_minimal()
      }
      
      myplot
    })
    
    #Output MA/volcano plot
    output[['plot']] <- renderPlot(DEplot())
    
    #Define output table when you click on gene with all genes or only DE
    output[['data']] <- renderTable({
      req(input[['plot_click']])
      results = DEresults()
      if (input[['allGenes']]){
        data <- results$DEtable
      }else{
        data <- results$DEtableSubset
      }
      data <- data %>% dplyr::mutate(`-log10pval` = -log10(pvalAdj))
      nearPoints(df = data, coordinfo = input[['plot_click']], threshold = 20, maxpoints = 10)
    }, digits = 4)
    
  })
}

DEplotPanelApp <- function(){
  shinyApp(
    ui = navbarPage("DE", tabPanel("", tabsetPanel(DEpanelUI('RNA'), MApanelUI('RNA')))),
    server = function(input, output, session){
      getPlotData.DE <- DEpanelServer('RNA')
      MApanelServer('RNA', getPlotData.DE)
    }
  )
}
