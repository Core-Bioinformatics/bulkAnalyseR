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
      if (input[['plotType']] == 'MA'){
        myplot <- ma_plot(
          genes.de.results = results$DEtable,
          pval.threshold = results$pvalThreshold, 
          lfc.threshold = results$lfcThreshold,
          add.labels.auto = input[["autoLabel"]],
          n.labels.auto = c(10, 10),
          add.labels.custom = length(input[["geneName"]]) > 0,
          genes.to.label = input[["geneName"]]
        )
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
    ui = navbarPage("DE", tabPanel("", tabsetPanel(DEpanelUI('RNA'), DEplotPanelUI('RNA')))),
    server = function(input, output, session){
      DEresults <- DEpanelServer('RNA')
      DEplotPanelServer('RNA', DEresults)
    }
  )
}
