MApanelUI <- function(id){
  ns <- NS(id)
  
  tabPanel(
    'Volcano and MA plots',
    shinyWidgets::dropdownButton(
      shinyWidgets::switchInput(
        inputId = ns('allGenes'),
        label = "Showing on click:", 
        labelWidth = "80px",
        onLabel = 'All genes',
        offLabel = 'Only DE genes',
        value = FALSE,
        onStatus = FALSE
      ),
      selectInput(ns('plotType'), 'Type of plot:', c('Volcano', 'MA')),
      selectInput(ns("geneName"), "Genes to highlight:", multiple = TRUE, choices = character(0)),
      status = "info",
      icon = icon("gear"), 
      tooltip = shinyWidgets::tooltipOptions(title = "Click to see inputs!")
    ),
    plotOutput(ns('plot'), click = ns('plot_click')),
    tableOutput(ns('data')) 
  )
}

MApanelServer <- function(id, getPlotData.DE){
  
  # check whether inputs (other than id) are reactive or not
  stopifnot({
    is.reactive(getPlotData.DE)
  })
  
  moduleServer(id, function(input, output, session){
    
    #Set up server-side search for gene names
    #updateSelectizeInput(session, "geneName", choices = getPlotData.DE()$all$gene_name, server = TRUE)
    
    #Define plot (MA or volcano)
    MAVolcanoPlot <- reactive({
      plotoutput = getPlotData.DE()
      plotdata = plotoutput$all
      lfcThreshold = plotoutput$logFC
      pValThreshold = plotoutput$pVal
      plotdata.DE = plotoutput$de
      plotdata.logFC = plotdata[((abs(plotdata$M) > lfcThreshold)), ]
      plotdata.pval = plotdata[(plotdata$adjustedpVal < pValThreshold), ]
      plotdata.mygene = plotdata[plotdata$gene_name %in% input[['geneName']], ]
      max.M = max(abs(plotdata$M))
      A <- NULL; M <- NULL; `-log10(adjustedpVal)` <- NULL
      #MA
      if (input[['plotType']] == 'MA'){
        #MA
        myplot <- ggplot(plotdata, aes(x = A, y = M)) +
          geom_point(color = 'black', alpha = 0.1) +
          ylim(-max.M, max.M) +
          geom_point(data = plotdata.DE, color = 'red', alpha=0.5) +
          geom_point(data = plotdata.mygene, color = 'green', alpha = 1, size = 2) +
          ggrepel::geom_text_repel(data = plotdata.mygene, label = plotdata.mygene$gene_name) + 
          theme_minimal()
      }else{   
        #volcano
        myplot <- ggplot(plotdata, aes(x = M, y = `-log10(adjustedpVal)`)) +
          geom_point(color = 'black', alpha = 0.1) +
          xlim(-max.M, max.M) +
          geom_point(data = plotdata.logFC, color = 'orange', alpha = 0.5) +
          geom_point(data = plotdata.logFC, color = 'blue', alpha = 0.5) +
          geom_point(data = plotdata.DE, color = 'red', alpha = 1) +
          geom_point(data = plotdata.mygene,color = 'green', alpha = 1, size = 2) +
          ggrepel::geom_text_repel(data = plotdata.mygene, label = input[['geneName']]) +
          theme_minimal() +
          geom_hline(yintercept = -log10(pValThreshold), color = 'gray')+
          geom_vline(xintercept = lfcThreshold, color = 'gray')+
          geom_vline(xintercept = -lfcThreshold, color='gray')
      }
      myplot
    })
    
    #Output MA/volcano plot
    output[['plot']] <- renderPlot(MAVolcanoPlot())
    
    #Define output table when you click on gene with all genes or only DE
    output[['data']] <- renderTable({
      req(input[['plot_click']])
      plotoutput = getPlotData.DE()
      data = plotoutput$all
      if (!(input[['allGenes']])){
        data = data[((abs(data$M) >  plotoutput$logFC) & (data$adjustedpVal < plotoutput$pVal)), ]
      }
      nearPoints(df = data, coordinfo = input[['plot_click']], threshold = 20, maxpoints = 10)
    }, digits = 4)
    
  })
}

MApanelApp <- function(){
  shinyApp(
    ui = navbarPage("DE", tabPanel("", tabsetPanel(DEpanelUI('RNA'), MApanelUI('RNA')))),
    server = function(input, output, session){
      getPlotData.DE <- DEpanelServer('RNA')
      MApanelServer('RNA', getPlotData.DE)
    }
  )
}
