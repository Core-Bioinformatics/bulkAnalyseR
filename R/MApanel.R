MApanelUI <- function(id){
  ns <- NS(id)
  
  tabPanel(
    'MA and volcano plots',
    sidebarLayout(
      
      # Sidebar panel for inputs ----
      sidebarPanel(
        
        #Should all genes or only DE be shown in table?
        shinyWidgets::switchInput(
          inputId = ns('allGenes'),
          label = "Genes in table", 
          labelWidth = "80px",
          onLabel = 'All genes',
          offLabel = 'Only DE genes',
          value = FALSE,
          onStatus = FALSE
        ),
        
        #choose type of plot for RNAseq
        selectInput(ns('plotType'), 'Type of plot:', c('MA' = 'MA', 'volcano' = 'volcano')),
        
        #Gene(s) to colour and label
        selectInput(ns("geneName"), "Genes to highlight:", multiple = TRUE, choices = character(0))
        
      ),
      
      #Main panel for displaying MA/volcano plot and table
      mainPanel(
        plotOutput(ns('plot'), click = ns('plot_click')),
        tableOutput(ns('data'))
      )
    )
  )
}

MApanelServer <- function(id,getPlotData.DE){
  
  # check whether inputs (other than id) are reactive or not
  stopifnot({
    is.reactive(getPlotData.DE)
  })
  
  # check whether inputs (other than id) are reactive or not
  stopifnot({TRUE})
  
  moduleServer(id, function(input, output, session){
    
    #Set up server-side search for gene names
    updateSelectizeInput(session, "geneName", choices = gene.df[gene.df$gene_id%in%RNAseqdata.normalised$gene_id,]$gene_name, server = TRUE)
    
    #Define plot (MA or volcano)
    MAVolcanoPlot <- reactive({
      plotoutput = getPlotData.DE()
      plotdata = plotoutput$all
      lfcThreshold=plotoutput$logFC
      pVal=plotoutput$pVal
      plotdata.DE = plotdata[((abs(plotdata$M) > lfcThreshold) & 
                                (plotdata$adjustedpVal < pVal)),]
      plotdata.logFC = plotdata[((abs(plotdata$M)>lfcThreshold)),]
      plotdata.pval = plotdata[(plotdata$adjustedpVal < pVal),]
      plotdata.mygene = plotdata[plotdata$gene_name %in% input[['geneName']],]
      max.M = max(abs(plotdata$M))
      A <- NULL; M <- NULL; `-log10(adjustedpVal)` <- NULL
      #MA
      if (input[['plotType']] == 'MA'){
        #MA
        myplot <- ggplot(plotdata, aes(x = A, y = M)) +
          geom_point(color = 'black', alpha = 0.1) +
          ylim(-max.M, max.M) +
          geom_point(data = plotdata.DE, color = 'red', alpha=0.5) +
          geom_point(data = plotdata.mygene, color = 'green', alpha=1, size=2) +
          ggrepel::geom_text_repel(data = plotdata.mygene, label = plotdata.mygene$gene_name) + 
          theme_bw()
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
          theme_bw() +
          geom_hline(yintercept = -log10(lfcThreshold), color = 'gray')+
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
        data = data[((abs(data$M) >  plotoutput$logFC) & (data$adjustedpVal < plotoutput$pVal)),]
      }
      nearPoints(df = data, coordinfo = input[['plot_click']], threshold = 20, maxpoints = 10)
    }, digits = 4)
    
  })
}

MApanelApp <- function(){
  shinyApp(
    ui = fluidPage(MApanelUI('RNA')),
    server = function(input, output, session){
      MApanelServer('RNA',reactive(gene.df[,'gene_id', drop = FALSE]))
    }
  )
}
