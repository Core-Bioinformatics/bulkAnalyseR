peaksPanelUI <- function(id, label, samples){
  ns <- NS(id)
  
  tabPanel(
    label,
    
    sidebarLayout(
      
      # Sidebar panel for inputs ----
      sidebarPanel(
        
        # Input: Selector for variable to plot against mpg ----
        selectInput(ns("variable1"), "Sample 1:", samples, selected = samples[1]),
        selectInput(ns("variable2"), "Sample 2:", samples, selected = samples[2]),
        
        #FC threshold line
        sliderInput(ns("lfcThreshold"), label = "logFC threshold",
                    min = 0, value = 0.5, max = 5,step = 0.5),
        textInput(ns("geneName"), "Gene name", value ='', placeholder = 'Mrgprb2')
        
      ),
      
      # Main panel for displaying outputs ----
      mainPanel(
        plotOutput(ns("plot"), click = ns("plot_click")),
        tableOutput(ns("data")),
      )
      
    )
  )
}

peaksPanelServer <- function(id, nearby.genes){
  # check whether inputs (other than id) are reactive or not
  stopifnot({
    !is.reactive(nearby.genes)
  })
  
  moduleServer(id, function(input, output, session){
    getPlotData <- reactive({
      
      l1 <- log2(nearby.genes[, input[["variable1"]]])
      l2 <- log2(nearby.genes[, input[["variable2"]]])
      
      plotdata <- data.frame(M = l1 - l2, A = (l1 + l2) / 2)
      rownames(plotdata) <- rownames(nearby.genes)
      
      nearby.genes.normalised.category1 <- nearby.genes[, input[["variable1"]]]
      nearby.genes.normalised.category2 <- nearby.genes[, input[["variable2"]]]
      nearby.genes.normalised.categories <- cbind(nearby.genes.normalised.category1, 
                                                  nearby.genes.normalised.category2)
      rownames(nearby.genes.normalised.categories) <- rownames(nearby.genes)
      plotdata$upstream_genes <- nearby.genes$upstream_gene_names
      plotdata$downstream_genes <- nearby.genes$downstream_gene_names
      
      plotdata
    })
    #define plot
    myplot <- reactive({
      plotdata <- getPlotData()
      plotdata.mygene.upstream <- 
        plotdata[(plotdata$upstream_genes != '') & 
                   input[['geneName']] != '' &
                   (as.vector(sapply(plotdata$upstream_genes, FUN = function(x){
                     grepl(input[['geneName']], x)
                     })) == TRUE), ]
      plotdata.mygene.downstream <- 
        plotdata[(plotdata$downstream_genes != '') & 
                   input[['geneName']] != '' &
                   (as.vector(sapply(plotdata$downstream_genes, FUN = function(x){
                     grepl(input[['geneName']], x)
                     })) == TRUE), ]
      max.M = max(abs(plotdata$M))
      
      A <- NULL; M <- NULL
      myplot <- ggplot(plotdata, aes(x = A, y = M)) +
        geom_point(color = 'black', alpha = 0.1) +
        ylim(-max.M, max.M) +
        geom_point(data = plotdata.mygene.upstream, color = 'green', alpha = 1, size = 2) +
        geom_point(data = plotdata.mygene.downstream, color = 'red', alpha = 1, size = 2) +
        geom_hline(yintercept = input[['lfcThreshold']], color = 'blue') +
        geom_hline(yintercept = -input[['lfcThreshold']], color = 'blue') +
        theme_bw()
      
      myplot
    })
    
    #output plot
    output[['plot']] <- renderPlot(myplot())
    
    #output click data
    output[['data']] <- renderTable({
      req(input[['plot_click']])
      data = getPlotData()
      nearPoints(df = data, coordinfo = input[['plot_click']], threshold = 20, maxpoints = 10)
    }, digits = 6)
  })
}

peaksPanelApp <- function(){
  shinyApp(
    ui = fluidPage(peaksPanelUI('chip', "ChIPseq")),
    server = function(input, output, session){
      peaksPanelServer('chip', ChIPseqdata)
    }
  )
}