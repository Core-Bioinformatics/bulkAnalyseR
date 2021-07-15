QCpanelUI <- function(id, metadata){
  ns <- NS(id)
  
  tabPanel(
    'Quality checks',
    tags$h1("Jaccard Similarity Index Heatmap"),
    sidebarLayout(
      sidebarPanel(
        checkboxGroupInput(ns('jaccard.annotations'), label = "Show annotations",
                           choices = colnames(metadata), selected = colnames(metadata)[ncol(metadata)]),
        sliderInput(ns('jaccard.n.abundant'), label = '# of (most abundant) genes',
                    min = 50, value = 500, max = 5000, step = 50, ticks = FALSE),
        checkboxInput(ns("jaccard.show.values"), label = "Show JSI values", value = FALSE)
      ),
      mainPanel(plotOutput(ns('jaccard')))
    ),
    # tags$h1("Principal Component Analysis"),
    # sidebarLayout(
    #   sidebarPanel(),
    #   mainPanel(plotOutput(ns('pca')))
    # )
  )
}

QCpanelServer <- function(id, expression.matrix, metadata){
  # check whether inputs (other than id) are reactive or not
  stopifnot({
    !is.reactive(expression.matrix)
    !is.reactive(metadata)
  })
  
  moduleServer(id, function(input, output, session){
    jaccard.plot <- reactive({
      myplot <- jaccard_heatmap(
        expression.matrix = expression.matrix,
        metadata = metadata,
        n.abundant = input[['jaccard.n.abundant']], 
        top.annotation.df = NULL, 
        top.annotation.colours = NULL,
        show.values = input[["jaccard.show.values"]]
      )
      myplot 
    })
    output[['jaccard']] <- renderPlot(jaccard.plot())
    
    # pca.plot <- reactive()
    # output[['pca']] <- renderPlot(pca.plot)
  })
}

QCpanelApp <- function(){
  shinyApp(
    ui = fluidPage(QCpanelUI('qc', metadata)),
    server = function(input, output, session){
      QCpanelServer('qc', expression.matrix, metadata)
    }
  )
}