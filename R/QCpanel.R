QCpanelUI <- function(id, metadata){
  ns <- NS(id)
  
  tabPanel(
    'Quality checks',
    tags$h1("Jaccard Similarity Index Heatmap"),
    shinyWidgets::dropdownButton(
      shinyjqui::orderInput(ns('jaccard.annotations'), label = "Show annotations",
                            items = colnames(metadata)[c(ncol(metadata), seq_len(ncol(metadata) - 1))][-2]),
      sliderInput(ns('jaccard.n.abundant'), label = '# of (most abundant) genes',
                  min = 50, value = 500, max = 5000, step = 50, ticks = FALSE),
      checkboxInput(ns("jaccard.show.values"), label = "Show JSI values", value = FALSE),
      
      status = "info",
      icon = icon("gear", verify_fa = FALSE), 
      tooltip = shinyWidgets::tooltipOptions(title = "Click to see inputs!")
    ),
    plotOutput(ns('jaccard')),
    
    tags$h1("Principal Component Analysis"),
    shinyWidgets::dropdownButton(
      radioButtons(ns('pca.annotation'), label = "Group by",
                   choices = colnames(metadata), selected = colnames(metadata)[ncol(metadata)]),
      sliderInput(ns('pca.n.abundant'), label = '# of (most abundant) genes',
                  min = 50, value = 500, max = 5000, step = 50, ticks = FALSE),
      checkboxInput(ns("pca.show.labels"), label = "Show sample labels", value = FALSE),
      
      status = "info",
      icon = icon("gear", verify_fa = FALSE), 
      tooltip = shinyWidgets::tooltipOptions(title = "Click to see inputs!")
    ),
    plotOutput(ns('pca'))
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
      meta <- lapply(metadata, function(x) factor(x, levels = unique(x))) %>% 
        as.data.frame() %>%
        dplyr::arrange(dplyr::across(input[['jaccard.annotations']]))
      
      myplot <- jaccard_heatmap(
        expression.matrix = expression.matrix[, meta[, 1]],
        metadata = meta,
        top.annotation.ids = match(input[['jaccard.annotations']], colnames(meta)),
        n.abundant = input[['jaccard.n.abundant']], 
        show.values = input[["jaccard.show.values"]]
      )
      myplot 
    })
    output[['jaccard']] <- renderPlot(jaccard.plot())
    
    pca.plot <- reactive({
      myplot <- plot_pca(
        expression.matrix = expression.matrix,
        metadata = metadata,
        annotation.id = match(input[['pca.annotation']], colnames(metadata)),
        n.abundant = input[['pca.n.abundant']],
        show.labels = input[['pca.show.labels']]
      )
      myplot
    })
    output[['pca']] <- renderPlot(pca.plot())
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