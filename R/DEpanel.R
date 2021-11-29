DEpanelUI <- function(id, metadata){
  ns <- NS(id)
  
  tabPanel(
    'Differential expression',
    sidebarLayout(
      
      # Sidebar panel for inputs ----
      sidebarPanel(
        
        # Input: Selector variables to compare
        selectInput(ns('variable1'), 'Condition 1:', unique(metadata[[ncol(metadata)]])),
        selectInput(ns('variable2'), 'Condition 2:', unique(metadata[[ncol(metadata)]]),
                    selected = unique(metadata[[ncol(metadata)]])[2]),
        
        #DE thresholds
        sliderInput(ns('lfcThreshold'), label = 'logFC threshold',
                    min = 0, value = 1, max = 5, step = 0.5),
        
        sliderInput(ns('pvalThreshold'), label = 'Adjusted p-value threshold',
                    min = 0, value = 0.05, max = 0.2, step = 0.005),
        
        #Only start DE when button is pressed
        actionButton(ns('goDE'), label = 'Start DE'),
        
        #download file name and button
        textInput(ns('fileName'),'File name for download', value ='DEset.csv', placeholder = 'DEset.csv'),
        downloadButton(ns('download'), 'Download Table'),
      ),
      
      #Main panel for displaying table of DE genes
      mainPanel(
        dataTableOutput(ns('data'))
      )
    )
  )
}


DEpanelServer <- function(id, expression.matrix, metadata, anno){
  # check whether inputs (other than id) are reactive or not
  stopifnot({
    !is.reactive(expression.matrix)
    !is.reactive(metadata)
    !is.reactive(anno)
  })
  
  moduleServer(id, function(input, output, session){
    
    DEresults <- eventReactive(input[["goDE"]], {
      condition.indices <- metadata[[ncol(metadata)]] %in% c(input[['variable1']], input[['variable2']])
      DEtable <- DEanalysis_edger(
        expression.matrix = expression.matrix[, condition.indices],
        condition = metadata[[ncol(metadata)]][condition.indices],
        var1 = input[['variable1']],
        var2 = input[['variable2']],
        anno = anno
      )
      
      DEtableSubset <- DEtable %>%
        dplyr::filter(abs(log2FC) > input[["lfcThreshold"]] & pvalAdj < input[["pvalThreshold"]])
      
      #the thresholds are returned here so that MA/volcano and table display 
      #don't use new thresholds without the button being used
      return(list('DEtable' = DEtable,
                  "DEtableSubset" = DEtableSubset,
                  'lfcThreshold' = input[["lfcThreshold"]], 
                  'pvalThreshold' = input[["pvalThreshold"]]))
    })
    
    #Define output table (only DE genes)
    output[['data']] <- renderDataTable(DEresults()$DEtableSubset)

    #DE data download
    output[['download']] <- downloadHandler(
      filename = function() {
        paste(input[['fileName']])
      },
      content = function(file) {
        utils::write.csv(x = DEresults()$DEtableSubset, file = file, row.names = FALSE)
      }
    )
    
    DEresults
    
  })
}

DEpanelApp <- function(){
  shinyApp(
    ui = fluidPage(DEpanelUI('RNA')),
    server = function(input, output, session){
      DEpanelServer('RNA')
    }
  )
}
