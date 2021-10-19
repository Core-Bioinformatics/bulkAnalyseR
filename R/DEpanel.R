DEpanelUI <- function(id, metadata){
  ns <- NS(id)
  
  tabPanel(
    'Differential expression',
    sidebarLayout(
      
      # Sidebar panel for inputs ----
      sidebarPanel(
        
        # Input: Selector variables to compare
        selectInput(ns('variable1'), 'Sample 1:', unique(metadata$condition)),
        selectInput(ns('variable2'), 'Sample 2:', unique(metadata$condition),
                    selected = unique(metadata$condition)[2]),
        
        #DE thresholds
        sliderInput(ns('lfcThreshold'), label = 'logFC threshold',
                    min = 0, value = 1, max = 5, step = 0.5),
        
        sliderInput(ns('pvalThreshold'), label = 'Adjusted p-value threshold',
                    min = 0, value = 0.05, max = 0.2, step = 0.005),
        
        #Only start DE when button is pressed
        actionButton(ns('goDE'), label = 'Start DE'),
        
        #download file name and button
        textInput(ns('fileName'),'File name for download', value ='DEset.csv', placeholder = 'DEset.csv'),
        downloadButton(ns('downloadData'), 'Download'),
      ),
      
      #Main panel for displaying table of DE genes
      mainPanel(
        dataTableOutput(ns('data'))
      )
    )
  )
}


DEpanelServer <- function(id, expression.matrix, metadata){
  # check whether inputs (other than id) are reactive or not
  stopifnot({
    !is.reactive(expression.matrix)
    !is.reactive(metadata)
  })
  
  moduleServer(id, function(input, output, session){
    
    #Run DE
    getPlotData <- eventReactive(input[["goDE"]], {
      
      #Extract columns for samples being compared
      expression.matrix.subset <- 
        expression.matrix[, metadata$condition %in% c(input[['variable1']], input[['variable2']])]
      
      #run DE with supplied thresholds
      reps1 <- sum(metadata$condition == input[['variable1']])
      reps2 <- sum(metadata$condition == input[['variable2']])
      design <- stats::model.matrix(~ 0 + as.factor(c(rep(input[['variable1']], reps1), 
                                                      rep(input[['variable2']], reps2))))
      edger <- edgeR::DGEList(counts = expression.matrix.subset, 
                              group = c(rep(0, reps1), rep(1, reps2)))
      edger <- edgeR::estimateDisp(edger, design)
      glm.fit = edgeR::glmFit(edger, design = design)
      glm.table <- edgeR::glmLRT(glm.fit, contrast = c(-1, 1))$table
      glm.table$adjustedp = stats::p.adjust(glm.table$PValue, method = 'BH')
      
      #create table for plotting
      plotdata = tibble::tibble(
        gene_id = rownames(expression.matrix.subset),
        gene_name = rownames(expression.matrix.subset),
        A = log2(rowSums(expression.matrix.subset) / ncol(expression.matrix.subset)),
        M = glm.table$logFC,
        pVal = glm.table$PValue,
        adjustedpVal = glm.table$adjustedp,
        `-log10(adjustedpVal)` = -log10(adjustedpVal)
      )
      
      #the thresholds are returned here so that MA/volcano and table display 
      #don't use new thresholds without the button being used
      return(list('all' = plotdata, 'logFC' = input[["lfcThreshold"]], 'pVal' = input[["pvalThreshold"]]))
    })
    
    #Define output table when you click on gene with all genes or only DE
    output[['data']] <- renderDataTable({
      plotoutput = getPlotData()
      data = plotoutput$all
      data = data[((abs(data$M) >  plotoutput$logFC) & (data$adjustedpVal < plotoutput$pVal)),]
      data
    })
    
    
    #DE data download
    output[['downloadData']] <- downloadHandler(
      filename = function() {
        paste(input[['fileName']])
      },
      content = function(file) {
        plotoutput = getPlotData()
        data <- plotoutput$all %>%
          dplyr::filter(abs(M) > plotoutput$logFC, adjustedpVal < plotoutput$pVal)
        utils::write.csv(x = data, file = file, row.names = FALSE)
      }
    )
    #this is passed on to MA/volcano and enrichment modules
    outputdf = reactive({list(all=getPlotData()$all,
                              de=getPlotData()$all[((abs(getPlotData()$all$M) > getPlotData()$logFC) & 
                                                      (getPlotData()$all$adjustedpVal < getPlotData()$pVal)),],
                              logFC=getPlotData()$logFC,
                              pVal=getPlotData()$pVal,
                              genelist = as.vector(getPlotData()$all$gene_id))})
    outputdf
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
