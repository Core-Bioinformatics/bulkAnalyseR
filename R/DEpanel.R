DEpanelUI <- function(id){
  ns <- NS(id)
  
  tabPanel(
    'Differential expression',
    sidebarLayout(
      
      # Sidebar panel for inputs ----
      sidebarPanel(
        
        # Input: Selector variables to compare
        selectInput(ns('variable1'), 'Sample 1:', 
                    c('control' = 'control', 'BCL11A' = 'BCL11A', 'CHD8' = 'CHD8')),
        selectInput(ns('variable2'), 'Sample 2:',
                    c('BCL11A' = 'BCL11A', 'control' = 'control', 'CHD8' = 'CHD8')),
        
        #DE thresholds
        sliderInput(ns('lfcThreshold'), label = 'logFC threshold',
                    min = 0, value = 1, max = 5, step = 0.5),
        
        sliderInput(ns('pvalThreshold'), label = 'Adjusted p-value threshold',
                    min = 0, value = 0.05, max = 0.2, step = 0.005),
        
        #Only start DE when button is pressed
        actionButton(ns('goDE'),label = 'Start DE'),
        
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


DEpanelServer <- function(id){
  # check whether inputs (other than id) are reactive or not
  stopifnot({TRUE})
  
  moduleServer(id, function(input, output, session){
    
    #Run DE
    getPlotData <- eventReactive(input$goDE,{
      #Add gene names
      RNAseqdata.normalised = merge(RNAseqdata.normalised, gene.df, by = c('gene_id'))
      rownames(RNAseqdata.normalised) = RNAseqdata.normalised$gene_id
      
      #Extract columns for samples being compared
      columns = list('control' = c('control_rep1','control_rep2','control_rep3'),
                     'BCL11A' = c('X11A_rep1','X11A_rep2','X11A_rep3'),
                     'CHD8' = c('CHD8_rep1','CHD8_rep2','CHD8_rep3'))
      RNAseqdata.normalised.category1 = RNAseqdata.normalised[, columns[[input[['variable1']]]]]
      RNAseqdata.normalised.category2 = RNAseqdata.normalised[, columns[[input[['variable2']]]]]
      RNAseqdata.normalised.categories = cbind(RNAseqdata.normalised.category1, RNAseqdata.normalised.category2)
      rownames(RNAseqdata.normalised.categories) = rownames(RNAseqdata.normalised)
      
      #run DE with supplied thresholds
      plotdata = data.frame(
        A = log2(rowSums(RNAseqdata.normalised.categories) / ncol(RNAseqdata.normalised.categories))
      )
      rownames(plotdata) = rownames(RNAseqdata.normalised)
      design <- stats::model.matrix(~ 0 + as.factor(c(rep(input[['variable1']], 3), 
                                                      rep(input[['variable2']], 3))))
      edger <- edgeR::DGEList(counts = RNAseqdata.normalised.categories, group = c(0, 0, 0, 1, 1, 1))
      edger <- edgeR::estimateDisp(edger, design)
      glm.fit = edgeR::glmFit(edger, design = design)
      glm.table <- edgeR::glmLRT(glm.fit, contrast=c(-1, 1))$table
      glm.table$adjustedp = stats::p.adjust(glm.table$PValue, method='BH')
      
      #create table for plotting
      plotdata$M = glm.table$logFC
      plotdata$pVal = glm.table$PValue
      plotdata$adjustedpVal = glm.table$adjustedp
      plotdata$gene_id = rownames(glm.table)
      plotdata$gene_name = RNAseqdata.normalised$gene_name
      plotdata$`-log10(adjustedpVal)` = -log10(plotdata$adjustedpVal)
      
      #the thresholds are returned here so that MA/volcano and table display 
      #don't use new thresholds without the button being used
      return(list('all'=plotdata,'logFC'=input[["lfcThreshold"]],'pVal'=input[["pvalThreshold"]]))
    })
    
    #Define output table when you click on gene with all genes or only DE
    output[['data']] <- renderDataTable({
      plotoutput = getPlotData()
      data = plotoutput$all
      data = data[((abs(data$M) >  plotoutput$logFC) & (data$adjustedpVal < plotoutput$pVal)),]
      data})

    
    #DE data download
    output[['downloadData']] <- downloadHandler(
      filename = function() {
        paste(input[['fileName']])
      },
      content = function(file) {
        utils::write.csv(x = getPlotData()$output[((abs(getPlotData()$all$M) > getPlotData()$logFC) & 
                                                     (getPlotData()$all$adjustedpVal < getPlotData()$pVal)),], 
                         file = file, 
                         row.names = FALSE)
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
