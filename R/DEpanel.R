DEpanelUI <- function(id){
  ns <- NS(id)
  
  tabPanel(
    'Differential expression',
    titlePanel('Interaction of BCL11A and CHD8 in triple negative breast cancer'),
    sidebarLayout(
      
      # Sidebar panel for inputs ----
      sidebarPanel(
        
        # Input: Selector variables to compare
        selectInput(ns('variable1'), 'Sample 1:', 
                    c('control' = 'control', 'BCL11A' = 'BCL11A', 'CHD8' = 'CHD8')),
        selectInput(ns('variable2'), 'Sample 2:',
                    c('BCL11A' = 'BCL11A', 'control' = 'control', 'CHD8' = 'CHD8')),
        
        #download file name and button
        textInput(ns('fileName'),'File name for download', value ='DEset.csv', placeholder = 'DEset.csv'),
        downloadButton(ns('downloadData'), 'Download'),
        
        #choose type of plot for RNAseq
        selectInput(ns('plotType'), 'Type of plot:', c('MA' = 'MA', 'volcano' = 'volcano')),
        
        #DE thresholds
        sliderInput(ns('lfcThreshold'), label = 'logFC threshold',
                    min = 0, value = 1, max = 5, step = 0.5),
        
        sliderInput(ns('pvalThreshold'), label = 'Adjusted p-value threshold',
                    min = 0, value = 0.05, max = 0.2, step = 0.005),
        
        #should all genes or only DE be shown in table?
        checkboxInput(ns('allGenes'), 'All genes', TRUE),
        
        #gene to colour
        textInput(ns('geneName'), 'Gene name', value ='', placeholder = 'Mrgprb2')
        
      ),
      
      #Main panel for displaying outputs
      mainPanel(
        plotOutput(ns('plot'), click = ns('plot_click')),
        tableOutput(ns('data'))
      )
    )
  )
}

DEpanelServer <- function(id){
  # check whether inputs (other than id) are reactive or not
  stopifnot({TRUE})
  
  moduleServer(id, function(input, output, session){
    getPlotData.DE <- reactive({
      
      # normalised counts
      # RNAseqdata.normalised <- readRDS('/servers/sutherland-scratch/ecw63/Teaching/shiny_WK/RNAseqdata.normalised.rds')
      # gene.df = readRDS('/servers/sutherland-scratch/ecw63/Teaching/shiny_WK/genedf.rds')
      
      #add gene names
      RNAseqdata.normalised = merge(RNAseqdata.normalised, gene.df, by = c('gene_id'))
      rownames(RNAseqdata.normalised) = RNAseqdata.normalised$gene_id
      
      #extract columns for samples being compared
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
      design <- model.matrix(~ 0 + as.factor(c(rep(input[['variable1']], 3), rep(input[['variable2']], 3))))
      edger <- DGEList(counts = RNAseqdata.normalised.categories,group = c(0, 0, 0, 1, 1, 1))
      edger <- estimateDisp(edger, design)
      glm.fit = glmFit(edger, design = design)
      glm.table <- glmLRT(glm.fit, contrast=c(-1, 1))$table
      glm.table$adjustedp = p.adjust(glm.table$PValue, method='BH')
      
      #create table for plotting
      plotdata$M = glm.table$logFC
      plotdata$pVal = glm.table$PValue
      plotdata$adjustedpVal = glm.table$adjustedp
      plotdata$gene_id = rownames(glm.table)
      plotdata$gene_name = RNAseqdata.normalised$gene_name
      plotdata$`-log10(adjustedpVal)` = -log10(plotdata$adjustedpVal)
      
      plotdata
    })
    
    #define plot (MA or volcano)
    plot_RNAseq <- reactive({
      plotdata = getPlotData.DE()
      plotdata.DE = plotdata[((abs(plotdata$M) > input[['lfcThreshold']]) & 
                                (plotdata$adjustedpVal < input[['pvalThreshold']])),]
      plotdata.logFC = plotdata[((abs(plotdata$M)>input[['lfcThreshold']])),]
      plotdata.pval = plotdata[(plotdata$adjustedpVal < input[['pvalThreshold']]),]
      plotdata.mygene = plotdata[plotdata$gene_name == input[['geneName']],]
      max.M = max(abs(plotdata$M))
      if (input[['plotType']] == 'MA'){
        #MA
        myplot <- ggplot(plotdata, aes(x = A, y = M)) +
          geom_point(color = 'black', alpha = 0.1) +
          ylim(-max.M, max.M) +
          geom_point(data = plotdata.DE, color = 'red', alpha=0.5) +
          geom_point(data = plotdata.mygene, color = 'green', alpha=1, size=2) +
          geom_text_repel(data = plotdata.mygene, label = input[['geneName']]) + 
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
          geom_text_repel(data = plotdata.mygene, label = input[['geneName']]) +
          theme_bw() +
          geom_hline(yintercept = -log10(input[['pvalThreshold']]), color = 'gray')+
          geom_vline(xintercept = input[['lfcThreshold']], color = 'gray')+
          geom_vline(xintercept = -input[['lfcThreshold']], color='gray')
      }
      myplot
    })
    
    #output MA/volcano plot
    output[['plot']] <- renderPlot(plot_RNAseq())
    
    #define output table when you click on gene with all genes or only DE
    output[['data']] <- renderTable({
      req(input[['plot_click']])
      data = getPlotData.DE()
      if (!(input[['allGenes']])){
        data = data[((abs(data$M) > input[['lfcThreshold']]) & (data$adjustedpVal < input[['pvalThreshold']])),]
      }
      nearPoints(df = data, coordinfo = input[['plot_click']], threshold = 20, maxpoints = 10)
    }, digits = 4)
    
    #DE data download
    output[['downloadData']] <- downloadHandler(
      filename = function() {
        paste(input[['fileName']])
      },
      content = function(file) {
        write.csv(x = getPlotData.DE()[((abs(getPlotData.DE()$M) > input[['lfcThreshold']]) & 
                                          (getPlotData.DE()$adjustedpVal < input[['pvalThreshold']])),], 
                  file = file, 
                  row.names = FALSE)
      }
    )
    
    getPlotData.DE
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