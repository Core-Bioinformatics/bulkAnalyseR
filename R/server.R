server <- shiny::shinyServer(function(input, output) {
  
  #RNAseq------------------------------------------------------------------------
  #------------------------------------------------------------------------------
  
  #Differential expression-------------------------------------------------------
  DEpanelServer("DE")
  
  #Enrichment--------------------------------------------------------------------
  
  #run enrichment
  getenrichmentData <- reactive({
    
    req(input$file1)
    tryCatch(
      {
        data <- read.csv(input$file1$datapath)
      },
      error = function(e) {
        # return a safeError if a parsing error occurs
        stop(safeError(e))
      }
    )
    sources = c(ifelse(input$GOBP,'GO:BP','a'),
                ifelse(input$GOMF,'GO:MF','a'),
                ifelse(input$GOCC,'GO:CC','a'),
                ifelse(input$KEGG,'KEGG','a'),
                ifelse(input$REAC,'REAC','a'),
                ifelse(input$TF,'TF','a'),
                ifelse(input$MIRNA,'MIRNA','a'),
                ifelse(input$CORUM,'CORUM','a'),
                ifelse(input$HP,'HP','a'),
                ifelse(input$HPA,'HPA','a'),
                ifelse(input$WP,'WP','a'))
    sources = sources[sources!='a']
    enrichment=gprofiler2::gost(query=data$gene_id,organism='mmusculus',correction_method = 'fdr',custom_bg = getPlotData.DE()$gene_id,sources=sources)
    enrichment$result
  })
  
  #plot enrichment data
  output$enrichmentPlot <- renderPlot({
    enrichment = getenrichmentData()
    enrichment = enrichment[,c("query","significant","p_value","term_size","query_size","intersection_size","precision","recall","term_id","source","term_name","effective_domain_size")]
    ggplot(enrichment,aes(x=source,y=p_value,color=source))+geom_point()+theme_bw()
  })
  
  #define clicking on enrichment data table
  output$enrichmentTable <-renderTable({
    req(input$plot_click2)
    enrichment = getenrichmentData()
    enrichment = enrichment[,c("query","significant","p_value","term_size","query_size","intersection_size","precision","recall","term_id","source","term_name","effective_domain_size")]    
    nearPoints(df=enrichment,coordinfo = input$plot_click2,maxpoints = 5)  })
  
  #download enrichment
  output$downloadEnrichment <- downloadHandler(
    filename = function() {
      paste(input$enrichmentfileName)
    },
    content = function(file) {
      write.csv(x = getenrichmentData()[,c("query","significant","p_value","term_size","query_size","intersection_size","precision","recall","term_id","source","term_name","effective_domain_size")], file, row.names = FALSE)
    }
  )
  
  #ChIP--------------------------------------------------------------------------
  #------------------------------------------------------------------------------
  
  getPlotData_chip <- reactive({
    # get data for plotting
    # ChIPseqdata = readRDS('/servers/sutherland-scratch/ecw63/Teaching/shiny_WK/ChIP_nearbygenes.rds')
    
    l1 = log2(ChIPseqdata[, input[["variable1_chip"]]])
    l2 = log2(ChIPseqdata[, input[["variable2_chip"]]])
    plotdata = data.frame(M=l1-l2,A=(l1+l2)/2)
    rownames(plotdata)=rownames(ChIPseqdata)
    ChIPseqdata.normalised.category1 = ChIPseqdata[,input[["variable1_chip"]]]
    ChIPseqdata.normalised.category2 = ChIPseqdata[,input[["variable2_chip"]]]
    ChIPseqdata.normalised.categories = cbind(ChIPseqdata.normalised.category1,ChIPseqdata.normalised.category2)
    rownames(ChIPseqdata.normalised.categories)=rownames(ChIPseqdata)
    plotdata$upstream_genes = ChIPseqdata$upstream_gene_names
    plotdata$downstream_genes = ChIPseqdata$downstream_gene_names
    return(plotdata)
  })
  #define plot
  plot_chip <- reactive({
    plotdata = getPlotData_chip()
    plotdata.mygene.upstream = plotdata[(plotdata$upstream_genes!='') & input$GeneName_chip!=''&(as.vector(sapply(plotdata$upstream_genes,FUN=function(x)grepl(input$GeneName_chip,x),simplify = TRUE))==TRUE),]
    plotdata.mygene.downstream = plotdata[(plotdata$downstream_genes!='') & input$GeneName_chip!=''&(as.vector(sapply(plotdata$downstream_genes,FUN=function(x)grepl(input$GeneName_chip,x),simplify = TRUE))==TRUE),]
    max.M = max(abs(plotdata$M))
    myplot<-ggplot(plotdata,aes(x=A,y=M))+
      geom_point(color='black',alpha=0.1)+
      ylim(-max.M,max.M) +
      geom_point(data=plotdata.mygene.upstream,color='green',alpha=1,size=2)+
      geom_point(data=plotdata.mygene.downstream,color='red',alpha=1,size=2) +
      geom_hline(yintercept = input$FC_threshold_chip,color='blue') +
      geom_hline(yintercept = -input$FC_threshold_chip,color='blue')+theme_bw()
    return(myplot)
  })
  
  #output plot
  output$maPlot_chip <- renderPlot({
    plot_chip()
  })
  
  #output click data
  output$maData_chip <- renderTable({
    req(input$plot_click_chip)
    data = getPlotData_chip()
    nearPoints(df = data, coordinfo = input$plot_click_chip,threshold=20,maxpoints = 10)
  },digits=6)
  
  #ATAC--------------------------------------------------------------------------
  #------------------------------------------------------------------------------
  
  #get data for plotting
  getPlotData_atac <- reactive({
    # ATACseqdata = readRDS('/servers/sutherland-scratch/ecw63/Teaching/shiny_WK/ATAC_nearbygenes.rds')
    
    l1 = log2(ATACseqdata[, input[["variable1_atac"]]])
    l2 = log2(ATACseqdata[, input[["variable2_atac"]]])
    plotdata = data.frame(M=l1-l2,A=(l1+l2)/2)
    rownames(plotdata)=rownames(ATACseqdata)
    ATACseqdata.normalised.category1 = ATACseqdata[,input[["variable1_atac"]]]
    ATACseqdata.normalised.category2 = ATACseqdata[,input[["variable2_atac"]]]
    ATACseqdata.normalised.categories = cbind(ATACseqdata.normalised.category1,ATACseqdata.normalised.category2)
    rownames(ATACseqdata.normalised.categories)=rownames(ATACseqdata)
    plotdata$upstream_genes = ATACseqdata$upstream_gene_names
    plotdata$downstream_genes = ATACseqdata$downstream_gene_names
    return(plotdata)
  })

  #define plot
  plot_atac <- reactive({
    plotdata = getPlotData_atac()
    plotdata.mygene.upstream = plotdata[(plotdata$upstream_genes!='') & input$GeneName_atac!=''&(as.vector(sapply(plotdata$upstream_genes,FUN=function(x)grepl(input$GeneName_atac,x),simplify = TRUE))==TRUE),]
    plotdata.mygene.downstream = plotdata[(plotdata$downstream_genes!='') & input$GeneName_atac!=''&(as.vector(sapply(plotdata$downstream_genes,FUN=function(x)grepl(input$GeneName_atac,x),simplify = TRUE))==TRUE),]
    max.M = max(abs(plotdata$M))
    myplot<-ggplot(plotdata,aes(x=A,y=M))+
      geom_point(color='black',alpha=0.1)+
      ylim(-max.M,max.M) +
      geom_point(data=plotdata.mygene.upstream,color='green',alpha=1,size=2)+
      geom_point(data=plotdata.mygene.downstream,color='red',alpha=1,size=2) +
      geom_hline(yintercept = input$FC_threshold_atac,color='blue') +
      geom_hline(yintercept = -input$FC_threshold_atac,color='blue')+theme_bw()
    return(myplot)
  })
  
  #output plot
  output$maPlot_atac <- renderPlot({
    plot_atac()
  })
  
  #output click data
  output$maData_atac <- renderTable({
    req(input$plot_click_atac)
    data = getPlotData_atac()
    nearPoints(df = data, coordinfo = input$plot_click_atac,threshold=20,maxpoints = 10)
  },digits=6)
  
  
})