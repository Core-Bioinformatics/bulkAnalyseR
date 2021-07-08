server <- shiny::shinyServer(function(input, output) {
  
  #RNAseq------------------------------------------------------------------------
  #------------------------------------------------------------------------------
  
  #Differential expression-------------------------------------------------------
  
  getPlotData.DE <- reactive({
    
    # normalised counts
    # RNAseqdata.normalised <- readRDS('/servers/sutherland-scratch/ecw63/Teaching/shiny_WK/RNAseqdata.normalised.rds')
    # gene.df = readRDS('/servers/sutherland-scratch/ecw63/Teaching/shiny_WK/genedf.rds')
    
    #add gene names
    RNAseqdata.normalised = merge(RNAseqdata.normalised,gene.df,by = c('gene_id'))
    rownames(RNAseqdata.normalised)=RNAseqdata.normalised$gene_id

    #extract columns for samples being compared
    columns = list('control'=c('control_rep1','control_rep2','control_rep3'),'BCL11A'=c('X11A_rep1','X11A_rep2','X11A_rep3'),'CHD8'=c('CHD8_rep1','CHD8_rep2','CHD8_rep3'))
    RNAseqdata.normalised.category1 = RNAseqdata.normalised[,columns[[input[["variable1"]]]]]
    RNAseqdata.normalised.category2 = RNAseqdata.normalised[,columns[[input[["variable2"]]]]]
    RNAseqdata.normalised.categories = cbind(RNAseqdata.normalised.category1,RNAseqdata.normalised.category2)
    rownames(RNAseqdata.normalised.categories)=rownames(RNAseqdata.normalised)
    
    #run DE with supplied thresholds
    plotdata = data.frame(A=log2(rowSums(RNAseqdata.normalised.categories)/ncol(RNAseqdata.normalised.categories)))
    rownames(plotdata)=rownames(RNAseqdata.normalised)
    design <- model.matrix(~ 0 + as.factor(c(rep(input[["variable1"]],3),rep(input[["variable2"]],3))))
    edger <- DGEList(counts = RNAseqdata.normalised.categories,group = c(0,0,0,1,1,1))
    edger <- estimateDisp(edger, design)
    glm.fit = glmFit(edger,design = design)
    glm.table <- glmLRT(glm.fit, contrast=c(-1,1))$table
    glm.table$adjustedp = p.adjust(glm.table$PValue,method='BH')
    
    #create table for plotting
    plotdata$M = glm.table$logFC
    plotdata$pVal = glm.table$PValue
    plotdata$adjustedpVal = glm.table$adjustedp
    plotdata$gene_id = rownames(glm.table)
    plotdata$gene_name = RNAseqdata.normalised$gene_name
    plotdata$`-log10(adjustedpVal)` = -log10(plotdata$adjustedpVal)
    return(plotdata)
  })
  
  #define plot (MA or volcano)
  plot_RNAseq <- reactive({
    plotdata = getPlotData.DE()
    plotdata.DE = plotdata[((abs(plotdata$M)>input[['FC_threshold']]) & (plotdata$adjustedpVal < input[['pval_threshold']])),]
    plotdata.logFC = plotdata[((abs(plotdata$M)>input[['FC_threshold']])),]
    plotdata.pval = plotdata[(plotdata$adjustedpVal < input[['pval_threshold']]),]
    plotdata.mygene = plotdata[plotdata$gene_name==input$GeneName,]
    max.M = max(abs(plotdata$M))
    if (input$plottype == 'MA'){
      #MA
      myplot<-ggplot(plotdata,aes(x=A,y=M))+
        geom_point(color='black',alpha=0.1)+
        ylim(-max.M,max.M)+
        geom_point(data=plotdata.DE,color='red',alpha=0.5)+
        geom_point(data=plotdata.mygene,color='green',alpha=1,size=2)+
        geom_text_repel(data=plotdata.mygene,label=input$GeneName)+theme_bw()
      return(myplot)
    }
    else {   
      #volcano
      myplot<-ggplot(plotdata,aes(x=M,y=`-log10(adjustedpVal)`))+
        geom_point(color='black',alpha=0.1)+
        xlim(-max.M,max.M)+
        geom_point(data=plotdata.logFC,color='orange',alpha=0.5)+
        geom_point(data=plotdata.logFC,color='blue',alpha=0.5)+
        geom_point(data=plotdata.DE,color='red',alpha=1)+
        geom_point(data=plotdata.mygene,color='green',alpha=1,size=2)+
        geom_text_repel(data=plotdata.mygene,label=input$GeneName)+
        theme_bw()+
        geom_hline(yintercept=-log10(input[['pval_threshold']]),color='gray')+
        geom_vline(xintercept = input[['FC_threshold']],color='gray')+
        geom_vline(xintercept = -input[['FC_threshold']],color='gray')
      return(myplot)}
  })
  
  #output MA/volcano plot
  output$maPlot <- renderPlot({
    plot_RNAseq()
  })
  
  #define output table when you click on gene with all genes or only DE
  output$maData.DE <- renderTable({
    req(input$plot_click)
    data = getPlotData.DE()
    if (!(input$all_genes)){
      data = data[((abs(data$M)>input[['FC_threshold']]) & (data$adjustedpVal < input[['pval_threshold']])),]}
    nearPoints(df = data, coordinfo = input$plot_click,threshold=20,maxpoints = 10)
  },digits=4)
  
  #DE data download
  output$downloadData <- downloadHandler(
    filename = function() {
      paste(input$fileName)
    },
    content = function(file) {
      write.csv(x = getPlotData.DE()[((abs(getPlotData.DE()$M)>input[['FC_threshold']]) & (getPlotData.DE()$adjustedpVal < input[['pval_threshold']])),], file, row.names = FALSE)
    }
  )
  
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