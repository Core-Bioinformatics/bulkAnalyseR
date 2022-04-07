## ----options, include = FALSE-------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "##>"
)

## ----bulkAnalyseR, include=FALSE----------------------------------------------
library("bulkAnalyseR")

## ----load Yang 2019 mRNAseq---------------------------------------------------
exp.yang <- read.csv(
  system.file("extdata", "expression_matrix_preprocessed.csv", package = "bulkAnalyseR"),
  row.names = 1) %>% as.matrix
head(exp.yang)

## ----metadata-----------------------------------------------------------------
meta <- data.frame(
  srr = colnames(exp.yang), 
  timepoint = rep(c("0h", "12h", "36h"), each = 2)
)

## ----convert type, include = FALSE, eval = FALSE------------------------------
#  meta$srr = as.character(meta$srr)
#  meta$timepoint = as.character(meta$timepoint)

## ----str gene expression------------------------------------------------------
str(exp.yang)
str(meta)

## ----gene table coord---------------------------------------------------------
gene.coord.table <- read.csv(
  url('https://raw.githubusercontent.com/Core-Bioinformatics/bulkAnalyseR/master/exampledata/Yang2019_ChIP/gene_coords_GRCm38.p6.csv'),
  row.names = 1)
str(gene.coord.table)

## ----chip coord table---------------------------------------------------------
chip.coord.table <- read.csv(
  url('https://raw.githubusercontent.com/Core-Bioinformatics/bulkAnalyseR/master/exampledata/Yang2019_ChIP/ChIP_peak_coords.csv'),
  row.names = 1)
str(chip.coord.table)

## -----------------------------------------------------------------------------
cis.integration <- tibble::tibble(
  reference.expression.matrix = 'exp.yang',
  reference.org.db = 'org.Mm.eg.db',
  reference.coord = 'gene.coord.table',
  comparison.coord = 'chip.coord.table',
  reference.table.name = 'mRNAseq',
  comparison.table.name = 'ChIPseq'
)

## ----generate cis app, eval=FALSE---------------------------------------------
#  generateShinyApp(
#    expression.matrix = exp.yang,
#    metadata = meta,
#    modality = "RNA",
#    shiny.dir = "shiny_Yang2019_CisIntegration",
#    app.title = "Shiny app for visualisation of three timepoints from the Yang 2019 data",
#    organism = "mmusculus",
#    org.db = "org.Mm.eg.db",
#    cis.integration = cis.integration
#  )
#  shiny::runApp('shiny_Yang2019_CisIntegration')

## ----CisIntegration, echo = FALSE, out.width = "80%"--------------------------
knitr::include_graphics("figures/CisIntegrationExample.png") 

## ----yang multiple modality, eval= FALSE--------------------------------------
#  exp.chip <- read.csv(
#    url('https://raw.githubusercontent.com/Core-Bioinformatics/bulkAnalyseR/master/exampledata/Yang2019_ChIP/ChIP_expression_matrix_preprocessed.csv'),
#    row.names = 1) %>% as.matrix
#  meta.chip = data.frame(
#    id = colnames(exp.chip),
#    timepoint = c('0h','12h','36h')
#  )
#  cis.integration.2 <- tibble::tibble(
#    reference.expression.matrix = c('exp.yang','exp.chip'),
#    reference.org.db = c('org.Mm.eg.db','NULL'),
#    reference.coord = c('gene.coord.table','chip.coord.table'),
#    comparison.coord = c('chip.coord.table','gene.coord.table'),
#    reference.table.name = c('mRNAseq','ChIPseq'),
#    comparison.table.name = c('ChIPseq','mRNAseq')
#  )
#  generateShinyApp(
#    expression.matrix = list(exp.yang,exp.chip),
#    metadata = list(meta,meta.chip),
#    modality = c('RNA','ChIP'),
#    shiny.dir = "shiny_Yang2019_CisIntegration2",
#    app.title = "Shiny app for visualisation of three timepoints from the Yang 2019 data",
#    organism = list("mmusculus",NA),
#    org.db = list("org.Mm.eg.db",NA),
#    cis.integration = cis.integration.2
#  )
#  shiny::runApp('shiny_Yang2019_CisIntegration2')

## ----li data------------------------------------------------------------------
exp.mirna <- read.csv(
  url('https://raw.githubusercontent.com/Core-Bioinformatics/bulkAnalyseR/master/exampledata/Li2021_miRNA_mRNA/expression_matrix_miRNA_preprocessed.csv'),
  row.names = 1) %>% as.matrix
str(exp.mirna)
exp.mrna <- read.csv(
  url('https://raw.githubusercontent.com/Core-Bioinformatics/bulkAnalyseR/master/exampledata/Li2021_miRNA_mRNA/expression_matrix_mRNA_preprocessed.csv'),
  row.names = 1) %>% as.matrix
str(exp.mrna)
meta.trans = data.frame(id = paste0(rep(c('control_','IDD_'),each = 3),1:3),
                        rep = rep(1:3,2),
                        type = rep(c('control','IDD'),each = 3))
meta.trans

## ----li trans app, eval=FALSE-------------------------------------------------
#  generateShinyApp(
#    shiny.dir = 'shiny_Li_2021',
#    app.title = 'Li 2021 Trans Regulatory Example',
#    modality=list('mRNA','miRNA'),
#    metadata = meta.trans,
#    expression.matrix = list(exp.mrna,exp.mirna),
#    org.db = list('org.Hs.eg.db',NA),
#    organism=list('hsapiens',NA),
#    trans.integration = tibble::tibble(
#      reference.expression.matrix='exp.mrna',
#      reference.org.db='org.Hs.eg.db',
#      comparison.expression.matrix='exp.mirna',
#      comparison.org.db='NULL',
#      reference.table.name='mRNA',
#      comparison.table.name='miRNA'
#    )
#  )
#  shiny::runApp('shiny_Li_2021')

## ----TransIntegration, echo = FALSE, out.width = "80%"------------------------
knitr::include_graphics("figures/TransIntegrationExample.png") 

## ---- messages = FALSE, eval = FALSE------------------------------------------
#  mirtarbase.comparison.table <- preprocess_miRTarBase(organism.code = 'mmu', org.db = 'org.Mm.eg.db')

## ---- messages=FALSE, eval=FALSE----------------------------------------------
#  mirtarbase.comparison.table <- preprocess_miRTarBase(
#    organism.code = 'mmu',
#    org.db = 'org.Mm.eg.db',
#    print.support.types = TRUE,
#    print.validation.methods = TRUE
#  )

## ---- eval=FALSE--------------------------------------------------------------
#  custom.integration <- tibble::tibble(
#    reference.expression.matrix = 'exp.yang',
#    reference.org.db = 'org.Mm.eg.db',
#    comparison.table = 'mirtarbase.comparison.table',
#    reference.table.name = 'RNA',
#    comparison.table.name = 'miRTarBase'
#  )
#  generateShinyApp(
#    expression.matrix = exp.yang,
#    metadata = meta,
#    modality = "RNA",
#    shiny.dir = "shiny_Yang2019_CustomIntegration",
#    app.title = "Shiny app for visualisation of three timepoints from the Yang 2019 data",
#    organism = "mmusculus",
#    org.db = "org.Mm.eg.db",
#    custom.integration = custom.integration
#  )
#  shiny::runApp('shiny_Yang2019_CustomIntegration')

## ----CustomIntegration, echo = FALSE, out.width = "80%"-----------------------
knitr::include_graphics("figures/CustomIntegrationExample.png") 

## ----EnrichmentIntegration, echo = FALSE, out.width = "80%"-------------------
knitr::include_graphics("figures/EnrichmentIntegrationExample.png") 

## ----sessionInfo--------------------------------------------------------------
sessionInfo()

