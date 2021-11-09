generateShinyApp <- function(
  expression.matrix,
  metadata,
  shiny.dir = "shiny_bulkAnalyser",
  app.title = "Visualisation of RNA-Seq data",
  organism = "hsapiens",
  org.db = "org.Hs.eg.db",
  theme = "flatly",
  panels.default = c("QC", "DE", "DEplot", "Enrichment"),
  panels.extra = tibble::tibble(
    UIfun = NULL, 
    UIvars = NULL, 
    serverFun = NULL, 
    serverVars = NULL
  ),
  data.extra = c(),
  packages.extra = c()
){
  generateAppFile(
    shiny.dir = shiny.dir,
    app.title = app.title,
    organism = organism,
    org.db = org.db,
    theme = theme,
    panels.default = panels.default,
    panels.extra = panels.extra,
    packages.extra = packages.extra
  )
  generateDataFiles(
    expression.matrix = expression.matrix,
    metadata = metadata,
    shiny.dir = shiny.dir,
    data.extra = data.extra
  )
  invisible(shiny.dir)
}

generateDataFiles <- function(
  expression.matrix,
  metadata,
  shiny.dir,
  data.extra
){
  if(any(c("expression_matrix", "metadata") %in% data.extra)){
    stop("expression_matrix and metadata are reserved names, please rename your extra objects")
  }
  save(expression.matrix, file = paste0(shiny.dir, "/expression_matrix.rda"))
  save(metadata, file = paste0(shiny.dir, "/metadata.rda"))
  lapply(data.extra, function(name){
    save(get(name), file = paste0(shiny.dir, "/", name, ".rda"))
  })
}

generateAppFile <- function(
  shiny.dir,
  app.title,
  organism,
  org.db,
  theme,
  panels.default,
  panels.extra,
  packages.extra
){
  lines.out <- c()
  
  packages.to.load <- c("shiny", "dplyr", "ggplot2", "bulkAnalyseR", packages.extra)
  code.load.packages <- paste0("library(", packages.to.load, ")")
  lines.out <- c(lines.out, code.load.packages, "")
  
  shiny.dir <- normalizePath(shiny.dir)
  code.source.objects <- c(
    paste0("r.files <- list.files(path = '", shiny.dir, "', pattern = '\\.R$')"),
    "r.files <- setdiff(r.files, 'app.R')",
    "for(fl in r.files) source(fl)",
    "rda.files <- list.files(pattern = '\\.rda$')",
    "for(fl in rda.files) load(fl)"
  )
  lines.out <- c(lines.out, code.source.objects, "")
  
  code.ui <- c(
    "ui <- navbarPage(",
    glue::glue("'{app.title}',"), 
    glue::glue("theme = shinythemes::shinytheme('{theme}'),"),
    "tabPanel('RNAseq',",
    "tabsetPanel("
  )
  if("QC" %in% panels.default) code.ui <- c(code.ui, "QCpanelUI('QC', metadata),")
  if("DE" %in% panels.default){
    code.ui <- c(code.ui, "DEpanelUI('DE', metadata),")
    if("DEplot" %in% panels.default) code.ui <- c(code.ui, "DEplotPanelUI('DEplot'),")
    if("Enrichment" %in% panels.default) code.ui <- c(code.ui, "enrichmentPanelUI('Enrichment'),")
  }
  for(i in seq_len(nrow(panels.extra))){
    code.ui <- c(code.ui, glue::glue("{panels.extra$UIfun[i]}({panels.extra$UIvars[i]}),"))
  }
  code.ui <- c(code.ui, ")),", ")")
  lines.out <- c(lines.out, code.ui, "")
  
  code.server <- c("server <- function(input, output, session){")
  if("QC" %in% panels.default){
    code.server <- c(code.server, "QCpanelServer('QC', expression.matrix, metadata)")
  }
  if("DE" %in% panels.default){
    code.server <- c(
      code.server, 
      glue::glue("DEresults <- DEpanelServer('DE', expression.matrix, metadata, '{org.db}')")
    )
    if("DEplot" %in% panels.default){
      code.server <- c(code.server, "DEplotPanelServer('DEplot', DEresults)")
    }
    if("Enrichment" %in% panels.default){
      code.server <- c(
        code.server, 
        glue::glue("enrichmentPanelServer('Enrichment', DEresults, organism = '{organism}')")
      )
    }
  }
  for(i in seq_len(nrow(panels.extra))){
    code.server <- c(code.server, glue::glue("{panels.extra$serverFun[i]}({panels.extra$serverVars[i]})"))
  }
  code.server <- c(code.server, "}")
  lines.out <- c(lines.out, code.server, "")
  
  lines.out <- c(lines.out, "shinyApp(ui, server)")
  
  lines.out <- gsub("\\\\", "\\\\\\\\", lines.out)
  
  write(lines.out, paste0(shiny.dir, "/app.R"))
  
}