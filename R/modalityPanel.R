#' Generate an app panel for a modality
#' @description These are the UI and server components of a modality
#' panel of the shiny app. Different modalities can be included by
#' specifying their inputs in \code{\link{generateShinyApp}}.
#' @inheritParams DEpanel
#' @return The UI and Server components of the shiny module, that can be used
#' within the UI and Server definitions of a shiny app.
#' @name modalityPanel
NULL

#' @rdname modalityPanel
#' @export
modalityPanelUI <- function(id, metadata, organism, panels.default){
  ns <- NS(id)
  
  tabsetPanel(
    landingPanelUI(ns('Landing'), show = "Landing" %in% panels.default),
    sampleSelectPanelUI(ns('SampleSelect'), show = "SampleSelect" %in% panels.default),
    QCpanelUI(ns('QC'), metadata, show = "QC" %in% panels.default),
    DEpanelUI(ns('DE'), metadata, show = "DE" %in% panels.default),
    DEplotPanelUI(
      ns('DEplot'), 
      show = ("DE" %in% panels.default) & ("DEplot" %in% panels.default)
    ),
    DEsummaryPanelUI(
      ns('DEsummary'), metadata, 
      show = ("DE" %in% panels.default) & ("DEsummary" %in% panels.default)
    ),
    enrichmentPanelUI(
      ns('Enrichment'), 
      show = ("DE" %in% panels.default) & ("Enrichment" %in% panels.default) & !is.null(organism)
    ),
    GRNpanelUI(
      ns('GRN'), 
      metadata, 
      show = ("DE" %in% panels.default) & ("GRN" %in% panels.default)
    ),
    patternPanelUI(ns('Patterns'), metadata, show = "Patterns" %in% panels.default),
    crossPanelUI(ns('Cross'), metadata, show = "Cross" %in% panels.default),
    
  )
}

#' @rdname modalityPanel
#' @export
modalityPanelServer <- function(id, expression.matrix, metadata, anno, organism, panels.default){
  ns <- NS(id)
  # check whether inputs (other than id) are reactive or not
  stopifnot({
    !is.reactive(expression.matrix)
    !is.reactive(metadata)
    !is.reactive(anno)
    !is.reactive(organism)
    !is.reactive(panels.default)
  })
  
  moduleServer(id, function(input, output, session){
    if("SampleSelect" %in% panels.default){
      filteredInputs <- sampleSelectPanelServer('SampleSelect', expression.matrix, metadata)
      expression.matrix <- reactive(filteredInputs()[['expression.matrix']])
      metadata <- reactive(filteredInputs()[['metadata']])
    }else{
      expression.matrix <- reactiveVal(expression.matrix)
      metadata <- reactiveVal(metadata)
    }
    if("QC" %in% panels.default){
      QCpanelServer('QC', expression.matrix, metadata, anno)
    }
    if("DE" %in% panels.default){
      DEresults <- DEpanelServer('DE', expression.matrix, metadata, anno)
      if("DEplot" %in% panels.default){
        DEplotPanelServer('DEplot', DEresults, anno)
      }
      if("DEsummary" %in% panels.default){
        DEsummaryPanelServer('DEsummary', expression.matrix, metadata, DEresults, anno)
      }
      if("Enrichment" %in% panels.default & !is.null(organism)){
        enrichmentPanelServer('Enrichment', DEresults, organism = organism)
      }
      if("GRN" %in% panels.default){
        GRNpanelServer('GRN', expression.matrix, metadata, DEresults, anno)
      }
    }
    if("Patterns" %in% panels.default){
      patternPanelServer('Patterns', expression.matrix, metadata, anno)
    }
    if("Cross" %in% panels.default){
      crossPanelServer('Cross', expression.matrix, metadata, anno)
    }
  })
}

modalityPanelApp <- function(){
  panels.default <- c("Landing", "SampleSelect", "QC", "DE", "DEplot", "DEsummary",
                      "Enrichment", "GRN", "Patterns", "Cross")
  shinyApp(
    ui = navbarPage(
      'Shiny app for the Yang 2019 data',
      theme = shinythemes::shinytheme('flatly'),
      tabPanel('RNA', modalityPanelUI('RNA', metadata, "mmusculus", panels.default))
    ),
    server = function(input, output, session){
      modalityPanelServer('RNA', expression.matrix.preproc, metadata, anno, "mmusculus", panels.default)
    }
  )
}