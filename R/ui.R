library(shiny)
library(ggplot2)
library(edgeR)
library(gprofiler2)
library(ggrepel)
library(shinythemes)
library(shinydashboard)


# Define UI for miles per gallon app ----
ui <- navbarPage("Interaction of BCL11A and CHD8 in triple negative breast cancer", 
                 theme = shinytheme("flatly"),
                 #RNAseq---------------------------------------------------------
                 tabPanel("RNAseq",
                          tabsetPanel(
                            # tab for DE
                            DEpanelUI("DE"),
                            
                            # tab for enrichment
                            enrichmentPanelUI("Enrichment")
                          )),
                 #ChIPseq--------------------------------------------------------
                 tabPanel("ChIPseq",
                          # App title ----
                          titlePanel("Interaction of BCL11A and CHD8 in triple negative breast cancer"),
                          
                          sidebarLayout(
                            
                            # Sidebar panel for inputs ----
                            sidebarPanel(
                              
                              # Input: Selector for variable to plot against mpg ----
                              selectInput("variable1_chip", "Sample 1:",
                                          c("control BCL11A IP" = "control_11AIP",
                                            "control CHD8 IP" = "control_CHD8IP",
                                            "BCL11A KD BCL11A IP" = "11AKD_11AIP",
                                            "BCL11A KD CHD8 IP" = "11AKD_CHD8IP",
                                            "CHD8 KD BCL11A IP" = "CHD8KD_11AIP",
                                            "CHD8 KD CHD8 IP" = "CHD8KD_CHD8AIP")),
                              selectInput("variable2_chip", "Sample 2:",
                                          c("control CHD8 IP" = "control_CHD8IP",
                                            "control BCL11A IP" = "control_11AIP",
                                            "BCL11A KD BCL11A IP" = "11AKD_11AIP",
                                            "BCL11A KD CHD8 IP" = "11AKD_CHD8IP",
                                            "CHD8 KD BCL11A IP" = "CHD8KD_11AIP",
                                            "CHD8 KD CHD8 IP" = "CHD8KD_CHD8AIP")),
                              
                              #FC threshold line
                              sliderInput(
                                "FC_threshold_chip", label = "logFC threshold",
                                min = 0, value = 0.5, max = 5,step = 0.5
                              ),
                              textInput("GeneName_chip","Gene name",value ='',placeholder = 'Mrgprb2')
                              
                            ),
                            
                            # Main panel for displaying outputs ----
                            mainPanel(
                              plotOutput("maPlot_chip", click = "plot_click_chip"),
                              tableOutput("maData_chip"),
                            )
                            
                          )),
                 
                 #ATACseq--------------------------------------------------------
                 tabPanel("ATACseq",
                          
                          # App title ----
                          titlePanel("Interaction of BCL11A and CHD8 in triple negative breast cancer"),
                          
                          # Sidebar layout with input and output definitions ----
                          sidebarLayout(
                            
                            # Sidebar panel for inputs ----
                            sidebarPanel(
                              
                              # Input: Selector for variable to plot against mpg ----
                              selectInput("variable1_atac", "Sample 1:",
                                          c("control" = "control",
                                            "BCL11A" = "BCL11A",
                                            "CHD8" = "CHD8")),
                              selectInput("variable2_atac", "Sample 2:",
                                          c("BCL11A" = "BCL11A",
                                            "control" = "control",
                                            "CHD8" = "CHD8")),
                              
                              #FC threshold line
                              sliderInput(
                                "FC_threshold_atac", label = "logFC threshold",
                                min = 0, value = 0.5, max = 5,step = 0.5
                              ),
                              textInput("GeneName_atac","Gene name",value ='',placeholder = 'Mrgprb2')
                              
                            ),
                            
                            # Main panel for displaying outputs ----
                            mainPanel(
                              plotOutput("maPlot_atac", click = "plot_click_atac"),
                              tableOutput("maData_atac"),
                            )
                            
                          ))
                 
)