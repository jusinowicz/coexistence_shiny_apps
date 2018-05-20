
library(shiny)
library(shinyIncubator)
library(rhandsontable)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  # Application title
  titlePanel("Two competitors (Leslie-Gower)"),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      matrixInput("ia",
                  "Initial abundances",
                  data.frame(matrix(c(1,1),nrow=2))),
     
      numericInput("ngens","Number of generations",20,min=10,max=1000,step=10),
     
      matrixInput("rs",
                  "Intrinsic growth rates",
                  data.frame(matrix(c(1.5,1.5),nrow=2))),
     
      matrixInput("alphamat",
                  "Competition coefficients matrix (α)",
                  data.frame(matrix(c(1,0.8,0.8,1),nrow=2))),
     
      matrixInput("surv",
                  "Survival rates (between 0 and 1)",
                  data.frame(matrix(c(0.9,0.9),nrow=2))),
     
      numericInput("cor","Correlation in reproduction",0,min=-1,max=1,step=0.05),
     
      matrixInput("var",
                  "Variance of reproduction",
                  data.frame(matrix(c(0.0,0.0),nrow=2))),
      # Button
      downloadButton("downloadData", "Download"),
     
      checkboxInput("draw_pop","Plot total population growth",value=F)
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      plotOutput("pop_space_lgf"),
      #plotOutput("pop_lgf"),
      conditionalPanel("input.draw_pop",plotOutput("pop_lgf"))
    )
  )
))