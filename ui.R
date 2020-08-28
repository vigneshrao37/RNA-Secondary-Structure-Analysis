
library(reticulate)
#py_run_file("/Users/vigneshrao/secondary/parser.py")
library(shiny)
shinyUI(
  pageWithSidebar(
    
    headerPanel("Analysis of Secondary Structure Loops of RNAs"),
    
    sidebarPanel(
      
      h2("File Upload for Analysis"),
      fileInput("file","Upload File to be Parsed: "),
      radioButtons("sep","Seperator",choices = c(Comma = ",",Period = ".",Tilde = "~")),
      tags$img(src="dna.png.jpeg",height = 150,width = 150)
      
      ),
      
    
    mainPanel(
      h2("File Visual: "),
      h5("How Input File Looks like: "),
      tableOutput("input_file"),
      br(),
      br(),
      br(),
      h2("Data Visual: "),
      textOutput("mychoice")
      
      )
  )
  
)