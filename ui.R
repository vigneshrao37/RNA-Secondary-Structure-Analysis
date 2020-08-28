library(shiny)
python_source = ("parser.py")
shinyUI(
  pageWithSidebar(
    
    headerPanel("Analysis of Secondary Structure Loops of RNAs"),
    
    sidebarPanel(
      
      h2("File Upload for Analysis"),
      fileInput("file","Upload dot file"),
      radioButtons("sep","Seperator",choices = c(Comma = ",",Period = ".")),

      
      ),
      
    
    mainPanel(
      h2("Main"),
      h5("Data Visual: "),
      textOutput("mychoice"),
      tableOutput("input_file")
      )
  )
  
)