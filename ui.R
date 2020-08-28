library(shiny)
shinyUI(
  pageWithSidebar(
    
    headerPanel("Analysis of Secondary Structure Loops of RNAs"),
    
    sidebarPanel(
      h2("File Upload for Analysis"),
      textInput("sequence","Enter dot file: ", ""),
      textInput("choice","What Analysis would you want to perform","")
      
      
      ),
      
    
    mainPanel(
      h2("Main"),
      h5("Data Visual: ")
      )
  )
  
)