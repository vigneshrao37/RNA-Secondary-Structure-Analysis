library(shiny)
shinyServer(
  
 
  function(input,output) {
    output$mysequence <- renderText(input$sequence)
    output$mychoice <- renderText(input$choice)
    
    
    
  }
)