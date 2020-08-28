library(shiny)
shinyServer(
  
 
  function(input,output) {
    output$input_file <- renderTable({
      
      file_to_read = input$file
      if(is.null(file_to_read)){
        return()
      }
      read.table(file_to_read$datapath, sep = input$sep )
    })
    
  }
)