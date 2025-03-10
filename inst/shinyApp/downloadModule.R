downloadUI <- function(id) {
  ns <- NS(id)
  hidden(downloadButton(ns('downloadData'), 'Download', style = "color: #fff;
  font-weight: bold;
  background-color: #27ae60;
  border-color: #fff;
  padding: 5px 14px 5px 14px;
  margin: 5px 5px 5px 5px;" ))
}

downloadServer <- function(id, gs){
  moduleServer(
    id,
    function(input, output, session){
      #Controls the Download button
      
      observe({
        if(is.null(gs())){
          hide("downloadData")
        } else {
          show("downloadData")
        }
      })
      
      output$downloadData <- downloadHandler(
        filename = function() {
          paste("gsva_es-", Sys.Date(), ".csv", sep="")
        },
        content = function(file) {
            if("ExpressionSet" %in% class(gs())) 
            {
              expressionSetObs <- exprs(gs())
              dataFrameObs <- as.data.frame(expressionSetObs)
              write.csv(dataFrameObs, file)
            }
            else
            {
              dataFrameObs <- as.data.frame(gs())
              write.csv(dataFrameObs, file)
            } 
        }
      )
    }
  )
}