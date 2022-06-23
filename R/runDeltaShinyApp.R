#' Run Delta Shiny App
#' 
#' Launch the Delta Shiny App to compute the Delta Measurement of Agreement for Nominal, from a user friendly interface.
#' The app also allows to download a report of results in .pdf or .tex (LaTeX) format. The latter format is compressed in a .zip file.
#' @export
#' @return No return value. A shiny app is launched.
runDeltaShinyApp <- function() {

  appDir <- system.file("shiny-examples", "deltaShiny", package = "DeltaMAN")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `DeltaMAN`.", call. = FALSE)
  }
  
  suppressWarnings({
    shiny::runApp(appDir, display.mode = "normal")
  })
  
}