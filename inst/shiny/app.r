#### WASP GUI ####
#---- Load libraries ----
library(shiny)
library(patchwork)
library(BSgenome.Hsapiens.UCSC.hg38)
library(chromVAR)
library(JASPAR2020)
library(TFBSTools)

#---- Run app ----
shiny::shinyApp(
  ui = app_ui(),
  server = app_server()
)
