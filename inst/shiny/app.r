#### WASP GUI ####
rm(list = ls(all.names = TRUE))
#---- Load libraries ----
library(shiny)
library(patchwork)
library(BSgenome.Hsapiens.UCSC.hg38)
library(chromVAR)
library(JASPAR2020)
library(TFBSTools)

#---- Run app ----
shiny::shinyApp(
  ui = WASP::app_ui(),
  server = WASP::app_server()
)
