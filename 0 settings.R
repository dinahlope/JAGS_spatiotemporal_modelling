#' JAGS model specification for spatiotemporal epidemiological modelling
#' @author Dinah Lope
#' @version 1.0
#' @created 20/06/2020
#' @license GPLv3
#' 
# includes all the the setup of
# packages and files

options(stringsAsFactors = FALSE,
        java.parameters = "-Xmx6G",
        scipen = 999)

library(rvest) #web scraping
#library(gdata)
library(dplyr) #data manipulations
library(magrittr) #pipe operator
#library(foreign) 
library(openxlsx) #reading and writing excel files
#library(readxl)
# library(read.table) #reading csv files
library(tidyr) #data organisation
library(glue) #concatenating strings
library(tidyselect) #detecting strings
# library(readr)
library(pdftools) #extracting data from pdf
library(stringr) #string operations
#library(stringi) #string processing
#library(xlsx)
#library(purrr)
#library(NSM3)
library(dineq) #compute weighted gini and theil indices
library(ineq) #plotting lorezn curve
library(DiagrammeR) #flowchart
library(leaflet) #creates interactive maps
library(htmlwidgets) #creates R bindings to JavaScript libraries
library(htmltools) #generating html and html outputs
library(lubridate)
library(dygraphs)
library(ggplot2)
library(rgeos)
library(maptools)
library(broom)
library(leaflet)
library(rgdal)
library(sp)
library(raster)
# library(xts)

