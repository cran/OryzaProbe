## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval = nzchar(Sys.getenv("MY_ENV_VARIABLE"))
)

## ----prep---------------------------------------------------------------------
#  #library(GEOquery)
#  #expr <- getGEO("GSE75471")[[1]]
#  library(OryzaProbe)
#  expr <- test_file
#  head(test_file)

## ----ID conversion without merging multiple probes----------------------------
#  library(OryzaProbe)
#  formatted_expr <- probeConvert(exprMatrix = expr)
#  head(formatted_expr)

## ----ID conversion with merging multiple probes-------------------------------
#  formatted_expr1 <- probeConvert(exprMatrix = test_file, probeMerge = T)
#  head(formatted_expr1)

