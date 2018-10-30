---
title: Making a CDF package
author: Martin Holub
date: Sept 04, 2018
output: html_document
---
# Making a CDF package

To obtain a valid binary CDF file  from either binary or ASCII cdf file:

```R
#' Get a valid CDF annotation file from ASCII cdf file
#'
#' @details
#' File can be obtain e.g. from BrainArray website
.make_binary_cdf_from_file <- function(annodb, annodata_path, chiptype){

  assertthat::assert_that(file.exists(annodb))

  out_path <- file.path(annodata_path, paste0(chiptype, ".cdf"))
  if (file.exists(out_path)){
    print(paste("Using file at:", out_path))
  } else {
    print("Converting ASCII CDF to binary....")
    affxparser::convertCdf(annodb, out_path, force = TRUE)
  }

  return(tools::file_path_sans_ext(basename(out_path)))
}
```

To pull out a CDF from existing package:

```R
#' Get a valid CDF annotation file from existing R package
#'
#' @details
#' package must be installed (e.g. from bioconductor)
.get_cdf_from_package <- function(annodb, annodata_path, rawdata_path){
  suppressPackageStartupMessages(library(annodb, character.only = T))
  fnames <- list.files(rawdata_path, full.names = T)
  targetf <- file.path(annodata_path, paste0(gsub("cdf$", "",annodb), ".cdf"))

  if (!(file.exists(targetf))){
    print("Pulling out CDF from package....")
    pathname <- aroma.affymetrix::env2Cdf(annodb, fnames[[1]])
    system(paste0("mv -v ", pathname, " ", annodata_path))
  } else{
    print(paste0("Using existing cdf at: ",targetf," ...."))
    pathname <- targetf
  }

  return(tools::file_path_sans_ext(basename(pathname)))
}
```

To make a CDF package from binary CDF:

```R
#' @title make_cdf_package: Create R package from binary CDF
#'
#' @description ....
#'
#' @param cdf_name path to cdf file
#' @param species_name caption for the species
#' @param verbose bool, controls vebrosity
#'
#' @return TRUE, invisibly
#'
#' @importFrom makecdfenv make.cdf.package
#'
#' @export
make_cdf_package <- function(cdf_file, species_name, verbose = TRUE){
  require(makecdfenv)

  filename <- cdf_file
  cdf_name <- tools::file_path_sans_ext(basename(cdf_file))
  if (!grepl(".cdf$", filename)){ # append extension if needed
    filename <- paste0(filename,".cdf")
  }

  cdf_path <- dirname(filename)
  makecdfenv::make.cdf.package(basename(filename), cdf_name, cdf_path, cdf_path,
                               author = "Nebion", unlink = TRUE, verbose = verbose,
                               species = species_name)
  invisible(TRUE)
}
```
