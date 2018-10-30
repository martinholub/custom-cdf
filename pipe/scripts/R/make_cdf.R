#!/usr/bin/Rscript --vanilla
## Rscript --vanilla /home/mholub/git/git-dev/CustomCDF/pipe/scripts/make_cdf.R -i /home/mholub/tmp/CDF/pipe/data/results/AraGene-1_1-st-v1_results.txt -o /home/mholub/tmp/CDF/pipe/data/cdf/aragene11stv1cdf.cdf -c 1190 -r 1190 -s "Arabidopsis Thaliana"
################################ DEFINE FUNCS #################################

#' @title flat2cdf: Creates a CDF package from mapping file
#'
#' @description ....
#'
#' @details
#' Sections where the name begins with "Unit" define the probes that
#' are a member of the unit (probe set). Each probe set is divided into subsections
#' termed "Blocks" which are referred to as "groups" in the Files SDK documentation.
#' A group will list the probes as its members. For expression probe sets there
#' is 1 and only 1 group per probe set.
#'
#' @param file_flat, str, path to a file with header and columns
#' 'probe_id, probe_x, probe_y, probe_seq, group_id, unit_id, strand, n_hits'
#' @param unit_col, int, position of unit column
#' @param group_col, int, position of group column
#' @param rows, int, number of rows on the array
#' @param cols, int, number of cols on the array
#' @param xynames, two-element vector indicating names of columns that carry x,y coords.
#' @param verbose, verbosity level
#'
#' @return list(cdfList = cdf_data, cdfHeader = hdr), invisibly
#'
#' @importFrom affxparser writeCdf
#'
#' @references https://www.affymetrix.com/support/developer/powertools/changelog/gcos-agcc/cdf.html
#' @references http://dept.stat.lsa.umich.edu/~kshedden/Courses/Stat545/Notes/AffxFileFormats/cdf.html
#' @references http://www.aroma-project.org/howtos/create_CDF_from_scratch/
#' @export
flat2cdf <- function(file_flat, cdf_name, rows = 1190, cols = 1190, verbose = 10,
                     xynames = c("probe_x","probe_y"), group_col = 5, unit_col = 6,
                     ...) {
  require(affxparser)

  # Get classes of data in flat file
  col_class <- rep("character", 8)
  fcon <- file(file_flat, open = "r")
  header <- unlist(strsplit(readLines(fcon, n = 1), "\t"))
  close(fcon)
  #col_class[c(1, match(xynames, header))] <- "integer"
  col_class[match(c(xynames, "n_hits"), header)] <- "integer"

  # Read Data from flat (txt) file
  if (verbose) cat("Reading TXT file ...")
  file_flat <- read.table(file_flat, header = TRUE, colClasses = col_class,
                          stringsAsFactors = FALSE, comment.char = "")
  if (verbose) cat(" Done.\n")

  # Split file by units
  if (verbose) cat("Splitting TXT file indices into units ...")
  unit_idxs <- split(seq_len(nrow(file_flat)), file_flat[, unit_col])
  if (verbose) cat(" Done.\n")

  # Make data on unit and eventually group level
  cdf_data <- vector("list",length(unit_idxs))

  if (verbose) cat("Creating structure for",length(unit_idxs),"units (dot=250):\n")
  for(i in  1:length(unit_idxs)) {
    curr_unit <- file_flat[ unit_idxs[[i]], ]
    ## Further split unit to groups, no splitting if unit=group
    groups_split <- split(curr_unit, factor(curr_unit[, group_col]))

    ## Make groups data for current unit
    groups_data <- vector("list",length(groups_split))
    for(j in 1:length(groups_split)) {
      grp_count <- nrow(groups_split[[j]])
      ## See references for format specification
      grp_dir <- unique(groups_split[[j]][["strand"]])
      groups_data[[j]] <- list(x = groups_split[[j]][,xynames[1]], y = groups_split[[j]][, xynames[2]],
                               pbase = rep("A",grp_count), tbase = rep("T",grp_count),
                               atom = 0:(grp_count-1), indexpos = 0:(grp_count-1),
                               groupdirection = grp_dir, natoms = grp_count, ncellsperatom = 1)
    }

    # Name of new "probe" is name of metafeature (e.g. gene) with sufix.
    names(groups_data) <- paste0(names(groups_split), "_at")
    unit_dir <- if (unique(curr_unit[["strand"]]) == "sense") 1 else 2

    ## See references for format specification
    cdf_data[[i]] <- list(unittype = 1, unitdirection = unit_dir, groups = groups_data,
                          natoms = nrow(curr_unit), ncells = nrow(curr_unit),
                          ncellsperatom = 1, unitnumber = i)
    if (verbose) { if(i %% 250==0) cat("."); if(i %% 5000==0) cat("(",i,")\n",sep="") }
  }

  # Cleanup
  names(cdf_data) <- names(unit_idxs)
  rm(file_flat, groups_data, groups_split, curr_unit, unit_idxs); gc()
  cat("\n")

  # Write data to binary CDF file format
  filename <- cdf_name
  cdf_name <- tools::file_path_sans_ext(basename(cdf_name))
  if (!grepl(".cdf$", filename)){ # append extension if needed
    filename <- paste0(filename,".cdf")
  }

  # Is there any way how to get the c(rows,cols)? E.g. from old CDF?
  hdr <- list(probesets = length(cdf_data), qcprobesets = 0, reference="",
              cdf_name = cdf_name, filename = filename,
              nqcunits = 0, nunits = length(cdf_data), rows = rows, cols = cols,
              refseq ="", nrows = rows, ncols = cols)

  # Warning: The API for this function is likely to be changed in future versions
  affxparser::writeCdf(hdr$filename, cdfheader = hdr, cdf = cdf_data,
                       cdfqc = NULL, overwrite = TRUE, verbose = verbose)

  invisible(list(cdfList = cdf_data, cdfHeader = hdr))
}

#' @title run test install of a package
#'
#' @param file_db, str path to file?
#' @param do_remove, bool, remove the package after installing?
#'
.do_test_install <- function(file_db, do_remove = TRUE){
  print(paste0("Running test install of ", file_db, " ..."))
  install.packages(file_db, repos=NULL, type="source")
  if (do_remove){
    remove.packages(basename(file_db))
  }
  print("Success.")
}

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

  # will fail for `example`, you can comment it out for this particular case
  .do_test_install(filename)
  invisible(TRUE)
}

#################################### DO STUFF ##################################
library(methods)
# Parse arguments
## Either use script: file.R directive
# file_flat <- snakemake@input[[1]]
# file_cdf <- snakemake@output[[1]]
# rows <- snakemake@params[[1]]
# cols <- snakemake@params[[2]]
# species_name <- snakemake@params[[3]]

## Or pass args via command line, python style
### https://www.r-bloggers.com/passing-arguments-to-an-r-script-from-command-lines/
### Later make input, output positional
suppressPackageStartupMessages(require(optparse))
option_list <- list(
  optparse::make_option(c("-i", "--input"), type="character", default=NULL,
                        help="dataset file name", metavar="character"),
  optparse::make_option(c("-o", "--output"), type="character", default="out.cdf",
                        help="output file name [default= %default]", metavar="character"),
  optparse::make_option(c("-r", "--rows"), type="integer", default=1190,
                        help="# of rows on chip [default= %default]", metavar="integer"),
  optparse::make_option(c("-c", "--cols"), type="integer", default=1190,
                        help="# of cols on chip [default= %default]", metavar="integer"),
  optparse::make_option(c("-s", "--species"), type="character", default="Homo Sapiens",
                        help="Species name [default= %default]", metavar="character"))

opt_parser <- optparse::OptionParser(option_list = option_list)
opts <- optparse::parse_args(opt_parser)

## Or use debug mode DEBUG
# opts = list(input = normalizePath("~/tmp/CDF/pipe/data/results/AraGene-1_1-st-v1_results.txt"),
#             output= normalizePath("~/tmp/CDF/pipe/data/cdf/aragene11stv1cdf.cdf"),
#             rows = 1190, cols = 1190, species = "Arabidopsis Thaliana")

file_flat <- opts$input
file_cdf <- opts$output
rows <- opts$rows
cols <- opts$cols
species_name <- opts$species

# Call functions
flat2cdf(file_flat, file_cdf, rows, cols)
make_cdf_package(file_cdf, species_name)

# TEST -------------------------------------------------------------------------
# require(makecdfenv)
# fnames <- affy::list.celfiles("CDF/DR/DR-99999/CEL", full.names = TRUE)
# affy::cleancdfname(affy::whatcdf(fnames[[1]]))
# # # testcdf <- makecdfenv::make.cdf.env("CDF/pipe/data/cdf/aragene11stv1.cdf")
# # # testcdf2 <- makecdfenv::make.cdf.env("CDF/AT/aragene11st_At_TAIRG_22.0.0/aragene11st_At_TAIRG.cdf")
# # # testcdf3 <- makecdfenv::make.cdf.env("CDF/AT/aragene11st_At_ENTREZG_22.0.0/aragene11st_At_ENTREZG.cdf")
# install.packages("/home/mholub/tmp/CDF/pipe_baT/data/cdf/zebgene11stv1cdf", repos = NULL, source = TRUE)
# affy_data <- affy::ReadAffy(filenames = fnames, cdfname = "zebgene11stv1cdf")
# eset_test <-affy::rma(affy_data)

# INSPECT
# my_cdf = affxparser::readCdf("~/tmp/CDF/MM/pipe/data/cdf/mogene21stv1cdf.cdf")
# my_cdf$ENSMUSG00000000001$groups$ENSMUSG00000000001_at$x
