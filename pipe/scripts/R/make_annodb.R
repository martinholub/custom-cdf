## This script will make annotation database for usage in you microarraay
## data analysis starting from probe<->reference alignement information.

#################################### DO STUFF ##################################
library(methods)
# Test Class -------------------------------------------------------------------
# setClass("TestSnake",
#          representation = list(input = "character", output = "character",
#                                params = "character", log = "character", config = "list"))
# snakemake = new("TestSnake",
#                 input = "~/tmp/CDF/MM/pipe_agi/data/results/MM_8x60K_results.txt",
#                 output = "~/tmp/CDF/MM/pipe_agi/data/cdf/mm8x60kdb.txt",
#                 params = "MM_8x60K",
#                 log = "",#"~/tmp/CDF/MM/pipe_agi/logs/cdf/mm8x60kdb.log",
#                 config = list(platform = "Agilent",
#                               species = "Mus Musculus",
#                               cdf = list(dset = "mmusculus_gene_ensembl",
#                                          biomart = "ensembl",
#                                          host = "www.ensembl.org",
#                                          in_id = "ensembl_gene_id",
#                                          out_id = "",
#                                          attrs = c("entrezgene", "ensembl_transcript_id",
#                                                    "description", "external_gene_name",
#                                                    "go_id", "embl"),
#                                          schema = "",
#                                          tax_id = 10090,
#                                          org_db = "mouse.db0")))
# Parse arguments -------------------------------------------------------------
## Either use script: file.R directive
file_flat <- snakemake@input[[1]]
file_db <- snakemake@output[[1]]
chip_name <- snakemake@params$chip
snake_dir <- snakemake@params$snake_dir
logfile <- snakemake@log[[1]]
platform <- tolower(snakemake@config[["platform"]])
species_name <- tools::toTitleCase(snakemake@config[["species"]])
dset <- snakemake@config$cdf[["dset"]]
biomart <- snakemake@config$cdf[["biomart"]]
host <- snakemake@config$cdf[["host"]]
in_id <- snakemake@config$cdf[["in_id"]]
out_id <- snakemake@config$cdf[["out_id"]]
attrs <- snakemake@config$cdf[["attrs"]]
schema <- snakemake@config$cdf[["schema"]]
tax_id <- toString(snakemake@config$cdf[["tax_id"]])
org_db <- snakemake@config$cdf[["org_db"]]
# version <- snakemake@config$version

species_parts <- unlist(strsplit(species_name, " "))
species_parts[2] <- tolower(species_parts[2]) #lowercase needed 
has_log <- FALSE

# Setup Redirect
if (length(logfile)>0) {
  print(paste("Logging to ", logfile))
  fh_logfile <- file(logfile, open = "wt")
  sink(fh_logfile)
  has_log <- TRUE
}

## Or pass args via command line, python style ---------------------------------
### https://www.r-bloggers.com/passing-arguments-to-an-r-script-from-command-lines/
# suppressPackageStartupMessages(require(optparse))
# option_list <- list(
#   optparse::make_option(c("-i", "--input"), type="character", default=NULL,
#                         help="dataset file name", metavar="character"),
#   optparse::make_option(c("-o", "--output"), type="character", default="outdb.db",
#                         help="output file name [default= %default]", metavar="character"),
#   optparse::make_option(c("-s", "--species"), type="character", default="Homo Sapiens",
#                         help="Species name [default= %default]", metavar="character"),
#   optparse::make_option(c("-n", "--chip-name"), type="character", default="Agilent Chip",
#                         help="Chip name [default= %default]", metavar="character"))
# opt_parser <- optparse::OptionParser(option_list = option_list)
# opts <- optparse::parse_args(opt_parser)

## Or use debug mode DEBUG
# opts = list(input = normalizePath("~/tmp/CDF/pipe_agi/data/results/Agilent_Oligo_V4_results.txt"),
#             output= normalizePath("~/tmp/CDF/pipe/data/cdf/agilentoligov4db.db"),
#             species = "Arabidopsis Thaliana", chipname = "Agilent_Oligo_V4")

# file_flat <- opts$input
# file_db <- opts$output
# species_name <- opts$species
# chip_name <- opts$chipname

################################################################################
# Fetch functions --------------------------------------------------------------
source(file.path(snake_dir, "scripts/R", "annodb.R"))
# Debug ------------------------------------------------------------------------
# save.image("debug_image.RData")
# stop(paste0("DEBUG: ", getwd()))

# Convert between naming schemes in_id <-> out_id ------------------------------
# Done only if maping desired by config
if (!is.null(out_id)){
  if (!(in_id == out_id | out_id == "")){
    # Set up mart
    mart <- biomaRt::useMart(biomart=biomart, host=host, dataset = dset)
    file_flat <- .remap_gene_symbols(file_flat, mart, in_id, out_id)
    in_id <- out_id
  }
}

# Write Gene annotation to separate file ---------------------------------------
#.write_gene_annotation(file_flat, mart, in_id, attrs, file_db)

# Create Also chip DB if possible ----------------------------------------------
# save.image("debug_image.RData")
# stop(paste0("DEBUG: ", getwd()))
if (!is.null(schema)  &&
  # Conversion `.remap_gene_symbols` was run ... I think this is not necessary
  # !is.null(out_id) && !(in_id == out_id | out_id == "") &&  
  !grepl("NOCHIP|NOSCHEMA", schema, ignore.case = T) &&
  grepl("CHIP", schema)){
  # We have a known schema available to use for database crreation.
  .make_db_from_schema(file_db, file_flat, schema, chip_name)
} else if (!is.null(org_db) && org_db != "") {
  # TODO: Consider tryCatch fallback on generation of org_db
  .make_chipdb_from_org(file_db, file_flat, org_db, tax_id, species_parts)
  
} else {
  org_db <- .make_organism_package(file_db, tax_id, species_parts)
  .make_chipdb_from_org(file_db, file_flat, org_db, tax_id, species_parts)
}

# Finish Redirect --------------------------------------------------------------
if (has_log) {
  sink()
  if (isOpen(fh_logfile)) close(fh_logfile)
}
