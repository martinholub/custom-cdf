########################## FUNCTIONS ###########################################
# Convert between naming schemes in_id <-> out_id V2----------------------------
#' @details 
#' For "Mart" expects sth like in_id: "ensmebl_gene_id", out_id: "entrezgene"
#' For "str" expects sth like  in_id: "ENSEMBL", out_id: "ENTREZID"
.remap_gene_symbols <- function(file_flat, mart, in_id, out_id, dropna = FALSE){

  align_data <- read.table(file_flat,
                           colClasses = c("character", "character"),
                           col.names = c("probe", in_id))
  attribs <- c(in_id, out_id)
  ensembl_ids <- align_data[, 2]
  
  if (class(mart) == "Mart") {
    map_result <- biomaRt::getBM(attributes = attribs, mart = mart,
                                 values = ensembl_ids, filters = in_id,
                                 uniqueRows = FALSE)
    
    # Expand to include duplicate occurences
    map_result_ <- map_result[match(align_data[[in_id]], map_result[[in_id]]), ]
    map_result <- as.data.frame(cbind(align_data$probe, map_result_[[out_id]]))
    
  }else{ # str giving name of package that points to OrgDB
    
    orgdb <- eval(parse(text=paste0(mart, ".db")))                         
    remapped_genes <- mapIds(orgdb, align_data[[in_id]], out_id, in_id, multiVals = "first")
    map_result <- as.data.frame(cbind(align_data$probe, unname(remapped_genes)))
    map_result[[in_id]] <- remapped_genes
  }
  
  if (dropna) {
    align_data_remapped <- subset(map_result, !is.na(map_result[, 2]))
  } else {
    align_data_remapped <- map_result
  }
  
  # Write results to file -----------------------------------------------------
  suffix <- paste0("_", substr(attribs[2], 1, 5), ".")
  file_flat <- file.path(dirname(file_flat),
                         paste0(tools::file_path_sans_ext(basename(file_flat)),
                                suffix, tools::file_ext(file_flat)))
  write.table(align_data_remapped, file = file_flat,
              quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  print(paste("Written file to:", file_flat))
  invisible(file_flat)
}

# Convert between naming schemes in_id <-> out_id ------------------------------
.remap_gene_symbols_old <- function(file_flat, mart, in_id, out_id){
  print(paste0("Remaping feature IDs from ", in_id, " to ", out_id))
  # Read data ---------------
  align_data <- read.table(file_flat,
                           colClasses = c("character", "character"),
                           col.names = c("probe", in_id))
  ensembl_ids <- (align_data[, 2])
  
  attribs <- c(in_id, out_id)
  map_result <- biomaRt::getBM(attributes = attribs, mart = mart,
                               values = ensembl_ids, filters = in_id,
                               uniqueRows = FALSE)
  
  # Approach #1 - SUPERSEDED ---------------------------------------------------
  # drop NAs and nonuiqnue mapings
  # map_result <- map_result[!(apply(map_result, 1, function(x) any(is.na(x)))), ]
  
  # Order by both IDs
  # this will alow us to drop such genes that have higher entrezgene id
  # this is done because the lower number is standard, wheeras the upper is even impossible to google, so probably quite obscure
  # map_result <- map_result[order(map_result[, 1], map_result[, 2]), ]
  
  # Remove all genes that have more out_id matches for a single in_id  
  # idxer <- (duplicated(map_result[, 1])) # |
  # duplicated(map_result[, 1], fromLast = TRUE)) #| # also first occurence
  # (duplicated(map_result[,2]) | # also duplicates in other column
  # duplicated(map_result[, 2], fromLast = TRUE))
  
  # map_result <- map_result[!idxer, ]
  ## Output to file
  # entrez_data <- merge(align_data, map_result, all = FALSE, by = in_id)[, c(2,3)]
  # by = intersect(names(align_data), names(map_result))
  
  # Approach #2 ----------------------------------------------------------------
  # Expand to include duplicate occurences, keep nans
  map_result_ <- map_result[match(align_data[[in_id]], map_result[[in_id]]), ]
  entrez_data <- as.data.frame(cbind(align_data$probe, map_result_[[out_id]]),
                               col.names = c("probe", out_id))
  
  # Write results to file -----------------------------------------------------
  suffix <- paste0("_", substr(attribs[2], 1, 5), ".")
  file_flat <- file.path(dirname(file_flat),
                         paste0(tools::file_path_sans_ext(basename(file_flat)),
                                suffix, tools::file_ext(file_flat)))
  write.table(entrez_data, file = file_flat,
              quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  print(paste("Written file to:", file_flat))
  invisible(file_flat)
}

## Test ------------------------------------------------------------------------
.do_test_install <- function(file_db, do_remove = TRUE){
  print(paste0("Running test install of ", file_db, " ..."))
  install.packages(file_db, repos=NULL, type="source")
  if (do_remove){
    remove.packages(basename(file_db))
  }
  print("Success.")
}

# Create Annotation Database ---------------------------------------------------
#' Create A chip package from alignments in Entrez Coordinates
#' 
#' Prerequisite for this function is an existing dbschema. This will be present in
#' `AnnotationForge::available.dbschemas()`. The function operates with entrez ("eg") 
#' as base gene id. This results to NaNs if there is no known mapping between in_id
#' and out_id (TO BE VERIFIED!!). To make this more transparent, gene id mapping is 
#' done separately. in `.remap_gene_symbols`.
#' 
#' If schema does not exist, then one should create own organism package (or use
#' an exisitng one), e.g. with `AnnotationForge::makeOrgPackageFromNCBI` or
#' `AnnotationForge::makeOrgPackage`. See the `anno_db_notebook.Rmd` for some 
#' details. You may have some succes using `OrganismDbi::makeOrganismDbFromBiomart`,
#' remains to be investigated. This organism database can then be used with
#' `AnnotationForge::makeChipPackage`. Some of them are available from Bioconductor.
#' 
.make_db_from_schema <- function(file_db, file_flat, schema, chip_name){
  
  dbprefix <- tools::file_path_sans_ext(basename(file_db))
  dbdir <- dirname(file_db)
  # Clean up
  unlink(file_db, recursive = TRUE)
  unlink(file.path(dbdir, paste0(dbprefix, ".sqlite")))
  
  AnnotationForge::makeDBPackage(schema = schema,
                                 affy = FALSE, prefix = dbprefix,
                                 fileName = file_flat,
                                 baseMapType = "eg",
                                 outputDir = dbdir,
                                 version = "0.0.1", manufacturer = "Agilent",
                                 chipName = chip_name,
                                 author = "MH <mholub@nebion.com>")
  
  file_db <- file.path(dbdir, dbprefix) # paste0(dbprefix, ".db")
  print(paste0("Created ChipDB in ", file_db))
  .do_test_install(file_db)
}


# Make Organism Package from Biomart --------------------------------------------
#' Make OrganismDb from Biomart
#' 
#' @description 
#' Among all the possible ways how to obtain an annotation of an organism that
#' is not among your installed packages, this one is by far the easiest to use.
#' It wraps approaches explored in `anno_db_notebook.Rmd`. See also examples,
#' to get an idea what are the important functions here. As the annotation of ana o
#' 
#' @details 
#' Note that it would be also possible to create an OrgDB package using
#' `OrganismDbi::makeOrganismPackage` which could be used for subsequent runs.
#' This is not implemented yet.
#' 
#' @note 
#' Can consider saving the package in some common location.
#' 
#' @param tax_id, int, taxonomy_id
#' @param file_db tba
#' @param species_parts, tba
#' 
#' @return orgdb_name, name of enviroment with the OrgDB in it.
#' 
#' @seealso 
#' OrganismDbi::createOrganismPackage.R, OrganismDbi::makeOrganismDbFromBiomart,
#' OrganismDbi::makeOrganismDbFromTxDb
#' 
.make_organism_package <- function(file_db, tax_id, species_parts){
  dbdir <- dirname(file_db)# pass as config

    # orgdb_name <- OrganismDbi:::.taxIdToOrgDbName(tax_id)
  orgdb_name <- AnnotationForge:::.generateOrgDbName(species_parts[1], species_parts[2])
  orgdb <- OrganismDbi:::.taxIdToOrgDb(tax_id)
  
  db_path_sqlite <- file.path(dbdir, paste0(orgdb_name, ".sqlite"))
  file.copy(orgdb$conn@dbname, db_path_sqlite)
  pkg_name <- paste0(orgdb_name, ".db")
  assign(pkg_name, orgdb, .GlobalEnv)
  
  ## This Approach is totaly valid, but requires change of source code ---------
  ## Pull request was submitted
  mdata <- AnnotationDbi::metadata(orgdb)
  seed <- new("AnnDbPkgSeed",
              Package= pkg_name,
              Version="0.0.1",
              Author="MH <mholub@nebion.com>",
              Maintainer="MH <mholub@nebion.com>",
              PkgTemplate="NOSCHEMA.DB",
              AnnObjPrefix=orgdb_name,
              organism = mdata$value[mdata$name == "ORGANISM"],
              species = mdata$value[mdata$name == "SPECIES"],
              biocViews = "annotation",
              manufacturerUrl = "none",
              manufacturer = "none",
              chipName = "none")

  AnnotationForge::makeAnnDbPkg(seed, db_path_sqlite, dest_dir = dbdir,
                                no.man = TRUE)
  .do_test_install(file.path(dbdir, pkg_name), do_remove = FALSE)
  return(pkg_name)
  
  ## This approach is invalid, will fail due to cyclic dependency --------------
  # gene_keytype <- AnnotationDbi::chooseCentralOrgPkgSymbol(orgdb)
  # graph_data <- list(join1 = setNames(object=c('GOID', 'GO'),
  #                                    nm=c('GO.db', orgdb_name)),
  #                   join2 = setNames(object=c(gene_keytype, 'GENEID'),
  #                                    nm=c(orgdb_name, txdb_name)))
  
  # pkg_path <- OrganismDbi::makeOrganismPackage(
  #   "vapackage", graph_data, 
  #   destDir = "~/tmp/CDF",
  #   organism = "Vigna angularis", 
  #   version = "0.0.1",
  #   author = "MH <mholub@nebion.com>",
  #   maintainer = "MH <mholub@nebion.com>")
  
  ## This is not needed --------------------------------------------------------
  # txdb <- GenomicFeatures::makeTxDbFromBiomart(
  #   biomart="plants_mart",
  #   dataset="vangularis_eg_gene",
  #   transcript_ids=NULL,
  #   filter="",
  #   id_prefix="ensembl_",
  #   host="plants.ensembl.org")
  # txdb_name <- GenomicFeatures::makePackageName(txdb)
  # assign(txdb_name, txdb)
  # assign(orgdb_name, orgdb, .GlobalEnv)
  # return(orgdb_name)
}

# Make Chip Package from Ensembl -----------------------------------------------
#' Create a chip database using created or available orgnism package.
#' #' 
#' @seealso `.makeDBPackage`, `anno_db_notebook.Rmd`
#' 
.make_chipdb_from_org <- function(file_db, file_flat, org_db, tax_id, species_parts){
  dbprefix <- tools::file_path_sans_ext(basename(file_db))
  dbdir <- dirname(file_db)# pass as config
  probe2gene <-read.table(file_flat, sep = "\t", col.names = list("probes", "genes"))
  accessions <- 
  
  pkg_path <- AnnotationForge::makeChipPackage(
    prefix = dbprefix,
    probeFrame = probe2gene,
    orgPkgName = org_db,
    outputDir = dbdir,
    version = "0.0.1",
    tax_id = tax_id,
    genus = species_parts[1],
    species = species_parts[2],
    author = "MH <mholub@nebion.com>",
    maintainer = "MH <mholub@nebion.com>"
    # optionalAccessionsFrame = dframes$accessions
  )
  
  file_db <- file.path(dbdir, paste0(dbprefix, ".db"))
  print(paste0("Created ChipDB in ", file_db))
  .do_test_install(file_db)
  # invisible(TRUE)
}

################################################################################
#' Write annotation for genes to separate file
#' 
#' This can be used if your ChipDB package does not provide this information,
#' but ideally it would.
.write_gene_annotation <- function(file_flat, mart, in_id, attrs, file_db){
  align_data <- read.table(file_flat,
                           colClasses = c("character", "character"),
                           col.names = c("probe", in_id))
  in_ids <- unique(align_data[, 2])
  if (!(in_id %in% attrs)){
    attrs <- c(in_id, attrs)
  }
  # Add platform specifc attributes ----------------------------------------------
  do_add_platform_attributes <- FALSE
  if (do_add_platform_attributes){
    additional_attrs <- biomaRt::listAttributes(mart, "feature_page")$name
    additional_attrs <- additional_attrs[grep(platform, additional_attrs)]
    if (length(additional_attrs)>0){
      attrs <- c(attrs, additional_attrs)
    }
  }
  # Combine data -----------------------------------------------------------------
  # `embl` and `go_id` if required are often mutiple for each gene_id, this makes the
  # table bloated. Preferentially do not require them in config.
  map_result <- biomaRt::getBM(attributes = unique(attrs), mart = mart,
                               values = in_ids, filters = in_id,
                               uniqueRows = FALSE)
  combi_data <- merge(align_data, map_result, all = FALSE, by = in_id)
  
  # Write results to file --------------------------------------------------------
  # This file is currently not used as creation of ChipDb package is ususally 
  # succesfull (see func `cv_rep_probes` in pipeline).
  write.table(combi_data, file = file_db,
              quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  print(paste0("Written annotation dataframe to ", file_db))
}