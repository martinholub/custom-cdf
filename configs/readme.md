---
title: "Custom CDFs Configuration"
author: Martin Holub
date: Oct 10 , 2018
output:  pdf_document  
highlight: zenburn
toc: true
---

# Introduction

This manual describes how one can control the execution of `CustomCDF` pipeline with config. In fact, config should be the only part of the pipeline that needs to be adjusted for any respective experiment.

# Common Options
## Input

### Probes
Files can be defined as local or can be obtained from remote location. The procedure is
automated and tries to achieve consistent naming and formating.

Files can be usually fetched from:

1. [Agilent](https://earray.chem.agilent.com/earray/)
2. [Affymetrix](http://www.affymetrix.com/analysis/downloads/data/) (or product page, e.g. [here](https://www.affymetrix.com/support/technical/byproduct.affx?product=arabtiling))

**Appears in:**

- dir,url,file: `setup_workflow` (and propagates to `normalize_probenames`, `unique_probes`, `unzip_files`, `bowtie_align`)

```yaml
probes:
  dir: data/probes
  url: # this is a local file
  # agilent probes often have "\" in name, may have to remove manually
  file: FASTA_028005_D_Fasta_20171030.txt
# or
probes:
  dir: data/probes
  # affymetrix zips come with weird names, set "file" to something you like
  url: http://www.affymetrix.com/analysis/downloads/data/Maize.probe_fasta.zip
  file: MaizeGenomeArray.probe.fa
```

### Reference
Reference, Annotation and Variants are obtained from Ensembl ([Biomed](https://www.ensembl.org/info/data/ftp/index.html) or [Plants](http://plants.ensembl.org/info/website/ftp/index.html)). Reference is preferentially transcriptome (with `is_transcriptome=True`) and a single file.

Arbitrary combination of chromosomes can be concatenated (see function `parse_chomosomes` in `setup_workflow` rule).

**Appears in:**

- dir, file: `unzip_files`
- dir, url, file, chromosomes: `setup_workflow`
- merged: `merge_ref`
- is_transcriptome: `snp_filter`, `extract_data`

```yaml
reference:
  dir: data/reference
  url: # this is a local file
  file: Mus_musculus_GRCm38_p5.transcripts.fa
  merged: # nothing done
  is_transcriptome: True
  # or
reference:
  dir: data/reference
  url: # this is a local file
  file: Danio_rerio.GRCz10.dna.chromosome.*.fa
  chromosomes: [{1..25}, MT] # allows arbitrary chromosomes to be selected as reference
  merged: Danio_rerio.GRCz10.dna.all.fa
```

### Index
To create index for alignment, [bowtie1](http://bowtie-bio.sourceforge.net/manual.shtml) is used. (*Suffix is fixed and is set here to allow definition of output for Snakemake.*)

**Appears in:**

- dir, prefix, suffix: `bowtie_build`

```yaml
index:
  dir: data/index
  prefix: MM_GRCm38
  suffix: ["1", "2", "3", "4", "rev.1", "rev.2"] # dont change
```
## Annotation
Here you will often have choice between `<species>.<version>.chr.gtf` and `<species>.<version>.gtf`. Prefer the latter as it includes also [nonchromosomal scaffolds](https://www.ensembl.org/info/genome/genebuild/chromosomes_scaffolds_contigs.html). In general, when downloading files from Ensembl, the `file` should corespond to what appears in `url` just with `.gz` removed.  

**Appears in:**

- dir, url, file: `setup_workflow`
- dir, file: `unzip_files`, `snp_transcriptome`, `count_features`, `extract_data`

``` yaml
annotation:
  dir: data/annotation
  # This file will be fetched from online location
  url: ftp://ftp.ensembl.org/pub/release-93/gtf/mus_musculus/Mus_musculus.GRCm38.93.gtf.gz
  file: Mus_musculus.GRCm38.93.gtf
```
## Output
### Generic

**Appears in:**

- platform, version, species: `make_annodb`/`make_cdf`

```yaml
project: "MM Transcriptomic CDF" # arbitrary
version: "0.1" # arbitrary
species: "Mus musculus" # Mind capitalization: "Genus species"
platform: "Agilent" # controls cdf or DB creation
```

### Alignment
Special care should be given to selecting the value for `prefix` as this one should be both informative and concise. It appears in multiple places in the pipeline to define names of intermediate files and also of the resulting annotation file (if no `cdf: name` provided). For affymetrix, the name **must** comply with the name appearing in the name of the probe in `FASTA` file e.g. `probe:<PREFIX>:<probe_id>-<cluster_id>;<probe_x>:<probe_y>;`.

By default, no mismatches (`max_mismatches`) are allowed, but this can be relaxed in design of custom CDF.

**Appears in:**

- dir, prefix: `bowtie_align`, `samtools_sort`, `snp_filter`, `samtools_index`, `bam_to_sam`
- prefix (additonaly): `count_features`, `extract_data`, `make_cdf`, `make_annodb`
- max_mismatches: `bowtie_align`

``` yaml
alignment:
  dir: data/alignment
  prefix: MM_8x60K # this is used a lot and must comply with Affy probe naming in FASTA
  max_mismatches: 0 # no mismatch allowed
```

### Variants (Input + Output)
The pipeline allows to filter out matches that have known SNV location in central part of the probe sequence. Central part is defined as `[0 + probe_clip, len(probe) - probe_clip)` interval along the sequence. Filtering can be turned off with `filter=False`.

**Appears in:**

- dir, url, file: `setup_workflow`
- dir, file: `unzip_files`, `snp_transcriptome`, `snp_filter`
- filter, probe_clip: `snp_filter`

``` yaml
variants:
  dir: data/variants
  url: ftp://ftp.ensembl.org/pub/release-93/variation/vcf/mus_musculus/mus_musculus.vcf.gz
  file: mus_musculus.vcf
  filter: False
  probe_clip: 5
```
### Counts

Counting is done on `feature` level and counts are aggregated on `metafeature` level. Approach differs significantly for genomic and transcriptomic alignment. In particular, for *transcriptomic* alignment we do both counting and results extraction ourselves. For *genomic* alignment, the counting is done with `FeatureCounts` and only results are extracted by us (in consistence both in form and approach with transcriptomic results extraction).

The important settings pertaining to these steps are:
- `feature`: On which level do we count hits to the reference as reported in SAM file by bowtie_align? For transcriptomic alignment (`reference: is_transcriptome = True`), always `transcript`. For genomic alignment, can be any of `exon`, `CDS`, `transcript`, `gene`.
- `metafeature`: On which level do we extract the counts? For transcriptomic and genomic alignment, this can be either `gene` (gene-level CDF ) or `transcript` (transcript-level CDF).
- `do_add_transcript_version`: If available, should the the transcript version be appended to transcript name (used when  `metafeature = "transcript"`).
- `controls`: On affymetrix chips, some probes come under generic names and do not comply to naming convention. This list defines which will be dropped.


Please consult rules `count_features` and `extract_data` for more information.

**Appears in:**

- dir, feature, metafeautre: `count_features`, `extract_data`
- do_add_transcript_version, controls: `extract_data`

``` yaml
counts:
  dir: data/counts
  feature: "transcript" #  ("transcript", "exon")
  metafeature: "gene" # ("gene", "transcript")
  do_add_transcript_version: True # e.g. <transcript`>.<version> == "ENST000000.1"
  controls: ["affx", "antigenomic", "bac_spike", "bgp-", "intron", "exon",
            "polya_spike", "unmapped", "FLmRNA"]
```
### Results

**Appears in:**

- dir: `extract_data`


``` yaml
results:
  dir:
    data/results
```
### Report
Report is constructed with `report.Rmd` file using information pulled out from logs that are saved in `./logs` directory. This approach is quite fragile but simple and should be ok, if steps of pipeline not changed much. If new step is included, also the `report.Rmd` file should be extended.

**Appears in:**

- dir: `workflow_graph`, `report`

``` yaml
report:
  dir:
    data/report
```

### CDF

We can impose some criteria on generation of CDF (affy) or ChipDB (agilent) files.

*The following two paragraphs need to be reviewed by molecular biologist. It is quote important to get the terminology correct else it becomes practically impossible to understand correctly or communicate.*

The `sense` of an array is either `ST` (usually for affymetrix) or `AT` (usually for agilent). For `ST`, probes that align (and that represent valid cDNA sequences that would match) to these probes should be `antisense` (SAM flag `16`, appearing as **reverse complemented (RC)** in SAM file). For `AT`, probes probes aligning to the refernce should be taken as-is, corresponding to SAM flag `0` (`fw` tag for forward). Optionally, you can choose `sense` to be `both`, which will keep probes aligning in both directions.

This is convoluted and the best way how to find out what is the correct option is to inspect SAM alignments file. The SAM tags (`0` or `16`) will usually be heavily skewed to ne of them, indicating the default choice for an array. Note that you could still inspect probes that align in opposite direction then the default one (minority), by setting the opposite `sense`. This would allow one to investigate binding cDNA translated from mRNA transcribed from the reverse (noncoding) strand, which can give some insight on regulatory activity.

On additional note, this becomes even more convoluted as old affymetrix arrays are by-default `AT` compared to the new ones that are mostly `ST`.

The `multimatch_level` controls how we handle probes that have some known multi-matches. Setting it to `gene` and `transcript` will discard all probes that match to multiple genes or transcripts respectively. Leaving it empty will keep all multimatching probes (For genomic alignment, set it to `gene` and control summarization level via *featureCounts*). `whitelist` should be empty or path to a file with one column, giving names of features (genes or transcripts as per `multimatch_level`)  that should not be considered when filtering probes for multimatches. These can be for example nonchromosmal genes/transcripts but also any other.

The `min_probeset_size` is a design choice. 3 is value used by BrainArray people, 1 should be always used with agilent that has just one probe per probeset.  

**Appears in:**

- dir: `make_cdf`, `make_annodb`
- sense, multimatch_level, min_probeset_size, whitelist: `extract_data`

```yaml
cdf:
  name: maizecustom # dashes and underscores will be removed, `cdf` will be appended
  dir: "data/cdf"
  sense: "AT" # ("AT", "ST", "both")
  multimatch_level: "gene" # ("transcript", "gene", "")
  min_probeset_size: 1 # 1 for agilent, 3 for affy
  whitelist: "" # path to file
```

## Agilent

### ChipDB

In order to produce a valid ChipDB file, annotation for an organism must be available. This can be either installed `OrgDb` package (available from [Bioconductor](https://bioconductor.org/packages/3.7/data/annotation/)),
or it can be represented as a `chipdbschema` (available ones can be listed with `AnnotationForge::available.chipdbschemas()`). If none of these exists, new organism package can be generated. The necessary information for this is extracted from `AnnotationHub` using [NCBI taxonomyID](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi) `tax_id` (see below for example how to obtain list of organisms at AnnotationHub).

If you do not find the organism even at `AnnotationHub`, then I suggest that you use some other organism and accept having all-NA entries in `duplicates_annotation.txt`. You can do this e.g  by setting `org_db` to `org.Hs.eg.db`.

Conversion of gene names using biomaRt is optional and should be run only if your are not using Ensembl reference. Please consult [this manual](http://www.ensembl.info/2015/06/01/biomart-or-how-to-access-the-ensembl-data-from-r/) for information how to talk to Biomart from R (see also example below). `attrs` describes names of attributes that are to be fetched from biomaRt as annotations for genes, an needs to be supplied only if the conversion step is required.

**Appears in:**

- dset, biomart, host, in_id, out_id, attrs, tax_id, org_db, schema: `make_annodb`

``` yaml
cdf: # continued from above, # agilent specific
tax_id: "10090"
org_db: "org.Mmusculus.eg.db"
schema:  # use if org_db NOT installed, but schema available `AnnotationForge::available.chipdbschemas()`
# Optional, required only if `out_id` specified
dset: mmusculus_gene_ensembl
biomart: ensembl
host: www.ensembl.org
in_id: ensembl_gene_id
out_id:  # entrezgene #  no conversion done if empty
attrs: #["entrezgene","ensembl_transcript_id","description","external_gene_name"]

# or for plants:
cdf:
  dset: osativa_eg_gene
  biomart: plants_mart
  host: plants.ensembl.org
  [...]

```

#### Example of using biomaRt
```R
host <- "plants.ensembl.org"
marts <- biomaRt::listEnsembl(host = host)
dsets <- biomaRt::listDatasets(biomaRt::useMart(biomart=marts$biomart[1],host=host))
mart <- biomaRt::useMart(biomart=marts$biomart[1], host=host, dataset = dsets$dataset[1])
attrs <- biomaRt::listAttributes(mart, page = "feature_page")$name
```

#### Example of listing organisms available at AnnotationHub
```R
ah <- AnnotationHub::AnnotationHub()
ah <- subset(ah, ah$rdataclass == "OrgDb")
mc <- S4Vectors::mcols(ah)[, "taxonomyid", drop = FALSE]
AHID <- rownames(mc[mc$taxonomyid == taxid, , drop = FALSE])
if (!length(AHID))
    stop("no OrgDb package found for taxid ", taxid)
```

## Affymetrix

### Chip

The information in this section (`rows` and `cols`) is used in the creation of affymetrix CDF file. In particular the dimensions of the array must be obtained e.g. by looking at some already existing package or in the header of platform annotation file (usually deposited at GEO or Array express).

**Appears in:**

- rows, cols: `make_cdf`
- probe_length: `count_features`

``` yaml
chip:
  rows: 1190 # affy only
  cols: 1190 # affy only
  probe_length: 60 # # used ONLY in report
```

## Illumina
Not implemented, but should be a straightforward extension of Agilent.


# Contributing

Config should be the only part of the pipeline that needs to be adjusted for any respective experiment. Keep this in mind when making changes.
