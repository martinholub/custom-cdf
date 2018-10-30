#!/bin/bash

# Parse Arguments --------------------------------------------------------------
# reference: https://stackoverflow.com/a/14203146
POSITIONAL=()
while [[ $# -gt 0 ]]; do
  key="$1"

  case $key in
      -v|--vcf)
      vcf="$2"
      shift # past argument
      shift # past value
      ;;
      -b|--bam)
      bam="$2"
      shift # past argument
      shift # past value
      ;;
      -o|--output)
      bam_noSNP="$2"
      shift # past argument
      shift # past value
      ;;
      -c|--clip)
      clip="$2"
      shift # past argument
      shift # past value
      ;;
      *)    # unknown option
      POSITIONAL+=("$1") # save it in an array for later
      shift # past argument
      ;;
  esac
done
set -- "${POSITIONAL[@]}" # restore positional parameters

echo vcf = "${vcf}"
echo bam = "${bam}"
echo output = "${bam_noSNP}"
echo clip = "${clip}"

# tmp_file="/tmp/tmp.bed.gz"
tmp_file="$(mktemp tmp.XXXXXXXXXX.bed.gz)"
# tmp_bam="/tmp/tmp.bam"
tmp_bam="$(mktemp tmp.XXXXXXXXXX.bam)"

# Exctract SNVs to BED file ----------------------------------------------------
echo Creating BED file of SNVs ...
# Convert to bed (--snvs -> only SNV), drop some SNVs (named tmp*), keep some colums, ...
# sort (MUST to save memory), zip and index (tmp_file is bed.gz file)
# example: yes n | gunzip --stdout ../variants/Danio_rerio.vcf.gz | convert2bed -i vcf --snvs - | awk -F"\t" -v OFS="\t" '$4 !~ /^tmp/ {print $1,$2,$3,$4 }' | sort -k 1,1 -k2,2n | bgzip -c > tmp.bed.gz &&  tabix -b2 -e3 -s4 -p bed tmp.bed.gz
yes n | gunzip --stdout "${vcf}" | convert2bed -i vcf --snvs - | \
awk -F"\t" -v OFS="\t" '$4 !~ /^tmp/ {print $1,$2,$3,$4 }' | \
sort -k 1,1 -k2,2n | \
bgzip -c > "${tmp_file}" &&  tabix -b2 -e3 -s4 -p bed "${tmp_file}"

# Filter alignments based on SNP loci ------------------------------------------
# Convert alignements to bed, sort to save memory, clip to central part, ...
# subtract ailgnements that have SNV locus in central part, ...
# undo clipping, keep only such original bam entries that satisify these conditions

# -v... reverse
# -F 1.0 ... minimum overlap required as fraction of B
# TODO(mholub): be smarter about adding all columns (ranges with awk?)
# examples:
# bedtools bamtobed -i ZebGene-1_1-st-v1.bam | sort -k 1,1 -k2,2n | awk -F"\t" -v halfshift=5 -v OFS="\t" '{ print $1, $2+halfshift, $3-halfshift, $4, $5, $6}' > test.bed
# bedtools intersect -v -F 1.0 -sorted -a test.bed -b tmp.bed.gz > test_is.bed
# awk -F"\t" -v halfshift=5 -v OFS="\t" '{ print $1, $2-halfshift, $3+halfshift, $4, $5, $6}' test_is.bed | bedtools intersect -F 1.0 -a ZebGene-1_1-st-v1.bam -b stdin > test_is_clean.bam

# see also: http://bedtools.readthedocs.io/en/latest/content/tools/slop.html
# Could also write list <probe> <SNV_ID>

if [[ $bam == *.bam ]]; then
  echo Filtering BAM file "${bam}"...
  # bedtools bamtobed -i "${bam}" | awk -F"\t" -v OFS="\t" '{sub("\\.[0-9]+", "",$1); print $0}'| sort -k 1,1 -k2,2n | \
  bedtools bamtobed -i "${bam}" | sort -k 1,1 -k2,2n | \
  awk -F"\t" -v halfshift="${clip}" -v OFS="\t" '{ print $1, $2+halfshift, $3-halfshift, $4, $5, $6}' | \
  bedtools intersect -v -F 1.0 -sorted -a stdin -b "${tmp_file}" | \
  awk -F"\t" -v halfshift="${clip}" -v OFS="\t" '{ print $1, $2-halfshift, $3+halfshift, $4, $5, $6}' | \
  bedtools intersect -F 1.0 -a "${bam}" -b stdin > "${tmp_bam}"

  # Change SAM header to inform about filtering ----------------------------------
  echo Adjusting File Header ...
  # Add @CO field to header, indicating that the file was changed
  samtools view -H "${tmp_bam}" | echo -e  "$(cat -)\n@CO\tfilter=snvs" | \
  samtools reheader --no-PG - "${tmp_bam}" > "${bam_noSNP}"
else
  echo Filtering GTF file "${bam}" ... NOT IMPLEMENTED!
fi

# Housekeeping -----------------------------------------------------------------
echo Written file "${bam_noSNP}"
echo Cleaning up...
# Remove tmp file and its index
rm "${tmp_file}"*
rm "${tmp_bam}"
echo Done.
