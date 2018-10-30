#!/bin/bash

# This scipts removes duplicate entries from list of probes on a chip formated as
# fasta file. The duplicates are often due to assignement of one probe to multiple genomic
# loci. Ideally, such duplication would not appear in this FASTA file, but once
# the authors decided to do it, we need to fix it.
#
# The script takes three arguments
# $1, fasta file with probe information including name and sequences
# $2, fasta file to write output to, the file will be additionaly gzipped
# $3, temporary file where counts, names and seqs of duplicates are written
#     this can be useful if we want to do sth with tese probes later
#
# Caution, this works only for particular format of probe name! This is why they
# should be normalized prior to running this script.

function probe_validator(){
    prb_file=$1
    out_file=$2
    dups_file=$3
    # MUST be kept harmonized with `normalized_probenames.py`
    prb_name_sub=(":[0-9]+-[0-9]+;" ";")

    echo "--- Filtering duplicate probes from input fasta ---"
    # Pull out all repeated probes
    awk -v r1=${prb_name_sub[0]} -v r2=${prb_name_sub[1]} 'BEGIN{{RS="(\\n|^)>";FS=" ";OFS="\t"}}NR>1{{gsub(r1,r2,$1); print $1,$NF}}' $prb_file |
    sort | uniq -d -c | sed  's/^[[:space:]]*\([0-9]\+\)[[:space:]]/\1\t/' > "${dups_file}"
    f1=$(mktemp)
    cut -f3 $dups_file > $f1

    # Retain only probes that do not appear multiple times
    echo "$(awk 'BEGIN{{RS="(\\n|^)>"}}NR>1{{sub("\\n","\t"); gsub("\\n",""); print ">"$0}}' $prb_file)" | grep -v -Fwf $f1 |
    awk -F"\t" -v OFS="\\n" '{{print $1,$2}}' > "${out_file}"

    # Add 1 ocurence for each nonunique probe
    # This approach is valid, because we care  only about probe identity (its ID) and its sequence.
    echo "$(awk 'BEGIN{{RS="(\\n|^)>"}}NR>1{{sub("\\n","\t"); gsub("\\n",""); print ">"$0}}' $prb_file)" | grep -Fwf $f1 |
    awk -F"\t" -v OFS="\\n" '!seen[$2]++ {{print $1,$2}}' >> "${out_file}"
    rm ${f1}

    # The above works mostly like:
    # awk 'BEGIN{{RS=">"}}NR>1{{sub("\\n","\\t"); gsub("\\n",""); print RS$0}}' $prb_file |
    # awk -F"\\t" '!seen[$2]++' | awk -F"\\t" -v OFS="\\n" '{{print $1,$2}}' > $output_file
    # but will not break if the probe fasta will have identical sequences at differetn places on chip
    echo "--- Filtering finished ---"
}

# probe_validator $prb_file $out_file $dupl_file
probe_validator "$@"

# Report on number of lines

echo -n "Total number of probes: " && num_lines=$(wc -l $1 | cut -f1 -d" ") &&
echo $(($num_lines/2)) && echo -n "Number of unique probes: " &&
num_lines=$(wc -l $2 | cut -f1 -d" ") && echo $(($num_lines/2))

# Gzip the probes
gzip --keep --force $2
