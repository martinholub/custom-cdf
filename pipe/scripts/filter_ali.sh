%%bash
#!/bin/bash

# This filters alignments to have all the expected lenght (should be always the case for bowtie)
# and drops probes aligning to both strands
# aditionaly, include probes that droped for nonuniqueness?

# Inputs ----------------------------------------------------------------------
# SAM file with alignments
al=$1
# Lenght of a read
len=$2
# Probes as fasta
prbs=$3
# Params -----------------------------------------------------------------------
# Line where alignment records start from
start_line=$(grep -m1 -n -E "^[^@]" $al | cut -f1 -d:)
# No of column that has sequences (tolerate shorter length in the beginning). Should be 10
# col_num=$(awk -v sl="$start_line" -v len="$len" -F"\t" 'NR==sl {for(i=1;i<=NF;i++){if($i~/[ACTG]{len-5,}/) {print i;exit;}}}')
col_num=10

# Get Sequences ----------------------------------------------------------------
# Sequences of all alignments
seq=$(awk -v sl="$start_line" -v cn="$col_num" -F"\t" 'NR>=sl {print $cn}' $al | sort)
# They should be all the same length
echo "--- Alignment length distribution ---"
echo "$seq" | awk '{ print length($0); }' | sort | uniq -c

# Get only good alignments
al_ok=$(awk -v sl="$start_line" -v cn="$col_num" -v len="$len" -F"\t" 'NR>=sl && length($cn)==len {print $0}' $al)
# Get only good sequences
# seq=$(awk -v sl="$start_line" -v cn="$col_num" -F"\t" 'NR>=sl {print $cn}' $al_ok | sort)

echo "--- Bitwise flags ---"
# Should be just 0 and 16
echo "$al_ok" | awk -F"\t" '{print $2}' | sort | uniq -c
# Forward Strand Alignments
fw_seqs=$(echo "$al_ok" | awk -v cn="$col_num" -F"\t" '$2==0 {print $cn}')
# Reverse Strand Alignments
rc_seqs=$(echo "$al_ok" | awk -v cn="$col_num" -F"\t" '$2==16 {print $cn}')
echo "--- Total number of alignments (fw+rc) ---"
echo $(echo "$fw_seqs" | wc -l)"+"$(echo "$rc_seqs" | wc -l) | bc


# Filter seqs aligning to both strands -----------------------------------------
f1=$(mktemp)
f2=$(mktemp)
echo "$fw_seqs" | sort | uniq > "$f1"
# Reads aligned to reverse strand were reverse complemented. Get back their originals.
echo "$rc_seqs" | rev | tr ATCG TAGC | sort | uniq > "$f2"
# Some probes align to both strands
fwrc_seqs=$(comm -12 $f1 $f2)
echo "--- Probes aligned to both strands (to be dropped) ---"
echo "$fwrc_seqs" | wc -l
rm ${f1} ${f2}

# Drop promiscous probes. Awk may allow column specificity but not needed for SAM format.
f3=$(mktemp)
echo "$fwrc_seqs" > "$f3"
al_fil=$(echo "$al_ok" | grep -v -Fwf $f3)
rm ${f3}
echo "--- Number of alignments after fw/rc filtering ---"
echo "$al_fil" | wc -l

# DONE
# Augment by seqs that you have dropped due to duplication ---------------------
dup_seqs=$(awk 'NR%2==0' $prbs | sort | uniq -d -c )

## TODO(mholub): pickup here by augmenting al_fil by repeating given line as many
## times as you have duplicates and changing the first field accordingly!!
## will need to pick up NR%2==1 probes for that as well

# Paste header from input file and combine with filtered probes
al_out=$(printf "$(awk -v sl="$start_line" -F"\t" 'NR<sl' $al)\n%s" "$al_fil")
# Overwrite file that came is input
outf=$(echo "$al" | sed 's/\.sam$/_filt\.sam/')
echo "$al_out" > "$outf"
