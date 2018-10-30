# TODO:  At the end, verify 0 or 1 based indexing (BED has zero)

## This is an archive of stuff that will be deleted soon.
# rule mask_reference:
#     input:
#         ref = ancient(join(ref_dir, config["reference"]["file"])),
#         bed = expand("{var_dir}/{var_file}.gz")
#     output:
#         expand("{ref_dir}/{prefix}.hm.fa.gz", ref_dir = ref_dir, prefix = prefix)
#     shell:
#         """
#         bedtools maskfasta -fi {input.ref} -bed {input.bed} -fo {output}
        # """
        # # will fail due to memory requirements :() --- good now
        # # bedtools maskfasta -fi data/reference/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz -bed data/variants/arabidopsis_thaliana.vcf.gz -fo data/reference/AT_tair10_toplevel.hm.fa
        # # other option is `faidx mygenome.fa --bed BED --out outgenome.fa -v` from https://github.com/mdshw5/pyfaidx

# rule bcftools_call:
#     input:
#         ref = join(ref_dir, config["reference"]["file"] + ".gz"),
#         bai = expand("{ali_dir}/{{align_prefix}}.bam.bai", ali_dir = ali_dir),
#         bam = expand("{ali_dir}/{{align_prefix}}.bam", ali_dir = ali_dir)
#     output:
#         "data/variants/{align_prefix}.vcf"
#     params:
#         pmr = config["prior_mutation_rate"]
#     threads: 2
#     log: # Log files must contain exactly the same wildcards as the output
#         "logs/bcftools/{align_prefix}.log"
#     shell:
#         "(samtools mpileup -g -f {input.ref} {input.bam} | "
#         "bcftools call -mv -P '{params.pmr}' -o {output} --threads {threads} "
#         "- ) &> {log}"
#
# rule snp_filter:
#     input:
#         #...
#     output:
#         #...
#     shell:
#         """
#         bcftools query -f '%INFO/TSA' data/variants/arabidopsis_thaliana.vcf.gz | head
#         bcftools view --include 'INFO/TSA="SNV"' --known --output-type z --output-file data/variants/test.vcf.gz --threads 2 data/variants/arabidopsis_thaliana.vcf.gz
#         bcftools annotate -a data/annotation/Arabidopsis_thaliana.TAIR10.39_sorted.gtf.gz -include 'INFO/TSA="SNV"' --columns CHROM,-,FROM,TO,-,-,-,- -x '^ID,INFO/TSA'-o data/variants.test.bed -O v data/variants/arabidopsis_thaliana.vcf.gz
#
#         bedtools sort -i data/annotation/Arabidopsis_thaliana.TAIR10.39.gtf.gz > data/annotation/Arabidopsis_thaliana.TAIR10.39_sorted.gtf
#         bgzip data/annotation/Arabidopsis_thaliana.TAIR10.39.gtf
#
#         tabix -p gff --force --print-header data/annotation/Arabidopsis_thaliana.TAIR10.39_sorted.gtf.gz
#         """
#
#         ## This resulting can be used for subsequent allignements!
#         ## Problem is that bowtie! does not support ambiguity in genome reference!!
#         """
#         # Convert to BED format
#         convert2bed -i vcf --snvs < data/variants/arabidopsis_thaliana.vcf > data/variants/test.vcf.bed
#         # Keep only necessary columns, drop tmp_.* IDs (what do they stand for anyway?)
#         awk -F"\t" -v OFS="\t" '$4 !~ /^tmp/ {print $1,$2,$3,$4 }' data/variants/test.vcf.bed > data/variants/test_filter.vcf.bed
#         # bgzip and index
#         bgzip data/variants/test_filter.vcf.bed
#         tabix -b2 -e3 -s4 -p bed data/variants/test_filter.vcf.bed.gz
#
#         # Mask genomic loci with known (or potentionally select) SNPs
#         bedtools maskfasta -fi data/reference/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa -bed data/variants/test_filter.vcf.bed.gz -fo data/reference/AraGene-1_1-st-v1.snp.fa
#         # Bgzip and index
#         bgzip data/reference/AraGene-1_1-st-v1.snp.fa
#         samtools faidx data/reference/AraGene-1_1-st-v1.snp.fa
#         """
#         # This approach is nice, but probabdly does not allow us to specify loci of mismatch wrt to probe!
#         """
#         # Subtract SNP locations from genome
#         bedtools intersect -v -a data/annotation/Arabidopsis_thaliana.TAIR10.39_sorted.gtf.gz -b data/variants/test_filter.vcf.bed.gz  > data/annotation/Arabidopsis_thaliana.TAIR10.39.sorted.noSNP.gtf
#         ## this is the same
#         # bedtools subtract -A -a data/annotation/Arabidopsis_thaliana.TAIR10.39_sorted.gtf.gz -b data/variants/test_filter.vcf.bed.gz  > data/annotation/Arabidopsis_thaliana.TAIR10.39.sorted.noSNP2.gtf
#         bgzip data/annotation/Arabidopsis_thaliana.TAIR10.39.sorted.noSNP.gtf
#         tabix data/annotation/Arabidopsis_thaliana.TAIR10.39.sorted.noSNP.gtf.gz
#         """
#         # Ultimatelly, I will have to filter the reads, probably with awk
#         """
#         bedtools bamtobed -i data/alignement/AraGene-1_1-st-v1.bam > data/alignement/AraGene-1_1-st-v1.bed
#         # sort, bgzip, index
#         ## Best shot is probably altering the colums start/stop, intersecting and altering back
#         ### Later may do sth like '{ out=$1, $2+halfshift, $3-halfshift; for(i=4;i<=NF;i++){out=out,$i}; print out}' but need to get it right
#         awk -F"\t" -v halfshift="7" -v OFS="\t" '{ print $1, $2+halfshift, $3-halfshift, $4, $5, $6}' data/alignement/AraGene-1_1-st-v1.bed > data/alignement/AraGene-1_1-st-v1.center.bed
#         bedtools intersect -v -a data/alignement/AraGene-1_1-st-v1.center.bed -b data/variants/test_filter.vcf.bed.gz > data/alignement/AraGene-1_1-st-v1.center.noSNP.bed
#         awk -F"\t" -v halfshift="7" -v OFS="\t" '{ print $1, $2-halfshift, $3+halfshift, $4, $5, $6}' data/alignement/AraGene-1_1-st-v1.center.noSNP.bed > data/alignement/AraGene-1_1-st-v1.noSNP.bed
#         bedtools bedtobam -i data/alignement/AraGene-1_1-st-v1.noSNP.bed -g data/reference/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz > data/alignement/AraGene-1_1-st-v1.noSNP.bam
#         # make genome file
#         bowtie-inspect --summary data/index/AT_tair10_toplevel | awk -F"\t" -v OFS="\t" 'NR>3 { sub(" .*","",$2); print $2, $NF}' > data/reference/AT_tair10_toplevel.genome
#         bedtools bedtobam -i data/alignement/AraGene-1_1-st-v1.noSNP.bed -g data/reference/AT_tair10_toplevel.genome > data/alignement/AraGene-1_1-st-v1.noSNP.bam
#         # sort, index
#         # with this you should be able to proceeed into feature counts
#         # unfortunately not, it does not have sequence :D
#         # but probably can do:
#         bedtools intersect -F 1.0 -a data/alignement/AraGene-1_1-st-v1.bam -b data/alignement/AraGene-1_1-st-v1.noSNP.bed > data/alignement/AraGene-1_1-st-v1.noSNP.bam
#         index samtools index -b data/alignement/AraGene-1_1-st-v1.noSNP.bam data/alignement/AraGene-1_1-st-v1.noSNP.bam.bai
#         """
#
# include: "rules/bwa_map.smk"
# include: "rules/report.smk"

# Results similar as with featureCounts but less convenient and slower
# rule count_features_htseq:
#     input:
#         anno = ancient(join(ann_dir, config["annotation"]["file"] + ".gz")),
#         # bam = rules.snp_filter.output,
#         # bai = rules.samtools_index.output,
#         sam = rules.clean_alignments.output.filt
#     output:
#         samout = expand("{cnt_dir}/{align_prefix}_htseq.sam",
#                         cnt_dir = cnt_dir, align_prefix = align_prefix)
#     threads: 2
#     params:
#         feature_type = "transcript"
#     log:
#         "logs/counts/{}_htseq.log".format(align_prefix)
#     shell:
#         "htseq-count {input.sam} {input.anno} "
#         "-f sam -t {params.feature_type} -i gene_id "
#         "--additional-attr transcript_id " # or transcript_id
#         "--nonunique all --secondary-alignments score "
#         "--supplementary-alignments score --mode union --stranded yes " # ?? strandednes
#         "-o {output.samout} &> {log}" # prints counts to log -.-
