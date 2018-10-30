# This is a stub and is not used elswhere.
# Was used in development to copare bowtie and blast alignments
rule blast_align:
    """
    Align probes to genomic reference using BLAST

    This is useful to debug implementation fo Bowtie alignment. With current settings
    the alignments are identical for the two. Note that SAM output of BLAST is
    buggy. For documentation of command line options, see
    https://www.ncbi.nlm.nih.gov/books/NBK279684/#_appendices_Options_for_the_commandline_a_

    -word_size ... length of initial exact match
    -perc-identity ... 100 means no mismatch
    -max_target_seqs,-evalue ... control how hard blast tries to find alignments
    -dust ... whether to mask low complexity parts of genome with DUST algorithm

    Notes
      First you need to prepare a db: `makeblastdb -dbtype nucl -in {input.index}`
    """
    input:
        index = rules.merge_ref.output,
        probes = rules.unique_probes.output.uniq
    output:
        sam =expand("{ali_dir}/{align_prefix}_blast.sam",
                    ali_dir = ali_dir, align_prefix = align_prefix)
    log:
        "logs/blast_align/{}.log".format(align_prefix)
    threads: 2
    params:
        strategy = join(ali_dir, align_prefix + "strategy")
    shell: # or use -db {input.index} instead of -subject
        "blastn -task 'blastn-short' -query {input.probes} -db {input.index} "
        "-word_size=25 -ungapped -out {output} " # "-outfmt '17 SQ SR' "
        "-outfmt '6 qseqid sstrand qlen sstart qseq sseq length nident mismatch sseqid stitle qcovs qcovhsp qcovus' "
        "-dust=no -perc_identity=100 -max_target_seqs=1000000 "
        "-evalue=10 -strand=both "
        "-num_threads {threads} &> {log}"

rule tab_to_sam:
    input:
        # not needed, using bowtie

rule clean_alignments:
    """
    Filters alignments to default length and single strand

    Occasionaly, few spurious alignments with length different from the probe length
    appear (only for BLAST). Furthermore, multiple probe sequences will be matched
    to both forward and reverse strand. This rule remedies both issues by
    dropping such alignments.

    The rule preserves original and appends '_all' tag to its filename.
    """
     # TODO: Dont use this rule
    # Length of alignmnets is ok with bowtie, strandedness should be handled by counter/mapper

    input:
        sam = rules.bowtie_align.output.sam,
        # prb = join(prb_dir, config["probes"]["file"]),
        prb = rules.unique_probes.output.uniq_gz
    output:
        # orig= expand("{ali_dir}/{align_prefix}_all.sam",
        #             ali_dir = ali_dir, align_prefix = align_prefix),
        filt= expand("{ali_dir}/{align_prefix}_filt.sam",
                    ali_dir = ali_dir, align_prefix = align_prefix)
    params:
        script = join(snake_dir, "scripts/filter_ali.sh"),
        length = config["chip"]["probe_length"]
    log:
        "logs/probes/{}_clean.log".format(align_prefix)
    run:
        # shell("cp {input.sam} {output.orig}"),
        shell("chmod u+x {params.script}"),
        # shell("bash {params.script} {input.sam} {params.length} {input.prb} &> {log}")
