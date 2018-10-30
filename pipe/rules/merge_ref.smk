from os.path import join

rule merge_ref:
    """
    Merges multiple reference FASTA files into single file

    It is simpler to work with single reference file. However, occasionally we may
    want to align only to select subset of chromosomes or selct part of genome.
    In such case, the select files are placed in `reference` directory
    and this rule merges them into a single file.

    The file is not zipped is at is needed later.
    """
    input:
        rules.unzip_data.output.ref
    output:
        join(ref_dir, config["reference"]["merged"])
    shell:
        "cat {input} > {output}"
