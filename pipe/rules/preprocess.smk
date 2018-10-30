from os.path import join
import glob

rule unzip_data:
    """
    Unzip input reference and probe data

    Some tools may require input to be unzipped. Bowtie is one of them and thus
    its inputs are unzipped here. `ancient` flag prevents snakemake for checking
    for changes in file's timestamp. This is useful in developement, but should be removed
    in production.
    """
    input:
        ancient(expand("{path}.gz", path = paths)),
        # ancient(join(ref_dir, "{ref}.gz"))
        ancient(glob.glob(join(ref_dir, config["reference"]["file"] + ".gz")))
    output:
        other = temp(expand("{path}", path = paths)),
        ref = temp([re.sub("\.gz$","",x) for x in \
                    glob.glob(join(ref_dir, config["reference"]["file"] + ".gz"))])
    shell:
        """
        for f in {input}
        do
            # mkdir --parents $f
            gunzip --keep $f
        done
        """

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

rule unique_probes:
    """
    Remove duplicate occurences of probes in chip FASTA file

    This rule removes duplicate entries from list of probes on a chip formated as
    fasta file. The duplicates are often due to assignement of one probe to multiple genomic
    loci (ideally, such duplication would not appear in probes FASTA file).

    Notes:
      https://www.biostars.org/p/143617/
    """
    input:
        ancient(join(prb_dir, config["probes"]["file"])),
    priority: 5 # to make execute before align
    output:
        uniq = join(prb_dir, rm_ext(config["probes"]["file"]) + "_uniq.fa"),
        uniq_gz = join(prb_dir, rm_ext(config["probes"]["file"]) + "_uniq.fa.gz"),
        dups = join(prb_dir, rm_ext(config["probes"]["file"]) + ".dups.txt"),
    log:
        "logs/probes/{}.log".format(align_prefix)
    params:
        script = join(snake_dir, "scripts/filter_probes.sh"),
    # 1. Braces not used for variable access have to be escaped by repeating them
    # 2. Escape special characters in quotes, e.g. "\\n"
    run:
        shell("chmod u+x {params.script}"),
        shell("bash {params.script} {input} {output.uniq} {output.dups} &> {log}")
