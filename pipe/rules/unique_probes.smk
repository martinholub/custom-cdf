from os.path import join

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
