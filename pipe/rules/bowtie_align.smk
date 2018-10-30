# neednt to refer to config anymore

rule bwa_map_sub:
    input:
        "data/genome.fa",
        lambda wildcards: config["samples"][wildcards.sample]
    output:
        temp("mapped_reads/{sample}.bam") # dispose not useful to unclutter
    params:
        # -R STR : read group header line such as '@RG\tID:foo\tSM:bar'
        rg = "@RG\tID:{sample}\tSM:{sample}"
    log:
        "logs/bwa_map/{sample}.log"
    benchmark:
        "benchmarks/{sample}.bwa_map.txt" # if wildcards, must be same as output
    threads: 2 # also use snakemake --cores
    shell:
        "(bwa mem -t {threads} {input} | " # -R '{params.rg}'
        "samtools view -Sb - > {output}) 2> {log} 1>&2"
