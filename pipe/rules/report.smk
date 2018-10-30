
rule report_sub:# need to name the rule differently than in parent Snakemake file
    input:
        T1 = "calls/all.vcf",
        T2 = expand("benchmarks/{sample}.bwa_map.txt", sample=config["samples"])
    output:
        "report.html"
    script:
        # The script path is always relative to the Snakefile (in contrast to the input and output file paths, which are relative to the working directory)
        "../scripts/report.py"
