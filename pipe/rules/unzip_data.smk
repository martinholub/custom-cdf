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
