import subprocess
import tempfile
import os
import timeit
import sys
import re
from cdf_utils import rm_ext


class TranscriptomicSNP(object):
    """ Generate a VCF file with transcriptomic coordinates

    Gets transcriptomic snp loci in format <transc> <SNPpos> <additional_info...>.

    # available:
    # <chrom> <transc> <start> <end>
    # <chrom> <SNVpos>

    """
    def __init__(self, anno_path, vcf_path, out_path, do_add_transcript_version = True):
        self.anno_path = anno_path
        self.vcf_path = vcf_path
        self.out_path = out_path
        self.do_add_transcript_version = do_add_transcript_version

    @property
    def anno_path(self):
        """Path to annotation GTF file"""
        return self._anno_path

    @anno_path.setter
    def anno_path(self, value):
        value = value.rstrip("/")
        value = os.path.expanduser(value)
        assert os.path.isfile(value), "File {} doesn't exist".format(value)
        self._anno_path = value

    @property
    def vcf_path(self):
        """Path to annotation VCF file"""
        return self._vcf_path

    @vcf_path.setter
    def vcf_path(self, value):
        value = value.rstrip("/")
        value = os.path.expanduser(value)
        assert os.path.isfile(value), "File {} doesn't exist".format(value)
        self._vcf_path = value

    @property
    def out_path(self):
        return self._out_path

    @out_path.setter
    def out_path(self, value):
        value = rm_ext(value, "gz")
        value = value.rstrip("/")
        value = os.path.expanduser(value)
        self._out_path = value

    def get_transcript_snvs(self):
        """Generate a VCF file with transcriptomic coordinates

        The generated file is passed to `snp_filter`. We make some assumptions on
        the GTF and VCF file, if anything goes wrong, double check if they apply
        also on your input.

        We use only ready-made tools. The motivation was compactness, saving time
        on development and increased performance. The last metric does not seem to be met,
        as the code execution is rarther slow on big files.

        Notes:
            - Reimplementation in Python would allow profiling and optimizing performance.
            - I use `shell=True` and pass single string as command as otherwise not working
            (despite it should). This can be looked into.

        Refernces:
            https://docs.python.org/3/library/subprocess.html
            https://bedops.readthedocs.io/en/latest/content/reference/file-management/conversion/convert2bed.html#convert2bed
            http://bedtools.readthedocs.io/en/latest/content/tools/intersect.html
        """

        # Extract transcripts from GTF file, sort
        _, transcripts_gtf = tempfile.mkstemp()
        command = 'grep -P "\ttranscript\t" {} | sort -k1,1 -k4,4n > {}'.format(self.anno_path,transcripts_gtf)
        _ = subprocess.run(command, check = True, shell = True)

        # Get only SNVs from VCF, to save some memory/time
        _, snv_vcf = tempfile.mkstemp()
        command = 'grep -P "(^#|TSA=SNV)" {} > {}'.format(self.vcf_path, snv_vcf)
        subprocess.run(command, check = True, shell = True)

        # Get such transcripts, that have an overlap with a known SNV
        # Exploits fact that both coordinates are chromosome-based
        # writes information from both files next to each other
        _, transcripts_vcf = tempfile.mkstemp()
        command = "bedtools intersect -sorted -F 1 -wo -a {} -b {} > {}"
        command = command.format(transcripts_gtf, snv_vcf, transcripts_vcf)
        subprocess.run(command, check = True, shell = True)

        # Pull out annotation of transcripts as obtained in previous step
        _, transcripts = tempfile.mkstemp()
        if self.do_add_transcript_version:
            command = "cut -f9 {} | sed -E 's/.*\stranscript_id \"([^;]+)\"; "
            command += "transcript_version \"([^;]+)\";.*/\\1.\\2/' > {}"
        else:
            command = "cut -f9 {} | sed -E 's/.*\stranscript_id \"([^;]+)\".*/\\1/' > {}"
        command = command.format(transcripts_vcf, transcripts)
        subprocess.run(command, check = True, shell = True)

        # Combine transcript name with the SNV information
        # Converts the coordinates to transcriptomic ones

        # `$11-$4+1` computes coordinate of SNV along transcript,where
        # $4 start of transcript and $11 postion of SNV both in chromosme-based coords.
        # transcripts file carries information on transcript
        # columns $12..$17 carry information on the SNV itself
        command = "awk -F\"\t\" -v OFS=\"\t\" '{{print $11-$4+1,$12,$13,$14,$15,$16,$17}}' {} | paste {} - > {}"
        command = command.format(transcripts_vcf, transcripts, self.out_path)
        subprocess.run(command, shell = True, check = True)

        # Zip and index
        command = "bgzip --force {0} && tabix -p vcf {0}.gz".format(self.out_path)
        subprocess.run(command, shell = True, check = True)

        # Cleanup
        for f in [transcripts_gtf, snv_vcf, transcripts, transcripts_vcf]:
            os.remove(f)

def main(snakemake):
    """Execution Wrapper

    Note: Call pdb as `import pdb; pdb.Pdb(stdout = sys.__stdout__).set_trace()
        https://docs.python.org/3/library/pdb.html#pdb.Pdb
    """
    if not re.match(".*\.gtf$", snakemake.input[0]): # deal with unpredictable order
        snakemake.input = snakemake.input[::-1]

    # Preparation
    logfile = snakemake.log[0]
    if logfile:
        orig_stdout = sys.stdout
        lf = open (logfile, mode = "w")
        sys.stdout = lf

    start_time = timeit.default_timer()

    # Execution
    tsnp= TranscriptomicSNP(anno_path = snakemake.input[0],
                            vcf_path = snakemake.input[1],
                            out_path = snakemake.output[0],
                            do_add_transcript_version = snakemake.config["counts"]["do_add_transcript_version"])
    tsnp.get_transcript_snvs()

    # Cleanup
    end_time = timeit.default_timer()
    print("Finished in {:.3f} s".format(end_time - start_time))

    if logfile:
        sys_stdout = orig_stdout
        lf.close()

if __name__ == '__main__':
    main(snakemake)
