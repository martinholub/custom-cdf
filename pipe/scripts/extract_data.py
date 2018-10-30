import sys
import timeit
import re
sys.path.insert(0, "./pipe/scripts") # or from pipe.scripts.Module import Class
from CdfGenomicCounter import CdfGenomicCounter
from CdfTranscriptomicCounter import FeatureHitCounter
from CdfResultsExtractor import CdfResultsExtractor
from CdfMultimatchFilter import CdfMultimatchFilter

VERBOSITY = 3 # TODO: Do not hardcode verbosity

def main(snakemake):
    """ Execution wrapper for Feature Hit Mapping
    """
    if not re.match(".*\.sam$", snakemake.input[0]): # deal with unpredictable order
        snakemake.input = snakemake.input[::-1]

    logfile = snakemake.log[0]
    if logfile:
        orig_stdout = sys.stdout
        lf = open (logfile, mode = "w")
        sys.stdout = lf

    start_time = timeit.default_timer()

    cfm = CdfMultimatchFilter(
        fpath = snakemake.config["cdf"]["whitelist"],
        multimatch_level = snakemake.config["cdf"]["multimatch_level"]
    )

    if snakemake.config["reference"]["is_transcriptome"]:
        fhc = FeatureHitCounter(
            ali_path = snakemake.input[0],
            anno_path= snakemake.input[1],
            sense = snakemake.config["cdf"]["sense"],
            feature = snakemake.config["counts"]["feature"],
            metafeature = snakemake.config["counts"]["metafeature"],
            min_probeset_size = snakemake.config["cdf"]["min_probeset_size"],
            control_probes = snakemake.config["counts"]["controls"],
            multimatch_filter = cfm,
            do_add_transcript_version = snakemake.config["counts"]["do_add_transcript_version"],
            verbose = VERBOSITY
        )
        fhc.make_lookup()
        fhc.get_alignments()
        fhc.count_feature_hits()

    else:
        cgc = CdfGenomicCounter(
            cnt_file = snakemake.input[0], # sys.argv[1]
            hits_file = snakemake.input[1],
            metafeature=snakemake.config["counts"]["metafeature"],
            min_probeset_size = snakemake.config["cdf"]["min_probeset_size"],
            sense = snakemake.config["cdf"]["sense"].lower(),
            control_probes = snakemake.config["counts"]["controls"],
            multimatch_filter = cfm,
            verbose = VERBOSITY
        )
        cgc.get_counts()
        cgc.count_metafeature_hits()
        cgc.lookup = {}
        fhc = cgc

    cre = CdfResultsExtractor(
        chip_type= snakemake.params[0],
        platform = snakemake.config["platform"],
        res_path = snakemake.output[0],
        hits = fhc.hits,
        seq_lookup = fhc.seq_lookup,
        lookup = fhc.lookup
    )
    cre.save_results()

    end_time = timeit.default_timer()
    print("Finished in {:.3f} s".format(end_time - start_time))

    if logfile:
        sys_stdout = orig_stdout
        lf.close()

if __name__ == '__main__':
    main(snakemake)
