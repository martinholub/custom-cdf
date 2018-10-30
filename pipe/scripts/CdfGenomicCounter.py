from collections import OrderedDict, Counter
import re
from itertools import islice, chain
import os
import json
import sys
from time import strftime, localtime
import numpy as np
import pandas as pd
from tqdm import tqdm as TQDM
from copy import deepcopy
from cdf_utils import CdfGeneMapperException, rev_complement_seq, uniq_list

version = strftime('%Y-%m-%d %H:%M:%S', localtime(os.path.getmtime(sys.argv[0])))

class CdfGenomicCounter(object):
    """Class of an CDF feature hit counter

    * Logging is handled by snakemake by redirecting sys.stdout output to file.
    """
    def __init__(self, cnt_file, hits_file, sense, metafeature, min_probeset_size,
                control_probes, multimatch_filter, verbose = 0):
        # paths
        self.cnt_file = cnt_file
        self.hits_file = hits_file
        self.metafeature = metafeature
        # params
        self.sense = sense
        self.min_probeset_size = min_probeset_size
        # default init
        self.verbose = verbose
        # Objects, # TODO: Add property
        self.multimatch_filter = multimatch_filter


    @property
    def cnt_file(self):
        """Path to count SAM file"""
        return self._cnt_file

    @cnt_file.setter
    def cnt_file(self, value):
        value = value.rstrip("/")
        value = os.path.expanduser(value)
        assert os.path.isfile(value), "File {} doesn't exist".format(value)
        self._cnt_file = value

    @property
    def hits_file(self):
        """Path to hits as reported by featureCounts"""
        return self._hits_file

    @hits_file.setter
    def hits_file(self, value):
        value = value.rstrip("/")
        value = os.path.expanduser(value)
        assert os.path.isfile(value), "File {} doesn't exist".format(value)
        self._hits_file = value

    @property
    def sense(self):
        return self._sense

    @sense.setter
    def sense(self, value):
        assert isinstance(value, (str, ))
        value = value.lower()
        assert value in ("at","st","both")
        self._sense = value

    @property
    def metafeature(self):
        return self._metafeature

    @metafeature.setter
    def metafeature(self, value):
        assert isinstance(value, (str, ))
        value = value.lower()
        assert value in ("gene", "transcript")
        self._metafeature = value

    @property
    def min_probeset_size(self):
        return self._min_probeset_size

    @min_probeset_size.setter
    def min_probeset_size(self, value):
        assert isinstance(value, (int, ))
        self._min_probeset_size = value

    @property
    def control_probes(self):
        return self._control_probes

    @control_probes.setter
    def control_probes(self, value):
        if not isinstance(value, (list,)): value = [value]
        assert all(isinstance(x, (str, )) for x in value)
        self._control_probes = set(value)

    @property
    def verbose(self):
        return self._verbose

    @verbose.setter
    def verbose(self, value):
        assert isinstance(value, (bool, int, ))
        value = int(value)
        self._verbose = value

    def get_counts(self):
        """Gets counts as reported by featureCounts in SAM file

        The alignments are stored in form {transcript: {probes: [probes], strand: [strands]}}.
        Assumptions are made on the columns and information therein as SAM is standardized.

        Parameters
        --------------
        control_probes: list of strings
            Begginings of probenames that are controls
        cnt_file: str, path to counts SAM file produced by featureCounts
        hits_file: str, additional annotation for counts from featureCounts
        sense: str, {"at", "st", "both"}
            sense of the array
        verbose: bool, verbosity level

        Returns
        --------------
        blacklist: set
            set of alignments that should be dropped from output
        probes: set
            set of all alignments seen in the SAM file
        seq_lookup: dict
            {alignment: sequence} 1-to-1 mapping
        alignments: dict
            {transcript: {probes: [probes], strand: [strands]}} mapping
        """
        alignments = OrderedDict()
        strand_map = {  "16": "antisense", "0": "sense",
                        "-": "antisense", "+": "sense"} # SAM flag correspondence
        all_probes = []
        targets = []
        other_strand_probes = []
        seq_lookup = {}
        control_probes = {}
        control_probes_list = list(self.control_probes)

        unassigned_count = 0
        num_headlines = 0
        filtered_out_count = 0
        other_strand_count = 0

        target_hits = pd.read_table(self.hits_file, skiprows = 1, usecols = [0, 4, 6],
                                    skip_blank_lines = False, header = 0, index_col = 0,
                                    names = ["ID", "strand", "counts"],
                                    dtype = {"ID":str, "strand":str, "counts":int})

        target_hits = target_hits.to_dict(orient = "index")
        target_hits = {k: (next(iter(set(filter(bool, v["strand"].split(";"))))),
                            v["counts"]) for k,v in target_hits.items()}

        with open(self.cnt_file, 'r') as f:

            for line in f:
                # or equivalently line.startswith("#")
                if re.match("^@[HCSP][QDGO]", line):
                    num_headlines += 1
                    continue  # header line

                try:
                    status = re.search(r"\tXS:Z:(.+?)(\t|\n)", line).group(1)
                except AttributeError:
                        print("XS:Z:<status> field not found on line {}".format(i))
                        status = "Unknown"

                if status == "Assigned":
                    target = re.search(r"\tXT:Z:(.+?)(\t|\n)", line).group(1)
                else: # This probe is not matched to any metafeature
                    unassigned_count += 1
                    continue

                if "," in target: continue # multimapping probe

                line = line.split("\t")
                probe = line[0] # sam is standardized, can make assumptions
                probe_strand = strand_map[line[1]]

                seq = line[9]

                # Some probes are only control probes and come under generic names
                # This makes them hard to pass on, so we drop them here
                # Happens for v2.0+ of GeneChips
                # if probe in self.control_probes:
                # This is more robust but potentially slow
                if any(filter(bool, map(lambda x: re.search(":"+x, probe, re.IGNORECASE),
                                        control_probes_list))):
                    if probe not in control_probes.keys():
                        control_probes[probe] = [seq]
                    else:
                        control_probes[probe].extend([seq])
                    continue

                strand = strand_map[target_hits[target][0]] # assumption

                # TODO(mholub): Evaluate if equivalent with approach in Transcriptomic Counting
                # -s 2 -> probe_strand == "antisense", only `16` alignments are retained
                # -s 1 -> probe_strand == "sense", only `0` alignments are retained
                # At first sight it seems that the present approach is correct.
                if self.sense == "st" and probe_strand == "antisense":
                        seq = rev_complement_seq(seq)
                elif self.sense == "at" and probe_strand == "sense":
                        seq = rev_complement_seq(seq)
                elif self.sense == "both":
                    if probe_strand == "antisense":
                        seq = rev_complement_seq(seq)
                else:
                    other_strand_count += 1
                    other_strand_probes.append(probe)

                all_probes.append(probe)
                targets.extend(target)

                if probe not in seq_lookup.keys():
                    seq_lookup.update({probe: seq})
                else: # probes appearing multiple times have same seqs
                    if not seq_lookup[probe] == seq:
                        message = "Multiple sequences for probe '{}''".format(probe)
                        raise CdfGeneMapperException(message) # or just assert

                if target not in alignments.keys():
                    alignments.update({target: {"probes": [probe],
                                                "strands":[strand]}})
                else:
                    alignments[target]["probes"].append(probe)
                    alignments[target]["strands"].append(strand)


        # Antisense Probes (_at) for Sense Array and Vice Versa
        # See config readme for discussion on strandedness
        self.blacklist = set(other_strand_probes)

        # List of unique probes appearing in the SAM file
        # unless --no-unal, this can include also probes without alignment
        self.probes = set(all_probes)
        # All the unique targets we have some matching probe for
        self.targets = set(targets)
        # {probe: seq} lookup
        self.seq_lookup = seq_lookup
        # alignments
        self.alignments = alignments

        if self.verbose > 0:
            # Count how many hits we have per transcript
            template = "Saw {} targets ({} unique) hit by in total {} probes ({} unique)"
            print(template.format(  len(targets),
                                    len(alignments.keys()), # keys are always unique
                                    len(all_probes),
                                    len(self.probes)))

            # Count occurences of probe instancess
            probe_counter = Counter(all_probes)
            template = "Number of hits from multimapping probes: {}"
            print(template.format(len([k for k,v in probe_counter.items() if v>1])))

            template = "Number of all probes: {} ({} unique)"
            print(template.format(len(all_probes), len(set(all_probes))))

            print("Number of unassigned reads (probes): {}.".format(unassigned_count))
            print("Number of probes from other strand: {}.".format(other_strand_count))
            print("Number of filtered reads (probes): {}.".format(filtered_out_count))

            control_probes_seqs = list(chain.from_iterable(control_probes.values()))
            template = "{} unique control probes (of {} types) with {} hits to reference."
            print(template.format(  len(set(control_probes_seqs)),
                                    len(control_probes), len(control_probes_seqs)))
            template = "Types of contol probes: {}"
            print(template.format(list(control_probes.keys())))

        if self.verbose > 1:
            print("Example of data: ")
            print(dict(islice(alignments.items(), 1)))

            print("The most multimapping probes: ")
            print(probe_counter.most_common(3))

    def count_metafeature_hits(self):
        """ count read hits to metafeatures

        `metafeature` is controlled by featureCounts. Usually this will be metafeature,
        but one can also opt to summarize on level of `features` with option `-f`.
        Here we call it simple `metafeature`, whatever it is.

        Both are tracked simultaneously in the script as this allows to
        filter multimatches based on either gene or transcript assignments
        (`multimatch_level` in config).

        Parameters
        ----------------
        sense: str, {"at", "st", "both"}
            sense of the array
        verbose: bool, verbosity level
        min_probeset_size: int,
            minimum number of unique probe hits per metafeature for it to be retained
        metafeature: str,  {"gene", "transcript"}
        lookup: dict
                {transcript: (gene, strand, chromosome)} for all transcripts
        blacklist: set
            set of alignments that should be dropped from output
        probes: set
            set of all alignments seen in the SAM file
        seq_lookup: dict
            {alignment: sequence} 1-to-1 mapping
        alignments: dict
            {transcript: {probes: [probes], strand: [strands]}} mapping
        multimatch_filter: CdfMultimatchFilter
            object implementing `read_whitelist` and `get_multimatching` methods

        Returns
        -------------------
        genes: set, set of metafeatures corresponding to all mapped features
        hits: dict
            {metafeature: {"probes" = [probes], "nhits": len([probes])}

        """
        # B - Make filter for probes matching more (meta)features --------------
        metafeature_hits = deepcopy(self.alignments)
        if len(self.multimatch_filter.multimatch_level) > 0:
            self.multimatch_filter.hits = metafeature_hits
            multimatching = self.multimatch_filter.get_multimatching()
        else:
            multimatching = set()
        self.blacklist.update(multimatching) # update blacklist set with multimatching

        if self.verbose > 0: # Debug level messages
            template = "There are {} multimatching probes."
            print(template.format(len(multimatching)))

            template = "{} unique probes were blacklisted ({} will be dicarded from results)"
            print(template.format(len(self.blacklist), len(self.blacklist.intersection(self.probes))))

        # C - walk through meta(features) and apply filtering criteria ------------------
        base_hits = deepcopy(metafeature_hits)

        msg = "Filtering {}-specific probesets ...".format(self.metafeature)
        hits = {}
        tqdm_params = { 'unit': 'blocks',
                        'unit_scale': True,
                        'leave': False,
                        'miniters': 1,
                        'total': len(base_hits),
                        'desc': msg}
        with TQDM(**tqdm_params) as tqdm: # progbar
            for base, val in base_hits.items():
                tqdm.update(1)
                try:
                    # self.blacklist must be set
                    filt = [False if p in self.blacklist else True for p in val["probes"]]
                    probes = np.array(val["probes"])[filt]
                    strands = np.array(val["strands"])[filt]
                except Exception as e:
                    message = "No probes found for {} '{}''".format(self.metafeature, base)
                    raise CdfGeneMapperException(message) from e

                if len(set(strands)) > 1:
                    if self.sense != "both":
                        message = "Multiple strands for {} '{}''".format(self.metafeature, base)
                        raise CdfGeneMapperException(message)
                    else:
                        # Or do strand specific probesets; not used unless sense == "both"
                        assert len(set(strands)) == 2
                        fw, rc = sorted(set(strands))
                        probesf = probes[strands == fw]
                        probesr = probes[strands == rc]

                        for pr, st in zip([probesf, probesr], [fw, rc]):
                            if len(set(pr)) >= self.min_probeset_size: # threshold
                                hits[base+":"+st] = {"probes" : uniq_list(pr),
                                                    "nhits": len(pr)}
                else:
                    if len(set(probes)) >= self.min_probeset_size: # threshold
                        hits[base+":"+strands[0]] = {"probes": uniq_list(probes),
                                                    "nhits": len(probes)}

                # if i > 500: break # DEBUG

        # metafeature (probeset) hits by probes in format
        # {"probeset:strand" : {probes: set(probes), nhits: len(probes)}}
        self.hits = hits # TODO: add @property

        if self.verbose > 0:

            nhits_all = sum([float(v["nhits"]) for v in hits.values()])
            print("Sum of all hits = {:d}".format(int(nhits_all)))

            template = "{} assigned probes ({} unique) out of {} hits (from {} unique probes)"
            all_good = list(chain.from_iterable([v["probes"] for v in hits.values()]))
            sense_probes = list(chain.from_iterable(map(lambda v: v['probes'], self.alignments.values())))
            print(template.format(  len(all_good), len(set(all_good)),
                                    len(sense_probes),len(self.probes)))

        if self.verbose > 1:
            print("Example of a metafeature w/ fewest hits:")
            print(dict(islice(sorted(hits.items(), key =lambda x: x[1]["nhits"]), 3)))

            print("Metafeature w/ the most hits: ")
            print(dict(islice(sorted(hits.items(), key =lambda x: x[1]["nhits"], reverse = True), 1)))
