from collections import OrderedDict, Counter
import re
from itertools import islice, chain
import os
import json
import sys
from time import strftime, localtime
import numpy as np
from tqdm import tqdm as TQDM
import random
from copy import deepcopy
from cdf_utils import CdfGeneMapperException, rev_complement_seq, uniq_list

version = strftime('%Y-%m-%d %H:%M:%S', localtime(os.path.getmtime(sys.argv[0])))
# this parameter should be kept fixed as its functionality is mostly covered
# by changing `sense` in *config*
IS_OLD_CHIP = False

class FeatureHitCounter(object):
    """Class of an CDF feature hit counter

    * Logging is handled by snakemake by redirecting sys.stdout output to file.

    TODO: Combine common functionality with CdfGenomicCounter (!!!)
    """
    def __init__(self, anno_path, ali_path, sense, feature, metafeature,
                min_probeset_size, control_probes, multimatch_filter,
                do_add_transcript_version = True, verbose = 0):
        # paths
        self.anno_path = anno_path
        self.ali_path = ali_path
        # params
        self.sense = sense
        self.min_probeset_size = min_probeset_size
        self.feature = feature
        self.metafeature = metafeature
        self.control_probes = control_probes
        self.do_add_transcript_version = do_add_transcript_version
        # default init
        self.verbose = verbose
        self.lookup = {}
        self.seq_lookup = {}
        self.alignments = {}
        self.probes = set()
        self.transcripts = set()
        self.probes = set()
        # Objects, # TODO: Add property
        self.multimatch_filter = multimatch_filter

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
    def ali_path(self):
        """Path to SAM alignments"""
        return self._ali_path

    @ali_path.setter
    def ali_path(self, value):
        value = value.rstrip("/")
        value = os.path.expanduser(value)
        assert os.path.isfile(value), "File {} doesn't exist".format(value)
        self._ali_path = value

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
    def control_probes(self):
        return self._control_probes

    @control_probes.setter
    def control_probes(self, value):
        if not isinstance(value, (list,)): value = [value]
        value = list(filter(bool, value)) # drop None and ""
        assert all(isinstance(x, (str, )) for x in value) # passes on empty list
        self._control_probes = set(value)

    @property
    def feature(self):
        return self._feature

    @feature.setter
    def feature(self, value):
        assert isinstance(value, (str, ))
        value = value.lower()
        assert value in ("transcript")
        self._feature = value

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
    def verbose(self):
        return self._verbose

    @verbose.setter
    def verbose(self, value):
        assert isinstance(value, (bool, int, ))
        value = int(value)
        self._verbose = value

    @property
    def do_add_transcript_version(self):
        return self._do_add_transcript_version

    @do_add_transcript_version.setter
    def do_add_transcript_version(self, value):
        assert isinstance(value, (bool, int, ))
        value = bool(value)
        self._do_add_transcript_version = value

    @property
    def probes(self):
        return self._probes

    @probes.setter
    def probes(self, value):
        assert isinstance(value, (set, ))
        self._probes = value

    @property
    def genes(self):
        return self._genes

    @genes.setter
    def genes(self, value):
        assert isinstance(value, (set, ))
        self._genes = value

    @property
    def transcripts(self):
        return self._transcripts

    @transcripts.setter
    def transcripts(self, value):
        assert isinstance(value, (set, ))
        self._transcripts = value

    @property
    def alignments(self):
        return self._alignments

    @alignments.setter
    def alignments(self, value):
        assert isinstance(value, (dict, ))
        self._alignments = value

    @property
    def seq_lookup(self):
        return self._seq_lookup

    @seq_lookup.setter
    def seq_lookup(self, value):
        assert isinstance(value, (dict, ))
        self._seq_lookup = value

    @property
    def lookup(self):
        return self._lookup

    @lookup.setter
    def lookup(self, value):
        assert isinstance(value, (dict, ))
        self._lookup = value

    def make_lookup(self):
        """Constructs a lookup dictionary {transcript: (gene, strand, chromosome)}

        You decide if lookup includes only some chrs, nonchromosmal loci, scaffolds by
        selection of the annotation file. The lookup makes sense on level
        transcript <-> gene because we are aligning to cDNA (i.e. transcripts).

        Parameters
        ----------------
        anno_path: str, path to annotation file
        do_add_transcript_version: bool, add version number to transcript name?
        verbose: bool, verbosity level

        Returns
        ----------------
        lookup: dict
            {transcript: (gene, strand, chromosome)} for all transcripts
        """
        save_path = os.path.join(os.path.dirname(self.anno_path), "known_isoforms.json")
        # Lookup found
        if os.path.isfile(save_path):
            print("Found lookup at {}. Loading.".format(save_path))
            with open(save_path, 'r') as f:
                lookup = json.load(f)

        # Create lookup
        else:
            lookup = OrderedDict()

            with open(self.anno_path, 'r') as f:
                for line in f:
                    if line.startswith("#"): continue # a header line

                    if re.search("\ttranscript\t", line): # this line decribes transcript

                        # Extract information on transcripts
                        pattern = re.compile("^gene_id .*")
                        line = line.split("\t")
                        annot = next(filter(pattern.match, line)).split(";")
                        annot = [x.strip() for x in annot]

                        pattern = re.compile("gene_id \"(.*)\"")
                        gene_id = next(filter(pattern.match, annot))
                        gene_id = re.match(pattern, gene_id).group(1)

                        pattern = re.compile("transcript_id \"(.*)\"")
                        transc_id = next(filter(pattern.match, annot))
                        transc_id = re.match(pattern, transc_id).group(1)

                        # could also pick e.g. exon_id

                        chrom = line[0] # chromosome # or name

                        # Extract strandedness of the transcript
                        pattern = re.compile("^(\+|-)$")
                        try:
                            strand = next(filter(pattern.match, line))
                        except StopIteration as e:
                            strand = "." # SAM spec for empty field

                        # If transcriptomic reference comes with transcript versioning, pick it
                        # ALTERNATIVE: drop it here, as well as (potentially) in alignments
                        if self.do_add_transcript_version:
                            try:
                                pattern = re.compile("transcript_version \"(.*)\"")
                                transc_v = next(filter(pattern.match, annot))
                                transc_v = re.match(pattern, transc_v).group(1)
                                transc_id = "{}.{}".format(transc_id, transc_v)
                            except:
                                pass # no information on version

                        if not transc_id in lookup.keys():
                            lookup.update({transc_id : (gene_id, strand, chrom)})

                        else: # should not happen
                            msg = "Transcript {} already in dict.".format(transc_id)
                            raise CdfGeneMapperException(msg)
                            # lookup[transc_id].append((gene_id, strand, chrom))

            with open(save_path, 'w') as f:
                json.dump(lookup, f, sort_keys = True, indent = 4)
                print("Saved lookup dict to {}".format(save_path))

        if self.verbose > 0: # mainly debug level information
            template = "{} transcripts for {} genes"
            print(template.format(len(lookup.keys()),
                                 len(set(map(lambda x: x[0], lookup.values())))))
        if self.verbose > 1:
            print("Example of data:")
            print(dict(islice(lookup.items(), 3)))

        self.lookup = lookup

    def get_alignments(self):
        """Gets alignments from SAM file

        The alignments are stored in form {transcript: {probes: [probes], strand: [strands]}}.
        Assumptions are made on the columns of input and information therein as
        SAM is standardized.

        Parameters
        --------------
        control_probes: list of strings
            Begginings of probenames that are controls
            # TODO: pass patterns instead of begginings of strings
        ali_path: str, path to alignment SAM file
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


        References:
            http://bowtie-bio.sourceforge.net/manual.shtml#sam-bowtie-output
        """
        alignments = OrderedDict()

        # Some old chips appear to have reverse namig convention
        # keep False, # this can be replaced by `setting {cdf: {sense: "AT"}}`
        # Only concern is: Do you need IS_OLD_CHIP and "ST" to have also strandedness
        # correct in results? This needs to be checked, but will not infulence CDF/ChipDB,
        # only appear as inconsistency.
        if IS_OLD_CHIP:
            strand_map = {"0": "rc", "16": "fw"}
        else:
            strand_map = {"16": "rc", "0": "fw"} # SAM flag correspondence

        sense_probes = []
        transcripts = []
        seq_lookup = {}
        antisense_probes = []
        control_probes = {}
        control_probes_list = list(self.control_probes)

        with open(self.ali_path, 'r') as f:
            for line in f:
                if re.match("^@[HCSP][QDGO]", line): continue  # header line

                line = line.split("\t")
                probe = line[0] # sam is standardized, can make assumptions
                 # strandedness of ALIGNMENT (controlled by --nofw, --norc, independent from transcript strand)
                strand = strand_map[line[1]]
                transcript = line[2]
                seq = line[9]

                # Some probes are only control probes and come under generic names
                # This makes them hard to pass on, so we drop them here
                # Happens for v2.0+ of GeneChips

                # if probe in self.control_probes:
                # This is much more robust but may be slow
                if any(filter(bool, map(lambda x: re.search(":"+x, probe, re.IGNORECASE),
                                        control_probes_list))):
                    if probe not in control_probes.keys():
                        control_probes[probe] = [seq]
                    else:
                        control_probes[probe].extend([seq])
                    continue

                if strand == "fw":
                    # probes aligned to cDNA w/out RC == targets antisense
                    sense_probes.append(probe)
                    if self.sense == "st": continue
                elif strand == "rc":
                    # get sequence as appearing in probes FASTA file, bacues bowtie does RC on these
                    seq = rev_complement_seq(seq)
                    antisense_probes.append(probe)
                    if self.sense == "at": continue
                else:
                    message = "Unknown strandendes for probe {} ({})".format(probe, strand)
                    raise CdfGeneMapperException(message)

                transcripts.append(transcript)

                if probe not in seq_lookup.keys():
                    seq_lookup.update({probe: seq})
                else: # probes appearing multiple times have same seqs
                    if not seq_lookup[probe] == seq:
                        message = "Multiple sequences for probe '{}''".format(probe)
                        raise CdfGeneMapperException(message) # or just assert

                if transcript not in alignments.keys():
                    alignments.update({transcript: {"probes": [probe],
                                                    "strands":[strand]}})
                else:
                    alignments[transcript]["probes"].append(probe)
                    alignments[transcript]["strands"].append(strand)


        # Antisense Probes (_at) for Sense Array and Vice Versa
        # See config readme for discussion on strandedness
        if self.sense == "both":
            all_probes = sense_probes + antisense_probes
            self.blacklist = set()
        # Per BA protocol, we drop probes that match to both strands for single-stranded cdfs
        elif self.sense == "st":
            all_probes = antisense_probes
            self.blacklist = set(sense_probes)
        elif self.sense == "at":
            all_probes = sense_probes
            self.blacklist = set(antisense_probes)

        # List of unique probes appearing in the SAM file
        # unless --no-unal, this can include also probes without alignment
        self.probes = set(all_probes)
        # All the unique transcripts we have some matching probe for
        self.transcripts = set(transcripts)
        # {probe: seq} lookup
        self.seq_lookup = seq_lookup
        # alignments
        self.alignments = alignments

        if self.verbose > 0:
            # Count how many hits we have per transcript
            template = "Saw {} transcripts ({} unique) hit by in total {} probes ({} unique)"
            print(template.format(  len(transcripts),
                                    len(alignments.keys()), # keys are always unique
                                    len(all_probes),
                                    len(self.probes)))
            # check with
            # grep "^probe" ZebGene-1_1-st-v1.sam | awk -FS="\t" '{print $1}' | cut -f1 | sort | uniq | wc -l

            # Count occurences of probe instancess
            probe_counter = Counter(all_probes)
            template = "Number of hits from multimapping probes: {}"
            print(template.format(len([k for k,v in probe_counter.items() if v>1])))

            template = "Number of antisense probes: {} ({} unique)"
            print(template.format(len(antisense_probes), len(set(antisense_probes))))

            template = "Number of sense probes: {} ({} unique)"
            print(template.format(len(sense_probes), len(set(sense_probes))))

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

    def count_feature_hits(self):
        """ count read hits to (meta)features

        `meta` is metafeature and can be either 'gene' or 'transcript'. 'feature'
        is always 'transcript'. This allows us to create transcript specific
        CDF files by using `transcript` as metafeature.

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
        feature: str, {"gene", "transcript"}
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
        # A - Map transcript hits to features --------------------------------------
        meta_hits = OrderedDict() # used for filtering
        feature_hits = OrderedDict()
        unknown_features = {}
        metas = []
        features = []
        nonchrom_metas = []
        alignments_deepcopy = deepcopy(self.alignments)

        for feature, val in alignments_deepcopy.items(): # meta == gene
            try:
                meta, strand, chrom = self.lookup[feature]
            except KeyError as e:
                unknown_features.update({feature: val})
                continue

            if not re.match("^[0-9]{1,2}", chrom):
                # Used to disallow discarding probes that have multiple genes they map to
                # but the other match is on nonchromosomal gene.
                nonchrom_metas.append(meta)

            metas.append(meta)
            features.append(feature)

            ali_dict = {"probes": val["probes"],
                        "strands": val["strands"]}

            # Probes assignment to features
            if feature not in feature_hits.keys():
                feature_hits[feature]= deepcopy(ali_dict)
                # Could pick up some stats here
                # No. of transcripts, hits per trancript, ...
            else:
                msg = "Duplicate keys in alignment for {}".format(feature)
                raise CdfGeneMapperException(msg) # should not occur

            # Metafeature, often gene
            if meta not in meta_hits.keys():
                meta_hits[meta] = deepcopy(ali_dict)
            else:
                meta_hits[meta]["probes"] += val["probes"]
                meta_hits[meta]["strands"] += val["strands"]

        self.genes = set(metas)
        assert set(features).issubset(self.transcripts)
        # B - Make filter for probes matching more (meta)features ----------------------
        # Features to be whitelisted must be included in whitelist_file which is
        # read on init of the multi_match object

        # here you can also whitelist matches to nonchromosmal metafeatures
        # self.multimatch_filter.whitelist.update(set(nonchrom_metas))

        if self.multimatch_filter.multimatch_level == "gene":
            self.multimatch_filter.hits = meta_hits
            multimatching = self.multimatch_filter.get_multimatching()
        elif self.multimatch_filter.multimatch_level == "transcript":
            self.multimatch_filter.hits = feature_hits
            multimatching = self.multimatch_filter.get_multimatching()
        else:
            multimatching = set()

        self.blacklist.update(multimatching) # update blacklist set with multimatching

        if self.verbose > 0: # Debug level messages
            feat_with_probes = {k for k,v in feature_hits.items() if v["probes"]}
            template = "{} {}s (out of which {} unique) have some probes assigned"
            print(template.format(  len(feat_with_probes), self.feature,
                                    len(self.transcripts)))

            feat_with_multistrand = \
                {k for k,v in feature_hits.items() if len(set(v["strands"])) > 1}
            template = "{} {}s both strands. Should be 0."
            print(template.format(len(feat_with_multistrand), self.feature))

            meta_with_probes = {k for k,v in meta_hits.items() if v["probes"]}
            template = "{} genes (out of which {} unique) have some probes assigned"
            print(template.format(len(meta_with_probes), len(self.genes)))

            meta_with_multistrand = \
                {k for k,v in meta_hits.items() if len(set(v["strands"])) > 1}
            template = "{} genes seen on both strands. Should be 0."
            print(template.format(len(meta_with_multistrand)))

            try:
                template = "{} {}s are not known by annotation file (example: {})"
                print(  template.format(len(unknown_features), self.feature,
                        random.choice(list(unknown_features.keys()))))
            except IndexError as e:
                template="All {} from alignment were matched to genes."
                print(template.format(self.feature))

            template = "There are {} multimatching probes."
            print(template.format(len(multimatching)))

            template = "{} unique probes were blacklisted ({} will be dicarded from results)"
            print(template.format(len(self.blacklist), len(self.blacklist.intersection(self.probes))))

            template = "{} unique nonchrom. features"
            print(template.format(len(set(nonchrom_metas))))

        # C - walk through meta(features) and apply filtering criteria ------------------
        if self.metafeature == "transcript":
            base_hits = deepcopy(feature_hits)
        elif self.metafeature == "gene":
            base_hits = deepcopy(meta_hits)

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
                    # self.blacklist must be set for performance reasons
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

        # Feature (probeset) hits by probes in format (TODO: add @property)
        # {"probeset:strand" : {probes: set(probes), nhits: len(probes)}}

        # self.hits = dict(sorted(hits.items(), key = lambda kv: kv[0].lower()))
        self.hits = hits # do not sort(!) to keep order as per SAM, which is order along (c)DNA

        if self.verbose > 0:
            template = "{} metafeatures passing criteria (out of {} unique features)"
            if self.metafeature == "gene":
                print(template.format(len(hits), len(self.genes)))
            elif self.metafeature == "transcript":
                print(template.format(len(hits), len(self.transcripts)))

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
