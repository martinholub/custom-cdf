import os
import re
from copy import deepcopy

class CdfResultsExtractor(object):
    """Class to extract results from hits dictionary

    TODO:  Allow different inputs than dict (e.g. json)
    """
    def __init__(self, res_path, chip_type, platform, hits, seq_lookup, lookup = {}):
        self.chip_type = chip_type
        self.res_path = res_path
        self.platform = platform
        self.hits = hits
        self.seq_lookup = seq_lookup
        self.lookup = lookup

    @property
    def chip_type(self):
        return self._chip_type

    @chip_type.setter
    def chip_type(self, value):
        assert isinstance(value, (str, ))
        self._chip_type = value

    @property
    def res_path(self):
        return self._res_path

    @res_path.setter
    def res_path(self, value):
        value = value.rstrip("/")
        value = os.path.expanduser(value)
        parent = os.path.dirname(value)
        if not os.path.isdir(parent):
            os.mkdir(parent) # only this level, os.makedirs -> also parents
        assert os.path.isdir(parent), "Folder {} doesn't exist".format(parent)
        if os.path.isfile(value):
            print("File at {} will be overwritten.".format(value))
        self._res_path = value

    @property
    def platform(self):
        return self._platform

    @platform.setter
    def platform(self, value):
        if isinstance(value, (list, tuple, )): value = value[0]
        assert isinstance(value, (str, ))
        value = value.lower()
        is_good = False
        for tst in ["affy", "illu", "agi"]:
            is_good = is_good | value.startswith(tst)
        assert is_good, "platform {} not supported".format(value)

        self._platform = value

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

    def save_results(self):
        """ Save results of mapping and filtering

        File formated so as to correspond to R CDF constructors. Must be kept in
        correspondence with approach for extraction of genomic alignment results.

        "For expression probe sets there is 1 and only 1 group per probe set.",
        "The term unit is an internal term which means probe set."
        => group == unit == probeset(see References).

        Parameters
        ------------
        chip_type: str
            name of chip, as present in Affymetrix probes FASTA file
        hits: dict
            {metafeature: {"probes" = [probes], "nhits": len([probes])}
        platform: str, {"affy", "agilent"}
        seq_lookup: dict
            {"alignment": sequence} 1-to-1 mapping

        See also: `normalize_probenames`

        References
        ----------
            http://dept.stat.lsa.umich.edu/~kshedden/Courses/Stat545/Notes/AffxFileFormats/cdf.html
            https://www.affymetrix.com/support/developer/powertools/changelog/gcos-agcc/cdf.html

        """
        # AffySpecific code
        # `*.?` drops information on the probe name, which is fine
        probe_id_re = r"probe:" + re.escape(self.chip_type) + r":.*?([0-9]+)-([0-9]+)$"
        probe_coord_re = r"^([0-9]+):([0-9]+)"

        strand_map = {"fw": "sense", "rc": "antisense"}

        with open(self.res_path, 'w') as wf:
            if self.platform.startswith("affy"): # for affy you need header
                header_out = [  "probe_id", "probe_x", "probe_y", "probe_seq",
                                "group_id", "unit_id", "strand", "n_hits"]
                wf.write("\t".join(header_out) + "\n")

            hits_deepcopy = deepcopy(self.hits)
            for metafeaturestrand, val in hits_deepcopy.items():
                metafeature, strand = metafeaturestrand.split(":")

                # TODO: how to deal with genes that have + and - probeset?
                # would need to adapt group and unit ID, but currently not used

                # Not used, see docs of this functon
                # try: # get gene corresponding to the transcript
                #     unit, _ , _ = self.lookup[metafeature]
                # except KeyError as e: # metafeature is gene already
                #     unit = metafeature

                unit = metafeature

                # Probe/Target strand information: just for affy consistency

                if strand not in ("sense", "antisense"):
                    strand = strand_map[strand] # strand of probe is strand of target.
                hit_count = val["nhits"]

                for probe in val["probes"]:

                    if self.platform.startswith("affy"):
                        probe_seq = self.seq_lookup[probe]
                        probe_info = probe.split(";")
                        probe_id = re.search(probe_id_re, probe_info[0]).group(1)
                        probe_x, probe_y = re.findall(probe_coord_re, probe_info[1])[0]
                        wf.write("\t".join([probe_id, probe_x, probe_y, probe_seq,
                                            metafeature, unit, strand,
                                            str(hit_count)]) + "\n")

                    elif self.platform.startswith("agi"):
                        wf.write("\t".join([probe, metafeature]) + "\n")
