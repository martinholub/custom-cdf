import re
import sys
import tempfile
import shutil
import os
import timeit
from cdf_utils import CdfPipeException

class ProbeNameNormalizer(object):
    """ProbeNameNormalizer Object class

    Normalizes probes' names from FASTA headers to common format.

    For Affymetrix:
        probe:<PREFIX>:<probe_id>-<cluster_id>;<probe_x>:<probe_y>;
    For Agilent:
        <probe_id>
    """
    def __init__(self, fin, fout, chip, platform):
        self.fin = fin
        self.fout = fout
        self.chip = chip
        self.platform = platform
        # TODO: Add basic checks

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

    def normalize_probenames(self):
        """Normalizes probenames from probe.fa file to expected format

        This needs to be extended if new format is encountered. The remaining fields
        of FASTA header (everything after first space) are kept unchanged.

        The expected format of probenames is:
            AFFY: ">probe:<chip>:<probeid>-<transcriptclusterid>;<probeX>:<probeY>;"
            AGI:  "><probename>"

        TODO: Add reference for using temporary file.
        """
        # valid for affy
        probe_re = r">probe:" + re.escape(self.chip)
        case1 = probe_re + r":([0-9]+)(-[0-9]+)$" # probeid-transcriptcluster
        case2 = probe_re + r":([0-9]+)$" # probeid, occurs for control probes

        _ , ftemp = tempfile.mkstemp()
        i = 0 # lines_come_in_pairs_checker
        with open(self.fin, 'r') as rf, open(ftemp, 'w') as wf:
            for line in rf:
                if line.startswith(">"): # valid for all platforms
                    i += 1
                    if self.platform.startswith("affy"):
                        probe, *rest = line.split(" ")
                        probe, coords = probe.split(";")[:2]# drops empty string
                        if coords: # two occurences of ";" in probe
                            if re.match(case1, probe):
                                pass
                            elif re.match(case2, probe):
                                probe += "-00000000" # add mock transcriptcluster
                        else:
                            probe, chip, *info = probe.split(":")
                            coords = ":".join(info[1:]) # get coords

                            # Construct some sort of unique probe ID
                            ## This keeps the original name but is lengthy
                            probe_id = re.sub("_[ast]{2}$", "", info[0]) # remove sense identifier
                            probe_id += "___{}___".format(chip) # add clearly separated chape name
                            # construct id from (x,y) coordinates
                            probe_id += "".join([x.rjust(4,"0") for x in coords.split(":")])

                            ## This uses coordinates, to come up with sth like ID and is just a number
                            # probe_id = "".join([x.rjust(4,"0") for x in coords.split(":")])
                            probe = ":".join([probe, chip, probe_id]) + "-00000000"

                        probe = ";".join([probe, coords]) + ";"
                        wf.write(" ".join([probe, *rest]))

                    elif self.platform.startswith("agi"):
                        wf.write(line) # for agilent, don't do anything
                    else:
                        msg = "Probe Fasta normalization not implemented for {} platform".format(self.platform)
                        raise NotImplementedError(msg)
                else:
                    # Note: old affy chips may have oposite notion of 'sense'
                    # One way how to address it would be to re_complement the seq here.
                    # It is probably better idea to do it later and more easily.
                    i = 0 # reset counter
                    wf.write(line)

                if i > 1:
                    msg = "More than 1 sequence for seen on line: {}".format(line)
                    raise CdfPipeException(message = msg)

        # Copy over and clean up
        shutil.copyfile(ftemp, self.fout)
        os.remove(ftemp)

def main(snakemake):
    """Execution wrapper"""

    logfile = snakemake.log[0]
    if logfile:
        orig_stdout = sys.stdout
        lf = open (logfile, mode = "w")
        sys.stdout = lf

    start_time = timeit.default_timer()

    pnn = ProbeNameNormalizer(  snakemake.input[0],
                                snakemake.output[0],
                                snakemake.params[0],
                                snakemake.config["platform"])
    pnn.normalize_probenames()
    end_time = timeit.default_timer()
    print("Finished in {:.3f} s".format(end_time - start_time))

    if logfile:
        sys_stdout = orig_stdout
        lf.close()

if __name__ == "__main__":
    main(snakemake)
