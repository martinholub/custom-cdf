# some utility funcs and classes
from os.path import join
import re
from collections import OrderedDict

class CdfPipeException(Exception):
    def __init__(self, platform=None, message="CDF Pipeline Exception."):
        if platform:
            ex_message = "Platform {} not implemented.".format(platform)
        else:
            ex_message = message
        super().__init__(ex_message)

def uniq_list(array):
    """ Returns unique elements in list in original order

    set() does not preserve order. In python 3.5, OrderedDict is implemented in C,
    and thus acceptable alternative. For python 3.6+, it can be replaced with regular
    dict as it became ordered.

    References
      https://stackoverflow.com/a/39835527/6877443
    """
    return list(OrderedDict.fromkeys(array))


def normalize_cdf_name(chip_name, sense):
    """Regularize cdf name

    sense: str, any of "AT", "ST", "both"
    """
    cdf_name = chip_name.lower()
    cdf_name = re.sub("-(st|at)-", sense.lower(), cdf_name)
    cdf_name = cdf_name.replace("-", "_").replace("_", "") # must not have dashes
    return cdf_name

def output_cdf_selector(platform, cdf_dir, cdf_name):
    """Select output cdf or annodb based on platform"""

    if platform.startswith("affy"):
        cdf = join(cdf_dir, cdf_name + "cdf.cdf")
    elif platform.startswith(("agi", "illu")):
        cdf = join(cdf_dir, cdf_name + "db.db")
    else:
        raise CdfPipeException(platform)
    return cdf

def rm_ext(fname, ext = None):
    """Removes extension from filename"""
    if not ext:
        ext = re.match(".*\.(.*)$", fname).group(1)
    else:
        ext = ext.strip(".")
    return re.sub("\." + ext + "$","", fname)

def extract_script_selector(platform, snake_dir):
    """Select results exctraction script based on platform

    This is mainly for backwards compatibility when different scripts where used
    for different platforms.
    """
    if platform.startswith(("affy", "agi")):
        script = join(snake_dir, "scripts/extract_results_all.py")
    else:
        raise CdfPipeException(platform)
    return script

def do_merge_reference(config_ref):
    """Should the reference be merged?"""
    try:
        fname = config_ref["file"]
    except (KeyError, AttributeError) as e:
        raise e

    try:
        mname = config_ref["merged"]
        if mname == "": mname = None
    except (KeyError, AttributeError) as e:
        mname = None

    if mname is None or fname == mname:
        do_merge = False
    elif "*" in fname:
        do_merge = True
    else:
        do_merge = False
        # raise CdfPipeException(message = "Malformed reference filename ({})".format(fname))
    return do_merge

def select_merge_filename(dirname, config_ref):
    """Select filename for merged references"""
    if do_merge_reference(config_ref):
        fname = join(dirname, config_ref["merged"])
    else:
        fname = join(dirname, config_ref["file"])
    return fname


def rev_complement_seq(seq):
    """ reverse complements a sequence

    For ANY probe to hybridise it has to be the RC of cDNA.
    Does not work with any special characters in the sequence.

    #TODO: Assure that not confilciting with approach in filter_ali.sh
    """
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    # alphabet = set(["C","T","A","G"])
    assert all([i in complement.keys() for i in set(seq)])
    revseq = reversed(seq)
    revcompseq = "".join(complement.get(base) for base in revseq)
    return revcompseq

class CdfGeneMapperException(Exception):
    """Custom Exception class"""
    def __init__(self, message=None):
        def_message="CDF Gene Mapper Exception."
        if message:
            ex_message = def_message +" " + message
        else:
            ex_message = def_message
        super().__init__(ex_message)
