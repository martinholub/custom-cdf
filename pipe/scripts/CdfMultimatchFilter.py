import os
from collections import Counter
from itertools import chain

class CdfMultimatchFilter(object):
    def __init__(self, multimatch_level, fpath = None, hits = {}, whitelist = set()):
        self.fpath = fpath
        self.hits = hits
        self.whitelist = whitelist
        self.multimatch_level = multimatch_level
        self.read_whitelist()

    @property
    def fpath(self):
        """Path to filter file"""
        return self._fpath

    @fpath.setter
    def fpath(self, value):
        if value not in (None, ""):
            print("Feature MM whitelist provided: {}".format(value))
            value = value.rstrip("/")
            value = os.path.expanduser(value)
            assert os.path.isfile(value), "File {} doesn't exist".format(value)
        else:
            value = None
        self._fpath = value

    @property
    def hits(self):
        return self._hits

    @hits.setter
    def hits(self, value):
        assert isinstance(value, (dict, ))
        self._hits = value

    @property
    def whitelist(self):
        return self._whitelist

    @whitelist.setter
    def whitelist(self, value):
        if isinstance(value, (list, )): value = set(value)
        assert isinstance(value, (set, ))
        self._whitelist = value

    @property
    def multimatch_level(self):
        return self._multimatch_level

    @multimatch_level.setter
    def multimatch_level(self, value):
        if value not in (None, ""):
            assert isinstance(value, (str, ))
            value = value.lower()
            try:
                assert value.startswith(("ge", "tr"))
                if value.startswith("ge"):
                    value = "gene"
                if value.startswith("tr"):
                    value = "transcript"
            except AssertionError as e:
                value = None
        else:
            value = None
        self._multimatch_level = value

    def read_whitelist(self):
        """ Reads features that should not be considered for MM filter from file

        File is expected to be one-column listing of features

        Parameters:
        ---------------
        fpath: str, ""
            path to one-column file listing features not be considered for multimap filtering

        Returns
        -------------
        whitelist: set
            Set of features not to be considered as targets for multimapping probes.
        """
        content = set()

        if self.fpath:
            print("Reading whitelisted features for MM from {}".format(self.fpath))
            with open(self.fpath, 'r') as ff:
                for line in ff:

                    if not line.strip() or line.startswith("#"): # empty line or comment
                        continue
                    line = line.rstrip() # Remove new line
                    line = line.strip('"').strip("'") # remove quotes if any
                    content.update([line])

            print("{} whitelisted features loaded".format(len(content)))

        self.whitelist.update(content)

    def get_multimatching(self):
        """ Gets probes that match to multiple features

        Feature is a key in hits dictoniary. Some features can be excluded
        from this criterion by including them in whitelist

        Parameters
        ------------
        hits: dict
            Hits of probes per features

        Returns
        -------------
        multimatching: set
            set of probes having more than one distinct feature as target
        """
        probe_filter = {k: list(set(v["probes"])) for k,v in self.hits.items() \
                        if k not in self.whitelist}

        probe_filter = Counter(chain.from_iterable(probe_filter.values()))
        # TODO(mholub): Allow more than 0 multimatches 
        multimatching = set(k for k,v in probe_filter.items() if v > 1) # is set anyway

        return multimatching
