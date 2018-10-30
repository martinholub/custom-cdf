import shutil
import sys
import os
import requests
import ftplib
import zipfile
import gzip
import numpy as np
import re
from tqdm import tqdm as TQDM

# TODO: it may be nice to have files as objects

class CdfSetupException(Exception):
    def __init__(self, url=None, message="CDF Setup Exception."):
        if url:
            ex_message = "Failed on url: {}".format(url)
        else:
            ex_message = message
        super().__init__(ex_message)

class CdfPipeInitializer(object):
    def __init__(self, cfg, workdir = "."):
        self.cfg = cfg
        self.workdir = workdir
        assert isinstance(cfg, (dict, ))

    # MH: Learning how to use properties ---------------------------------------
    # Likely overkill here, but useful elsewhere
    # https://docs.python.org/3/library/functions.html#property
    # https://www.programiz.com/python-programming/property
    @property
    def is_debug(self):
        """bool, if we are in debug mode"""
        return self._is_debug

    @is_debug.setter
    def is_debug(self, value):
        self._is_debug = bool(value)

    @property
    def cfg(self):
        """dict, snakemake config"""
        return self._cfg

    @cfg.setter
    def cfg(self, value):
        # potentially add checks on fields in dict
        assert isinstance(value, (dict, ))
        self._cfg = value
    # --------------------------------------------------------------------------

    def create_dir_with_check(self, dirname):
        """Create directory from 'dir' field in config

        Parameters:
            dirname: str, path to dir to be created
        """
        if isinstance(dirname, list): dirname = dirname[0]

        if not os.path.isdir(dirname):
            try:
                os.makedirs(dirname, exist_ok = False)
            except OSError as e:
                print("Directory already exists. You should not be seeing this.")
                raise e

    def check_filename(self, fname, ext = None, chrs = []):
        """Normalize filename

        Parameters:
            fname: str, filename
            ext: str, extension with or without leading dot
            chrs: list, list of strings defining chromosmes to use
        """
        # Default extension is .gz
        if not ext:
            ext = ".gz"
        if not ext[0]==".":
            ext = "." + ext

        # Append extnesion to filename. It is expected that it arrives without it.
        # TODO: check if you do not need non-greedines with ".*(\..*?)$"
        curr_ext = re.match(".*(\..*)$", fname).group(1)
        if ext != curr_ext and curr_ext != "zip":
            fname += ext

        # Expand to one filename per chromosome fasta file
        if re.match("\*", fname):
            if chrs:
                if not type(chrs) is list: chrs = [chrs]
                fname = self.expand_filenames(fname, chrs)
            else:
                msg = "Wildcard in filename, but no 'chromosomes' defined. Aborting."
                raise CdfSetupException(msg)

        return fname

    def parse_chromosomes(self, chrs):
        """Parses list of chromosomes to download

        Expects the config entry 'chromosomes' to hold entry that can look like:
        [{1..3}, 5, MT], a list defining which chromosmes to expand to.
        """
        chrs_list = []
        if not isinstance(chrs, list): chrs = [chrs]
        for i in chrs:
            if re.match("\.{2,}", i): # two dot pattern indicates range
                lims = i.split("..")
                lims = [int(float(x.strip(".")) + y) for x,y in zip(lims, [0,1])]
                assert len(lims) == 2
                chrs_list.extend([str(x) for x in np.arange(*lims)])
            else:
                chrs_list.extend(i)
        return chrs_list

    def expand_filenames(self, fname, chrs):
        """Expand filename on '*' with values from 'chromosomes'

        Returns
        -----------
        list
            file names expanded on asterisk as a list of strings

        Examples
        ----------------
        reference:
          file: ref.*.fa
          chromsomes: [{1..3}, 8, MT]
        # will produce [ref.1.fa, ref.2.fa, ref.8.fa, ref.MT.fa]

        Raises
        -----------
        AssertionError
            When `chrs` are not list or `fname` a string
        """

        chrs = self.parse_chromosomes(chrs)
        assert isinstance(chrs, list), "`chrs` must be parsed to a list of strings."
        assert isinstance(fname, str), "`fname` must be a single string."
        assert re.match("\*",fname), "`fname` must include an asterisk `*` to expand on."
        return [re.sub("\*", x, fname) for x in chrs]

    def _is_gzipfile(self, fname):
        """Check if file is valid gzip file"""
        with open(fname, "rb") as f:
            return zipfile.binascii.hexlify(f.read(2)) == b"1f8b"

    def validate_file(self, fname):
        """
        Check if file is a valid GZIP file and convert if necessary
        """

        assert os.path.isfile(fname), "File {} does not exist.".format(fname)
        assert os.path.getsize(fname) > 0, "File {} has size 0.".format(fname)

        if not self._is_gzipfile(fname) and zipfile.is_zipfile(fname):
            try:
                fparts = list(filter(bool, fname.split("/")))
                print("Repacking {} to valid GZIP archive...".format(fparts[-1]))
                fdir = os.path.join(*fparts[:-1])

                with zipfile.ZipFile(fname, 'r') as zf:
                    unzip_f_new = fparts[-1].replace(".gz","")
                    unzip_f = [x for x in zf.namelist() if any(a in x.lower() for a in ("probe", "fasta"))]
                    if len(unzip_f) == 0:
                        unzip_f = unzip_f_new # Try the anticipated name
                    else:
                        if len(unzip_f) > 1:
                            print("Warning: Multiple files ({}) in archive may be valid FASTA.".format(unzip_f))
                        unzip_f = unzip_f[0]

                    zfpath = zf.extract(unzip_f, fdir)

                os.remove(fname)
                with open(zfpath, 'rb') as f_in:
                    with gzip.open(fname, 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
                os.remove(zfpath)

            except Exception as e:
                print("File {} didn't pass check.".format(fname))
                raise e

    def download_file_batch(self, baseurl, fnames, dirname):
        """Download files specified as list 'fnames' from 'baseurl'"""

        if not isinstance(fnames, list): fnames = [fnames]

        for f in fnames:
            f_url = os.path.join(baseurl, f)
            self.download_file(f_url, dirname, f)

    def _ftp_callback(self, block, file_, progbar):
        """Write to file from ftp in blocks and update TDQM progress bar"""
        file_.write(block)
        progbar.update(len(block))

    def _download_checker(self, fname):
        """Check if we need to download anything
        """
        do_download = False
        try: # overwrite bad trials
            # 10000 is arbitrary, else would have to peek for size of remote file
            do_download = os.path.getsize(fname) < 10000
        except FileNotFoundError:
            do_download = True
        return do_download

    def download_file(self, url, dirname, fname = None):
        """Download single file from url to specified folder"""

        if isinstance(dirname, list): dirname = dirname[0]
        if isinstance(url, list): url = url[0]

        if not fname:
            fname = list(filter(bool, url.split("/")))[-1]
        assert len(fname) > 0
        assert isinstance(fname, str)

        fname_ = fname
        fname = os.path.join(dirname, fname)
        if self._download_checker(fname):
            if url.startswith("http"):
                # Download file from HTTP location
                # TODO: test w/ HTTPS
                with requests.Session() as sess:
                    r =  sess.get(url, stream = True)
                    with open(fname, "wb") as f:
                        print("Downloading {}".format(fname_))
                        tqdm_params = { 'unit': 'blocks',
                                        'unit_scale': True,
                                        'leave': False,
                                        'miniters': 1,
                                        'total': int(r.headers["Content-Length"])}
                        with TQDM(**tqdm_params) as tqdm: # progress bar
                            for chunk in r.iter_content(chunk_size = 512):
                                f.write(chunk)
                                tqdm.update(len(chunk))
                                # shutil.copyfileobj(r.raw, f)

            elif url.startswith("ftp"):
                # Download file from FTP server
                url = url.replace("ftp://", "") # re.sub("^ftp://","",url)
                fparts = list(filter(bool, url.split("/")))
                ftp = ftplib.FTP(fparts[0]) # base remote address
                ftp.login()
                ftp.cwd(os.path.join(*fparts[1:-1]))

                assert fparts[-1] in ftp.nlst(), "File {} not in remote location.".format(fparts[-1])

                cmd = "RETR {}".format(fparts[-1])
                # Perhaps can condition only here
                # if ftp.size(fparts[-1]) > os.path.getsize(fname):
                with open(fname, "wb") as f:
                    print("Downloading {}".format(fparts[-1]))
                    tqdm_params = { 'unit': 'blocks',
                                    'unit_scale': True,
                                    'leave': False,
                                    'miniters': 1,
                                    'total': ftp.size(fparts[-1])}
                    with TQDM(**tqdm_params) as tqdm: # progress bar
                        ftp.retrbinary( cmd, lambda block: self._ftp_callback(block, f, tqdm),
                                        blocksize = 512)
                ftp.quit()
        else:
            print("File {} exists. Not overwritten.".format(fname))

        return fname

    def down_me_maybe(self):
        """Download files as specified in configuration file

        To download file, the required fields are 'url', 'file', 'dir'. Optionally,
        the 'file' can contain wildcard '*'. It will be expanded on the wildcard
        to all chromsomes as specified in 'chromosomes'.
        See [expand_filenames](/expand_filenames).

        # Example:
        reference:
          file: ref.fa
          url: ftp.ensembl.org/current_fasta/ref.fa
          dir: data/reference
        """

        for k,v in self.cfg.items():
            if isinstance(v, (dict, )):
                if all(x in v.keys() for x in ["url", "dir"]):
                    # if "filter" in v.keys():
                    #     if not v["filter"]: continue # dont download vcf
                    fname = None
                    if "file" in v.keys():
                        fname = v["file"]
                        if "chromosomes" in v.keys():
                            fname = self.check_filename(fname, v["chromosomes"])
                        else:
                            fname = self.check_filename(fname)


                    if v["dir"]: self.create_dir_with_check(v["dir"])
                    if v["url"]:
                        fname = self.download_file(v["url"], v["dir"], fname)
                        self.validate_file(fname)


if __name__ == '__main__':
    # snakemake object is defined by Snakemake workflow.
    cpi = CdfPipeInitializer(snakemake.config, os.getcwd())
    cpi.down_me_maybe()
