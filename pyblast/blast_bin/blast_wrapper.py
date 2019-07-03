import pprint
import ftplib
import subprocess
import os
import re
import shutil
from glob import glob


class BlastWrapper(object):
    """A really shallow wrapper for running blast commands"""

    def __init__(self, path=None):
        if path is None:
            path = self.find_global_installations()
        if path is None:
            paths = self.find_bin_installations()
            if len(paths) > 0:
                path = paths[0]
                print('Found previous installation at "{}"'.format(path))
        self.path = path
        if path is not None:
            self._add_path_to_environment(path)

    @staticmethod
    def find_bin_installations():
        here = os.path.dirname(os.path.abspath(__file__))
        glob_pattern = os.path.join(here, "bin", "ncbi*/bin")
        dirs = sorted(list(glob(glob_pattern)))
        dirs = [d for d in dirs if os.path.isdir(os.path.abspath(d))]
        return dirs

    @staticmethod
    def find_global_installations():
        install_path = shutil.which("makeblastdb")
        if install_path is None:
            return None
        return os.path.split(install_path)[0]

    @staticmethod
    def _add_path_to_environment(path):
        """Add the blast installations located in 'bin' to the environment temporarily"""
        real_path = os.path.abspath(path)
        os.environ["PATH"] += ":" + real_path

    def is_installed(self):
        return shutil.which("makeblastdb") is not None

    # def install(self, email, platform, dir):
    #     return self.install_blast_using_ftp(email, platform, dir)

    def install(self, email, platform):
        """
        Installs blast to bin directory.

        :param email: user's email
        :type email: basestring
        :param platform: platform (mac, linux, etc.)
        :type platform: basestring
        :return: None
        :rtype: None
        """
        dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "bin")
        if not os.path.isdir(dir):
            os.makedirs(dir)
        path = self.install_blast_using_ftp(email, platform, dir)
        self.path = os.path.join(path, "bin")

    @staticmethod
    def _platform_formats():
        version_pattern = "\d+\.\d+\.\d+\+"
        blast_pattern = "(ncbi-blast-{ver})-{platform}{ext}"
        tarball_ext = "\.tar\.gz"
        tarball_formats = {
            "mac": blast_pattern.format(
                ver=version_pattern, platform=".*?macosx.*?", ext=tarball_ext
            ),
            "linux(pentium)": blast_pattern.format(
                ver=version_pattern, platform="ia32-linux", ext=tarball_ext
            ),
            "linux(X64 chip)": blast_pattern.format(
                ver=version_pattern, platform="x64-linux", ext=tarball_ext
            ),
        }
        return tarball_formats

    @classmethod
    def _valid_platforms(cls):
        return list(cls._platform_formats().keys())

    @classmethod
    def _get_format(cls, platform):
        tarball_formats = cls._platform_formats()
        fmt = tarball_formats.get(platform, None)
        if fmt is None:
            raise Exception(
                'Cannot find platform "{}". Please select from {}'.format(
                    platform, list(tarball_formats.keys())
                )
            )
        return fmt

    @classmethod
    def install_blast_using_ftp(
        cls,
        email,
        platform,
        dir,
        user="anonymous",
        domain="ftp.ncbi.nlm.nih.gov",
        cwd="blast/executables/blast+/LATEST",
    ):
        # print config
        if not os.path.isdir(os.path.abspath(dir)):
            raise NotADirectoryError('Directory "{}" is not a directory'.format(dir))

        fmt = cls._get_format(platform)
        print("ftp config:")
        pp = pprint.PrettyPrinter(indent=4)
        pp.pprint(
            {
                "email": email,
                "platform": platform,
                "fileformat": fmt,
                "dir": dir,
                "user": user,
                "domain": domain,
                "cwd": cwd,
            }
        )
        print()

        # login
        ftp = ftplib.FTP(domain)
        ftp.login(user, email)

        # get files
        data = []
        ftp.cwd(cwd)
        ftp.dir(data.append)

        # find filename with format
        filename = None
        blastver = None
        for line in data:
            m = re.search(fmt, line)
            if m:
                filename = m.group()
                blastver = m.group(1)
                break
        if filename is None:
            raise Exception(
                "Unable to find blast binary. Data retrieved:\n{}".format(data)
            )

        input_file = os.path.join(dir, filename)

        # download tarball
        if not os.path.isfile(input_file):
            if not os.path.isdir(blastver):
                print("Downloading {}".format(filename))
                with open(input_file, "wb") as handle:
                    ftp.retrbinary("RETR {}".format(filename), handle.write)
                    ftp.quit()
                print("Download complete")

        # unzip and remove tarball
        if os.path.isfile(input_file):
            subprocess.call(["tar", "zxvpf", input_file, "-C", dir])
            os.remove(input_file)
        print("Blast installed to {}".format(input_file))
        return input_file

    def ask_to_install(self):
        response = input(
            "Blast not installed and so script cannot run. Would you like to install it now?\n"
            + "Note that this will require an internet connection and your email to download blast "
            "from NCBI: (y|n)"
        )
        if response == "y":
            email = input("email address: ")
            platform = input("platform ({}): ".format(BlastWrapper._valid_platforms()))
            self.install(email, platform)
        else:
            raise Exception(
                "Blast not installed. Please run 'pyblast install [EMAIL] [PLATFORM]' from the commandline"
                " For help, type 'pyblast install' in the commandline."
            )
