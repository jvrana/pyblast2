import ftplib
import os
import pprint
import re
import shutil
import subprocess
from glob import glob
from os.path import abspath
from os.path import isdir
from os.path import isfile
from os.path import join
from os.path import split
from typing import Optional

import fire
from termcolor import colored
from termcolor import cprint

from pyblast.__bin__ import __bin__


def get_platform_formats():
    version_pattern = r"\d+\.\d+\.\d+\+"
    blast_pattern = "(ncbi-blast-{ver})-{platform}{ext}"
    tarball_ext = r"\.tar\.gz"
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


platform_formats = get_platform_formats()
valid_platforms = list(platform_formats.keys())


def get_format(platform):
    tarball_formats = platform_formats
    fmt = tarball_formats.get(platform, None)
    if fmt is None:
        raise Exception(
            'Cannot find platform "{}". Please select from {}'.format(
                platform, list(tarball_formats.keys())
            )
        )
    return fmt


def install_blast_using_ftp(
    email: str,
    platform: str,
    directory: str,
    user: str = "anonymous",
    domain: str = "ftp.ncbi.nlm.nih.gov",
    cwd: str = "blast/executables/blast+/LATEST",
    clean: bool = True,
):
    """Install Blast by downloading using ftp.

    :param email: users email
    :param platform: platform options
    :param directory: directory to install BLAST
    :param user: user for server
    :param domain: domain for NCBI server
    :param cwd: working directory for ftp server
    :param clean: if True (default), clean up file after extraction
    :return:
    """
    # check directory
    if not isdir(abspath(directory)):
        raise NotADirectoryError('Directory "{}" is not a directory'.format(directory))

    # get the platform format
    fmt = get_format(platform)

    # print config
    print("ftp config:")
    pp = pprint.PrettyPrinter(indent=4)
    pp.pprint(
        {
            "email": email,
            "platform": platform,
            "fileformat": fmt,
            "dir": directory,
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
        raise Exception("Unable to find blast binary. Data retrieved:\n{}".format(data))
    input_file = os.path.join(directory, filename)

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
        subprocess.call(["tar", "zxvpf", input_file, "-C", directory])
        os.remove(input_file)

    cprint("Blast installed to {}".format(input_file), "green")
    if clean:
        os.remove(input_file)
    return input_file


def _is_valid_bin(path):
    return isfile(join(path, "makeblastdb"))


def find_local_installations():
    glob_pattern = join(__bin__, "ncbi*/bin")
    dirs = sorted(list(glob(glob_pattern)))
    dirs = [d for d in dirs if isdir(abspath(d)) if _is_valid_bin(d)]
    return dirs


def get_install_path():
    return shutil.which("makeblastdb")


def is_installed():
    if get_install_path():
        return True
    return False


def is_locally_installed():
    if find_local_installations():
        return True
    return False


def find_global_installations():
    install_path = get_install_path()
    if install_path is None:
        return []
    return split(install_path)[0]


def status():
    """Checks status of blast.

    :param path: installation path
    :type path: basestring
    """
    path = get_install_path()
    if not path:
        print("BLAST is not installed.")
    else:
        print("BLAST is installed at {}".format(path))


def ask_for_email(email=None):
    if email is None:
        email = input("email address: ")
    return email


def ask_for_platform(platform=None):
    if platform is None:
        platform = input("platform ({}): ".format(valid_platforms))
    return platform


def force_install(directory, email=None, platform=None):
    email = ask_for_email(email)
    platform = ask_for_platform(platform)
    install_blast_using_ftp(email, platform, directory)


def ask_to_install(directory=None, email=None, platform=None):
    """Installs blast.

    :param email: users email
    :type email: basestring
    :param platform: platform (mac, linux, etc.)
    :type platform: basestring
    :param directory: installation path
    :type directory: basestring
    :return: None
    :rtype: None
    """
    if directory is None:
        directory = __bin__
    if not isdir(directory):
        raise NotADirectoryError("{} is not a directory".format(directory))

    if is_installed():
        warning = colored("Warning: BLAST is already installed at ", "red")
        warning += colored(get_install_path(), "blue")
        print(warning)
    else:
        cprint(
            "Blast not installed and so script cannot run. If you have BLAST installed, "
            "be sure to add it to your $PATH. Otherwise you can install it",
            "red",
        )

        global_installations = find_global_installations()
        bin_installations = find_local_installations()
        existing_installations = global_installations + bin_installations
        if existing_installations:
            warning = "Warning: Found existing installations at\n{}".format(
                "\n  ".join(existing_installations)
            )
            cprint(warning, "yellow")
            cprint(
                "You can use these insallations by adding them to your $PATH.", "yellow"
            )

    cprint("Install BLAST at {}?".format(directory), "green")
    cprint(
        "Note that this will require an internet connection and your email to "
        "download blast from NCBI.",
        "green",
    )
    answer = input(colored("Install? (y|n): ", "green"))
    if answer == "y":
        cprint("Installing BLAST", "green")
    else:
        cprint(
            "User canceled.\n"
            + "\nFor help, type 'pyblast install -- --help' in the commandline.",
            "red",
        )
        return

    force_install(directory, email=email, platform=platform)


def install(
    directory: Optional[str] = None,
    email: str = None,
    platform: str = None,
    interactive: bool = True,
    force: bool = False,
    use_bin: bool = True,
):
    if not force and is_installed():
        print("Install path: ", end="")
        cprint(get_install_path(), "blue")
        cprint(
            "BLAST already installed. Use '--force' argument to force installation",
            "red",
        )
    if directory is None:
        if use_bin:
            directory = __bin__
        if not isdir(__bin__):
            assert isdir(os.path.dirname(__bin__))
            os.makedirs(__bin__)
    if interactive:
        ask_to_install(directory, email, platform)
    else:
        cprint("Installing BLAST at {}".format(directory), "green")
        force_install(directory, email, platform)
    cprint("Be sure to add the installation to your $PATH!", "green")


def uninstall_local():
    """Uninstalls any local installations located in `pyblast/bin`"""
    if isdir(__bin__):
        cprint("Removing\n", "green")
        dirs = glob(join(__file__, "*"))
        cprint("\n".join(dirs), "blue")
        shutil.rmtree(__bin__)


def main():
    return fire.Fire(
        {"install": install, "status": status, "uninstall_local": uninstall_local}
    )


def entrypoint():
    return main()


if __name__ == "__main__":
    main()
