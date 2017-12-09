"""
Script for installing BLAST
"""

import argparse
import ftplib
import json
import os
import pprint
import re
import shutil
import subprocess

PATHS = "_paths.json"
DIR_PATH = os.path.dirname(os.path.realpath(__file__))
BIN_DIR = os.path.join(DIR_PATH, 'bin')
INSTALL_PATHS = os.path.abspath(os.path.join(DIR_PATH, PATHS))


def which(program):
    """Return path of executable"""
    import os
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    return None


def remove_bin_installations():
    """Remove the _paths.json and any installations in the bin folder"""
    if os.path.isdir(BIN_DIR):
        shutil.rmtree(BIN_DIR)
    if os.path.isfile(PATHS):
        os.remove(PATHS)


def open_install_paths():
    """Open paths located in _paths.txt"""
    if not os.path.isfile(INSTALL_PATHS):
        initialize_files()
    with open(INSTALL_PATHS, 'r') as in_file:
        return json.load(in_file)


def initialize_files():
    if not os.path.isdir(BIN_DIR):
        os.mkdir(BIN_DIR)
    if os.path.isfile(PATHS):
        try:
            open_install_paths()
        except json.decoder.JSONDecodeError as e:
            os.remove(PATHS)
    else:
        with open(INSTALL_PATHS, 'w') as paths:
            json.dump([], paths)


def append_to_install_paths(new_path):
    """Append a new bin path to the list of paths"""
    paths = open_install_paths()
    if new_path not in paths:
        paths.append(new_path)
    with open(INSTALL_PATHS, 'w') as out_file:
        json.dump(paths, out_file)


def add_paths_to_environment():
    """Add the blast installations located in 'bin' to the environment temporarily"""
    for path in open_install_paths():
        real_path = os.path.abspath(os.path.join(path, 'bin'))
        os.environ['PATH'] += ':' + real_path


def has_executable():
    """Whether blast is installed and executable"""
    b = which('makeblastdb')
    print(b)
    return b is not None


def check_installation():
    """Add the path located in blast_bin/_paths.txt to the environment in an attempt to run blast"""
    if not has_executable():
        add_paths_to_environment()
    if not has_executable():
        error_message = "BLAST executables not found in path. Be sure BLAST is correctly installed."
        help_message = "Please run 'install_pyblast <youremail> <yourplatform' in your terminal." \
                       "Run 'install_pyblast -h' for help."
        raise Exception(error_message + "\n" + "*" * 50 + "\n" + help_message)


def get_formats():
    """Return blast installation platform to installation filename patterns"""
    version_pattern = "\d+\.\d+\.\d+\+"
    blast_pattern = '(ncbi-blast-{ver})-{platform}{ext}'
    tarball_ext = '\.tar\.gz'
    tarball_formats = {
        "mac": blast_pattern.format(ver=version_pattern, platform='.*?macosx.*?', ext=tarball_ext),
        "linux(pentium)": blast_pattern.format(ver=version_pattern, platform='ia32-linux', ext=tarball_ext),
        "linux(X64 chip)": blast_pattern.format(ver=version_pattern, platform='x64-linux', ext=tarball_ext)
    }
    return tarball_formats


def get_blast_formats(platform):
    """Return installation filename from a platform name"""
    tarball_formats = get_formats()
    return tarball_formats[platform]


def install(user_email, platform, force=False):
    # setup config
    initialize_files()
    config = dict(
        cwd='blast/executables/blast+/LATEST',
        domain="ftp.ncbi.nlm.nih.gov",
        user="anonymous",
        email=user_email,
        platform=platform,
        filename=None,
        dir=BIN_DIR
    )
    if not os.path.isdir(config['dir']):
        os.mkdir(config['dir'])
    config['fmt'] = get_blast_formats(config['platform'])
    # detect blast executable
    path = shutil.which('makeblastdb')
    if path is None or force:
        print(" ********** Installing BLAST ********** ")
        install_blast_using_ftp(config)
    else:
        print("BLAST installed at " + path)


def install_blast_using_ftp(config):
    # print config
    print("ftp config:")
    pp = pprint.PrettyPrinter(indent=4)
    pp.pprint(config)
    print()
    # login
    ftp = ftplib.FTP(config['domain'])
    ftp.login(config['user'], config['email'])
    # get files
    data = []
    ftp.cwd(config['cwd'])
    ftp.dir(data.append)
    # find filename with format
    for l in data:
        f = config['fmt']
        m = re.search(f, l)
        if m:
            config['filename'] = m.group()
            config['blastver'] = m.group(1)
            config['in'] = os.path.join(config['dir'], config['filename'])
            config['out'] = os.path.join(config['dir'], config['blastver'])
            config['exec'] = os.path.join(config['out'])
            break
    if config['filename'] is None:
        raise Exception("Unable to find blast binary. Data retrieved:\n{}".format(data))
    print("filename: {}".format(config['filename']))
    # download tarball
    if not os.path.isfile(config['in']):
        if not os.path.isdir(config['blastver']):
            with open(config['in'], 'wb') as handle:
                print("downloading {0}".format(config['filename']))
                ftp.retrbinary('RETR {}'.format(config['filename']), handle.write)
                ftp.quit()
            print("download complete")

    # unzip and remove tarball
    if os.path.isfile(config['in']):
        subprocess.call(['tar', 'zxvpf', config['in'], '-C', config['dir']])
        os.remove(config['in'])

    # add to path
    append_to_install_paths(config['out'])


def main():
    parser = argparse.ArgumentParser(description="Install BLAST from ncbi")
    parser.add_argument("user_email", type=str, help="A user email is required for BLAST download from ncbi.")
    parser.add_argument("platform", type=str, help="Choose your platform. Choose from {}".format(
        ', '.join(["'{}'".format(x) for x in get_formats().keys()])))
    parser.add_argument("--f", help="Forces install")
    args = parser.parse_args()
    force = False
    if args.f:
        force = True
    install(args.user_email, args.platform, force=force)


if __name__ == "__main__":
    main()
