import argparse
import ftplib
import os
import pprint
import re
import shutil
import subprocess

dir_path = os.path.dirname(os.path.realpath(__file__))


class PathManager(object):
    def __init__(self, path):
        self.pathfile = os.path.abspath(path)
        self.pathdir = os.path.dirname(self.pathfile)

    def append_path(self, path):
        paths = self.lines()
        paths.append(path)
        paths = list(set(paths))
        with open(self.pathfile, 'w') as f:
            f.writelines(paths)

    def lines(self):
        lines = []
        with open(self.pathfile, 'r') as f:
            lines = f.readlines()
        lines = [line.strip() for line in lines]
        lines = list(set(lines))
        if '' in lines:
            lines.remove('')
        return lines

    def paths(self):
        lines = self.lines()
        paths = [os.path.join(self.pathdir, line) for line in lines]
        return list(set(paths))

    def append_paths_to_env(self):
        for path in self.paths():
            os.environ['PATH'] += ':'+os.path.join(path, 'bin')


def install_blast(user_email, force=False):
    def get_blast_format(platform):
        global tarball_formats
        version_pattern = "\d+\.\d+\.\d+\+"
        blast_pattern = '(ncbi-blast-{ver})-{platform}{ext}'
        tarball_ext = '\.tar\.gz'
        tarball_formats = {
            "mac"            : blast_pattern.format(ver=version_pattern, platform='.*?macosx.*?', ext=tarball_ext),
            "linux(pentium)" : blast_pattern.format(ver=version_pattern, platform='ia32-linux', ext=tarball_ext),
            "linux(X64 chip)": blast_pattern.format(ver=version_pattern, platform='x64-linux', ext=tarball_ext)
        }
        return tarball_formats[platform]

    # setup config
    config = dict(
            cwd='blast/executables/blast+/LATEST',
            domain="ftp.ncbi.nlm.nih.gov",
            user="anonymous",
            email=user_email,
            platform="mac",
            filename=None,
            dir=dir_path,
    )
    config['fmt'] = get_blast_format(config['platform'])
    # detect blast executable
    path = shutil.which('makeblastdb')
    if path is None or force:
        print(" ********** Installing BLAST ********** ")
        install_blast_using_ftp(config)
    else:
        print("BLAST installed at "+path)


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
    pm = PathManager(os.path.abspath(os.path.join(dir_path, "_paths.txt")))
    pm.append_path(os.path.join(config['blastver']))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Install BLAST from ncbi")
    parser.add_argument("user_email", type=str, help="A user email is required for BLAST download from ncbi.")
    parser.add_argument("--f", help="Forces install")
    args = parser.parse_args()
    force = False
    if args.f:
        force = True
    install_blast(args.user_email, force=force)
