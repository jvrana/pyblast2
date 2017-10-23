import ftplib
import re
import pprint
import subprocess
import os

dir_path = os.path.dirname(os.path.realpath(__file__))

version_pattern = "\d+\.\d+\.\d+\+"
blast_pattern = '(ncbi-blast-{ver})-{platform}{ext}'
tarball_ext = '\.tar\.gz'
tarball_formats = {
    "mac": blast_pattern.format(ver='\d+\.\d+\.\d+\+', platform='.*?macosx.*?', ext=tarball_ext),
    "linux(pentium)": blast_pattern.format(ver='\d+\.\d+\.\d+\+', platform='ia32-linux', ext=tarball_ext),
    "linux(X64 chip)": blast_pattern.format(ver='\d+\.\d+\.\d+\+', platform='x64-linux', ext=tarball_ext)
}

config = dict(
cwd = 'blast/executables/blast+/LATEST',
domain = "ftp.ncbi.nlm.nih.gov",
user = "anonymous",
email = "justin.vrana@gmail.com",
platform = "mac",
filename = None,
dir = dir_path,
)
config['fmt'] = tarball_formats[config['platform']]

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

for l in data:
    f = config['fmt']
    m = re.search(f, l)
    if m:
        config['filename'] = m.group()
        config['blastver'] = m.group(1)
        config['in'] = os.path.join(config['dir'], config['filename'])
        config['out'] = os.path.join(config['dir'], config['blastver'])

if config['filename'] is None:
    raise Exception("Unable to find blast binary. Data retrieved:\n{}".format(data))

print("filename: {}".format(config['filename']))

if not os.path.isfile(config['in']):
    if not os.path.isdir(config['in'].split('.')[0]):
        with open(config['in'], 'wb') as handle:
            print("downloading {0}".format(config['filename']))
            ftp.retrbinary('RETR {}'.format(config['filename']), handle.write)
            ftp.quit()

print("download complete")

# unzip tarball
pp.pprint(config)
subprocess.call(['tar', 'zxvpf', config['in'], '-C', config['dir']])