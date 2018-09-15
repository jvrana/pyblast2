import fire
from pyblast.blast_bin import BlastWrapper


def install(email, platform, path=None):
    """
    Installs blast

    :param email: users email
    :type email: basestring
    :param platform: platform (mac, linux, etc.)
    :type platform: basestring
    :param path: installation path
    :type path: basestring
    :return: None
    :rtype: None
    """
    b = BlastWrapper(path)
    return b.install(email, platform)


def main():
    return fire.Fire({
        "install": install
    })


if __name__ == '__main__':
    main()