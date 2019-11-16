import fire

from pyblast.blast_bin import BlastWrapper


def status(path=None):
    """Checks status of blast.

    :param path: installation path
    :type path: basestring
    """
    b = BlastWrapper(path)
    path = b.is_installed()
    if not path:
        print("FALSE: Blast is not installed.")
    else:
        print("TRUE: Blast is installed at {}".format(path))


def install(email=None, platform=None, path=None):
    """Installs blast.

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
    if email is None and platform is None:
        b.ask_to_install()
    else:
        return b.install(email, platform)


def main():
    return fire.Fire({"install": install, "status": status})


def entrypoint():
    return main()


if __name__ == "__main__":
    main()
