"""
Script for installing BLAST
"""

import argparse

from blast_bin.install_manager import get_formats, install, initialize_files


def main():
    parser = argparse.ArgumentParser(description="Install BLAST from ncbi")
    parser.add_argument("user_email", type=str, help="A user email is required for BLAST download from ncbi.")
    parser.add_argument("platform", type=str, help="Choose your platform. Choose from {}".format(get_formats().keys()))
    parser.add_argument("--f", help="Forces install")
    args = parser.parse_args()
    force = False
    if args.f:
        force = True
    install(args.user_email, args.platform, force=force)


if __name__ == "__main__":
    main()
