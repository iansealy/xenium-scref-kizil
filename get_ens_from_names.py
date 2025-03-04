#!/usr/bin/env python3

"""Get Ensembl IDs from gene names.

The script takes multiple tab-delimited files of Ensembl IDs and their corresponding
gene names. It then takes a list of gene names, which are converted to Ensembl IDs.
"""

import argparse
import typing
from collections import defaultdict

__version__ = "1.0.0"


def main(arg: argparse.Namespace) -> None:
    """Get Ensembl IDs from gene names.

    Wrapper function used when file is run as a script.
    """
    get_ens_from_names(arg.files, arg.names)


def get_ens_from_names(files: list[typing.TextIO], names: typing.TextIO) -> None:
    """Get Ensembl IDs from gene names."""
    name_ens_dict: defaultdict[str, str] = defaultdict()
    for file in files:
        for line in file:
            ensembl_id, name = line.strip().split("\t")
            if name in name_ens_dict:
                print(f"{name} maps to {name_ens_dict[name]} and {ensembl_id}")
            else:
                name_ens_dict[name] = ensembl_id
    for name in names:
        if name not in name_ens_dict:
            print(f"Missing {name}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "files",
        metavar="FILE",
        nargs="+",
        type=argparse.FileType("r"),
        help="TSV files linking Ensembl IDs to gene names",
    )
    parser.add_argument(
        "--names",
        "-n",
        dest="names",
        metavar="FILE",
        type=argparse.FileType("r"),
        required=True,
        help="the file of gene names",
    )
    parser.add_argument(
        "--version",
        "-v",
        action="version",
        version="%(prog)s " + __version__,
    )

    args = parser.parse_args()

    main(args)
