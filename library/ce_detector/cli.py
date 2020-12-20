"""Console script for ce_ddtertor."""

import argparse
import sys


def main():
    """Console script for ce_ddtertor."""
    parser = argparse.ArgumentParser(description='__doc__')

    parser.add_argument('_', nargs='*')

    args = parser.parse_args()

    print("Arguments: " + str(args._))
    print("Replace this message by putting your code into "
          "ce_ddtertor.cli.main")

    return 0


if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover
