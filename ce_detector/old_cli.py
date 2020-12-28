# !/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Console script for ce_detector.
@author: YangyangLi
@contact:li002252@umn.edu
@license: MIT Licence
@file: re_cli.py.py
@time: 2020/12/21 6:01 PM
"""
import sys

from annotator import Annotator
from detector import JunctionDetector
from utils import Timer, get_logger, get_parser


def main():
    """ function to integrate the annotation usage
    """

    log1 = get_logger('junction_detector', create_file=False)
    log2 = get_logger('annotate_junctions', create_file=False)

    args = get_parser()

    # detect junction reads
    detector = JunctionDetector(args.input, args.reference,
                                args.quality)
    with Timer() as t:
        junctionmap = detector.run(log1)

    log1.info(f'FINISHED FIND JUNCTION CONSUMING {t.elapsed:.2f}s')
    # annotate junction reads

    annotator = Annotator(junctionmap, args.gffdb)
    annotator.run(log2, args.out)


if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover
