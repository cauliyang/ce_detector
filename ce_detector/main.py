#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: YangyangLi
@contact:li002252@umn.edu
@version: 0.0.1
@license: MIT Licence
@file: main.py
@time: 2021/1/29 7:30 AM
"""
from .annotator import Annotator
from .detector import JunctionDetector
from .scanner import Scanner


def main(chrom, ann_chrom, bam, reference, quality, gffdb, cutoff, verbose):
    """
    :param chrom:
    :type chrom:
    :param ann_chrom:
    :type ann_chrom:
    :param bam:
    :type bam:
    :param reference:
    :type reference:
    :param quality:
    :type quality:
    :param gffdb:
    :type gffdb:
    :param cutoff:
    :type cutoff:
    :param out:
    :type out:
    :param verbose:
    :type verbose:
    :return:
    :rtype:
    """
    detector = JunctionDetector(
        bam,
        reference,
        quality,
    )
    junctionmap = detector.run(chrom=chrom, ann_chrom=ann_chrom, verbose=verbose)

    annotator = Annotator(gffdb)
    scanner = Scanner(cutoff=cutoff)

    if not junctionmap.is_empty():
        junctionmap = annotator.run(junctionmap=junctionmap, verbose=verbose)
        junctionmap = scanner.run(junctionmap, verbose=verbose)

    return junctionmap

    # scanner.write2file(verbose=verbose)
