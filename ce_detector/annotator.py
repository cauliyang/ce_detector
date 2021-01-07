#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""class for annotating junction reads
"""
from __future__ import annotations

from collections import defaultdict
from typing import Any

import gffutils
import numpy as np

from .utils import get_yaml
from .utils import timethis

# HardCode the information of chromosome because its name of two ref are not identical
CHROMS = get_yaml()["chr2hg38"]


class Annotator:
    """annotate junction reads

    :param junctionmap:  instance return :class:`ce_detector.detector.JunctionMap`
    :type junctionmap: instance
    :param database: database of annotation files
    :type database: Any
    :param output: filename of annotated junction reads. Defaults to None
    :type output: TestIo
    """

    def __init__(self, junctionmap, database: Any, output=None):
        """Constructor of Annotator"""
        self.junctionMap = junctionmap
        self.database = gffutils.FeatureDB(database)
        self.output = output

    @staticmethod
    def detect_property(start, end, junction_list):
        """detect type of slice, number of skipped donors and number of skipped acceptors

        type of slice including D A DA N NDA

        :param start: start of junction read
        :type start: int
        :param end: end of junction read
        :type end: int
        :param junction_list: gene list of junction reads
        :type junction_list: numpy.array
        :return: type of slice, number of skipped donors, number of skipped acceptors
        """
        known_donors = junction_list[:, 0]
        known_acceptors = junction_list[:, 1]

        donors_skipped = ((start + 1 < known_donors) & (known_donors < end)).sum()

        acceptors_skipped = (
            (start + 1 < known_acceptors) & (known_acceptors < end)
        ).sum()

        if [start + 1, end] in junction_list.tolist():

            reads_type = "DA"

        elif (start + 1 in known_donors) and (end in known_acceptors):

            reads_type = "NDA"

        elif (start + 1 in known_donors) and (end not in known_acceptors):

            reads_type = "D"

        elif (start + 1 not in junction_list[:, 0]) and (end in junction_list[:, 1]):

            reads_type = "A"

        else:

            reads_type = "N"

        return reads_type, donors_skipped, acceptors_skipped

    def annotate_junction(self, read, result, db):
        """annotate junction reads and write results to file

        :param read: junction read return :class:`ce_detector.detector.Read`
        :type read: instance
        :param result: gene list used for annotation
        :type result: defaultdict[Any, Any]
        :param db: database of annotation file
        :type db: instance of file
        """
        chrom = read.chrom
        start, end = read.start, read.end
        region = f"{CHROMS[chrom]}:{start}-{end}"  # change chromosome

        gene_list = []

        for gene in db.features_of_type(("gene"), limit=region):

            gene_list.append(gene.id)
            # new gene
            if gene.id not in result:
                # transcript

                for transcript in db.children(
                    gene,
                    level=1,
                    featuretype=("primary_transcript", "transcript", "mRNA"),
                ):  # mRNA may drop out

                    # find all position of introns for every gene known junctions
                    for junction in db.interfeatures(
                        db.children(
                            transcript,
                            level=1,
                            featuretype="exon",
                            order_by="start",
                        ),
                        new_featuretype="intron",
                    ):
                        result[gene.id].append([junction.start, junction.end])

                if (
                    len(result[gene.id]) >= 1
                ):  # drop gene with only one exon in terms of number of intron

                    result[gene.id] = np.unique(result[gene.id], axis=0)

                else:
                    result.pop(gene.id)
                    gene_list.pop()

        # annotate junctions reads

        for gene in gene_list:
            junction_list = result[gene]

            reads_type, donors_skipped, acceptors_skipped = self.detect_property(
                start,
                end,
                junction_list,
            )
            #         read.type = reads_type
            #     output.write(f'{reads_information}\t{reads_type}\t{gene}\n')

            read.information.append(
                [reads_type, donors_skipped, acceptors_skipped, gene],
            )
            if self.output:
                self.output.write(
                    f"{read}\t{reads_type}\t{donors_skipped}\t{acceptors_skipped}\t{gene}\n",
                )

    @timethis(name="Junction Annotator", message="FINISHED")
    def run(self):
        """main function used to annotate junction reads

        pick all genes covered by one junction read and annotate all of them:
        type of slice, number of skipped donors and number of skipped acceptors

        """

        result = defaultdict(list)

        # logger.info(f'Begin to annotate junctions!')

        for _index, read in enumerate(self.junctionMap):
            #         if index < 10:
            self.annotate_junction(read, result, self.database)
