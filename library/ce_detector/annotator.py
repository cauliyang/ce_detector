#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""class for annotating junction reads
"""
from __future__ import annotations

from collections import defaultdict

import gffutils
import numpy as np

# HardCode the information of chromosome because its name of two ref are not identical
CHROMS = {
    'chr1': 'NC_000001.11',
    'chr2': 'NC_000002.12',
    'chr3': 'NC_000003.12',
    'chr4': 'NC_000004.12',
    'chr5': 'NC_000005.10',
    'chr6': 'NC_000006.12',
    'chr7': 'NC_000007.14',
    'chr8': 'NC_000008.11',
    'chr9': 'NC_000009.12',
    'chr10': 'NC_000010.11',
    'chr11': 'NC_000011.10',
    'chr12': 'NC_000012.12',
    'chr13': 'NC_000013.11',
    'chr14': 'NC_000014.9',
    'chr15': 'NC_000015.10',
    'chr16': 'NC_000016.10',
    'chr17': 'NC_000017.11',
    'chr18': 'NC_000018.10',
    'chr19': 'NC_000019.10',
    'chr20': 'NC_000020.11',
    'chr21': 'NC_000021.9',
    'chr22': 'NC_000022.11',
    'chrX': 'NC_000023.11',
    'chrY': 'NC_000024.10',
    'chrM': 'NC_012920.1'
}


class Annotator:
    """annotate junction reads"""

    def __init__(self, junctionmap, database):
        """
        :param junctionmap:  instance from JunctionMap
        :type junctionmap: instance
        :param database: filename of database of annotation files
        :type database: str
        """
        self.junctionMap = junctionmap
        self.database = gffutils.FeatureDB(database)

    @staticmethod
    def detect_property(start, end, junction_list):
        """ detect type of slice, number of skipped donors and number of skipped acceptors

        type of slice including D A DA N NDA

        :param start: start of junction read
        :type start: int
        :param end: end of junction read
        :type end: int
        :param junction_list: gene list of junction reads
        :type junction_list: list
        :return: type of slice, number of skipped donors, number of skipped acceptors
        """
        known_donors = junction_list[:, 0]
        known_acceptors = junction_list[:, 1]

        donors_skipped = ((start + 1 < known_donors) & (known_donors < end)).sum()

        acceptors_skipped = ((start + 1 < known_acceptors) &
                             (known_acceptors < end)).sum()

        if [start + 1, end] in junction_list.tolist():

            reads_type = 'DA'

        elif (start + 1 in known_donors) and (end in known_acceptors):

            reads_type = 'NDA'

        elif (start + 1 in known_donors) and (end not in known_acceptors):

            reads_type = 'D'

        elif (start +
              1 not in junction_list[:, 0]) and (end in junction_list[:, 1]):

            reads_type = 'A'

        else:

            reads_type = 'N'

        return reads_type, donors_skipped, acceptors_skipped

    def annotate_junction(self, read, result: 'list', db, output=None):
        """annotate junction reads and write results to file

        :param read: junction read
        :type read: instance
        :param result: gene list used for annotation
        :type result: list
        :param db: database of annotation file
        :type db: instance of file
        :param output:  file handler of output. Defaults to None
        :type output: TestIo
        """
        chrom = read.chrom
        start, end = read.start, read.end
        region = f'{CHROMS[chrom]}:{start}-{end}'  # change chromosome

        gene_list = []

        for gene in db.features_of_type(('gene'), limit=region):

            gene_list.append(gene.id)
            # new gene
            if gene.id not in result:
                # transcript

                for transcript in db.children(
                        gene,
                        level=1,
                        featuretype=('primary_transcript', 'transcript',
                                     'mRNA')):  # mRNA may drop out

                    # find all position of introns for every gene known junctions
                    for junction in db.interfeatures(
                            db.children(transcript,
                                        level=1,
                                        featuretype='exon',
                                        order_by='start'),
                            new_featuretype='intron'):
                        result[gene.id].append([junction.start, junction.end])

                if len(result[gene.id]) >= 1:  # drop gene with only one exon

                    result[gene.id] = np.unique(result[gene.id], axis=0)

                else:
                    result.pop(gene.id)
                    gene_list.pop()

        # annotate junctions reads

        for gene in gene_list:

            junction_list = result[gene]

            reads_type, donors_skipped, acceptors_skipped = self.detect_property(
                start, end, junction_list)
            #         read.type = reads_type
            #     output.write(f'{reads_information}\t{reads_type}\t{gene}\n')
            if output:
                output.write(
                    f'{read}\t{reads_type}\t{donors_skipped}\t{acceptors_skipped}\t{gene}\n'
                )

    def run(self, logger, output='ohoannnotate_junction.bed'):
        """ main function used to annotate junction reads

        pick all genes covered by one junction read and annotate all of them:
        type of slice, number of skipped donors and number of skipped acceptors

        :param logger: logging handler
        :type logger: instance
        :param output: filename of outpput. Defaults to None
        :type output: str
        """

        result = defaultdict(list)

        logger.info(f'Begin to annotate junctions!')

        with open(output, 'w') as f:
            for index, read in enumerate(self.junctionMap):
                #         if index < 10:
                self.annotate_junction(read, result, self.database, f)