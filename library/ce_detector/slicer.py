#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
"""
from __future__ import annotations

from collections import defaultdict

import gffutils
import numpy as np
import pysam as ps

from utils import get_logger, get_parser, Timer

# TODO: Change annotate to a class

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

POSITIVE_SITE, NIGATIVE_SITE = [('GT', 'AG'), ('AT', 'AC'), ('GC', 'AG')], \
                               [('CT', 'AC'), ('GT', 'AT'), ('CT', 'GC')]


class Read:
    """ build a read class for storing information of junction reads
    """
    __slots__ = [
        'chrom', 'start', 'end', 'idn', 'score', 'strand', 'anchor',
        'acceptor', 'identifiers'
    ]

    def __init__(self, chrom, start, end, idn, score, strand, anchor,
                 acceptor):
        """
        :param chrom: chromosome of genome
        :type chrom: str
        :param start: start position of junction read
        :type start: int
        :param end: end position of junction read
        :type end: int
        :param idn: index of junction read
        :type idn: int
        :param score: support of junction read
        :type score: int
        :param strand: direction of junction read (-|+)
        :type strand: str
        :param anchor: anchor of junction read
        :type anchor: str
        :param acceptor: acceptor of junction read
        :type acceptor: str
        """
        self.chrom, self.start, self.end = chrom, start, end
        self.idn, self.score, self.strand = idn, score, strand
        self.anchor, self.acceptor = anchor, acceptor
        self.identifiers = f'{self.chrom}_{self.start}_{self.end}'

        # self.gene, self.type, self.donorSkipped, self.acceptorSkipped = [None
        #                                                                  ] * 4

    def __repr__(self):
        return '\t'.join(
            map(str, [
                self.chrom, self.start, self.end, self.idn, self.score,
                self.strand, f'{self.acceptor}-{self.acceptor}'
            ]))


class JunctionMap:
    """ build a class to store information of all junction reads
    """

    def __init__(self):
        self.junctionList = {}

    def add_read(self, read):
        """ add read to  junctionlist

        :param read: instance from Read
        :type read: instance
        """
        self.junctionList[read.identifiers] = read

    def get_read(self, identifiers):
        """ get read from junctionlist according to identifiers

        :param identifiers: identifier for every read: chrom_start_end
        :type identifiers: int
        :return: instance from Read
        :rtype: instance
        """
        return self.junctionList[identifiers]

    def __contains__(self, identifiers):
        """ check if read is in junctionlist according to identifiers

        :param identifiers: identifier for every read: chrom_start_end
        :type identifiers: str
        :return: whether Read is in JunctionMap
        :rtype: bool
        """
        return identifiers in self.junctionList

    def __iter__(self):
        """ iterate every read in junctionlist
        """
        return iter(self.junctionList.values())

    def write2file(self, output, header=None):
        """ write all reads in junctionlist to file

        :param output: file name of output
        :type output: str
        :param header: header of output
        :type header: str
        """
        with open(output, 'w') as f:
            f.write(f'{header}\n')
            for read in self:
                f.write(f'{read}\n')


#                 print(f'{read[0]}\n')

class JunctionDetector:
    """class for detecting junction reads and record position
    """

    def __init__(self,
                 bam_file=None,
                 output=None,
                 reference=None,
                 quality=None):
        """
        :param bam_file: bam file
        :type bam_file: str
        :param output: filename of output
        :type output: str
        :param reference: filename of genome reference
        :type reference: str
        :param quality: quality for filtering junction reads
        :type quality: int
        """
        self.bam_file, self.output = bam_file, output

        self.reference = reference

        self.quality = quality

    @staticmethod
    def check_strand(anchor: str, acceptor: str) -> str:
        """ check type of strand

        :param anchor: anchor of read
        :type anchor: str
        :param acceptor: acceptor of read
        :type acceptor: str
        :return: type of strand (-|+)
        :rtype: str
        """
        strand = 'N'

        if (anchor, acceptor) in POSITIVE_SITE:

            strand = '+'

        elif (anchor, acceptor) in NIGATIVE_SITE:

            strand = '-'

        return strand

    def worker(self, bam_file, reference, chrom, quality, idn, junctionMap):
        """ find junction reads and annotate slice site

        :param bam_file: handle of bam_file
        :type bam_file: instance
        :param reference: handle of reference
        :type reference: instance
        :param chrom: chromosome
        :type chrom: int
        :param quality: quality for filtering reads
        :type quality: int
        :param idn: identifier of reads
        :type idn: int
        :param junctionMap: instance from junctionMap
        :type junctionMap: instance
        :return: instance from junctionMap
        :rtype: instance
        """
        # detect junction reads
        junction_regions = bam_file.find_introns([
            r for r in bam_file.fetch(contig=chrom)
            if r.mapping_quality > quality  # chrom
        ])

        # annotate slice sites
        for ((start, end), score) in junction_regions.items():
            junction_bases = reference.fetch(reference=CHROMS[chrom],  # change chrom
                                             start=start,
                                             end=end)
            anchor, acceptor = junction_bases[:2].upper(), junction_bases[-2:].upper()

            strand = self.check_strand(anchor, acceptor)
            read = Read(chrom, start, end, idn + 1, score, strand, anchor,
                        acceptor)
            junctionMap.add_read(read)
            #             line = f'{chrom}\t{start}\t{end}\t{idn+1}\t{score}\t{strand}\t{anchor}-{acceptor}'

            idn += 1
        return junctionMap

    def run(self, log):
        """
        :param log:
        :type log:
        :return:
        :rtype:
        """

        ID = 0

        #         output = open(self.output, 'w') # change

        reference = ps.FastaFile(self.reference)

        bam = ps.AlignmentFile(self.bam_file)

        junctionMap = JunctionMap()

        # write junction reads information
        log.info(
            f'Junction Detector Start.\nParameters:\nReference: {self.reference}\nQuality Threshold: {self.quality}\nOutput: {self.output}\n'
        )
        for chrom in CHROMS.keys():  # change CHROMS
            self.worker(bam, reference, chrom, self.quality, ID, junctionMap)

            log.info(f'{chrom} FINSHED. NO PROBLEM FOUND')

        if self.output:
            header = f'chrom\tstart\tend\tidn\tscore\tstrand\tsplice_site'
            junctionMap.write2file(self.output, header)

        return junctionMap


class Annotator:
    """annotate junction reads"""

    def __init__(self, junctionMap, database):
        """Constructor for Annotator"""
        self.junctionMap = junctionMap
        self.database = database

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
        :param output: filename of output. Defaults to None
        :type output: str
        """
        chrom = read.chrom
        start, end = read.start, read.end
        region = f'{CHROMS[chrom]}:{start}-{end}'

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
        db = gffutils.FeatureDB(self.database)

        result = defaultdict(list)

        logger.info(f'Begin to annotate junctions!')

        with open(output, 'w') as f:
            for index, read in enumerate(self.junctionMap):
                #         if index < 10:
                self.annotate_junction(read, result, self.database, f)


def main():
    """ function to integrate the annotation usage
    """

    log = get_logger('junction_detector', create_file=False)
    log2 = get_logger('annotate_junctions', create_file=False)

    args = get_parser()

    # detect junction reads
    detector = JunctionDetector(args.input, args.out, args.reference,
                                args.quality)
    with Timer() as t:
        junctionMap = detector.run(log)

    log.info(f'FINISHED FIND JUNCTION CONSUMING {t.elapsed:.2f}s')
    # annotate junction reads

    annotator = Annotator(junctionMap, args.gffdb)
    annotator.run()


if __name__ == '__main__':
    main()
