#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""class for detecting junction reads

"""
import pysam as ps

from .utils import get_yaml
from .utils import timethis

# change keys to be consist with chromosome's values
CHROMS = get_yaml()["chr2hg38"]

POSITIVE_SITE, NIGATIVE_SITE = [("GT", "AG"), ("AT", "AC"), ("GC", "AG")], [
    ("CT", "AC"),
    ("GT", "AT"),
    ("CT", "GC"),
]


class Read:
    """build a read class for storing information of every junction read

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

    __slots__ = (
        "chrom",
        "start",
        "end",
        "idn",
        "score",
        "strand",
        "anchor",
        "acceptor",
        "information",
    )

    def __init__(self, chrom, start, end, idn, score, strand, anchor, acceptor):
        self.chrom, self.start, self.end = chrom, start, end
        self.idn, self.score, self.strand = idn, score, strand
        self.anchor, self.acceptor = anchor, acceptor
        self.information = []

    @property
    def identifiers(self):
        return f"{self.chrom}_{self.start}_{self.end}"

    def __repr__(self):
        return (
            f"Read({self.chrom!r}, {self.start!r}, {self.end!r}, {self.idn!r}, "
            f"{self.score!r}, {self.strand!r}, {self.anchor!r}, {self.acceptor!r})"
        )

    def __str__(self):
        return "\t".join(
            map(
                str,
                [
                    self.chrom,
                    self.start,
                    self.end,
                    self.idn,
                    self.score,
                    self.strand,
                    f"{self.anchor}-{self.acceptor}",
                ],
            ),
        )


class JunctionMap:
    """build a class to store information of all junction reads"""

    def __init__(self):
        self.junctionList = {}

    def add_read(self, read):
        """add read to  junctionlist

        :param read: instance from Read
        :type read: instance
        """
        self.junctionList[read.identifiers] = read

    def get_read(self, identifiers):
        """get read from junctionlist according to identifiers

        :param identifiers: identifier for every read: chrom_start_end
        :type identifiers: int
        :return: instance from Read
        :rtype: instance
        """
        return self.junctionList[identifiers]

    def __contains__(self, identifiers):
        """check if read is in junctionlist according to identifiers

        :param identifiers: identifier for every read: chrom_start_end
        :type identifiers: str
        :return: whether Read is in JunctionMap
        :rtype: bool
        """
        return identifiers in self.junctionList

    def __iter__(self):
        """iterate every read in junctionlist"""
        return iter(self.junctionList.values())

    def write2file(self, output, header=None):
        """write all reads in junctionlist to file

        :param output: file name of output
        :type output: str or TextIo
        :param header: header of output
        :type header: str
        """
        try:
            output.write(f"{header}\n")
        except AttributeError:
            output = open(output, "w")
            output.write(f"{header}\n")
        finally:
            for read in self:
                output.write(f"{read}\n")
            output.close()


class JunctionDetector:
    """class for detecting junction reads and record position

    :param bam_file: bam file
    :type bam_file: str
    :param output: filename of output
    :type output: str
    :param reference: filename of genome reference
    :type reference: str
    :param quality: quality for filtering junction reads
    :type quality: int
    """

    def __init__(self, bam_file, reference, quality=0, output=None):

        self.bam_file, self.output = bam_file, output

        self.reference = reference

        self.quality = quality

    @staticmethod
    def check_strand(anchor, acceptor):
        """check type of strand

        :param anchor: anchor of read
        :type anchor: str
        :param acceptor: acceptor of read
        :type acceptor: str
        :return: type of strand (-|+)
        :rtype: str
        """
        strand = "N"

        if (anchor, acceptor) in POSITIVE_SITE:

            strand = "+"

        elif (anchor, acceptor) in NIGATIVE_SITE:

            strand = "-"

        return strand

    def worker(self, bam_file, reference, chrom, quality, idn, junctionmap):
        """find junction reads and annotate slice site

        :param bam_file: handle of bam_file
        :type bam_file: instance
        :param reference: handle of reference
        :type reference: instance
        :param chrom: chromosome
        :type chrom: str
        :param quality: quality for filtering reads
        :type quality: int
        :param idn: identifier of reads
        :type idn: int
        :param junctionmap: instance from junctionmap
        :type junctionmap: instance
        :return: instance from junctionmap
        :rtype: instance
        """
        # detect junction reads
        junction_regions = bam_file.find_introns(
            [
                r
                for r in bam_file.fetch(contig=chrom)  # chrom
                if r.mapping_quality > quality
            ],
        )

        # annotate slice sites
        for ((start, end), score) in junction_regions.items():
            junction_bases = reference.fetch(
                reference=CHROMS[chrom],
                start=start,
                end=end,  # change chrom
            )
            anchor, acceptor = junction_bases[:2].upper(), junction_bases[-2:].upper()

            strand = self.check_strand(anchor, acceptor)
            read = Read(chrom, start, end, idn + 1, score, strand, anchor, acceptor)
            junctionmap.add_read(read)
            #             line = f'{chrom}\t{start}\t{end}\t{idn+1}\t{score}\t{strand}\t{anchor}-{acceptor}'

            idn += 1

    @timethis(name="Junction detector", message="FINISHED")
    def run(self):
        """detect junction reads and annotate slice site, write results to file

        :return: instance from junctionmap
        :rtype: instance
        """

        idn = 0

        reference = ps.FastaFile(self.reference)

        bam = ps.AlignmentFile(self.bam_file)

        junctionmap = JunctionMap()

        # write junction reads information
        for chrom in CHROMS.keys():  # change CHROMS
            self.worker(bam, reference, chrom, self.quality, idn, junctionmap)

            # log.info(f'{chrom} FINSHED. NO PROBLEM FOUND')

        if self.output:
            junctionmap.write2file(self.output)  # write to file

        return junctionmap
