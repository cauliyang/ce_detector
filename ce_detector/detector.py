#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""class for detecting junction reads

"""
import re
import time

import numpy as np
import pysam as ps
from statsmodels.stats.multitest import fdrcorrection

from .utils import timethis


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
        "pvalue",
    )

    def __init__(self, chrom, start, end, idn, score, strand, anchor, acceptor, pvalue):
        self.chrom, self.start, self.end = chrom, start, end
        self.idn, self.score, self.strand = idn, score, strand
        self.anchor, self.acceptor = anchor, acceptor
        self.information = []
        self.pvalue = pvalue

    @property
    def identifiers(self):
        return f"{self.chrom}_{self.start}_{self.end}"

    def __repr__(self):
        return (
            fr"Read({self.chrom}, {self.start}, {self.end}, {self.idn}, "
            fr"{self.score}, {self.strand}, {self.anchor}, {self.acceptor})"
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

    def __init__(self, chrom):
        self.chrom = chrom
        self._junctionList = []
        self.result = None

    @classmethod
    def build(cls, chrom, data):
        new_instance = cls(chrom)
        new_instance.add_data(data)
        return new_instance

    def add_data(self, data):
        self._junctionList = data

    @property
    def data(self):
        return self._junctionList

    def is_empty(self):

        if len(self._junctionList):
            return False
        else:
            return True

    def size(self):
        return len(self._junctionList)

    def add_read(self, read):
        """add read to  junctionlist

        :param read: instance from Read
        :type read: instance
        """
        self._junctionList.append(read)

    def __contains__(self, read):
        """check if read is in junctionlist according to identifiers

        :type read: object

        :return: whether Read is in JunctionMap
        :rtype: bool
        """
        return read in self._junctionList

    def __iter__(self):
        """iterate every read in junctionlist"""
        for item in self._junctionList:
            yield item

    def __repr__(self):
        return f"ce_detector.detector.JunctionMap(chrom = {self.chrom})"

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

    SPLICE_SITE = dict(
        zip([("GT", "AG"), ("AT", "AC"), ("GC", "AG")], "+++")
    )  # positive site
    SPLICE_SITE.update(
        dict(zip([("CT", "AC"), ("GT", "AT"), ("CT", "GC")], "---"))
    )  # negative site

    PATTERN = re.compile(r"\d*?S*(\d+)M(\d+)N(\d+)M")

    def __init__(self, bam_file, reference, quality, output=None):

        self.bam, self.reference = ps.AlignmentFile(bam_file), ps.FastaFile(reference)

        self.output, self.quality = output, quality

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

        strand = JunctionDetector.SPLICE_SITE.get((anchor, acceptor), "N")

        return strand

    @staticmethod
    def fdr_correction(junctionmaps):
        chroms = junctionmaps.keys()
        junctionmaps_len = [jmap.size() for jmap in junctionmaps.values()]
        split_ind = np.cumsum(junctionmaps_len)[:-1]
        all_pvalues = np.array(
            [read.pvalue for jmap in junctionmaps.values() for read in jmap]
        )
        cond, fdr = fdrcorrection(all_pvalues)
        cond_pvalue = np.split(cond, split_ind)

        for ind, chrom in enumerate(chroms):
            data = np.array(junctionmaps[chrom].data)[cond_pvalue[ind]]
            junctionmaps[chrom] = JunctionMap.build(chrom, data.tolist())

        return junctionmaps

    def get_pvalue(self, chrom, start, end):

        jun_reads = self.bam.fetch(contig=chrom, start=start, end=end)
        max_ratio = float("-inf")
        max_info = tuple()

        for read in jun_reads:

            if "N" in read.cigarstring and read.mapping_quality > self.quality:

                block1, gap, block2 = map(
                    int, self.PATTERN.search(read.cigarstring).groups()
                )
                block_ratio = min(block1, block2) / max(block1, block2)
                if block_ratio > max_ratio:
                    max_ratio = block_ratio
                    max_info = (block1, gap, block2)

        R = min(max_info[0], max_info[-1])
        L = max_info[1]
        p_value = 1 - (1 - (1 / 4) ** R) ** (L - R + 1)

        return max_info, p_value

    def worker(self, bam_file, reference, chrom, ann_chrom, quality, idn, junctionmap):
        """find junction reads and annotate slice site

        :param ann_chrom:
        :type ann_chrom:
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
                reference=ann_chrom,
                start=start,
                end=end,  # change chrom
            )
            anchor, acceptor = junction_bases[:2].upper(), junction_bases[-2:].upper()
            strand = self.check_strand(anchor, acceptor)
            _, p_value = self.get_pvalue(chrom, start, end)
            read = Read(
                chrom, start, end, idn + 1, score, strand, anchor, acceptor, p_value
            )
            junctionmap.add_read(read)

            idn += 1

        return junctionmap

    @timethis(name="Junction detector", message=" ")
    def run(self, chrom, ann_chrom, logger, verbose=False):
        """detect junction reads and annotate slice site, write results to file

        :return: instance from junctionmap
        :rtype: instance
        """

        idn = 0

        junctionmap = JunctionMap(chrom)

        # write junction reads information
        if verbose:
            logger.info(f"Chrom {chrom} Beginning")
            start = time.time()

        junctionmap = self.worker(
            self.bam, self.reference, chrom, ann_chrom, self.quality, idn, junctionmap
        )

        if verbose:
            logger.info(f"Chrom {chrom} Finished {time.time() - start:.2f}s")

        return junctionmap
