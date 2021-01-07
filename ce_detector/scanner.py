#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" class for scanning cryptic exons based on previous _result ::
Junction detector and junction annotator
"""
from typing import Iterable

import numpy as np
import pandas as pd

from .utils import timethis


def check(axis) -> np.array:
    """check if cryptic exon has children

    Which means that check if the cryptic exon is split by others junction reads
    in terms of start and end position

    :param axis: an array, a junction read, included start and end
    :type axis: numpy.array
    :return: whether junction read can split cryptic exon
    :rtype: bool
    """

    return (axis[0] < 0) & (axis[1] > 0)


def assign_value(df_ce, ces, ns, ce_id, ns_id) -> None:
    """assign value of `child column` for every cryptic exons that contains junction reads with N type

    Given the start and end of cryptic exons and junction reads,
    Note: types of junction reads contains N, D, A, DA, NDA. For details:
    https://regtools.readthedocs.io/en/latest/commands/junctions-annotate/

    :param df_ce: pandas.DataFrame of cryptic exons
    :type df_ce: pandas.DataFrame
    :param ces: cryptic exons' pandas.DataFrame,
        which has gene ids that both cryptic exon and junction
        read own. It only contain gene id, start and end
    :type ces: pandas.DataFrame
    :param ns: pandas.DataFrame of junction reads whose type is N,
        as well as has same gene id as cryptic exons
    :type ns: pandas.DataFrame
    :param ce_id: index of gene id of cryptic exons which has junction reads as children
    :type ce_id: numpy.array
    :param ns_id: index of gene id of junction reads with N type, which have junction reads as children
    :type ns_id: numpy.array
    """
    children_info = f"{ns.iloc[ns_id, 0]}-{ns.iloc[ns_id, 1]},"  # start-end
    gene_id = ces.index[ce_id]
    length = ces.iloc[ce_id, 1] - ces.iloc[ce_id, 0]
    df_ce.loc[(gene_id, length), "children"] += children_info


def split_ce(df_ce, df_n) -> Iterable:
    """Iterator: check whether detected cryptic exons are split by other junction reads

    :param df_ce: pandas.DataFrame of cryptic exons return from :func:`find_ce`
    :type df_ce: pandas.DataFrame
    :param df_n: pandas.DataFrame of junction reads with N type
    :type df_n: pandas.DataFrame
    :return: param:`df_ce` with new column `children`
    :rtype: iterator
    """
    # set gene as index
    df_ce.set_index("gene", inplace=True)
    df_n.set_index("gene", inplace=True)
    # get genes contained in both cryptic exons and junction reads with N type
    gene_ind = set(df_ce.index) & set(df_n.index)

    ces = df_ce.loc[gene_ind, ["start", "end"]]
    ns = df_n.loc[gene_ind, ["start", "end"]]
    ces_values = ces.values
    ns_values = ns.values
    # set another index in case that one gene has more than one cryptic exons
    df_ce.set_index("length", append=True, inplace=True)

    if ces_values.size > 1:

        x = np.apply_along_axis(
            lambda row: row - ns_values,
            1,
            ces_values,
        )  # iter for row

        bool_re = np.apply_along_axis(check, 2, x)

        position = np.argwhere(
            bool_re,
        )  # return row id, n_df id[[ce_id,n_id],[ce_id,n_id]]

        if position.size >= 1:

            for ce_id, ns_id in position:
                assign_value(df_ce, ces, ns, ce_id, ns_id)

    else:

        x = ces_values - ns_values
        bool_x = (x[:, 0] < 0) & (x[:, 1] > 0)
        position = np.argwhere(bool_x)  # return [[n_id],[n_id]]

        if position.size >= 1:
            for ns_id in position.flatten():
                assign_value(df_ce, ces, ns, ce_id=0, ns_id=ns_id)

    yield df_ce


def find_ce(groups) -> Iterable:
    """parse _result getting from annotations in order to detect cryptic exons

    :param groups:  annotated junction reads are grouped by `strand` and `type`
    :type groups: Groupby object return from ``Dataframe.groupby``
    :return: pd.DataFrame of two strands
    :rtype: Iterable
    """
    # old column  and new column
    pick_col = (
        "chrom",
        "end_D",
        "start_A",
        "start",
        "end",
        "strand",
        "score_DA",
        "score_D",
        "score",
        "gene",
    )
    rename_col = (
        "chrom",
        "start",
        "end",
        "start_DA",
        "end_DA",
        "strand",
        "score_DA",
        "score_D",
        "score_A",
        "gene",
    )

    for strand in ("+", "-"):
        da = groups.get_group((strand, "DA"))
        d = groups.get_group((strand, "D"))
        a = groups.get_group((strand, "A"))
        n = groups.get_group((strand, "N"))

        yield from (
            da.merge(
                d,
                left_on=["chrom", "start", "gene"],
                right_on=["chrom", "start", "gene"],
                suffixes=("_DA", "_D"),
            )
            .merge(
                a,
                left_on=["chrom", "end_DA", "gene"],
                right_on=["chrom", "end", "gene"],
                suffixes=("", "_A"),
            )[lambda df: df["end_D"] < df["start_A"]]
            .pipe(lambda df: df.loc[:, pick_col])
            .rename(columns=dict(zip(pick_col, rename_col)))
            .assign(length=lambda df: df.end - df.start)
            .assign(children="")
            .pipe(split_ce, n)
        )


class Scanner:
    """class for scanning cryptic exons based on annotated junction reads

    :param cutoff: cutoff used to filter junction reads with relatively low score or depth
    :type cutoff: int
    :param output: filename of _result
    :type output: str
    """

    def __init__(self, cutoff, output):
        """Constructor for Scanner"""
        self.cutoff = cutoff
        self.output = output
        self._result = None

    def __repr__(self):
        return f"Scanner({self.cutoff!r})"

    @timethis(name="File Writer for Scanner", message="FINISHED")
    def write2file(self):
        """ start iterator and write _result to file"""
        pd.concat(self._result).to_csv(self.output, sep="\t", encoding="utf8")

    @timethis(name="Cryptic Exon Scanner", message="FINISHED")
    def run(self, junctionmap) -> Iterable:
        """run program to scan cryptic exons

        :param junctionmap: instance from :class:`ce_detector.detector.JunctionMap`
        :type junctionmap: instance
        :return: temporary result used to store cryptic exons
        :rtype: Iterable
        """
        self._result = (
            pd.DataFrame(
                [
                    [read.chrom, read.start, read.end, read.strand, read.score, *info]
                    for read in junctionmap
                    for info in read.information
                    if read.score >= self.cutoff
                ],
                columns=[
                    "chrom",
                    "start",
                    "end",
                    "strand",
                    "score",
                    "type",
                    "dk",
                    "ak",
                    "gene",
                ],
            )
            .pipe(lambda df: df.groupby(["strand", "type"]))
            .pipe(find_ce)
        )

        return self._result
