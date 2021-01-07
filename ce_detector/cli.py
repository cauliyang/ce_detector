#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: YangyangLi
@contact:li002252@umn.edu
@version: 0.0.1
@license: MIT Licence
@file: cli.py.py
@time: 2020/12/28 10:21 PM
"""
import click
import gffutils
from rich.traceback import install

from . import __version__
from .annotator import Annotator
from .detector import JunctionDetector
from .scanner import Scanner

install()

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"], max_content_width=150)


@click.group(
    context_settings=CONTEXT_SETTINGS,
    options_metavar="<options>",
    subcommand_metavar="<command>",
)
@click.version_option(__version__)
def cli():
    """program designed for detecting cryptic exons

    Needed file types:

    \b
    1. bam file
    2. genome reference
    3. annotation file
    """
    pass


@cli.command(
    "build",
    short_help="build database for annotation file",
    options_metavar="<options>",
)
@click.argument("gff", type=click.Path(exists=True))
@click.option(
    "--out-directory",
    "-o",
    type=click.Path(exists=True),
    default=".",
    show_default=True,
    metavar="<path>",
)
def build(gff, out_directory):
    """build database for annotation file

    build database of annotation file in order to use the database to annotate junctions reads later.
    this process is time-confusing so that you may need to prepare a cup of coffee!
    file of database named {prefix of annotation file}.db

    \f
    :param gff: the path of annotation file
    :type gff: str
    :param out_directory: the path of _result of database
    :type out_directory: str
    :return: {out directory}/{prefix of annotation file}.db
    """
    _ = gffutils.create_db(
        gff,
        out_directory,
        merge_strategy="create_unique",
        keep_order=True,
    )


@cli.command("detect", short_help="scan cryptic exons", options_metavar="<options>")
@click.option(
    "--bam",
    "-b",
    help="The bam file (Bam or Sam format)",
    required=True,
    type=click.Path(exists=True),
    metavar="<path>",
)
@click.option(
    "--reference",
    "-r",
    help="The reference fasta (.fna) file, which contains index file(*.fai)",
    required=True,
    type=click.Path(exists=True),
    metavar="<path>",
)
@click.option(
    "--quality",
    "-q",
    help="The threshold to filter low quality reads",
    default=0,
    type=click.INT,
    show_default=True,
    metavar="<int>",
)
@click.option(
    "--out",
    "-o",
    help="The output file of detected cryptic exons",
    default="cryptic_exons.bed",
    type=click.STRING,
    show_default=True,
    metavar="<path>",
)
@click.option(
    "--gffdb",
    "-db",
    help="The database of annotation file",
    type=click.Path(exists=True),
    required=True,
    metavar="<path>",
)
@click.option(
    "--cutoff",
    "-c",
    help="cutoff for filtering junction reads wth low depth during scanning cryptic exons",
    type=click.INT,
    default=1,
    show_default=True,
    metavar="<int>",
)
@click.option(
    "--out-ann",
    "-oa",
    help="The output file of annotated junction reads",
    default="annotated_junctions.bed",
    type=click.File(mode="w", encoding="utf-8"),
    show_default=True,
    metavar="<path>",
)
def detect(bam, reference, quality, gffdb, cutoff, out, out_ann):
    """detect junction reads and scan cryptic exons

    \b
    analysis following steps below:
    1. detect junction reads in terms of bam file
    2. annotate junction reads in terms of genome reference and annotation file
    3. scan cryptic exons according to its definition

    \f
    :param out_ann:
    :type out_ann:
    :param cutoff:
    :type cutoff:
    :param bam: bam file
    :type bam: str
    :param reference: genome reference file
    :type reference: str
    :param quality: quality used to filter low quality reads.Default:0
    :type quality: int
    :param gffdb: database file of annotation file
    :type gffdb: str
    :param out: file name of detected cryptic exons
    :type out: str
    """

    # click.secho((f'\nProgram Start.\nParameters:\nReference: {reference}'
    #              f'\nQuality Threshold: {quality}\nOutput: {out}\n '), fg='red', bold=True)

    detector = JunctionDetector(
        bam,
        reference,
        quality,
    )
    junctionmap = detector.run()

    annotator = Annotator(junctionmap, gffdb, out_ann)
    annotator.run()

    scanner = Scanner(cutoff=cutoff, output=out)
    scanner.run(annotator.junctionMap)
    scanner.write2file()


if __name__ == "__main__":
    cli()
