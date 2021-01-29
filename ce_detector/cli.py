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
from concurrent import futures

import click
import gffutils
import pandas as pd
from rich.traceback import install

from . import __version__
from .main import main
from .utils import get_yaml

install()

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"], max_content_width=150)


@click.group(
    context_settings=CONTEXT_SETTINGS,
    options_metavar="<options>",
    subcommand_metavar="<command>",
)
@click.version_option(__version__)
@click.option("--verbose", is_flag=True, default=False, help="Show the verbose mode")
@click.pass_context
def cli(ctx, verbose):
    """program designed for detecting cryptic exons

    Needed file types:

    \b
    1. bam file
    2. genome reference
    3. annotation file
    """
    ctx.ensure_object(dict)
    ctx.obj["verbose"] = verbose


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
@click.pass_context
def build(ctx, gff, out_directory):
    """build database for annotation file

    build database of annotation file in order to use the database to annotate junctions reads later.
    this process is time-confusing so that you may need to prepare a cup of coffee!
    file of database named {prefix of annotation file}.db

    \f
    :param ctx: click context used to pass parameters
    :type ctx: ``click.Context``
    :param gff: the path of annotation file
    :type gff: str
    :param out_directory: the path of _result of database
    :type out_directory: str
    :return: {out directory}/{prefix of annotation file}.db
    """
    verbose = ctx.obj["verbose"]
    if verbose:
        pass
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
    default=30,
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
@click.option("--parallel", is_flag=True, default=False, help="using parallel mode")
@click.pass_context
def detect(ctx, bam, reference, quality, gffdb, cutoff, out, parallel):
    """detect junction reads and scan cryptic exons

    \b
    analysis following steps below:
    1. detect junction reads in terms of bam file
    2. annotate junction reads in terms of genome reference and annotation file
    3. scan cryptic exons according to its definition

    \f
    :param parallel:
    :type parallel:
    :param ctx: click context used to pass parameters
    :type ctx: ``click.Context``
    :param cutoff: threshold for filtering junction reads with low quality
    :type cutoff: int
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

    verbose = ctx.obj["verbose"]

    # change keys to be consist with chromosome's values
    CHROMS = get_yaml()["chr2hg38"]

    worker = (
        futures.ProcessPoolExecutor()
        if parallel
        else futures.ThreadPoolExecutor(max_workers=24)
    )

    tasks = []
    junctionmaps = []
    with worker as executor:
        for chrom, ann_chrom in CHROMS.items():
            future = executor.submit(
                main, chrom, ann_chrom, bam, reference, quality, gffdb, cutoff, verbose
            )
            tasks.append(future)
        for future in futures.as_completed(tasks):
            junctionmaps.append(future.result())

    final_result = [
        jmap.result for jmap in junctionmaps if isinstance(jmap.result, pd.DataFrame)
    ]
    pd.concat(final_result).to_csv(out, sep="\t", encoding="utf8", index=False)


if __name__ == "__main__":
    cli()
