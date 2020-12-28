#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: YangyangLi
@contact:li002252@umn.edu
@version: 0.0.1
@license: MIT Licence
@file: re_cli.py.py
@time: 2020/12/28 10:21 PM
"""
import click
import gffutils


@click.group()
def cli():
    """ program designed for detecting cryptic exons

    Needed file types:

    \b
    1. bam file
    2. genome reference
    3. annotation file

    return: result of cryptic exons (BED)
    """
    # click.echo(f'This is ce detector')
    pass


@cli.command('build',short_help='build database for annotation file')
@click.argument('gff', type=click.Path(exists=True))
@click.option('--out-directory', '-o',
              type=click.Path(exists=True),
              default='.',
              show_default=True
              )
def build(gff, out_directory):
    """ build database for annotation file

    build database of annotation file in order to use the database to annotate junctions reads later.
    this process is time-confusing so that you may need to prepare a cup of coffee!
    file of database named {prefix of annotation file}.db

    \b
    :param gff: the path of annotation file
    :type gff: str
    :param out_directory: the path of result of database
    :type out_directory: str
    :return: {out directory}/{prefix of annotation file}.db
    """
    # _ = gffutils.create_db(gff, out_directory, merge_strategy='create_unique', keep_order=True)
    click.echo(f'gff={gff} out_directory={out_directory}')
    click.echo(f'Finished building database!')


@cli.command('detect', short_help='scan cryptic exons')
@click.option('--bam', '-b', help='The bam file (Bam or Sam format)', required=True, type=click.Path(exists=True))
@click.option('--reference',
              '-r',
              help='The reference fasta (.fna) file, which contains index file(*.fai)',
              required=True,
              type=click.Path(exists=True))
@click.option('--quality', '-q',
              help='The threshold to filter low quality reads',
              default=0,
              type=click.INT,
              show_default=True)
@click.option('--out', '-o',
              help='The output file of detected cryptic exons',
              default='cryptic_exons.bed',
              type=click.File('w', encoding='utf-8'),
              show_default=True)
@click.option('--gffdb', '-db',
              help='The database of annotation file',
              type=click.Path(exists=True),
              required=True
              )
def detect(bam, reference, quality, gffdb, out):
    """detect junction reads and scan cryptic exons

    \b
    analysis following steps below:
    1. detect junction reads in terms of bam file
    2. annotate junction reads in terms of genome reference and annotation file
    3. scan cryptic exons according to its definition

    \b
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
    click.echo(f'bam={bam} reference={reference} quality={quality} gffdb={gffdb} out={out}')


if __name__ == '__main__':
    cli()
