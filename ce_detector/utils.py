import argparse
import json
import logging
import time
from os.path import join, dirname

import importlib_resources


class Timer:
    """construct  Timer to show working time of tasks
    """

    def __init__(self, func=time.perf_counter):
        """init values

        Args:
            func : Defaults to time.perf_counter.
        """
        self.elapsed = 0.0

        self._func = func

        self._start = None

    def start(self):
        """start a task

        Raises:
            RuntimeError:[if task has started then raise error]
        """

        if self._start is not None:
            raise RuntimeError('Already started')

        self._start = self._func()

    def stop(self):
        """end a task

        Raises:
            RuntimeError[if task has not started then raise error]
        """
        if self._start is None:
            raise RuntimeError('Not started')

        end = self._func()

        self.elapsed += end - self._start

        self._start = None

    def reset(self):
        """reset the working time
        """
        self.elapsed = 0

    def running(self):
        """check if task is running

        Returns:
            bool
        """

        return self._start is not None

    def __enter__(self):
        """function used to address 'with text'
        """

        self.start()

        return self

    def __exit__(self, *args):
        """function used to address 'with text'
        """

        self.stop()


def get_logger(logger_name, create_file=False):
    """ set logger and output console

    :param logger_name: logger name
    :type logger_name: str
    :param create_file: whether creat file to store log
    :type create_file: bool
    :return: logger
    :rtype: instance
    """
    # create logger for prd_ci
    log = logging.getLogger(logger_name)

    log.setLevel(level=logging.INFO)

    # create formatter and add it to the handlers
    formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s')

    if create_file:
        # create file handler for logger.
        fh = logging.FileHandler(f'{logger_name}_log.txt')

        fh.setLevel(level=logging.DEBUG)

        fh.setFormatter(formatter)

    # relate console handler for logger.
    ch = logging.StreamHandler()
    ch.setLevel(level=logging.DEBUG)
    ch.setFormatter(formatter)

    # add handlers to logger.
    if create_file:
        log.addHandler(fh)

    log.addHandler(ch)

    return log


# TODO: add subcommand to create database for gtf file
def get_parser():
    """ get parameter from terminal
    """

    parser = argparse.ArgumentParser(prog='PROG', description='Program designed to find CE',
                                     formatter_class=argparse.MetavarTypeHelpFormatter)
    # get bam file
    parser.add_argument('-i',
                        '--input',
                        help='The input file (Bam or Sam format)',
                        required=True,
                        type=str)
    # get genome reference
    parser.add_argument(
        '-r',
        '--reference',
        help='The reference fasta (.fna) file, which contains index file(*.fai).',
        required=True,
        type=str)
    # get filename of junction reads
    parser.add_argument('-o',
                        '--out',
                        help='The output file.',
                        default='junctions.bed',
                        type=str)
    # get quality for filtering junction reads
    parser.add_argument(
        '-q',
        '--quality',
        help='The threshold to filter low quality reads;Default:0',
        default=0,
        type=int)
    # get filename of database of gtf file
    parser.add_argument('-d',
                        '--gffdb',
                        help='The annotated file',
                        required=True,
                        type=str)

    args = parser.parse_args()

    return args


def get_json():
    """ get information of chromosome stored in json file
    :return: chromosome values
    :rtype: dict
    """
    if __package__:
        path = join(importlib_resources.files(__package__).as_posix(), 'chromosome.json')
    else:
        path = join(dirname(__file__), 'chromosome.json')
    return json.load(open(path))
