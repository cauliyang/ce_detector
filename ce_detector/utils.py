# -*- coding: utf-8 -*-
import logging
import time
from functools import partial
from functools import wraps
from os.path import dirname
from os.path import join

import importlib_resources
import yaml
from rich.logging import RichHandler


class Timer:
    """construct  Timer to show working time of tasks"""

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
            raise RuntimeError("Already started")

        self._start = self._func()

    def stop(self):
        """end a task

        Raises:
            RuntimeError[if task has not started then raise error]
        """
        if self._start is None:
            raise RuntimeError("Not started")

        end = self._func()

        self.elapsed += end - self._start

        self._start = None

    def reset(self):
        """reset the working time"""
        self.elapsed = 0

    def running(self):
        """check if task is running

        Returns:
            bool
        """

        return self._start is not None

    def __enter__(self):
        """function used to address 'with text'"""

        self.start()

        return self

    def __exit__(self, *args):
        """function used to address 'with text'"""

        self.stop()


def rich_logger(logger_name, create_file=False, level=logging.INFO):
    """set logger and output console

    :param level: level for logging
    :param logger_name: logger name
    :type logger_name: str
    :param create_file: whether creat file to store log
    :type create_file: bool
    :return: logger
    :rtype: instance
    """
    # create logger for prd_ci
    log = logging.getLogger(logger_name)

    # create formatter and add it to the handlers
    template = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    formatter = logging.Formatter(template)

    logging.basicConfig(
        level=level,
        format=template,
        handlers=[
            RichHandler(markup=True, show_level=False, show_time=False),
        ],
    )

    if create_file:
        # create file handler for logger.
        fh = logging.FileHandler(f"{logger_name}_log.txt")

        fh.setLevel(level=logging.DEBUG)

        fh.setFormatter(formatter)

    # add handlers to logger.
    if create_file:
        log.addHandler(fh)

    return log


def get_logger(logger_name, create_file=False, level=logging.INFO):
    """set logger and output console

    :param level: level for logging
    :param logger_name: logger name
    :type logger_name: str
    :param create_file: whether creat file to store log
    :type create_file: bool
    :return: logger
    :rtype: instance
    """
    # create logger for prd_ci
    log = logging.getLogger(logger_name)

    log.setLevel(level=level)

    # create formatter and add it to the handlers
    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    )

    if create_file:
        # create file handler for logger.
        fh = logging.FileHandler(f"{logger_name}_log.txt")

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


def get_yaml():
    """get information of chromosome stored in yaml file
    :return: chromosome values
    :rtype: dict
    """
    if __package__:
        path = join(importlib_resources.files(__package__).as_posix(), "chromosome.yml")
    else:
        path = join(dirname(__file__), "chromosome.yml")
    return yaml.safe_load(open(path))


def timethis(
    func=None,
    level=logging.INFO,
    name=None,
    message=None,
    creat_file=False,
    rich=True,
):
    if func is None:
        return partial(
            timethis,
            level=level,
            name=name,
            message=message,
            creat_file=creat_file,
        )

    logname = name if name else func.__module__
    log = (
        rich_logger(logname, create_file=creat_file)
        if rich
        else get_logger(logname, create_file=creat_file)
    )
    logmsg = message if message else func.__name__

    @wraps(func)
    def wrapper(*args, **kwargs):
        with Timer() as t:
            temp = func(*args, **kwargs)

        log.info(f"[bold green]{logmsg} CONSUMING {t.elapsed:.2f}s")
        return temp

    return wrapper
