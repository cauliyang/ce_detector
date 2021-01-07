#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Tests for `ce_detector` package."""
import pytest


# from ce_detector import ce_detector


@pytest.fixture
def response():
    """Sample pytest fixture.
    ß
        See more at: http://doc.pytest.org/en/latest/fixture.html
    """
    # import requests
    # return requests.get('https://github.com/audreyr/cookiecutter-pypackage')


def test_content(response):
    """Sample pytest test function with the pytest fixture as an argument."""
    # from bs4 import BeautifulSoup
    # assert 'GitHub' in BeautifulSoup(response.content).title.string
