#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# from click.testing import CliRunner
import sys

sys.path.append("..")
# from ce_detector import cli

#
# class TestCli:
#     def test_build(self):
#         runner = CliRunner()
#         with runner.isolated_filesystem():
#             with open('test_annotation.gff', 'w') as f:
#                 f.write('This is annotation file')
#
#             _result = runner.invoke(cli, [
#                 'build',
#                 'test_annotation.gff',
#             ])
#             assert _result.exit_code == 0
#             assert 'database' in _result.output
#
#     def test_detect(self):
#         runner = CliRunner()
#
#         with runner.isolated_filesystem():
#             with open('test_bam.bam', 'w') as f:
#                 f.write('This is bam file')
#             with open('test_reference.fa', 'w') as f:
#                 f.write('This is genome reference file')
#             with open('test_gffdb.db', 'w') as f:
#                 f.write('This is gffdb file')
#
#             _result = runner.invoke(cli, [
#                 'detect',
#                 '--bam', 'test_bam.bam',
#                 '--reference', 'test_reference.fa',
#                 '--gffdb', 'test_gffdb.db',
#                 '--quality', '30',
#                 '--out', 'test_out.bed'
#             ])
#
#             assert _result.exit_code == 0
#             assert 'test_bam.bam' in _result.output
#             assert 'test_reference.fa' in _result.output
#             assert 'test_gffdb.db' in _result.output
#             assert 'quality=30' in _result.output
#             assert 'test_out.bed' in _result.output
