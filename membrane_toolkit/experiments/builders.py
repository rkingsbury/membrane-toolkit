# Copyright (c) Ryan Kingsbury
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

"""
Home for Builder classes that aggregate raw experimental data
"""

from maggma.builders import Builder


class ExcelBuilder(Builder):
    """
    Abstract class for builders that ingest Excel .xlsx files
    """


class TxtBuilder(Builder):
    """
    Abstract class for builders that ingest text-based files like .csv
    """
