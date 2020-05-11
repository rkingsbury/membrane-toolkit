# membrane-toolkit

[![testing](https://github.com/rkingsbury/membrane-toolkit/workflows/testing/badge.svg)](https://github.com/rkingsbury/membrane-toolkit/actions?query=workflow%3Atesting) [![codecov](https://codecov.io/gh/rkingsbury/membrane-toolkit/branch/master/graph/badge.svg?token=I7RP0QML6S)](https://codecov.io/gh/rkingsbury/membrane-toolkit) [![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

## What is membrane-toolkit?

membrane-toolkit is a suite of Python tools for membrane science. There are currently two primary ways to
use it:

### A Calculator

membrane-toolkit's `core` module packages complex equations from the
literature into rigorous, validated Python routines that make it fast and easy to perform state-of-the-science
analysis. Each core function has detailed documentation, including formulas, assumptions, and references, and
is thoroughly tested for correctness using published data.

Compared to manual calculations, using core functions to perform calculations saves time, improves quality,
increases reproducibility, and reduces the chances of error. In addition, every core function is available with
support for automatic unit conversion, eliminating another potential source of error.

Core functions can be used for standalone calculations, as part of a data analysis pipline, or as building blocks
for sophisticated transport models.

### An automatic Data Aggregator

membrane-toolkit's `pipeline` module makes it easy to automatically aggregate and process experimental
data. It leverages the [maggma](https://materialsproject.github.io/maggma/) framework to collect parse data from
Excel spreadsheets into versatible objects called Stores that can interface with numerous database providers,
including MongoDB, Pandas dataframes, Amazon storage, and more.

A typical data aggregation pipeline consists of the following steps:

1. Choose or customize a Drone, which defines the file format to use when collecting experimental data.
2. Collect data using the template, one file per experiment.
3. Run the Drone to aggregate all the data into a Store
4. View, summarize, and analyze the data in the Store
5. Update as you generate more data. membrane-toolkit's Drones are smart enough to know when files have
   been updated and to re-parse accordingly.

---

membrane-toolkit is written in [Python](http://docs.python-guide.org/en/latest/) and supports Python 3.+.

Please note that membrane-toolkit is released with a Contributor [Code of Conduct](https://github.com/rkingsbury/membrane-toolkit/blob/master/CODE_OF_CONDUCT.MD). By participating in this project you 
agree to abide by its terms.