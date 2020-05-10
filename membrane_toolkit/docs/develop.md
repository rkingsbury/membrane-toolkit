# Contributing to membrane-toolkit

We welcome contributions of all kinds, from reporting bugs or requesting features via a
[Github issue](https://github.com/rkingsbury/membrane-toolkit/issues) to contributing
documentation and examples to hacking the code and submitting a pull request. Please note
that membrane-toolkit is released with a Contributor
[Code of Conduct](https://github.com/rkingsbury/membrane-toolkit/blob/master/CODE_OF_CONDUCT.MD).
By participating in this project you agree to abide by its terms.




### Code Organization

Core - fully documented, thoroughly tested, standalone functions for key calculations,
    implemented with floating point math and available in unitized versions with Quantity support






### Coding Concepts

membrane-toolkit is designed to work closely with the family of codes maintained by
[The Materials Project](https://github.com/materialsproject/). In particular, 
we follow the major architectural decisions of the 
[pymatgen](https://github.com/materialsproject/pymatgen) maintainers regarding
coding style, test framework, Python compatbility, etc. 


### Key third-party libraries

Wherever possible, we 
import functionality from pymatgen or other codes rather than duplicate 
functionality. For example, we use pymatgen's `Composition` class to represent 
chemical information about different materials. If you are new to this family 
of codes, it will be especially helpful to familiarize yourself with the 
following methods and packages, in addition to standard scientific Python codes
like numpy and scipy:

 - `MSONable`, `loadfn` and `dumpfn` serialization routines from the
   [monty](https://github.com/materialsvirtuallab/monty) package
 - Units-aware computation using [pint](https://pint.readthedocs.io/en/0.9/)
   Quantity objects.
 - `Element`, `Specie`, and `Composition` classes from 
   [pymatgen](https://github.com/materialsproject/pymatgen)
 - String operations in [pymatgen.util.string](https://github.com/materialsproject/pymatgen/blob/2a813c172f3be38efc3205a102021eaba1da156f/pymatgen/util/string.py)
 - Type Hinting functionality in the 
   [typing](https://docs.python.org/3/library/typing.html) module
 - `Builder` class from [maggma](https://github.com/materialsproject/maggma)


## Coding Guidelines

To ensure the long-term sustainability of the membrane-toolkit codebase, we enforce very
strict quality control standards for all contributions. Specifically, all new
code must adhere to the following:

1. **Tests** are required for all new modules and methods. The only way to
   minimize regression and ensure confidence in membrane-toolkit's calculations is to 
   thoroughly test all code. If the maintainer cannot test your code, 
   the contribution will be rejected. Test should conform to the following schema:
        
        * Use [pytest](https://docs.pytest.org/en/latest/)
        * Place tests in files named `test_xxxx.py` in a `tests` subdirectory 
          of the directory where your code resides.
        * If your test requires external files, place those in `test_files` in
          the root membrane-toolkit directory.
        * Test all realistic combinations of input arguments 

2. **Documentation** required for all modules, classes and methods. In
   particular:
   
        * The method docstrings should explain ALL required and optional
          arguments and return values
        * Methods and docstrings should use 
          [TypeHinting](https://docs.python.org/3/library/typing.html)
        * Use [Google-style](https://google.github.io/styleguide/pyguide.html)
          docstrings
        * Include a `References:` section with citations where appropriate

3. **PyCodeStyle**, **flake8**, and **mypy** are used to lint (check) all code, and
   all code is expected to adhere to this style. To make linting less tedious, we
   use [black](https://github.com/psf/black) to automatically format all files.
   We highly recommend that you install black and pycodestyle in your IDE to 
   check your code style BEFORE submitting a pull request.

4. **Python 3**. We only support Python 3.7+.


If in doubt about any of the above, please refer to the source code in `membrane_toolkit.core` for
examples of what is expected.
