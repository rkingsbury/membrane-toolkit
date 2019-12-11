Coding Concepts
===============

pymemsci is designed to work closely with the family of codes maintained by
[The Materials Project](https://github.com/materialsproject/). In particular, 
we follow the major architectural decisions of the 
[pymatgen](https://github.com/materialsproject/pymatgen) maintainers regarding
coding style, test framework, Python compatbility, etc. 


Key third-party libraries
=========================
Wherever possible, we 
import functionality from pymatgen or other codes rather than duplicate 
functionality. For example, we use pymatgen's `Composition` class to represent 
chemical information about different materials. If you are new to this family 
of codes, it will be especially helpful to familiarize yourself with the 
following methods and packages, in addition to standard scientific Python codes
like numpy and scipy:

* `MSONable`, `loadfn` and `dumpfn` serialization routines from the
  [monty](https://github.com/materialsvirtuallab/monty) package
* Units-aware computation using [pint](https://pint.readthedocs.io/en/0.9/)
  Quantity objects.
* `Element`, `Specie`, and `Composition` classes from 
  [pymatgen](https://github.com/materialsproject/pymatgen)
* String operations in [pymatgen.util.string](https://github.com/materialsproject/pymatgen/blob/2a813c172f3be38efc3205a102021eaba1da156f/pymatgen/util/string.py)
* Type Hinting functionality in the 
  [typing](https://docs.python.org/3/library/typing.html) module


Coding Guidelines
=================

To ensure the long-term sustainability of the pymemsci codebase, we enforce very
strict quality control standards for all contributions. Specifically, all new
code must adhere to the following:

1. **Unittests** are required for all new modules and methods. The only way to
   minimize regression and ensure confidence in pymemsci's calculations is to 
   thoroughly test all code. If the maintainer cannot test your code, 
   the contribution will be rejected. Test should conform to the following schema:
        
        * Use [pytest](https://docs.pytest.org/en/latest/)
        * Place tests in files named `test_xxxx.py` in a `/tests` subdirectory 
          of the directory where your code resides.
        * If you test requires external files, place those in `test_files` in
          the root pymemsci directory.
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

3. **PyCodeStyle** is used to lint (check) all code, and all code is expected
   to adhere to its style. To make pycodestyle compliance less tedious, we
   use [black](https://github.com/psf/black) to automatically format all files.
   We highly recommend that you install black and pycodestyle in your IDE to 
   check your code style BEFORE submitting a pull request.

4. **Python 3**. We only support Python 3.5+.


For the above, if in doubt, please refer to the core classes in pymemsci for
examples of what is expected.