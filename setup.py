from setuptools import setup, find_packages

setup(
    name="membrane_toolkit",
    version="0.1.0",
    description="Python tools for membrane science",
    url="https://github.com/rkingsbury/membrane-toolkit",
    author="Ryan Kingsbury",
    author_email="RyanSKingsbury@alumni.unc.edu",
    license="BSD 3-clause license",
    packages=find_packages(),
    install_requires=["numpy>=1.14.3", "monty>=3.0.2", "maggma>=0.19", "pint>=0.11", "openpyxl", "pandas"],
    include_package_data=True,
    classifiers=[
        "Programming Language :: Python :: 3",
        "Development Status :: 3 - Alpha",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: BSD License",
        "Operating System :: OS Independent",
        "Topic :: Software Development :: Libraries :: Python Modules",
    ],
)
