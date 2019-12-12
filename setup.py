from setuptools import setup, find_packages

setup(
    name="pymemsci",
    version="0.1.0",
    description="Python membrane science",
    url="n/a",
    author="Ryan Kingsbury",
    author_email="RyanSKingsbury@alumni.unc.edu",
    license="TBD",
    packages=find_packages(),
    install_requires=["numpy>=1.14.3", "monty>=3.0.2"],
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Development Status :: 1 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Software Development :: Libraries :: Python Modules",
    ],
)
