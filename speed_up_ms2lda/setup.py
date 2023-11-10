#!/usr/bin/env python
import os

from setuptools import setup
from setuptools import find_packages

here = os.path.abspath(os.path.dirname(__file__))

version = {}
with open(os.path.join(here, "__version__.py")) as f:
    exec(f.read(), version)

with open("README.md") as readme_file:
    readme = readme_file.read()

setup(
    name="ms2lda_reloaded",
    version=version["__version__"],
    description="Fast gensim implementation of ms2lda",
    long_description=readme,
    author="Niek de Jonge",
    author_email="niek.dejonge@wur.nl",
    url="https://github.com/niekdejonge/MS2LDA_Reloaded",
    packages=find_packages(),
    include_package_data=True,
    license="Apache Software License 2.0",
    zip_safe=False,
    test_suite="tests",
    python_requires='>=3.7',
    install_requires=[
        "matchms>=0.6.1",
        "spec2vec",
        "gensim>=4.0.0",
        "tqdm",
        "matplotlib",
    ],
    extras_require={"dev": ["pytest",
                            "pytest-cov",
                            "prospector[with_pyroma]"],
    }
)
