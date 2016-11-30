import ez_setup
ez_setup.use_setuptools()

import setuptools
import os
import sys


# Grab the version from the application itself for consistency
import libgwas

# Use the README as the long description
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


setuptools.setup(name="libgwas",
    version=libgwas.__version__,
    author="Eric Torstenson",
    author_email="eric.s.torstenson@vanderbilt.edu",
    url="https://github.com/edwards-lab/libGWAS",
    download_url="https://github.com/edwards-lab/libGWAS/archive/v1.0.0.tar.gz",
    packages=["libgwas", "tests"],
    license="GPL",
    description=["GWAS Parser Library"],
    install_requires=["scipy","numpy"],
    long_description=read('README'),
    keywords=["GWAS","genetic analysis"],
    test_suite='tests',
    package_data={'tests/bedfiles/':['*'],
                  'doc':['*']},
    classifiers=[
        "Development Status :: 4 - Beta",
        "Topic :: Utilities",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Software Development :: Libraries",
        "Programming Language :: Python :: 2.7"
    ],
)
