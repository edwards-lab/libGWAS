from setuptools import setup


version = {}
with open("libgwas/version.py") as fp:
    exec(fp.read(), version)


import os


# Use the README as the long description
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


setup(
    name="libgwas",
    version=version["__version__"],
    author="Eric Torstenson",
    author_email="eric.s.torstenson@vanderbilt.edu",
    url="https://github.com/edwards-lab/libGWAS",
    download_url="https://github.com/edwards-lab/libGWAS/archive/v1.1.0.tar.gz",
    packages=["libgwas", "libgwas.tests"],
    license="GPL",
    description=["GWAS Parser Library"],
    install_requires=["scipy", "numpy", "pytabix", "bgen_reader"],
    long_description=read("README"),
    keywords=["GWAS", "genetic analysis"],
    test_suite="libgwas.tests",
    package_data={"libgwas/tests/bedfiles/": ["*"], "doc": ["*"]},
    classifiers=[
        "Development Status :: 4 - Beta",
        "Topic :: Utilities",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Software Development :: Libraries",
        "Programming Language :: Python :: 2.7",
    ],
)
