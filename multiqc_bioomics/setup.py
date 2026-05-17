"""
Setup file for the multiqc_bioomics MultiQC plugin.

Install with:
    pip install .
or in development mode:
    pip install -e .

MultiQC discovers plugins via the ``multiqc.modules`` entry point group.
"""

from setuptools import setup, find_packages

setup(
    name="multiqc_bioomics",
    version="0.1.0",
    author="BioMultiOmics Contributors",
    author_email="",
    description="MultiQC plugin for BioMultiOmics multi-omics reports",
    long_description=open("../README.md").read() if __import__("os").path.exists("../README.md") else "",
    long_description_content_type="text/markdown",
    url="https://github.com/diladeniz/multiomics",
    license="Apache-2.0",
    packages=find_packages(),
    python_requires=">=3.8",
    install_requires=[
        "multiqc>=1.21",
    ],
    entry_points={
        # MultiQC plugin discovery — the key must be ``multiqc.modules``
        "multiqc.modules": [
            "bioomics = multiqc_bioomics:MultiqcModule",
        ],
        # Register the search pattern so MultiQC can find our JSON files
        "multiqc.cli_options": [
            "bioomics = multiqc_bioomics:cli_options",
        ],
    },
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: Apache Software License",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)
