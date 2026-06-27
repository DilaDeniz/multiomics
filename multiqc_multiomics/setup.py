"""
Setup file for the multiqc_multiomics MultiQC plugin.

Install with:
    pip install .
or in development mode:
    pip install -e .

MultiQC discovers plugins via the ``multiqc.modules`` entry point group.
"""

from setuptools import setup, find_packages

setup(
    name="multiqc_multiomics",
    version="0.1.0",
    author="Multiomics Contributors",
    author_email="",
    description="MultiQC plugin for Multiomics multi-omics reports",
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
            "multiomics = multiqc_multiomics:MultiqcModule",
        ],
        # Register the search pattern so MultiQC can find our JSON files
        "multiqc.cli_options": [
            "multiomics = multiqc_multiomics:cli_options",
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
