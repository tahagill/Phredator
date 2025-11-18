"""
Phredator: Intelligent Next-Generation Sequencing Quality Control Tool
"""
from setuptools import setup, find_packages
from pathlib import Path

# Read the long description from README
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text(encoding="utf-8")

setup(
    name="phredator",
    version="1.0.2",
    author="Taha Ahmad",
    author_email="tahagill99@gmail.com",
    description="Intelligent NGS Quality Control with automated fix suggestions",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/tahagill/Phredator",
    project_urls={
        "Bug Tracker": "https://github.com/tahagill/Phredator/issues",
        "Documentation": "https://github.com/tahagill/Phredator#readme",
        "Source Code": "https://github.com/tahagill/Phredator",
    },
    packages=find_packages(exclude=["tests*", "test_*"]),
    include_package_data=True,
    package_data={
        'phredator': ['config/*.yaml'],
    },
    install_requires=[
        "numpy>=1.20.0",
        "PyYAML>=5.4.0",
    ],
    python_requires=">=3.8",
    keywords=[
        'bioinformatics', 'ngs', 'quality-control', 'fastqc', 'sequencing', 
        'genomics', 'rna-seq', 'chip-seq', 'whole-genome-sequencing',
        'next-generation-sequencing', 'qc', 'phred', 'multiqc',
        'data-analysis', 'computational-biology', 'pipeline'
    ],
    entry_points={
        'console_scripts': [
            'phredator=phredator.cli.cli:main',
        ],
    },
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Operating System :: OS Independent",
    ],
)