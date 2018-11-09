# ESI for "Applying molecular simulation to the analysis of lipid monolayer reflectometry"

[![arXiv](https://img.shields.io/badge/arXiv-1810.07616-orange.svg)](https://arxiv.org/abs/1810.07616) [![DOI](https://zenodo.org/badge/144010644.svg)](https://zenodo.org/badge/latestdoi/144010644)

[Andrew R. McCluskey](https://orcid.org/0000-0003-3381-5911)&ast;, [James Grant](https://orcid.org/0000-0003-1362-2055), [Andrew J. Smith](https://orcid.org/0000-0003-3745-7082), [Jonathan L. Rawle](https://orcid.org/0000-0001-8767-4443), [David J. Barlow](https://orcid.org/0000-0002-0094-5122), [M. Jayne Lawrence](https://orcid.org/0000-0003-4738-4841), [Stephen C. Parker](https://orcid.org/0000-0003-3804-0975), [Karen J. Edler](https://orcid.org/0000-0001-5822-0127)&ast;.

&ast;[a.r.mccluskey@bath.ac.uk](mailto:a.r.mccluskey@bath.ac.uk) & [k.edler@bath.ac.uk](mailto:k.edler@bath.ac.uk)

This is the electronic supplementary information (ESI) associated with the publication "Applying molecular simulation to the analysis of lipid monolayer reflectometry".
This ESI provides full details of the analyses performed in the work and access to an automated analysis workflow, through this we aim to provide better access to analysis reproduciblility.
The [Supplementary Information document](reports/si.pdf) can be found in the [reports](/reports) folder, alongside a preprint copy of the publication.
For more information about reproducible workflows in Python, check out [Tania Allard's talk from Pycon2018](http://bitsandchips.me/Talks/PyCon.html#/title).

## [Data]()

## Analysis

This ESI aims to provide a fully reproducible workflow to the data analysis presented within the paper.

Requirements:

- anaconda or miniconda python

The supplied Snakemake file, will reproduce all of the analysis, plot the figures, and build a preprint version of the paper (`reports/paper.pdf`) when run. Be aware that the analyses within this work are non-trivial and take many hours to run so **use caution** before re-running.

If you **still** want to re-run all of the analysis, please download the [experimental data](), place it in a directory named `data` before running the following commands:

```
conda create --name paper_env python

source activate paper_env

pip install --upgrade pip

pip install -r config/requirements.txt

snakemake clean # this will remove all of the output from previous runs

snakemake
```

## [Figures](/reports/figures)

PDF versions of the figures, can be found in the `reports/figures` directory.

## Acknowledgements

A. R. M. is grateful to the University of Bath and Diamond Light Source for co-funding a studentship (Studentship Number STU0149).

## Project Organization

    .
    ├── AUTHORS.md
    ├── LICENSE         # CC-BY-SA-4.0
    ├── README.md       # You are here
    ├── Snakefile       # Makefile to outline workflow
    ├── bin             # Some python scripts
    ├── output          # Files and data output by analysis scripts
    │   ├── simulation
    │   └── traditional
    ├── config          # requirements.txt file
    ├── notebooks       # Notebooks for analysis
    │   ├── simulation
    │   └── traditional
    ├── reports         # Paper and ESI
    │   └── figures
    └── models          # mol_vol.py custom model for refnx
