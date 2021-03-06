# ESI for "Assessing molecular simulation for the analysis of lipid monolayer reflectometry"

[![J. Phys. Comm.](https://img.shields.io/badge/manuscript%20DOI-10.1088%2F2399--6528%2Fab12a9-yellow.svg)](https://doi.org/10.1088/2399-6528/ab12a9)
[![DOI](https://zenodo.org/badge/156537189.svg)](https://zenodo.org/badge/latestdoi/156537189) [![arXiv](https://img.shields.io/badge/arXiv-1901.05514-orange.svg)](https://arxiv.org/abs/1901.05514) [![License](https://img.shields.io/github/license/arm61/sim_vs_trad.svg?color=lightgrey)](https://github.com/arm61/sim_vs_trad/blob/master/LICENSE)

[Andrew R. McCluskey](https://orcid.org/0000-0003-3381-5911)&ast;, [James Grant](https://orcid.org/0000-0003-1362-2055), [Andrew J. Smith](https://orcid.org/0000-0003-3745-7082), [Jonathan L. Rawle](https://orcid.org/0000-0001-8767-4443), [David J. Barlow](https://orcid.org/0000-0002-0094-5122), [M. Jayne Lawrence](https://orcid.org/0000-0003-4738-4841), [Stephen C. Parker](https://orcid.org/0000-0003-3804-0975), [Karen J. Edler](https://orcid.org/0000-0001-5822-0127)&ast;.

&ast;[a.r.mccluskey@bath.ac.uk](mailto:a.r.mccluskey@bath.ac.uk)/[andrew.mccluskey@diamond.ac.uk](mailto:andrew.mccluskey@diamond.ac.uk) & [k.edler@bath.ac.uk](mailto:k.edler@bath.ac.uk)

*Comparing different molecular dynamics potential models to improve reflectometry analysis.*

This is the electronic supplementary information (ESI) associated with the publication "Assessing molecular simulation for the analysis of lipid monolayer reflectometry".
This ESI provides full details of the analyses performed in the work and access to an automated analysis workflow, through this we aim to provide better access to analysis reproduciblility.
The [Supplementary Information document](reports/si.pdf) can be found in the [reports](/reports) folder, alongside a preprint copy of the publication.
For more information about reproducible workflows in Python, check out [Tania Allard's talk from Pycon2018](http://bitsandchips.me/Talks/PyCon.html#/title).

## [Data](https://researchdata.bath.ac.uk/id/eprint/586)

The reduced neutron reflectometry data and simulation inputs and trajectories can be obtained from the University of Bath Research Data Archive.

DOI: [10.15125/BATH-00586](https://doi.org/10.15125/BATH-00586)

## Analysis

This ESI aims to provide a fully reproducible workflow to the data analysis presented within the paper.

Requirements:

- anaconda or miniconda python
- [REVTeX](https://journals.aps.org/revtex)

The supplied Snakemake file, will reproduce all of the analysis, plot the figures, and build a preprint version of the paper (`reports/preprint.pdf`) when run. Be aware that the analyses within this work are non-trivial and take many hours to run so **use caution** before re-running.

If you **still** want to re-run all of the analysis, please download the [experimental data zip file](https://doi.org/10.15125/BATH-00586), and unzip it (in the `sim_vs_trad` directory) using the following command:

```
unzip sim_vs_trad_data.zip
```

Then run the following commands:


```
conda env create -f config/environment.yml

source activate sim_vs_trad

snakemake clean # this will remove all of the output from previous runs

snakemake
```

## [Figures](/reports/figures)

PDF versions of the figures, can be found in the `reports/figures` directory.

## Citations

**Paper**: A. R. McCluskey, J. Grant, A. J. Smith, J. L. Rawle, D. J. Barlow, M. J. Lawrence, S. C. Parker and K. J. Edler, *J. Phys. Comm.*, 2019, **3**(7), 075001.

**ESI**: A. R. McCluskey, J. Grant, A. J. Smith, J. L. Rawle, D. J. Barlow, M. J. Lawrence, S. C. Parker and K. J. Edler, sim_vs_trad (Version 1.0), http://doi.org/10.5281/zenodo.3254719.

**Data**: A. R. McCluskey, J. Grant, A. J. Smith, J. L. Rawle, D. J. Barlow, M. J. Lawrence, S. C. Parker and K. J. Edler, Dataset for ‘Assessing molecular simulation for the analysis of lipid monolayer reflectometry’, https://doi.org/10.15125/BATH-00586.

### BibTeX

```
@article{sim_vs_trad_paper, 
    title={Assessing molecular simulation for the analysis of lipid monolayer reflectometry}, 
    volume={3}, 
    DOI={10.1088/2399-6528/ab12a9}, 
    number={7}, 
    journal={Phys. Chem. Chem. Phys.}, 
    author={McCluskey, A. R. and Grant, J. and Smith, A. J. and Rawle, J. L. and Barlow, D. J. and Lawrence, M. J. and Parker, S. C. and Edler, K. J.}, 
    year={2019}, 
    pages={075001} 
}

@misc{sim_vs_trad_esi, 
    title={sim_vs_trad (Version 1.0)}, 
    url={http://doi.org/10.5281/zenodo.3254719}, 
    author={McCluskey, A. R. and Grant, J. and Smith, A. J. and Rawle, J. L. and Barlow, D. J. and Lawrence, M. J. and Parker, S. C. and Edler, K. J.}, 
    year={2019} 
}

@misc{sim_vs_trad_data, 
    title={Dataset for “Assessing molecular simulation for the analysis of lipid monolayer reflectometry”}, 
    url={https://doi.org/10.15125/BATH-00586}, 
    author={McCluskey, A. R. and Grant, J. and Smith, A. J. and Rawle, J. L. and Barlow, D. J. and Lawrence, M. J. and Parker, S. C. and Edler, K. J.}, 
    year={2019} 
}
```

## Acknowledgements

A. R. M. is grateful to the University of Bath and Diamond Light Source for co-funding a studentship (Studentship Number STU0149).

## Project Organization

    .
    ├── AUTHORS.md
    ├── LICENSE         # CC BY-SA-4.0
    ├── README.md       # You are here
    ├── Snakefile       # Makefile to outline workflow
    ├── output          # Files and data output by analysis scripts
    │   ├── simulation
    │   └── traditional
    ├── config          # requirements.txt file
    ├── reports         # Paper and ESI
    │   └── figures
    └── scripts       # Python scripts for analysis
        ├── simulation
        └── traditional
