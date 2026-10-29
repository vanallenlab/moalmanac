[![DOI](https://img.shields.io/badge/DOI-10.1038%2Fs43018--021--00243--3-blue)](https://doi.org/10.1038/s43018-021-00243-3)
[![License: CC BY 4.0](https://img.shields.io/badge/Database-CC%20BY%204.0-informational)](./LICENSE)
[![License: GPL v2](https://img.shields.io/badge/Code-GPL%20v2-informational)](./LICENSE-code)

Molecular Oncology Almanac is a clinical interpretation algorithm for cancer genomics to annotate and evaluate whole-exome and transcriptome genomic data from individual patient samples. Specifically, Molecular Oncology Almanac can:

- Identify mutations and genomic features related to therapeutic sensitivity and resistance and of prognostic relevance.
- Annotate somatic and germline variants based on their presence in [several datasources](datasources/).
- Sort and evaluate somatic mutations from single nucleotide variants, insertions and deletions, copy number alterations, and fusions based on clinical and biological relevance.
- Integrate data types to observe which genes have been altered in both the somatic and germline setting.
- Extract and evaluate germline mutations relevant to adult and hereditary cancers.
- Identify overlap between somatic variants observed from both DNA and RNA, or any other source of validation sequencing.
- Identify somatic and germline variants that may be related to microsatellite stability.
- Calculate coding mutational burden and compare your patient to TCGA.
- Identify genomic features that may be related to one another.
- Create portable web-based actionability reports, summarizing clinically relevant findings.

It is supported by our [underlying manually curated precision oncology knowledge base](https://github.com/vanallenlab/moalmanac).

You can view additional documentation, including [descriptions of inputs](docs/description-of-inputs.md) and [outputs](docs/description-of-outputs.md), within the [docs](docs/) folder of this repository.

## Getting Molecular Oncology Almanac

The codebase is available for download through this GitHub repository, [DockerHub](https://hub.docker.com/r/vanallenlab/moalmanac/), and [Terra](https://app.terra.bio/#workflows/vanallenlab/moalmanac). **Accessing Molecular Oncology Almanac through GitHub will require building some of the [datasources](datasources/) but they are also contained in the Docker container**.

## Installation

### Download

This repository along with datasources and dependencies are packaged with Docker and can be downloaded with the command:

 ```bash
docker pull vanallenlab/moalmanac:latest
```

Alternatively, this can be downloaded through GitHub by either using the website or terminal. To download on the website, navigate to the top of this page, click the green `Clone or download` button, and select `Download ZIP` to download this repository in a compressed format. To install using GitHub on terminal, type:

```bash
git clone https://github.com/vanallenlab/moalmanac.git
cd moalmanac
```

### Python dependencies

This repository supports Python 3.14. We recommend managing your environment with [Conda-forge](https://conda-forge.org)'s installer [Miniforge](https://conda-forge.org/download/).

Run the following from this repository's directory to create a virtual environment and install dependencies with Miniforge:

```bash
conda env create -f environment.yaml
conda activate moalmanac
```

The [environment.yaml](./environment.yaml) file specifies Python 3.14 and delegates package installation to pip using [pyproject.toml](./pyproject.toml). No separate `pip install` step is needed.

Or, if you prefer not to use conda, you can use a standard [virtual environment](https://docs.python.org/3/tutorial/venv.html):

```bash
python3.14 -m venv venv
source venv/bin/activate
```

If using VS Code, open the repository folder and select the interpreter via **Command Palette → Python: Select Interpreter**, then choose `moalmanac`, or what you named your virtual environment, from the list.

Then install dependencies using one of the following two options:

#### Reproducible install (recommended)

Install from the [requirements-lock.txt](./requirements-lock.txt) file to get the exact same package versions used in Docker and the verified test environment:

```bash
pip install -r requirements-lock.txt
```

#### Standard install

Install directly from [pyproject.toml](./pyproject.toml) with the latest compatible versions:

```bash
pip install .
```

A new requirements lock file can be generated with:

```bash
pip-compile --output-file=requirements-lock.txt pyproject.toml
```

## Usage

Usage documentation can be found within the [moalmanac/](moalmanac) directory of this repository.

## Contributing

Please follow [our contribution instructions](./CONTRIBUTING.md) if you are interested in contributing to this project.

## Citation

The Molecular Oncology Almanac database is available under CC BY 4.0, which requires appropriate credit. If you use the database, or find this tool or any code herein useful, please cite:

> [Reardon, B., Moore, N.D., Moore, N.S., *et al*. Integrating molecular profiles into clinical frameworks through the Molecular Oncology Almanac to prospectively guide precision oncology. *Nat Cancer* (2021). https://doi.org/10.1038/s43018-021-00243-3](https://www.nature.com/articles/s43018-021-00243-3)

Where relevant, please also cite the version of the database used with its Zenodo DOI.

- The concept DOI always points to the latest version: [https://doi.org/10.5281/zenodo.17089309](https://doi.org/10.5281/zenodo.17089309).
- Each release also has a version-specific DOI. We began minting Zenodo DOIs with our 2025-09-04 release; if you would like a DOI for a prior version, please [contact us](https://moalmanac.org/about#contact).

## License

- Database content: [Creative Commons Attribution 4.0 International (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/), see [LICENSE](./LICENSE).
- Code: [GNU General Public License, Version 2 (GPL v2)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html), see [LICENSE-code.md](./LICENSE-code).
- Third-party files in [datasources/](./datasources) are not covered by the CC BY license; see [datasources/README.md](./datasources/README.md).

## Disclaimer

For purposes of this Disclaimer, the “Service” means the Molecular Oncology Almanac database, software, API, browsers, and related materials.

THIS DATABASE AND SOFTWARE HAVE NOT BEEN VALIDATED FOR DIAGNOSTIC OR CLINICAL USE, AND DANA-FARBER CANCER INSTITUTE (DFCI) MAKES NO WARRANTY OF FITNESS FOR DIAGNOSTIC OR CLINICAL USE. USERS ARE SOLELY RESPONSIBLE FOR ENSURING THAT THEIR ACCESS TO AND USE OF THE DATABASE, SOFTWARE, API, BROWSERS, AND RELATED MATERIALS COMPLIES WITH ALL APPLICABLE LAWS, REGULATIONS, INSTITUTIONAL POLICIES, AND ETHICAL REQUIREMENTS, INCLUDING ANY REQUIREMENTS APPLICABLE TO RESEARCH INVOLVING HUMAN SUBJECTS OR HUMAN-DERIVED DATA. DFCI MAKES NO REPRESENTATIONS OR WARRANTIES OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING, WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NONINFRINGEMENT OR VALIDITY OF ANY INTELLECTUAL PROPERTY RIGHTS OR CLAIMS, WHETHER ISSUED OR PENDING, AND THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER OR NOT DISCOVERABLE.

In no event shall DFCI or its trustees, directors, officers, employees, students, affiliates, core faculty, associate faculty and contractors be liable for incidental, punitive, consequential or special damages, including economic damages or injury to persons or property or lost profits, regardless of whether the party was advised, had other reason to know or in fact knew of the possibility of the foregoing, regardless of fault, and regardless of legal theory or basis. To the extent liability cannot be disclaimed, DFCI’s liability shall be limited to direct, actual damages resulting from DFCI’s failure to perform its obligations under the applicable terms.

The Molecular Oncology Almanac database is made available under the [Creative Commons Attribution 4.0 International License (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/), which permits use, including commercial use, provided appropriate credit is given. The software is licensed under the [GNU General Public License, Version 2 (GPL v2)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html). Certain third-party code or materials may be subject to separate license terms. For additional terms governing access to and use of the Service, please see the [Terms of Use](https://docs.moalmanac.org/about/terms-of-use/).
