# Molecular Oncology Almanac

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

You can view additional documentation, including [descriptions of inputs](docs/description-of-inputs.md) and [outputs](docs/description-of-outputs.md), within the [docs](docs/) folder of this repository.

## Getting Molecular Oncology Almanac

The codebase is available for download through this GitHub repository, [DockerHub](https://hub.docker.com/r/vanallenlab/moalmanac/), and [Terra](https://portal.firecloud.org/#methods/vanallenlab/moalmanac/2). **Accessing Molecular Oncology Almanac through GitHub will require building some of the [datasources](datasources/) but they are also contained in the Docker container**.

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

Install directly from `pyproject.toml` with the latest compatible versions:

```bash
pip install .
```

A new requirements lock file can be generated with:

```bash
pip-compile --output-file=requirements-lock.txt pyproject.toml
```

## Usage

Usage documentation can be found within the [moalmanac/](moalmanac) directory of this repository.

## How to contribute

Please follow [our contribution instructions](docs/how-to-contribute.md) if you are interested in contributing to this project.

## Citation

If you find this tool or any code herein useful, please cite:
> [Reardon, B., Moore, N.D., Moore, N.S., *et al*. Integrating molecular profiles into clinical frameworks through the Molecular Oncology Almanac to prospectively guide precision oncology. *Nat Cancer* (2021). https://doi.org/10.1038/s43018-021-00243-3](https://www.nature.com/articles/s43018-021-00243-3)

## Disclaimer - For research use only

DIAGNOSTIC AND CLINICAL USE PROHIBITED. DANA-FARBER CANCER INSTITUTE (DFCI) and THE BROAD INSTITUTE (Broad) MAKE NO REPRESENTATIONS OR WARRANTIES OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING, WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NONINFRINGEMENT OR VALIDITY OF ANY INTELLECTUAL PROPERTY RIGHTS OR CLAIMS, WHETHER ISSUED OR PENDING, AND THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER OR NOT DISCOVERABLE.

In no event shall DFCI or Broad or their Trustees, Directors, Officers, Employees, Students, Affiliates, Core Faculty, Associate Faculty and Contractors, be liable for incidental, punitive, consequential or special damages, including economic damages or injury to persons or property or lost profits, regardless of whether the party was advised, had other reason to know or in fact knew of the possibility of the foregoing, regardless of fault, and regardless of legal theory or basis. You may not download or use any portion of this program for any non-research use not expressly authorized by DFCI or Broad. You further agree that the program shall not be used as the basis of a commercial product and that the program shall not be rewritten or otherwise adapted to circumvent the need for obtaining permission for use of the program other than as specified herein.
