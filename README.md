# Molecular Oncology Almanac

Molecular Oncology Almanac is a clinical interpretation algorithm for cancer genomics to annotate and evaluate whole-exome and transcriptome genomic data from individual patient samples. Specifically, Molecular Oncology Almanac can:
- Identify mutations and genomic features related to therapeutic sensitivity and resistance and of prognostic relevance.
- Annotate somatic and germline variants based on their presence in [several datasources](https://github.com/vanallenlab/moalmanac/tree/main/moalmanac/datasources).
- Sort and evaluate somatic mutations from single nucleotide variants, insertions and deletions, copy number alterations, and fusions based on clinical and biological relevance. 
- Integrate data types to observe which genes have been altered in both the somatic and germline setting.
- Extract and evaluate germline mutations relevant to adult and hereditary cancers.
- Identify overlap between somatic variants observed from both DNA and RNA, or any other source of validation sequencing.
- Identify somatic and germline variants that may be related to microsatellite stability.
- Calculate coding mutational burden and compare your patient to TCGA.
- Calculate contribution of known [COSMIC mutational signatures](https://cancer.sanger.ac.uk/signatures/signatures_v2/) with [deconstructsigs](https://github.com/raerose01/deconstructSigs).
- Identify genomic features that may be related to one another.
- Create portable web-based actionability reports, summarizing clinically relevant findings. 

You can view additional documentation, including [descriptions of inputs](docs/description-of-inputs.md) and [outputs](docs/description-of-outputs.md), within the [docs](docs/) folder of this repository.

## Getting Molecular Oncology Almanac
The codebase is available for download through this Github repository, [Dockerhub](https://hub.docker.com/r/vanallenlab/moalmanac/), and [Terra](https://portal.firecloud.org/#methods/vanallenlab/moalmanac/2). The method can also be run on Terra, without having to use Terra, by using [our portal](https://portal.moalmanac.org/). **Accessing Molecular Oncology Almanac through GitHub will require building some of the [datasources](moalmanac/datasources/) but they are also contained in the Docker container**.

### Installation
Molecular Oncology Almanac is a Python application using Python 3.6 but also utilizes R to run [deconstructSigs](https://github.com/raerose01/deconstructSigs) as a subprocess. This application, datasources, and all dependencies are packaged on Docker and can be downloaded with the command
 ```bash
docker pull vanallenlab/moalmanac
```

Alternatively, the package can be built from this Github repository. To download via Github,
```bash
git clone https://github.com/vanallenlab/moalmanac.git
```

We recommend using a [virtual environment](https://docs.python.org/3/tutorial/venv.html) and running Python with either [Anaconda](https://www.anaconda.com/download/) or  [Miniconda](https://conda.io/miniconda.html). After installing Anaconda or Miniconda, you can set up by running
```bash
conda create -n moalmanac python=3.6 -y
source activate moalmanac
pip install -r requirements.txt
```

You can install [deconstructSigs](https://github.com/raerose01/deconstructSigs) after [installing R](https://www.r-project.org/) with the following commands
```bash
Rscript -e 'install.packages("RCurl", repos = "http://cran.rstudio.com/")' \
    && Rscript -e 'source("http://bioconductor.org/biocLite.R"); biocLite("BSgenome"); biocLite("BSgenome.Hsapiens.UCSC.hg19"); biocLite("GenomeInfoDb")' \
    && Rscript -e 'install.packages("reshape2", repos = "http://cran.rstudio.com/")' \
    && Rscript -e 'install.packages("deconstructSigs", repos = "http://cran.rstudio.com/")'
```

## Usage
Molecular Oncology Almanac will run based on any combination of input data but does require a patient_id to label outputs. Additional settings can be set by modifying the [config.ini](moalmanac/config.ini) file and column names may be modified by editing the [colnames.ini](moalmanac/colnames.ini) file.

Required arguments:
```
    --patient_id            <string>    patient identifier
```

Optional arguments:
```
    --tumor_type            <string>    tumor ontology, default=Unknown
    --stage                 <string>    tumor stage, default=Unknown
    --snv_handle            <string>    handle for MAF file of somatic single nucleotide variants
    --indel_handle          <string>    handle for MAF file of somatic insertions and deletions
    --bases_covered_handle  <string>    handle for text file which contains the number of calcable somatic bases
    --called_cn_handle      <string>    handle for text file which contained genes and copy number calls, will be used over `--cnv_handle`
    --cnv_handle            <string>    handle for annotated seg file for somatic copy number
    --fusion_handle         <string>    handle for STAR fusion output, .final.abridged
    --germline_handle       <string>    handle for MAF file of germline single nucleotide variants and insertions and deletions
    --validation_handle     <string>    handle for MAF file of somatic single nucleotide variant called from validation sequencing
    --ms_status             <string>    microsatellite status as deemed by MSI sensor, MSI or MSS, default=Unknown
    --purity                <float>     tumor purity
    --ploidy                <float>     tumor ploidy
    --wgd                   <boolean>   specify the occurence of whole genome duplication
    --disable_matchmaking   <boolean>   remove patient-to-cell line matchmaking from report
    --description           <string>    description of patient
```

## Example
A file `run_example.py` was packaged with this application to run Molecular Oncology Almanac on data found in the folder `example_data`. The outputs produced are the same as those hosted in the folder `example_output`. From the application's folder, run
```bash
python run_example.py
```

## How to contribute
Please follow [our contribution instructions](docs/how-to-contribute.md) if you are interested in contributing to this project. 

## Citation
If you find this tool or any code herein useful, please cite:  
> [Reardon, B., Moore, N.D., Moore, N.S., *et al*. Integrating molecular profiles into clinical frameworks through the Molecular Oncology Almanac to prospectively guide precision oncology. *Nat Cancer* (2021). https://doi.org/10.1038/s43018-021-00243-3](https://www.nature.com/articles/s43018-021-00243-3)

## Disclaimer - For research use only
DIAGNOSTIC AND CLINICAL USE PROHIBITED. DANA-FARBER CANCER INSTITUTE (DFCI) and THE BROAD INSTITUTE (Broad) MAKE NO REPRESENTATIONS OR WARRANTIES OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING, WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NONINFRINGEMENT OR VALIDITY OF ANY INTELLECTUAL PROPERTY RIGHTS OR CLAIMS, WHETHER ISSUED OR PENDING, AND THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER OR NOT DISCOVERABLE.

In no event shall DFCI or Broad or their Trustees, Directors, Officers, Employees, Students, Affiliates, Core Faculty, Associate Faculty and Contractors, be liable for incidental, punitive, consequential or special damages, including economic damages or injury to persons or property or lost profits, regardless of whether the party was advised, had other reason to know or in fact knew of the possibility of the foregoing, regardless of fault, and regardless of legal theory or basis. You may not download or use any portion of this program for any non-research use not expressly authorized by DFCI or Broad. You further agree that the program shall not be used as the basis of a commercial product and that the program shall not be rewritten or otherwise adapted to circumvent the need for obtaining permission for use of the program other than as specified herein.
