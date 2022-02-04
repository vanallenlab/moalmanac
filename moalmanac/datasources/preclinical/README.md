# Directly leveraging cancer cell lines for clinical interpretation
Data processing of cancer cell lines for usage by Molecular Oncology Almanac is discussed in the [MOAlmanac paper Github repository](https://github.com/vanallenlab/moalmanac-paper). This repository contains the files directly utilizes by MOAlmanac to test for preclinical efficacy of relationships and to perform patient model matchmaking.

So that users do not have to follow the procedure of downloading, annotating, and evaluating raw data, you can download the utilized cell line somatic variants, copy number alterations, and fusions to this directory with `download-files.sh`. This should be run with the repository's virtual environment active as it utilizes the package [gdown](https://github.com/wkentaro/gdown). The files are hosted on a public Google Drive folder through the Broad Institute,

## Usage
```bash
bash download-files.sh
```
