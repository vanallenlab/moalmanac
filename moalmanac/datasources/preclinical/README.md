# Directly leveraging cancer cell lines for clinical interpretation
Data processing of cancer cell lines for usage by Molecular Oncology Almanac is discussed in the [MOAlmanac paper Github repository](https://github.com/vanallenlab/moalmanac-paper). This repository contains the files directly utilizes by MOAlmanac to test for preclinical efficacy of relationships and to perform patient model matchmaking.

So that users do not have to follow the procedure of downloading, annotating, and evaluating raw data, you can download the utilized cell line somatic variants, copy number alterations, and fusions to this directory with `download-files.sh`. This should be run with the repository's virtual environment active as it utilizes the package [gdown](https://github.com/wkentaro/gdown). The files are hosted on a public Google Drive folder through the Broad Institute and are,
- `ccle.copy-numbers.evaluated.txt` (100.5 MB, MD5: a57eee56310a7c96b4a2cfa53aa6a2de), copy number alterations annotated and evaluated with MOAlmanac
- `ccle.variants.evaluated.txt` (53.8 MB, MD5: 970123c95fdbcd3138eba8b6b71047ce), somatic variants annotated and evaluated with MOAlmanac
- `sanger.fusions.evaluated.txt` (527 KB, MD5: df830391473ea1e73fddd7e449216434), fusions with the strongest match per feature from the `gene1` and `gene2` files
- `sanger.fusions.gene1.evaluated.txt` (483 KB, MD5: 59398e768f48a369724ba619ee6dc177), fusions annotated and evaluated with MOAlmanac relative to gene 1
- `sanger.fusions.gene2.evaluated.txt` (477 KB, MD5: 357e374a78c4200fa8fcca66979a3ebb), fusions annotated and evaluated with MOAlmanac relative to gene 2


For more details, please refer to the [paper repository](https://github.com/vanallenlab/moalmanac-paper) and/or [the present study](https://www.nature.com/articles/s43018-021-00243-3).

## Usage
```bash
bash download-files.sh
```
