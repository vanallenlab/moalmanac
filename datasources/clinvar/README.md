# Databases: ClinVar
Molecular Oncology Almanac utilizes [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/) to annotate the clinical significance of germline variants. Germline alterations which are annotated as benign will not be reported in the actionable output from Molecular Oncology Almanac, but will still be reported in all other outputs. 

## About ClinVar
ClinVar is a freely accessible database which reports relationships between human genetic variation and phenotypes. It is hosted by the National Center for Biotechnology Information (NCBI) and funded by the National Institutes of Health (NIH). ClinVar has recently partnered with the [Clinical Genome Resource](https://www.clinicalgenome.org/).

## Usage: Downloading ClinVar for use
Releases of ClinVar [are available for download from their FTP server](ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar) and, as of this writing, Molecular Oncology Almanac utilizes [variant_summary.txt.gz](ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/), accessed on October 25th, 2018. 

The script `prepare_clinvar.py` is used to extract relevant columns for use
```bash
python prepare_clinvar.py --input variant_summary.txt --date 2023-03-09
```

## References
1. [Landrum MJ, Lee JM, Riley GR, et al. ClinVar: public archive of relationships among sequence variation and human phenotype. Nucleic Acids Res. 2014;42(Database issue):D980-5.](https://academic.oup.com/nar/article/42/D1/D980/1051029)
