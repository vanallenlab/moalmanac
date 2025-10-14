# Databases: Catalogue of Somatic Mutations in Cancer
The Molecular Oncology Almanac's heuristic utilizes the [Catalogue of Somatic Mutations In Cancer (COSMIC)](http://cancer.sanger.ac.uk/cosmic) to identify alterations that may be relevant to cancer. Specifically, our heuristic evaluates alterations at the gene level with respect to COSMIC.

## About COSMIC
COSMIC is developed and maintained by the [Wellcome Trust Sanger Institute](http://www.sanger.ac.uk/) for exploring the impact of somatic mutations in human cancer.  

COSMIC's is distributed under their license and terms, [described here](https://www.cosmickb.org/terms/). Users must register and download COSMIC's datasource themselves for use with the Molecular Oncology Almanac. 

## Usage: Downloading and formatting COSMIC
Data can be downloaded from [COSMIC's one click data downloads portal](http://cancer.sanger.ac.uk/cosmic/download). The Molecular Oncology Almanac leverages the `COSMIC Mutation Data` file, labeled as `CosmicMutantExport.tsv.gz`. Please download the whole file (4 GB) and uncompress the file.

The script `prepare_cosmic.py` is used to extract Gene names and Protein Changes for use. For example,
```bash
python prepare_cosmic.py --input CosmicMutantExport.tsv --version "v97" --gene_column_name "Gene name" --protein_column_name "Mutation AA"
```

This script will produce an output containing unique gene and protein change pairs for annotating with the Molecular Oncology Almanac. The produced output will be named `CosmicMutantExport_{version}.lite.txt`. 

## References
1. [Forbes SA, Beare D, Boutselakis H, et al. COSMIC: somatic cancer genetics at high-resolution. Nucleic Acids Res. 2017;45(D1):D777-D783.](https://academic.oup.com/nar/article/45/D1/D777/2605743)
