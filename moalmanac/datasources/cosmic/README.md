# Databases: Catalogue of Somatic Mutations in Cancer
The Molecular Oncology Almanac's heuristic utilizes the [Catalogue of Somatic Mutations In Cancer (COSMIC)](http://cancer.sanger.ac.uk/cosmic) to identify alterations that may be relevant to cancer. Specifically, our heuristic evaluates alterations at the gene level with respect to COSMIC.

## About COSMIC
COSMIC is developed and maintained by the [Wellcome Trust Sanger Institute](http://www.sanger.ac.uk/) for exploring the impact of somatic mutations in human cancer. At the time of this writing, COSMIC v85 contains nearly 30,000 genes and 6 million protein changes. 

## Usage: Downloading and formatting COSMIC
Data can be downloaded from [COSMIC's one click data downloads portal](http://cancer.sanger.ac.uk/cosmic/download). The Molecular Oncology Almanac leverages the `COSMIC Mutation Data` file, labeled as `CosmicMutantExport.tsv.gz`. Please download the whole file (4 GB) and uncompress the file.

The script `prepare_cosmic.py` is used to extract Genes and Protein Changes for use. Here, we have renamed `CosmicMutantExport.tsv` to be `CosmicMutantExport_v85.tsv` in order to designate the datasource version used in the present study.

```
python prepare_cosmic.py --cosmicMutantExport CosmicMutantExport_v85.tsv
```

This script will produce an output with the suffix `.lite.txt`, which will be used by MOAlmanac.

## References
1. [Forbes SA, Beare D, Boutselakis H, et al. COSMIC: somatic cancer genetics at high-resolution. Nucleic Acids Res. 2017;45(D1):D777-D783.](https://academic.oup.com/nar/article/45/D1/D777/2605743)
