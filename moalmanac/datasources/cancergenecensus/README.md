# Databases: Cancer Gene Census
The Molecular Oncology Almanac's heuristic utilizes the [Cancer Gene Census (CGC)](http://cancer.sanger.ac.uk/census) to identify alterations that may be relevant to cancer, even if they are not a hotspot or present in our action-alteration database. Specifically, our heuristic evaluates alterations at the gene level with respect to the CGC, both tier 1 and 2.

## About the Cancer Gene Census
The Cancer Gene Census is developed and maintained by the [Wellcome Trust Sanger Institute](http://www.sanger.ac.uk/) and is a project of [COSMIC](http://cancer.sanger.ac.uk/cosmic) to catalogue genes with mutations that are causally implicated in cancer.

## Usage: Downloading the Cancer Gene Census
Data can be downloaded for the CGC from either [COSMIC's data downloads page](http://cancer.sanger.ac.uk/cosmic/download) or directly exported as a csv or tsv from [the CGC's webpage](http://cancer.sanger.ac.uk/census). Once downloaded, gene names can be extracted with the script `extract_genes.py`. For example,
```bash
python extract_genes.py -input cancer_gene_census.tsv -output cancer_gene_census_v97.genes.tsv -gene_column_name "Gene Symbol"
```

The script `diff_versions.py` may be also be used to compare two versions of the Cancer Gene Census to identify genes added or removed between versions. For example,
```bash
python diff_versions.py -old_version cancer_gene_census_v85.tsv -n cancer_gene_census_v97.genes.tsv 
```

## References
1. [Futreal PA, Coin L, Marshall M, et al. A census of human cancer genes. Nat Rev Cancer. 2004;4(3):177-83.](https://www.ncbi.nlm.nih.gov/pubmed/14993899)
