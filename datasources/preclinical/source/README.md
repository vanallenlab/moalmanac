# Source data for cancer cell lines
Data from the following sources was acquired and _not_ committed to Github. Instead, this directory's README contains instructions for how each set of data was acquired and each subfolder's README contains a manifest of files downloaded. 

## Source: Cancer Cell Line Encyclopedia (CCLE)
Genomic and sample metadata for cell lines characterized by the CCLE was downloaded from [cBioPortal](https://www.cbioportal.org/study/summary?id=ccle_broad_2019). In addition, [supplementary tables](https://www.nature.com/articles/s41586-019-1186-3#Sec60) 1, 3, 6 and [source data](https://www.nature.com/articles/s41586-019-1186-3#Sec61) for figure 5 were also downloaded.

Files were placed into the [ccle-2019/](ccle-2019/) folder in this directory.

## Source: Dependency Map (DepMap)
`sample_info.csv` from DepMap Public 20Q3 was downloaded by navigating to their [Data Downloads](https://depmap.org/portal/download/), selecting `All Downloads` on the left side, and searching for the file name.

This file was placed into the [depmap/](depmap/) folder in this directory.

## Source: Genomics of Drug Sensitivity in Cancer (GDSC)
Fusion, sample metadata, and therapeutic response were downloaded from the GDSC. Fusions were downloaded and unzipped from the Sanger Institute's [Cell Model Passports downloads page](https://cellmodelpassports.sanger.ac.uk/downloads) by selecting `Fusions Data`, `View All Versions`, and `fusions_20191101.zip`. Sample metadata was downloaded from the same webpage by selecting `Model and Gene Annotation`, `View All Versions` under `Model List`, and selecting `model_list_20200204.csv`. Therapeutic screen data was downloaded from [the GDSC download page](https://www.cancerrxgene.org/downloads/bulk_download), specifically the `GDSC1-dataset`, `GDSC2-dataset`, and `IC50 data definitions` from the `Drug Screening - IC50s` section under `Screening Data`. 

Files were placed into the [gdsc/](gdsc/) folder in this directory.
