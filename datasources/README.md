# Datasources used by MOAlmanac
MOAlmanac leverages several datasources to annotate and evaluate molecular features for clinical and biological relevance. **If you are running this software after cloning it from Github, you will need to preprocess some before running this tool** because some files exceed the Github storage limit of 100 Mb. For instructions on how to build the relevant datasources, please refer to their respective folders in this directory.

Some files require download over Google. If you are unable to download through Google or use Docker to access the processed datasources, please contact us and we will work with you to figure something out.

| Name                                                    | Immediately ready for use from Github | Included in Docker |
|---------------------------------------------------------|---------------------------------------|--------------------|
| [COSMIC](cosmic/)                                       | :x:                                   | :x: |
| [ExAC](exac/)                                           | :x:                                   | :white_check_mark: |
| [American College of Medical Genetics v2 (ACMG)](acmg/) | :white_check_mark:                    | :white_check_mark: |
| [Cancer cell line data](preclinical/)                   | :white_check_mark:                    | :white_check_mark: |
| [Cancer Gene Census (CGC)](cancergenecensus/)           | :white_check_mark:                    | :white_check_mark: |
| [Cancer Hotspots](cancerhotspots/)                      | :white_check_mark:                    | :white_check_mark: |
| [ClinVar](clinvar/)                                     | :white_check_mark:                    | :white_check_mark: |
| [Genes associated with hereditary cancers](hereditary/) | :white_check_mark:                    | :white_check_mark: |
| [MSigDB](gsea_gene_sets/)                               | :white_check_mark:                    | :white_check_mark: |
| [Molecular Oncology Almanac](moalmanac/)                | :white_check_mark:                    | :white_check_mark: |
| [Oncotree](oncotree/)                                   | :white_check_mark:                    | :white_check_mark: |
| [TCGA mutational burden](lawrence/)                     | :white_check_mark:                    | :white_check_mark: |
