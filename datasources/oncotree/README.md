# Databases: OncoTree
The Molecular Oncology Almanac utilizes OncoTree to standardize the input of disease ontology in order to frame a given patient's mutational burden and match assertions from the almanac as closely as possible. 

## About Oncotree
OncoTree was developed by Memorial Sloan Kettering Cancer Center as a tool for the cBioPortal and distributed under a Creative Commons Attribution 4.0 International License. Please visit their [website](https://oncotree.mskcc.org/), [API](https://oncotree.mskcc.org/swagger-ui.html), or [GitHub](https://github.com/cBioPortal/oncotree) to learn more.

## Usage: Generating OncoTree for use
The following command can be run to generate a new OncoTree table for usage with the Molecular Oncology Almanac.

```bash
python get_oncotree.py -d 2023-03-09
```

The output will be of the following name, `oncotree.2023-03-09.txt`

## References
[Kundra, R. et al. OncoTree: A Cancer Classification System for Precision Oncology. JCO Clin Cancer Inform 5, 221-230 (2021).](https://ascopubs.org/doi/10.1200/CCI.20.00108)
