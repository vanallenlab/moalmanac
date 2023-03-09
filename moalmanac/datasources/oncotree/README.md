# Databases: OncoTree
The Molecular Oncology Almanac utilizes OncoTree to standardize the input of disease ontology in order to frame a given patient's mutational burden and match assertions from the almanac as closely as possible. 

## About Oncotree
OncoTree was developed by Memorial Sloan Kettering Cancer Center as a tool for the cBioPortal and distributed under a Creative Commons Attribution 4.0 International License. Please visit their [website][1], [swagger][2], or [GitHub][3] to learn more.

## Usage: Generating OncoTree for use
The following command can be run to generate a new OncoTree table for usage with the Molecular Oncology Almanac.

```bash
python get_oncotree.py -d 2023-03-09
```

The output will be of the following name, `oncotree.2023-03-09.txt`

## References
[1]: http://oncotree.mskcc.org/oncotree/#/home "Oncotree"
[2]: http://oncotree.mskcc.org/oncotree/swagger-ui.html "Oncotree Swagger (API)"
[3]: https://github.com/cBioPortal/oncotree "Oncotree Github"
