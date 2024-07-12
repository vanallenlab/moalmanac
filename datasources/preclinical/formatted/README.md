# Formatted cancer cell line data
Cancer cell line data from external sources must be formatted for annotation and use by the Molecular Oncology Almanac. 

The jupyter notebook `0.process-source-data.ipynb` was modified from [02.process-cell-line-molecular-features](https://github.com/vanallenlab/moalmanac-paper/blob/main/analyses/preclinical/02.process-cell-line-molecular-features.ipynb) in the [moalmanac paper Github repository](https://github.com/vanallenlab/moalmanac-paper). We copied `cell-line-names.formatted.txt` from the paper repository and used this notebook to generate the following files,
- samples to consider: `cell-lines.summary.txt`
- somatic variants: `cell-lines.somatic_variants.txt`
- called copy number alterations: `cell-lines.copy-numbers.txt`
- fusions: `cell-lines.fusions.txt`
- samples with therapy response: `sanger.gdsc.txt`

The jupyter notebook `1.map-almanac-to-gdsc.ipynb` was modified from [00.map-almanac-to-gdsc.ipynb](https://github.com/vanallenlab/moalmanac-paper/blob/main/analyses/preclinical/00.map-almanac-to-gdsc.ipynb) in the [moalmanac paper Github repository](https://github.com/vanallenlab/moalmanac-paper). This notebook was used to generate the following file,
- map between therapies: `almanac-gdsc-mappings.json`
