# Directly leveraging cancer cell lines for clinical interpretation
Molecular Oncology Almanac leverages cancer cell lines for clinical interpretation in two ways, 
- [To test if relationships between molecular features and therapies show efficacy in cancer cell lines](../../../docs/description-of-outputs.md#preclinical-efficacy)
- [To perform profile-to-cell line matchmaking](../../../docs/description-of-outputs.md#profile-to-cell-line-matchmaking)

The following files must be configured for use with MOAlmanac,
- `formatted/almanac.gdsc.mappings.json` (14 KB)
   - Mappings between therapies cataloged in MOAlmanac and therapies utilized in GDSC
- `formatted/cell-lines.summary.txt` (260 KB)
   - A table of cancer cell line names across data sets, which data type they have available, and if they are used in any analyses
- `cell-lines.copy-numbers.annotated.txt` (3.3 MB)
   - A table of copy number alterations observed in cancer cell lines, annotated with MOAlmanac
- `cell-lines.somatic-variants.annotated.txt` (3.8 MB)
   - A table of somatic variants observed in cancer cell lines, annotated with MOAlmanac
- `cell-lines.fusions.gene1.annotated.txt` (314 KB)
   - A table of fusions observed in cancer cell lines, annotated with MOAlmanac relative to gene 1
- `cell-lines.fusions.gene2.annotated.txt` (311 KB)
   - A table of fusions observed in cancer cell lines, annotated with MOAlmanac relative to gene 2
- `cell-lines.fusions.annotated.txt` (363 KB)
   - A table of fusions observed in cancer cell lines, using annotations of the most clinically and biologically relevant gene per fusion
- `cell-lines.pkl` (7.8 MB)
   - A file containing summary information of the above, used when creating reports

**These files should be reproduced if the underlying [moalmanac/](../moalmanac/) is updated**.

## Usage: downloading and formatting preclinical data
Please follow the following steps to configure raw data for use with MOAlmanac,
1. Follow instructions under `source/` to download data used in the present study
2. Follow instructions under `formatted/` to format samples and molecular features
3. Follow instructions under `annotated/` to annotate molecular features after formatting
4. Open and execute the notebook `generate-dictionary.ipynb` to produce `cell-lines.pkl`, for easy look up by MOAlmanac

Much of these files and scripts have been repurposed and modified from the [MOAlmanac paper Github repository](https://github.com/vanallenlab/moalmanac-paper).

## References
1. [Ghandi, M. et al. Next-generation characterization of the Cancer Cell Line Encyclopedia. Nature 569, 503–508 (2019).](https://www.nature.com/articles/s41586-019-1186-3)
2. [Yang, W. et al. Genomics of Drug Sensitivity in Cancer (GDSC): a resource for therapeutic biomarker discovery in cancer cells. Nucleic Acids Res. 41, D955–61 (2013).](https://academic.oup.com/nar/article/41/D1/D955/1059448)
3. [Sondka, Z. et al. The COSMIC Cancer Gene Census: describing genetic dysfunction across all human cancers. Nat. Rev. Cancer 18, 696–705 (2018).](https://www.nature.com/articles/s41568-018-0060-1)
4. [Wang, B. et al. Similarity network fusion for aggregating data types on a genomic scale. Nat. Methods 11, 333–337 (2014).](https://www.nature.com/articles/nmeth.2810)
5. [Reardon, B., Moore, N.D., Moore, N.S., *et al*. Integrating molecular profiles into clinical frameworks through the Molecular Oncology Almanac to prospectively guide precision oncology. *Nat Cancer* (2021).](https://www.nature.com/articles/s43018-021-00243-3)
6. [Reardon, B. & Van Allen, E. M. Molecular profile to cancer cell line matchmaking. Protocol Exchange.](https://protocolexchange.researchsquare.com/article/pex-1539/v1)
