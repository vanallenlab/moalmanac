# Databases: Cancer Hotspots
The Molecular Oncology Almanac utilizes [Cancer Hotspots](http://cancerhotspots.org/#/home) and [3d Hotspots](http://3dhotspots.org/3d/#/home) to evaluate alterations that may be hotspots in Cancer, even if they are not present in our action-alteration database. Specifically, our heuristic evaluates alterations at the gene and residue level with respect to Cancer Hotspots & 3d Hotspots. 

At the time of this writing the Molecular Oncology Almanac leverages v2 of Cancer Hotspots.

## About Cancer Hotspots & 3d Hotspots
Both Cancer Hotspots and 3d Hotspots were developed and are maintained by Memorial Sloan Kettering Cancer Center as resources for identifying statistically significant mutations in cancer. Cancer Hotspots were originally reported for [single residue mutation hotspots identified in 11,119 tumor samples in 2016](https://www.ncbi.nlm.nih.gov/pubmed/26619011) and was recently updated to include hotspots for [both single residue mutations and in-frame insertions/deletions identified in 24,592 tumor samples in 2017](https://www.ncbi.nlm.nih.gov/pubmed/29247016). This update brought with it the removal of 184 genes and addition of 152 and the unique number of residues increased to 3004, up from 1168. The 2016 study later was expanded by [clustering single nucleotide variants in 3d space on protein structures](https://www.ncbi.nlm.nih.gov/pubmed/28115009) to identify additional hotspots.

## Usage: Downloading Cancer Hotspots
Data from both the 2016 and 2017 releases of Cancer Hotspots are available for [download as .xls files directly from their webpage](http://cancerhotspots.org/#/download) and [likewise for 3d Hotspots](http://3dhotspots.org/3d/#/download). An [API](http://cancerhotspots.org/swagger-ui.html) is also available, but at the time of this writing the direct download was used for retrieval. 

## References
1. [Chang MT, Asthana S, Gao SP, et al. Identifying recurrent mutations in cancer reveals widespread lineage diversity and mutational specificity. Nat Biotechnol. 2016;34(2):155-63.](https://www.ncbi.nlm.nih.gov/pubmed/26619011)
2. [Chang MT, Shrestha bhattarai T, Schram AM, et al. Accelerating discovery of functional mutant alleles in cancer. Cancer Discov. 2017;](https://www.ncbi.nlm.nih.gov/pubmed/29247016)
3. [Gao J, Chang MT, Johnsen HC, et al. 3D clusters of somatic mutations in cancer reveal numerous rare mutations as functional targets. Genome Med. 2017;9(1):4.](https://www.ncbi.nlm.nih.gov/pubmed/28115009)
