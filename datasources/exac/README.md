# Databases: Exome Aggregation Consortium (ExAC)
The Molecular Oncology Almanac utilizes [ExAC](http://exac.broadinstitute.org/) to annotate population allele frequency of alterations and identify variants that may not be related to cancer. Specifically, our heuristic will sort against alterations that appear in ExAC. While alterations will stay in their same bin, they will appear lower than other variants that do not appear in ExAC.

## About the Exome Aggregation Consortium
The Exome Aggregation Consortium was developed and is maintained by [several principle investigators](http://exac.broadinstitute.org/about). It is a dataset of germline sequencing data from over 60,000 unrelated individuals. 

ExAC is [released openly and publicly](https://gnomad.broadinstitute.org/policies). 

## Usage: Downloading and formatting ExAC for use
All releases of ExAC [are available for download on their webpage](https://gnomad.broadinstitute.org/downloads#exac) and, as of this writing, The file of the extension `.sites.vep.vcf.gz`, with corresponding index `.sites.vep.vcf.gz.tbi`, should be downloaded, this may take some time since the file is on the order of ~ 4.6 GB in size. The [GATK](https://gatk.broadinstitute.org/hc/en-us) tool [VariantsToTable](https://gatk.broadinstitute.org/hc/en-us/articles/360036711531-VariantsToTable) is then used by the Molecular Oncology Almanac to extract variants and desired columns which have passed the VQSR sensitivity filter and hard filters, as annotated by ExAC. The script `build_exac.sh` is then used to turn the vcf into a format readable by the Molecular Oncology Almanac.

The following steps should be performed to prepare ExAC for use with MOAlmanac,
1. Download [ExAC](https://gnomad.broadinstitute.org/downloads#exac-variants) from gnomAD's webpage, titled "All chromosomes VCF" under the "Exomes" section. Download this and the TBI (the VCF index)
2. Download [GATK](https://gatk.broadinstitute.org/hc/en-us)
  - Their website only has downloads for GATK4 currently. GATK3.8 can be downloaded [from their archives](https://console.cloud.google.com/storage/browser/gatk-software/package-archive/gatk?pli=1)
  - `gs://gatk-software/package-archive/gatk/GenomeAnalysisTK-3.8-0-ge9d806836.tar.bz2`
3. Download [hg19 reference genome files](https://console.cloud.google.com/storage/browser/gcp-public-data--broad-references/hg19/v0;tab=objects?prefix=&forceOnObjectsSortingFiltering=false) from [gcp-public-data--broad-references](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811). In particular, 
  - `gs://gcp-public-data--broad-references/hg19/v0/Homo_sapiens_assembly19.fasta`
  - `gs://gcp-public-data--broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai`
  - `gs://gcp-public-data--broad-references/hg19/v0/Homo_sapiens_assembly19.dict`

Place these files in this directory and execute the `build_exac.sh` script, 

```
bash build_exac.sh ExAC.r1.sites.vep.vcf Homo_sapiens_assembly19.fasta GenomeAnalysisTK-3.8.jar
```

After converting ExAC to a tab delimited file, sites with multiple alternate alleles are expanded to make up their own rows, using the Python script `expand_exac.py`. There are approximately ~650,000 sites in ExAC with multiple alternate alleles, so running this script may take a couple of hours. It may be better to simply copy this output off of the docker file. 
```
python expand_exac.py --exac exac.lite-pass.1.4-r1.txt
```

This script will produce an output named `exac.expanded.r1.txt` that is about 923 MB in size.

If you do not have access to Google Cloud or Docker and are having trouble building this datasource, please reach out. We are happy to try our best to help figure something out.

## References
1. [Lek M, Karczewski KJ, Minikel EV, et al. Analysis of protein-coding genetic variation in 60,706 humans. Nature. 2016;536(7616):285-91.](https://www.nature.com/articles/nature19057)
2. [Mckenna A, Hanna M, Banks E, et al. The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. Genome Res. 2010;20(9):1297-303.](http://genome.cshlp.org/content/20/9/1297)
