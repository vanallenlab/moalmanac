# Databases: Exome Aggregation Consortium (ExAC)
The Molecular Oncology Almanac utilizes [ExAC](http://exac.broadinstitute.org/) to annotate population allele frequency of alterations and identify variants that may not be related to cancer. Specifically, our heuristic will sort against alterations that appear in ExAC. While alterations will stay in their same bin,they will appear lower than other variants that do not appear in ExAC.

## About the Exome Aggregation Consortium
The Exome Aggregation Consortium was developed and is maintained by [several principle investigators](http://exac.broadinstitute.org/about). It is a dataset of germline sequencing data from over 60,000 unrelated individuals. 

## Usage: Downloading and formatting ExAC for use
All releases of ExAC [are available for download on their webpage](ftp://ftp.broadinstitute.org/pub/ExAC_release/) and, as of this writing, The Molecular Oncology Almanac is built for release with [release 1](ftp://ftp.broadinstitute.org/pub/ExAC_release/release1/). The file of the extention `.sites.vep.vcf.gz`, with corresponding index `.sites.vep.vcf.gz.tbi`, should be downloaded, this may take some time since the file is on the order of ~ 4.6 GB in size.

The GATK tool [VariantsToTable](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_gatk_tools_walkers_variantutils_VariantsToTable.php) is then used by the Molecular Oncology Almanac to extract variants and desired columns which have passed the VQSR sensitivity filter and hard filters, as annotated by ExAC. The script `build_exac.sh` is then used to turn the vcf into a format readable by the Molecular Oncology Almanac. 

The following inputs are used:
1. ExAC's raw vcf
2. Reference genome fasta
3. GATK jar file

```
bash build_exac.sh ExAC.r1.sites.vep.vcf Homo_sapeisn_assembly19.fasta GenomeAnalysisTK-3.8.jar
```

After converting ExAC to a tab delimited file, sites with multiple alternate alleles are expanded to make up their own rows, using the Python script `expand_exac.py`. There are approximately ~650,000 sites in ExAC with multiple alternate alleles, so running this script may take a couple of hours. It may be better to simply copy this output off of the docker file. 
```
python expand_exac.py --exac exac.lite-pass.1.4-r1.txt
```

## References
1. [Lek M, Karczewski KJ, Minikel EV, et al. Analysis of protein-coding genetic variation in 60,706 humans. Nature. 2016;536(7616):285-91.](https://www.nature.com/articles/nature19057)
2. [Mckenna A, Hanna M, Banks E, et al. The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. Genome Res. 2010;20(9):1297-303.](http://genome.cshlp.org/content/20/9/1297)
