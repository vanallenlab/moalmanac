# Description of inputs for Molecular Oncology Almanac

Example inputs can be found in the [`example_data/`](/example_data/) folder, found in the root directory of this repository. 

# Table of Contents
The following describes required and optional arguments to run `moalmanac/moalmanac.py`.
* [Required arguments](#required-arguments)  
    - [Patient id](#patient-id)
    - [Config](#config)
    - [Databases](#databases)
* [Optional arguments](#optional-arguments)
    - [Tumor type](#tumor-type)
    - [Stage](#stage)
    - [Somatic single nucleotide variants](#somatic-single-nucleotide-variants)
    - [Somatic insertion and deletion variants](#somatic-insertion-and-deletion-variants)
    - [Bases covered](#bases-covered)
    - [Copy number alterations](#copy-number-alterations)
    - [Fusions](#fusions)
    - [Germline variants](#germline-variants)
    - [Somatic variants from validation sequencing](#somatic-variants-from-validation-sequencing)
    - [Microsatellite status](#microsatellite-status)
    - [Mutational signatures](#mutational-signatures)
    - [Purity](#purity)
    - [Ploidy](#ploidy)
    - [Whole genome doubling](#whole-genome-doubling)
    - [Disable matchmaking](#disable-matchmaking)
    - [Description](#description)
    - [Output directory](#output-directory)
    - [Preclinical databases](#preclinical-databases)

Alternatively, a simplified version of the interpretation algorithm can be run by using `moalmanac/simplified_input.py`. 
- [Simplified input arguments](#simplified-input)
    
# Required arguments
The following arguments are required to run Molecular Oncology Almanac.

## Patient id
`--patient_id` expects a single string value which is used for labeling outputs. 

## Config
`--config` expects a file path to the [config.ini](https://github.com/vanallenlab/moalmanac/blob/main/moalmanac/config.ini) file. 

This config file contains the following sections,
- `function_toggle` - allows several features of the MOAlmanac algorithm to be enabled or disabled
- `versions` - specifies the versions of the [MOAlmanac algorithm (interpreter)](https://github.com/vanallenlab/moalmanac/releases) and [database](https://github.com/vanallenlab/moalmanac-db/releases).    
- `exac` - specifies the allele frequency threshold used with [ExAC](https://github.com/vanallenlab/moalmanac/tree/main/datasources/exac) to specify if a variant is a common variant or not
- `fusion` - specifies minimum spanning fragments required for review by MOAlmanac, column names expected from inputs, and how "Fusion" should be written from input
- `mutations` - specifies the minimum coverage and allelic fraction that a variant needs for review by MOAlmanac
- `seg` - specifies the percentile to evaluate copy gain and loss variants from segmented copy number input files, as well as how amplification and deletion should be written as strings
- `signatures` - specifies the minimum contribution required to review COSMIC mutational signatures by mMOAlmanac
- `validation_sequencing` - Thresholds for minimum power to detect variants and minimum allelic fraction for annotation from validation sequencing. This is further described in the [Methods section](https://www.nature.com/articles/s43018-021-00243-3#Sec8) of our paper.
- `feature_types` - String labels for each biomarker type passed to the algorithm. These values will be included in `feature_type` column of outputs.

## Databases
`--dbs` expects a file path to the [annotation-databases.ini](../moalmanac/annotation-databases.ini) file.

This config file contains a single section `databases` that lists the following:
- `root` - path to `datasources/` directory
- `almanac_handle` - path within `root` that points to the `molecular-oncology-almanac.json` datasource file
- `cancerhotspots_handle` - path within `root` that points to the Cancer Hotspots datasource file
- `3dcancerhotspots_handle` - path within `root` that points to the Cancer Hotspots 3D datasource file
- `cgc_handle` - path within `root` that points to the Cancer Gene Census file
- `cosmic_handle` - path within `root` that points to the COSMIC datasource file
- `gsea_pathways_handle` - path within `root` that points to the GSEA pathways datasource file
- `gsea_modules_handle` - path within `root` that points to the GSEA modules datasource file
- `exac_handle` - path within `root` that points to the ExAC datasource file
- `acmg_handle` - path within `root` that points to the ACMG datasource file
- `clinvar_handle` - path within `root` that points to the ClinVar datasource file
- `hereditary_handle` - path within `root` that points to the genes related to hereditary cancers datasource file
- `oncotree_handle` - path within `root` that points to the Oncotree datasource file
- `lawrence_handle` - path within `root` that points to the Lawrence et al. TCGA mutational burden datasource file

For more information about each datasource, view the [datasources directory](../datasources/README.md)

# Optional arguments
Molecular Oncology Almanac will run successfully given any combination of the following arguments: 

## Tumor type
`--tumor_type` expects a string representing the report tumor type. MOAlmanac will attempt to map this string to either an [Oncotree term or code](https://github.com/vanallenlab/moalmanac/tree/main/moalmanac/datasources/oncotree). MOAlmanac will consider clinically relevant matches of the same tumor type prior to considering matches of another tumor type. 

## Stage
`--stage` also expects a string and is intended for use to input disease stage. This is not functionally used within the method and only is outputted for display in the produced actionability report.

## Somatic single nucleotide variants 
`--snv_handle` anticipates a tab delimited file which contains somatic single nucleotide variants (snvs). This file should follow the guideline's set by either TCGA or GDC [Mutation Annotation Format (MAF)](https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/). Insertions and deletions can also be included in this input, the two MAFs are simply concatenated together. 

### Example
|Hugo Symbol|NCBI_Build|Chromosome|Start_position|End_position|Reference_Allele|Tumor_Seq_Allele1|Tumor_Seq_Allele2|Variant_Classification|Tumor_Sample_Barcode|Matched_Norm_Sample_Barcode|Annotation_Transcript|Protein_Change|t_ref_count|t_alt_count|  
|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|
|BRAF|37|7|140453136|140453136|A|T|A|Missense_Mutation|ProfileA-Tumor|ProfileA-Normal|ENST00000288602.6|p.V600E|70|35|
|MSH2|37|2|47739466|47739466|G|A|G|Missense_Mutation|ProfileA-Tumor|ProfileA-Normal|ENST00000406134.1|p.D887N|50|25| 
|STAG2|37|X|123191810|123191810|A|T|A|Missense_Mutation|ProfileA-Tumor|ProfileA-Normal|ENST00000371160.1|p.F467I|60|20|

### Required fields
Required fields can be changed from their default expectations by editing the appropriate section of [colnames.ini](https://github.com/vanallenlab/moalmanac/blob/main/moalmanac/colnames.ini). Column names are **not** case-sensitive. 
- `Hugo_Symbol`, gene symbol associated with the variant
- `NCBI_Build`, reference genome used
- `Chromosome`, chromosome of the variant
- `Start_position`, genomic start position of the variant
- `End_position`, genomic end position of the variant
- `Reference_Allele`, reference allele at the genomic location
- `Variant_Classification`, consequence of variant: `Missense`, `Nonsense`, `Nonstop`, `Splice_Site`, `Frame_Shift_Ins`, `Frame_Shift_Del`, `In_Frame_Ins`, or `In_Frame_Del`
- `Tumor_Seq_Allele1`, alternate allele at the genomic location
- `Tumor_Seq_Allele2`, second allele at the genomic location (will be the same as `Reference_Allele` for SNVs)
- `Tumor_Sample_Barcode`, string associated with the tumor profile
- `Matched_Norm_Sample_Barcode`, string associated with the corresponding normal profile
- `Annotation_Transcript`, transcript associated with variant
- `t_ref_count`, number of reference alleles observed at genomic position
- `t_alt_count`, number of alternate alleles observed at genomic position

At least one of the following also must be included: 
- `HGVSp_Short`, protein change associated with the variant using the one-letter amino acid codes
- `Protein_Change`, protein change associated with the variant using the one-letter amino-acid codes 

## Somatic insertion and deletion variants 
`--indel_handle` anticipates a tab delimited file which contains somatic insertions and deletions (indels). This file should follow the guideline's set by either TCGA or GDC [Mutation Annotation Format (MAF)](https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/). Single nucleotide variants can also be included in this input, the two MAFs are simply concatenated together. 

### Example
|Hugo Symbol|NCBI_Build|Chromosome|Start_position|End_position|Reference_Allele|Tumor_Seq_Allele1|Tumor_Seq_Allele2|Variant_Classification|Tumor_Sample_Barcode|Matched_Norm_Sample_Barcode|Annotation_Transcript|Protein_Change|t_ref_count|t_alt_count|  
|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|
|PMPCA|37|9|139312448|139312449|-|G|-|Intron|ProfileA-Tumor|ProfileA-Normal|ENST00000371717.3||92|31|
|C10orf2|37|10|102748300|102748301|TC|TC|-|Frame_Shift_Del|ProfileA-Tumor|ProfileA-Normal|ENST00000370228.1|p.L112fs|294|28| 
|MEST|37|7|130138285|130138285|C|C|-|Frame_Shift_Del|ProfileA-Tumor|ProfileA-Normal|ENST00000223215.4|p.L168fs|60|20|

### Required fields
Required fields can be changed from their default expectations by editing the appropriate section of [colnames.ini](https://github.com/vanallenlab/moalmanac/blob/main/moalmanac/colnames.ini). Column names are **not** case-sensitive. 
- `Hugo_Symbol`, gene symbol associated with the variant
- `NCBI_Build`, reference genome used
- `Chromosome`, chromosome of the variant
- `Start_position`, genomic start position of the variant
- `End_position`, genomic end position of the variant
- `Reference_Allele`, reference allele at the genomic location
- `Variant_Classification`, consequence of variant: `Missense`, `Nonsense`, `Nonstop`, `Splice_Site`, `Frame_Shift_Ins`, `Frame_Shift_Del`, `In_Frame_Ins`, or `In_Frame_Del`
- `Tumor_Seq_Allele1`, alternate allele at the genomic location
- `Tumor_Seq_Allele2`, second allele at the genomic location (will be the same as `Reference_Allele` for SNVs)
- `Tumor_Sample_Barcode`, string associated with the tumor profile
- `Matched_Norm_Sample_Barcode`, string associated with the corresponding normal profile
- `Annotation_Transcript`, transcript associated with variant
- `t_ref_count`, number of reference alleles observed at genomic position
- `t_alt_count`, number of alternate alleles observed at genomic position

At least one of the following also must be included: 
- `HGVSp_Short`, protein change associated with the variant using the one-letter amino acid codes
- `Protein_Change`, protein change associated with the variant using the one-letter amino-acid codes 

## Bases covered
`--bases_covered_handle` anticipates a tab delimited file which contains a single integer representing the number of bases tested for somatic variants. This is used as the denominator to calculate tumor mutational burden (number of coding somatic variants / somatic bases tested). 

### Example
|29908096|
|-|

### Required fields
This input is looking for an integer value. 

## Called copy number alterations
`--called_cn_handle` anticipated a tab delimited file which contains one column for gene name and a second for the copy number call. For the latter, only the values `Amplification` and `Deletion` will be used by the Molecular Oncology Almanac.

### Example
|gene|call|
|-|-|
|TP53|Deletion|
|CDKN2A|Deletion|
|BRAF|Baseline|
|EGFR|Amplification|

The rows associated with _TP53_, _CDKN2A_, and _EGFR_ will be interpreted and scored by Molecular Oncology Almanac while _BRAF_ will be filtered.

### Required fields
Required fields can be changed from their default expectations by editing the appropriate section of [colnames.ini](https://github.com/vanallenlab/moalmanac/blob/main/moalmanac/colnames.ini). Column names are **not** case-sensitive. 
- `gene`, gene symbol associated with the copy number alteration
- `call`, copy number event of the gene. `Amplification` and `Deletion` are accepted and all other values will be filtered.

## Copy number alterations
`--cnv_handle` anticipates a tab delimited file which contains total copy number from a source such as GATK CNV or ReCapSeg, support for allele specific copy number is in progress. This file should have genes associated with segments. Amplifications are called from the top 2.5% of all unique segments and deletions called from the bottom 2.5% of all unique segments. 

### Example
|gene|segment_contig|segment_start|segment_end|sample|segment_mean|  
|-|-|-|-|-|-|
|BRAF|7|140035556|142013739|ProfileA-Tumor|1.250062303|
|CDKN2A|9|21818453|27173612|ProfileA-Tumor|0.822108092|
|BOC|3|112282632|113393977|ProfileA-Tumor|0.957205107|

### Required fields
Required fields can be changed from their default expectations by editing the appropriate section of [colnames.ini](https://github.com/vanallenlab/moalmanac/blob/main/moalmanac/colnames.ini). Column names are **not** case-sensitive. 
- `gene`, gene symbol associated with the copy number alteration
- `segment_contig`, chromosome of the copy number alteration
- `segment_start`, genomic location of the segment's start position
- `segment_end`, genomic location of the segment's end position
- `sample`, string associated with the tumor profile
- `segment_mean`, normalized segment mean 

## Fusions
`--fusion_handle` anticipates a tab delimited file which contains fusions, specifically in the format of STAR Fusion. 

### Example
|#FusionName|SpanningFragCount|LeftBreakpoint|RightBreakpoint|
|-|-|-|-|
|EML4--ALK|0|6:47471176|11:66563752|
|COL1A2--APBA3|6|9:35657873|21:46320255|
|POLR2A--AP2M1|12|17:7406801|3:183898675|

### Required fields
Required fields can be changed from their default expectations by editing the appropriate section of [colnames.ini](https://github.com/vanallenlab/moalmanac/blob/main/moalmanac/colnames.ini). Column names are **not** case-sensitive. 
- `#FusionName`, gene symbols associated with the fusion separated by `--`. Genes are labeled from 5' to 3'. 
- `SpanningFragCount`, counts of RNA-seq fragments supporting the fusion
- `LeftBreakpoint`, genomic position of the fusion's left breakpoint
- `RightBreakpoint`, genomic position of the fusion's right breakpoint

## Germline variants
`--germline_handle` anticipates a tab delimited file which contains germline variants (both snvs and indels) associated with a case profile. This file should follow the guideline's set by either TCGA or GDC [Mutation Annotation Format (MAF)](https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/). Column names are **not** case-sensitive.

### Example
|Hugo Symbol|NCBI_Build|Chromosome|Start_position|End_position|Reference_Allele|Tumor_Seq_Allele1|Tumor_Seq_Allele2|Variant_Classification|Tumor_Sample_Barcode|Matched_Norm_Sample_Barcode|Annotation_Transcript|Protein_Change|t_ref_count|t_alt_count|  
|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|
|BRAF|37|7|140453136|140453136|A|T|A|Missense_Mutation|ProfileA-Tumor|ProfileA-Normal|ENST00000288602.6|p.V600E|40|25|
|MSH2|37|2|47739466|47739466|G|A|G|Missense_Mutation|ProfileA-Tumor|ProfileA-Normal|ENST00000406134.1|p.D887N|40|30| 
|STAG2|37|X|123191810|123191810|A|T|A|Missense_Mutation|ProfileA-Tumor|ProfileA-Normal|ENST00000371160.1|p.F467I|80|26|

### Required fields
Required fields can be changed from their default expectations by editing the appropriate section of [colnames.ini](https://github.com/vanallenlab/moalmanac/blob/main/moalmanac/colnames.ini). Column names are **not** case-sensitive. 
- `Hugo_Symbol`, gene symbol associated with the variant
- `NCBI_Build`, reference genome used
- `Chromosome`, chromosome of the variant
- `Start_position`, genomic start position of the variant
- `End_position`, genomic end position of the variant
- `Reference_Allele`, reference allele at the genomic location
- `Variant_Classification`, consequence of variant: `Missense`, `Nonsense`, `Nonstop`, `Splice_Site`, `Frame_Shift_Ins`, `Frame_Shift_Del`, `In_Frame_Ins`, or `In_Frame_Del`
- `Tumor_Seq_Allele1`, alternate allele at the genomic location
- `Tumor_Seq_Allele2`, second allele at the genomic location (will be the same as `Reference_Allele` for SNVs)
- `Tumor_Sample_Barcode`, string associated with the tumor profile
- `Matched_Norm_Sample_Barcode`, string associated with the corresponding normal profile
- `Annotation_Transcript`, transcript associated with variant
- `t_ref_count`, number of reference alleles observed at genomic position
- `t_alt_count`, number of alternate alleles observed at genomic position

At least one of the following also must be included: 
- `HGVSp_Short`, protein change associated with the variant using the one-letter amino acid codes
- `Protein_Change`, protein change associated with the variant using the one-letter amino-acid codes 

## Somatic variants from validation sequencing
`--validation_handle` anticipates a tab delimited file which contains somatic variants (snvs and/or indels) from any form of validation or orthogonal sequencing on the tumor; such as re-sequencing the same tissue or somatic variants called from RNA. This file should follow the guideline's set by either TCGA or GDC [Mutation Annotation Format (MAF)](https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/). 

Variants from this file are only used for confirmation and are not used for discovery. Specifically, MOAlmanac will look for reported somatic variants in the validation sequencing and identify if any supporting reads are present and if there is sufficient power for detection. This is consistent with best practices recommended by [Yizhak et al. 2019](https://science.sciencemag.org/content/364/6444/eaaw0726.long).

### Example
|Hugo Symbol|NCBI_Build|Chromosome|Start_position|End_position|Reference_Allele|Tumor_Seq_Allele1|Tumor_Seq_Allele2|Variant_Classification|Tumor_Sample_Barcode|Matched_Norm_Sample_Barcode|Annotation_Transcript|Protein_Change|t_ref_count|t_alt_count|  
|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|
|BRAF|37|7|140453136|140453136|A|T|A|Missense_Mutation|ProfileA-Tumor|ProfileA-Normal|ENST00000288602.6|p.V600E|40|25|
|MSH2|37|2|47739466|47739466|G|A|G|Missense_Mutation|ProfileA-Tumor|ProfileA-Normal|ENST00000406134.1|p.D887N|40|30| 
|STAG2|37|X|123191810|123191810|A|T|A|Missense_Mutation|ProfileA-Tumor|ProfileA-Normal|ENST00000371160.1|p.F467I|80|26|

### Required fields
Required fields can be changed from their default expectations by editing the appropriate section of [colnames.ini](https://github.com/vanallenlab/moalmanac/blob/main/moalmanac/colnames.ini). Column names are **not** case-sensitive. 
- `Hugo_Symbol`, gene symbol associated with the variant
- `NCBI_Build`, reference genome used
- `Chromosome`, chromosome of the variant
- `Start_position`, genomic start position of the variant
- `End_position`, genomic end position of the variant
- `Reference_Allele`, reference allele at the genomic location
- `Variant_Classification`, consequence of variant: `Missense`, `Nonsense`, `Nonstop`, `Splice_Site`, `Frame_Shift_Ins`, `Frame_Shift_Del`, `In_Frame_Ins`, or `In_Frame_Del`
- `Tumor_Seq_Allele1`, alternate allele at the genomic location
- `Tumor_Seq_Allele2`, second allele at the genomic location (will be the same as `Reference_Allele` for SNVs)
- `Tumor_Sample_Barcode`, string associated with the tumor profile
- `Matched_Norm_Sample_Barcode`, string associated with the corresponding normal profile
- `Annotation_Transcript`, transcript associated with variant
- `t_ref_count`, number of reference alleles observed at genomic position
- `t_alt_count`, number of alternate alleles observed at genomic position

At least one of the following also must be included: 
- `HGVSp_Short`, protein change associated with the variant using the one-letter amino acid codes
- `Protein_Change`, protein change associated with the variant using the one-letter amino-acid codes

## Microsatellite status
`--ms_status` is a categorical input for microsatellite status, anticipating one of four options: 
- `--unknown`, when status is unknown
- `--mss` for microsatellite stable tumors, MSS
- `--msil` for microsatellite instability "low", MSI-L
- `--msih` for microsatellite instability "high", MSI-H 

Microsatellite status is reported in the clinical actionability report. 

## Mutational signatures
`--mutational_signatures` anticipates a tab delimited file which contains contributions to Single Base Substitution (SBS) Mutational Signatures from [COSMIC version 3.4](https://cancer.sanger.ac.uk/signatures/sbs/). The file should only contain signature contributions for the tumor sample being studied. We recommend generating SBS mutational signatures with [SigProfilerAssignment](https://github.com/AlexandrovLab/SigProfilerAssignment), and have prepared [a wrapper GitHub repository](https://github.com/vanallenlab/SigProfilerAssignment-wrapper) to run SigProfilerAssignment and format signature contributions as expected.  

### Example
| signature | contribution | 
|---|--------------|
| SBS1 | 0.03846154   |
| SBS2 | 0            |
| SBS3 | 0.8525641    |
| ... | ...          |
| SBS95 | 0            |

### Required fields,
The required fields for this file can be changed from their default expectations by editing the appropriate section of `colnames.ini`. Column names are **not** case sensitive.
- `signature`, labels for each of the 79 SBS mutational signatures included in COSMIC mutational signatures [version 3.4](https://cancer.sanger.ac.uk/signatures/sbs/) 
- `contribution`, a float value between 0 and 1 for the row's associated signature weight. This column's values should sum to 1. 

## Purity
`--purity` anticipates a float value between 0.0 and 1.0 for the reported tumor purity. This is just used for reporting in the clinical actionability report.

## Ploidy
`--purity` anticipates a float value for the reported tumor ploidy. This is just used for reporting in the clinical actionability report.

## Whole genome doubling
`--wgd` is a boolean argument that, when passed, is interpreted as the tumor harboring a whole-genome doubling event. If passed, MOAlmanac will match against relevant assertions. 

## Disable matchmaking
`--disable_matchmaking` removes patient-to-cell line matchmaking from being included in the actionability report. 

## Description
`--description` is a string input that is generally a free text field for users to enter any additional comments. 

## Output directory
`--output-directory` allows users to specify an output directory to write outputs to, the current working directory will be used if unspecified. 

## Preclinical databases
`--preclinical-dbs` expects a file path to the [preclinical-databases.ini](../moalmanac/preclinical-databases.ini) file.

This config file contains a single section `preclinical` that lists the following:
- `root` - path to `datasources/preclinical/` directory
- `almanac_gdsc_mappings` - path within `root` that points to the `formatted/almanac-gdsc-mappings.json` datasource file
- `summary` - path within `root` that points to the `formatted/cell-lines.summary.txt` datasource file
- `variants` - path within `root` that points to the `annotated/cell-lines.somatic-variants.annotated.txt` datasource file
- `copynumbers` - path within `root` that points to the `annotated/cell-lines.copy-numbers.annotated.txt` file
- `fusions` - path within `root` that points to the `annotated/cell-lines.fusions.annotated.txt` datasource file
- `fusions1` - path within `root` that points to the `annotated/cell-lines.fusions.annotated.gene1.txt` datasource file
- `fusions2` - path within `root` that points to the `annotated/cell-lines.fusions.annotated.gene2.txt` datasource file
- `gdsc` - path within `root` that points to the `formatted/sanger.gdsc.txt` datasource file
- `dictionarey` - path within `root` that points to the `cell-lines.pkl` datasource file

For more information about each datasource, view the [datasources/preclinical/ directory](../datasources/preclinical/README.md)

# Simplified input
`--input` is an argument only used with `simplified_input.py`. It accepts a tab delimited file with one genomic alteration per row based on MOAlmanac's [standardized feature columns](../docs/description-of-outputs.md#standardized-feature-columns). In short the following columns are expected,
1. `feature_type`, the data type of the molecular features and accepts `Somatic Variant`, `Germline Variant`, `Copy Number`, or `Rearrangement`. These strings can be customized in the `feature_types` section of [config.ini](config.ini).
2. `gene` or `feature`, the gene name of the genomic alteration.
3. `alteration_type`, classification or consequence of the genomic alteration
    - For somatic and germline variants: `Missense`, `Nonsense`, `Nonstop`, `Splice_Site`, `Frame_Shift_Ins`, `Frame_Shift_Del`, `In_Frame_Ins`, or `In_Frame_Del`
    - For copy number alterations: `Amplification` or `Deletion`
    - For rearrangements: `Fusion` or `Translocation`
4. `alteration`, specific genomic alteration,
    - For somatic and germline variants: the protein change with [1-letter amino acid codes](https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/), `p.HGVSp_Short`
    - For copy number alterations: Leave blank
    - For rearrangements: the full fusion separated by two dashes, `--`

For example,
| feature_type | feature | alteration_type | alteration |
|---|---|---|---|
| Somatic Variant | BRAF | Missense | p.V600E |
| Copy Number | CDK4 | Amplification |  |
| Rearrangement | COL1A1 | Fusion | COL1A1--CITED4 |
| Germline Variant | BRCA2 | Frameshift | p.S1982fs |

If you use this method, please cite our publication:
> [Reardon, B., Moore, N.D., Moore, N.S., *et al*. Integrating molecular profiles into clinical frameworks through the Molecular Oncology Almanac to prospectively guide precision oncology. *Nat Cancer* (2021). https://doi.org/10.1038/s43018-021-00243-3](https://www.nature.com/articles/s43018-021-00243-3)
