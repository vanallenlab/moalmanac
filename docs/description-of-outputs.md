# Description of outputs for Molecular Oncology Almanac

All outputs will be produced by Molecular Oncology Almanac, though some may not contain content if the required input file was not passed.

# Table of contents
* [Common output columns used](#common-output-columns-used)
  * [Standardized feature columns](#standardized-feature-columns)
  * [Feature type specific columns](#feature-type-specific-columns)
    * [Somatic and germline variants](#somatic-and-germline-variants)
    * [Fusions](#fusions)
  * [Datasource bins](#datasource-bins)
* [Produced outputs](#produced-outputs)
  * [Germline](#germline)
  * [Integrated summary](#integrated-summary)
  * [Mutational burden](#mutational-burden)
  * [Mutational signatures](#mutational-signatures)
    * [Trinucleotide context counts](#trinucleotide-context-counts)
    * [COSMIC signature (v2) weights](#cosmic-signature-v2-weights)
    * [Trinucleotide context counts image](#trinucleotide-context-counts-image)
    * [Trinucleotide context normalized counts image](#trinucleotide-context-normalized-counts-image)
  * [Report](#report)

# Common output columns used
Many of these columns are used by several outputs produced by Molecular Oncology Almanac. 

## Identifiers
Most, if not all, table outputs will contain these three columns used to describe the case profiles:
* `patient_id` (str) - the string passed by `--patient_id` to label the given molecular profile
* `tumor_sample_barcode` (str) - label present in MAF for tumor sample
* `normal_sample_barcode` (str) - label present in MAF for normal sample

## Standardized feature columns
Molecular Oncology Almanac standardizes primary descriptors for molecular features into four columns: feature_type, feature, alteration_type, and alteration. 

* `feature_type` describes the data type of the molecular features.
  * Potential values are specified under the section [feature_types] in [config.ini](/moalmanac/config.ini)
  * Examples include: Somatic Variant, Germline Variant, Copy Number, Rearrangement, Microsatellite Stability, Mutational Burden, Mutational Signature, Aneuploidy
* `feature` is the main descriptor used to describe a molecular event. The contents vary by feature type, for first-order genomic events this will be the gene name.
  * Somatic variants: gene name
  * Germline variants: gene name
  * Copy number alterations: gene name
  * Rearrangements: gene name, Molecular Oncology Almanac will process each partner in the fusion separately
  * Microsatellite stability: microsatellite stability status (MSI-High or MSI-Low)
  * Mutational burden: High Mutational Burden, if the mutational burden is deemed to be high
  * Mutational signatures: the specific COSMIC (v2) mutational signature, formatted as "COSMIC Signature (number)"
  * Aneuploidy: Whole-genome doubling, this will only be populated if the `--wgd` value is passed to Molecular Oncology Almanac.
* `alteration_type` is a descriptor to provide more granular detail on the molecular event.
  * Somatic variants: variant classification of the variant (Missense, Nonsense, etc.)
  * Germline variants: variant classification of the variant (Missense, Nonsense, etc.)
  * Copy number alterations: type of copy number event (Amplification or Deletion)
  * Rearrangements: type of rearrangement, Fusion unless otherwise specified
* `alteration` is a further description to state the specific event. 
  * Somatic variants: protein change of the variant
  * Germline variants: protein change of the variant
  * Rearrangements: full fusion, "Gene 1--Gene 2"
  * Mutational signatures: contribution or weight of the mutational signature (range: 0 - 1)
  
## Feature type specific columns
Some columns are used to describe specific feature types, they are as follows:

### Somatic and germline variants
The following columns are included to further describe nucleotide variants:
* `chromosome` (str) - the chromosome location of the variant
* `start_position` (int) - genomic location of the start of the genomic event
* `end_position` (int) - genomic location of the end of the genomic event. This will be the same as `start_position` for single nucleotide variants.
* `reference_allele` (string) - reference allele observed
* `observed_allele_1` (string) - observed allele observed, matches reference for single nucleotide variants.
* `observed_allele_2` (string) - observed allele observed, base that does not match reference for single nucleotide variants.
* `tumor_f` (float) - variant allelic fraction, calculated by dividing the count of a given base by the count of all bases present at a genomic site. 
* `total_coverage` (int) - total bases called at a genomic site
* `exac_af` (float) - Allele frequency of allele across all populations in ExAC
* `exac_common` (boolean) - Deemed common if the allele frequency in ExAC is greater than a minimum allele frequency specified in [config.ini](/moalmanac/config.ini). The default value is 0.001, or more common than 1 in 1000 alleles. 
* `clinvar` (string) - Clinical Significance of the genomic event from ClinVar

### Fusions
The following columns are included to further describe fusions
* `spanningfrags` (int) - "the number of RNA-seq fragments that encompass the fusion junction such that one read of the pair aligns toa  different gene than the other paired-end read of that fragment", quoted from the [Star Fusion wiki](https://github.com/STAR-Fusion/STAR-Fusion/wiki)
* `left_gene` (str) - the gene name of the left gene composing a fusion
* `left_chr` (str) - the chromosome location of `left_gene`
* `left_position` (int) - the genomic location that the fusion is observed at in `left_gene`
* `right_gene` (str) - the gene name of the right gene composing a fusion
* `right_chr` (str) - the chromosome location of `right_gene`
* `right_position` (int) - the genomic location that the fusion is observed at in `right_gene`

## Datasource bins
Molecular Oncology Almanac will match molecular features to [several datasources](/moalmanac/datasources) based on how closely a given molecular feature matches a catalogued feature. Specifically, we assign a numeric value to each somatic variant, germline variants, copy number alteration, and fusion for each appropriate datasource. 
* Molecular Oncology Almanac (`almanac_bin`):
  * `0`: the molecular feature's gene is not present in the database
  * `1`: the molecular feature's gene is present in the database
  * `2`: the molecular feature matches at least one entry by both feature type and gene
  * `3`: the molecular feature matches at least one entry by feature type, gene, and alteration_type
  * `4`: the molecular feature matches at least one entry by feature type, gene, alteration_type, and alteration
* Cancer Hotspots (`cancerhotspots_bin`):
  * `0`: the molecular feature's gene is not present in the database
  * `1`: the molecular feature's gene is present in the database
  * `2`: the molecular feature matches at least one entry by both gene and protein change
* 3D Hotspots (`cancerhotspots3D_bin`):
  * `0`: the molecular feature's gene is not present in the database
  * `1`: the molecular feature's gene is present in the database
  * `2`: the molecular feature matches at least one entry by both gene and protein change
* Cancer Gene Census (`cgc_bin`):
  * `0`: the molecular feature's gene is not present in the database
  * `1`: the molecular feature's gene is present in the database
* GSEA Pathways (`gsea_pathways_bin`):
  * `0`: the molecular feature's gene is not present in the database
  * `1`: the molecular feature's gene is present in the database 
* GSEA Modules (`gsea_modules_bin`):
  * `0`: the molecular feature's gene is not present in the database
  * `1`: the molecular feature's gene is present in the database
* COSMIC (`cosmic_bin`):
  * `0`: the molecular feature's gene is not present in the database
  * `1`: the molecular feature's gene is present in the database
  * `2`: the molecular feature matches at least one entry by both gene and protein change

# Produced outputs
The following outputs are produced by the Molecular Oncology Almanac. Each section lists the filename suffix and then a details the contents of the output.

## *.actionable.txt

The file of the suffix `.actionable.txt` is a tab-delimited file and is the primary output produced by Molecular Oncology Almanac. Here, any molecular feature of biological or clinical relevance will be listed along with details of the assertion and citation, for therapeutic sensitivity, resistance, and prognosis.

## Germline
Three germline-specific outputs are produced and populated if a germline MAF is passed (through `--germline_handle`) for a given molecular profile. These three outputs are based on genes present in the [American College of Medical Genetics v2](/moalmanac/datasources/acmg), cancer-related genes (appearing in [Molecular Oncology Almanac](/moalmanac/datasources/almanac), Cancer Hotspots [Molecular Oncology Almanac](/moalmanac/datasources/cancerhotspots), or [Cancer Gene Census](/moalmanac/datasources/cancergenecensus)), or [genes related to heritable cancers](/moalmanac/datasources/hereditary). 

In addition to 

## Integrated summary
Filename suffix: `.integrated.summary.txt`

All genes catalogued in Molecular Oncology Almanac are listed along with a list of alterations that appear in somatic variants, germline variants, copy number alterations, and fusions per gene. This is written as a tab delimited file.

## Mutational burden
Filename suffix: `.mutational_burden.txt`

The mutational burden of a given profile can be computed if somatic variants (`--snv_handle`) and bases covered (`--bases_covered_handle`) are passed to Molecular Oncology Almanac. Somatic insertions and deletions (`--indel_handle`) are also included if they are provided. The mutational burden (also referred to as nonsynonymous mutational burden or coding mutational burden) is calculated by dividing the number of nonsynonymous variants by the number of bases considered for variant calling; Molecular Oncology Almanac uses Missense, Nonsense, Nonstop, Splice Site, Frame Shift, and Insertion and Deletion variants for this calculation.

Molecular Oncology Almanac will attempt to map the provided tumor type (`--tumor_type`) to Oncotree and also an ontology present in TCGA in order to compare the mutational burden relative to TCGA. Columns present in this file are as follows:
- `patient-id` (str) - the string associated with the given molecular profile (`--patient_id`)
- `reported_tumor_type` (str) - the string provided by the user (`--tumor_type`)
- `oncotree_term` (str) - mapped oncotree term for the reported tumor type
- `oncotree_code` (str) - mapped oncotree code for the reported tumor type
- `bases_covered` (int) - number of somatic bases considered for variant calling (`--bases_covered_handle`)
- `n_nonsyn_mutations` (int) - number of somatic coding variants (from `--snv_handle` and `--indel_handle`)
- `coding_mutational_burden_per_megabase` (float) - calculated by dividing `n_nonsyn_mutations` by `bases_covered` and converting the denominator from bases to megabases
- `percentile_tcga` (float, 0-100) - the percentile of `coding_mutational_burden_per_megabase` relative to TCGA from [Lawrence et al. 2013](/moalmanac/datasources/lawrence)
- `percentile_tcga_tissuetype` (float, 0-100) - the percentile of `coding_mutational_burden_per_megabase` relative to a matched TCGA tissue type, if the tumor type matched
- `high_burden?` (boolean, True or False) - boolean value for evaluated mutational burden

Molecular Oncology Almanac designates high mutational burden under two circumstances: 
- Mutations per Mb > 10 
- At least a mutational burden of 80th percentile of TCGA tumor type, if matched, or TCGA generally, if not matched.

## Mutational signatures
Molecular Oncology Almanac runs [deconstructSigs](https://github.com/raerose01/deconstructSigs) as a subprocess based on the MAF file passed with the input argument `--snv_handle`, performing NMF against the 30 COSMIC v2 signatures. 

### Trinucleotide context counts
Filename suffix: `.sigs.context.txt`

Trinucleotide context counts of observed somatic variants for all 96 bins are listed in this tab delimited file.

### COSMIC signature (v2) weights
Filename suffix: `.sigs.cosmic.txt`

Weights for the 30 COSMIC (v2) mutational signatures are listed in this tab delimited file. Thresholds for a signature to be considered present or not present by Molecular Oncology Almanac are specified in [config.ini](/moalmanac/config.ini) under the `[signatures]` heading.

### Trinucleotide context counts image
Filename suffix: `.sigs.tricontext.counts.png`

Trinucleotide context raw counts of observed somatic variants for all 96 bins are visualized in this png file.

### Trinucleotide context normalized counts image
Filename suffix: `.sigs.tricontext.normalized.png`

Trinucleotide context normalized counts of observed somatic variants for all 96 bins are visualized in this png file.

## Report
Filename suffix: `.report.txt`

The report is a web-based file generated with [Flask](https://flask.palletsprojects.com/en/1.1.x/) and [Frozen-Flask](https://pythonhosted.org/Frozen-Flask/) used to display results from the `.actionable.txt` output. An example report generated from the data found in [/example_data/](/example_data) is available in [/example_output/](/example_output/example.report.html). If you ran Molecular Oncology Almanac _not_ through Terra or with [run_example.py](/moalmanac/run_example.py), you may need move the report to the outputs folder by running `mv build/index.html report.html`.

The report contains three primary sections: profile information, the actionability report, and comparisons to cancer cell lines.  

The profile information contains metadata provided on the studied molecular profile. Specifically the following items are included: 
- identifiers such as profile id, tumor barcode, and normal barcode
- disease information such as oncotree code, oncotree term, and stage
- metrics such as purity, ploidy, and microsatellite status

The actionability report contains four sections which list variants and features associated with therapeutic sensitivity, therapeutic resistance, prognosis, and biological relevance. This section also features three models, two which describe how the Molecular Oncology Almanac ranks molecular features and another which lists database versions used by the tool. Each molecular feature associated with a clinical assertion will receive its own row in the report along with the following information:
- Associated evidence for the clinical assertion (FDA-Approved, Guideline, Clinical trial, Clinical evidence, Preclinical evidence, Inferential evidence)
- How closely the molecular feature matches the catalogued molecular feature associated with the clinical assertion (Putatively Actionable, Investigate Actionability)
- The feature type associated with the molecular feature
- The molecular feature
- The associated therapy or prognosis, rationale behind the clinical assertion, and a hyperlink to the source

Additional equivalent within a provided ontology or stronger matches from another ontology are available for viewing as a pop up, too, by selecting the `[More Details]` button.

For molecular features associated with therapeutic sensitivity that have a therapy evaluated on cancer cell lines, a button `[Preclinical evidence]` will appear below the therapy and rationale which will open a modal to compare the sensitivity to the therapy of interest between mutant and wild type cell lines.

Molecular features which are biologically relevant are listed without clinical association. Molecular features will appear here if the associated gene is catalogued in the Molecular Oncology Almanac but under a different feature type, variants are associated with microsatellite stability, and all present COSMIC version 2 mutational signatures not associated with a clinical assertion are reported.

The last section of the report, comparison of molecular profile to cancer cell lines, displays results from Molecular Oncology Almanac's patient-to-cell line matchmaking module. **This will not appear in the report if `--disable_matchmaking` is passed as an argument**. The 5 most similar cancer cell lines to the provided profile are listed each listing the cell line name, sensitive therapies from GDSC, and clinically relevant features present. Users can click `[More details]` under each cell line's name for more details about a given cell line: aliases, sensitive therapies, clinically relevant molecular features, all somatic variants, copy number alterations, and fusions occuring in cancer gene census genes, and the 10 most sensitive therapies to the cancer cell line.
