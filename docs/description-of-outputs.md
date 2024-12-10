# Description of outputs for Molecular Oncology Almanac

All outputs will be produced by Molecular Oncology Almanac, though some may not contain content if the required input file was not passed.

# Table of contents
* [Common output columns used](#common-output-columns-used)
  * [Identifiers](#identifiers)
  * [Standardized feature columns](#standardized-feature-columns)
  * [Feature type specific columns](#feature-type-specific-columns)
    * [Somatic and germline variants](#somatic-and-germline-variants)
    * [Germline variants](#germline-variants)
    * [Fusions](#fusions)
    * [Somatic variants overlap](#somatic-variants-overlap)
  * [Datasource bins](#datasource-bins)
    * [Sorting somatic molecular features](#sorting-somatic-molecular-features)
    * [Clinical and biological relevance](#clinical-and-biological-relevance)
    * [Evidence of clinical assertions](#evidence-of-clinical-assertions)
      * [Associated evidence, Predictive Implication](#associated-evidence-predictive-implication)
      * [Therapeutic sensitivity](#therapeutic-sensitivity)
      * [Therapeutic resistance](#therapeutic-resistance)
      * [Disease prognosis](#disease-prognosis)
* [Produced outputs](#produced-outputs)
  * [Log](#log)
  * [Actionable](#actionable)
  * [Germline](#germline)
    * [American College of Medical Genetics](#american-college-of-medical-genetics)
    * [Somatic cancers](#somatic-cancers)
    * [Hereditary Cancers](#hereditary-cancers)
  * [Input metadata](#input-metadata)
  * [Integrated summary](#integrated-summary)
  * [Microsatellite Instability variants](#microsatellite-instability-variants)
  * [MOAlmanac execution json](#moalmanac-execution-json)
  * [Mutational burden](#mutational-burden)
  * [Preclinical efficacy](#preclinical-efficacy)
  * [Profile-to-cell line matchmaking](#profile-to-cell-line-matchmaking)
  * [Report](#report)
  * [Somatic molecular features](#somatic-molecular-features)
    * [Somatic filtered](#somatic-filtered)
    * [Somatic scored](#somatic-scored)
  * [Validation sequencing image](#validation-sequencing-image)

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
  * Mutational signatures: the specific COSMIC (v3.4) mutational signature, formatted as "COSMIC Signature (number)"
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
* `reference_allele` (str) - reference allele observed
* `observed_allele_1` (str) - observed allele observed, matches reference for single nucleotide variants.
* `observed_allele_2` (str) - observed allele observed, base that does not match reference for single nucleotide variants.
* `tumor_f` (float) - variant allelic fraction, calculated by dividing the count of a given base by the count of all bases present at a genomic site. 
* `total_coverage` (int) - total bases called at a genomic site
* `exac_af` (float) - Allele frequency of allele across all populations in ExAC
* `exac_common` (boolean) - Deemed common if the allele frequency in ExAC is greater than a minimum allele frequency specified in [config.ini](/moalmanac/config.ini). The default value is 0.001, or more common than 1 in 1000 alleles. 
* `clinvar` (str) - Clinical Significance of the genomic event from ClinVar

### Germline variants
In addition to the columns the population allele frequency (`exac_af`) and a boolean to designate a variant as common (`exac_common`), germline variants are thoroughly annotated with ExAC by subpopulation. The following columns are included for germline outputs:
* `exac_ac` (int) - allele counts of the observed allele, all subpopulations
* `exac_an` (int) - allele number of the all alleles at genomic location, all subpopulations
* `exac_afr_ac` (int) - allele counts of the observed allele, African / African American subpopulation
* `exac_afr_an` (int) - allele number of the all alleles at genomic location, African / African American subpopulation
* `exac_amr_ac` (int) - allele counts of the observed allele, Latino/Admixed American subpopulation
* `exac_amr_an` (int) - allele number of the all alleles at genomic location, Latino/Admixed American subpopulation
* `exac_eas_ac` (int) - allele counts of the observed allele, East Asian subpopulation
* `exac_eas_an` (int) - allele number of the all alleles at genomic location, East Asian subpopulation
* `exac_fin_ac` (int) - allele counts of the observed allele, European (Finish) subpopulation
* `exac_fin_an` (int) - allele number of the all alleles at genomic location, European (Finish) subpopulation
* `exac_nfe_ac` (int) - allele counts of the observed allele, European (non-Finish) subpopulation
* `exac_nfe_an` (int) - allele number of the all alleles at genomic location, European (non-Finish) subpopulation
* `exac_sas_ac` (int) - allele counts of the observed allele, South Asian subpopulation
* `exac_sas_an` (int) - allele number of the all alleles at genomic location, South Asian subpopulation
* `exac_oth_ac` (int) - allele counts of the observed allele, Other subpopulation
* `exac_oth_an` (int) - allele number of the all alleles at genomic location, Other subpopulation

### Fusions
The following columns are included to further describe fusions
* `spanningfrags` (int) - "the number of RNA-seq fragments that encompass the fusion junction such that one read of the pair aligns toa  different gene than the other paired-end read of that fragment", quoted from the [Star Fusion wiki](https://github.com/STAR-Fusion/STAR-Fusion/wiki)
* `left_gene` (str) - the gene name of the left gene composing a fusion
* `left_chr` (str) - the chromosome location of `left_gene`
* `left_position` (int) - the genomic location that the fusion is observed at in `left_gene`
* `right_gene` (str) - the gene name of the right gene composing a fusion
* `right_chr` (str) - the chromosome location of `right_gene`
* `right_position` (int) - the genomic location that the fusion is observed at in `right_gene`

### Somatic variants overlap
The following columns are included to inspect the relationship between somatic variants, germline variants, and validation sequencing:
* `number_germline_mutations_in_gene` (int) - the number of nonsynonymous variants in gene in provided germline MAF 
* `validation_total_coverage` (int) - total coverage at genomic location in validation sequencing
* `validation_tumor_f` (float) - allele frequency of observed base in validation sequencing
* `validation_detection_power` (float) - power calculation for likelihood of observing variant in validation sequencing, if it existed. See Methods from the publication for more information.

## Datasource bins
Molecular Oncology Almanac will match molecular features to [several datasources](/datasources) based on how closely a given molecular feature matches a catalogued feature. 

### Sorting somatic molecular features
The following datasources are used to sort somatic molecular features. Specifically, we assign a numeric value to each somatic variant, germline variants, copy number alteration, and fusion for each appropriate datasource. 

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

### Clinical and biological relevance
Each molecular feature will also receive a label in the `score_bin` column based on values for each datasource, listed above in ascending order. Molecular features which appear in the Molecular Oncology Almanac will receive labels describing the clinical actionability or biological relevance, depending on how closely the molecular feature matches a catalogued feature. Specifically, the following labels will be applied under the specified conditions:
* `Putatively Actionable`
  * Somatic and germline variants - Gene, variant classification, and protein change match a catalogued variant 
  * Copy number alterations - Gene and copy number direction match a catalogued event
  * Fusions - Both genes involved in a fusion event match a catalogued event
* `Investigate Actionability`
  * Somatic and germline variants - Gene and feature type match a catalogued variant but not variant classification or a specific protein change
  * Copy number alterations - Gene and feature type match a catalogued copy number alteration but not direction
  * Fusions - One gene fusion partner is catalogued as a fusion in Molecular Oncology Almanac but not both
* `Biologically Relevance`
  * The gene(s) associated with the molecular feature is present in Molecular Oncology Almanac but under a different feature type
  
The following second-order molecular features are evaluated in `score_bin` as follows:
* High mutational burden is labeled as `Investigate Actionability` 
* MSI-High is labeled as `Investigate Actionability`
* Whole-genome doubling is labeled as `Investigate Actionability`
* Mutational signatures catalogued by Molecular Oncology Almanac are labeled as `Investigate Actionability` and otherwise labeled as `Biologically Relevant`
* Variants associated with microsatellite instability are listed as "Supporting variants" as `Biologically Relevant`

### Evidence of clinical assertions
If a molecular feature matched as `Putatively Actionable`, or `Investigate Actionability` in Molecular Oncology Almanac, the molecular feature will be associated with clinical evidence. A molecular feature will be matched independently on catalogued events associated with therapeutic sensitivity, therapeutic resistance, and disease prognosis.

#### Associated evidence, Predictive Implication
All catalogued events in Molecular Oncology Almanac are cited and have associated evidence. Evidence tiers are as follows:

|Evidence|Description|
|--------|-----------|
|FDA-Approved|The Food and Drug Association (FDA) recognizes an association between the alteration and recommended clinical action.|
|Guideline|This relationship is catalogued as a guideline for standard of care treatment.|
|Clinical trial|The alteration is or has been used as an eligibility criterion for a clinical trial.|
|Clinical evidence|The relationship is reported in a clinical study that did not directly involve a clinical trial.|
|Preclinical evidence|This relationship is reported in a study involving mice, cell lines, or patient derived models.|
|Inferential evidence|The relationship is inferred as a result of mathematical modeling or an association between molecular features.|

#### Therapeutic sensitivity
Based on the score of a moleculear feature in `almanac_bin`, Molecular Oncology Almanac will attempt to match to catalogued assertions associated with therapeutic sensitivity. The following columns will be populated, if possible:
* `sensitivity_predictive_implication` (str) - [evidence tier of associated assertion](#associated-evidence-predictive-implication) 
* `sensitive_score_bin` (str) - the same as `score_bin`, only for catalogued molecular features associated with therapeutic sensitivity
* `sensitive_therapy_name` (str) - the name of the therapy associated with therapeutic sensitivity
* `sensitive_therapy_type` (str) - the name of the therapy type associated with therapeutic sensitivity
* `sensitive_description` (str) - the description catalogued in Molecular Oncology Almanac for the assertion
* `sensitive_citation` (str) - the citation of the evidence
* `sensitive_url` (str) - URL associated with the evidence

#### Therapeutic resistance
Based on the score of a moleculear feature in `almanac_bin`, Molecular Oncology Almanac will attempt to match to catalogued assertions associated with therapeutic resistance. The following columns will be populated, if possible:
* `resistance_predictive_implication` (str) - [evidence tier of associated assertion](#associated-evidence-predictive-implication)
* `resistance_score_bin` (str) - the same as `score_bin`, only for catalogued molecular features associated with therapeutic resistance
* `resistance_therapy_name` (str) - the name of the therapy associated with therapeutic resistance
* `resistance_therapy_type` (str) - the name of the therapy type associated with therapeutic resistance
* `resistance_description` (str) - the description catalogued in Molecular Oncology Almanac for the assertion
* `resistance_citation` (str) - the citation of the evidence
* `resistance_url` (str) - URL associated with the evidence

#### Disease prognosis
Based on the score of a moleculear feature in `almanac_bin`, Molecular Oncology Almanac will attempt to match to catalogued assertions associated with disease prognosis. The following columns will be populated, if possible:
* `prognostic_predictive_implication` (str) - [evidence tier of associated assertion](#associated-evidence-predictive-implication)
* `prognostic_score_bin` (str) - the same as `score_bin`, only for catalogued molecular features associated with disease prognosis
* `favorable_prognosis` (str) - boolean value for `1` being favorable prognosis and `0` corresponding to unfavorable prognosis
* `prognostic_description` (str) - the description catalogued in Molecular Oncology Almanac for the assertion
* `prognostic_citation` (str) - the citation of the evidence
* `prognostic_url` (str) - URL associated with the evidence

# Produced outputs
The following outputs are produced by the Molecular Oncology Almanac. Each section lists the filename suffix and then a details the contents of the output.

## Log
Filename suffix: `.log`

A timestamped log of inputs provided, configuration variables set, and what happens step-by-step as moalmanac.py is running.

## Actionable
Filename suffix: `.actionable.txt`

All molecular features associated with clinical or biological relevance will appear in this tab delimited output. Columns included are:
- [Standardized features](#standardized-feature-columns)
- [Sorting somatic molecular features](#sorting-somatic-molecular-features)
- [Clinical and biological relevance](#clinical-and-biological-relevance)
- [Evidence of clinical assertions](#evidence-of-clinical-assertions)
- [Somatic variants overlap](#somatic-variants-overlap)
- [Identifiers](#identifiers)

## Germline
Three germline-specific outputs are produced and populated if a germline MAF is passed (through `--germline_handle`) for a given molecular profile. These three outputs are based on genes present in the [American College of Medical Genetics v2](/datasources/acmg), cancer-related genes (appearing in [Molecular Oncology Almanac](/datasources/almanac), Cancer Hotspots [Molecular Oncology Almanac](/datasources/cancerhotspots), or [Cancer Gene Census](/datasources/cancergenecensus)), or [genes related to heritable cancers](/datasources/hereditary). 

### American College of Medical Genetics
Filename suffix: `.germline.acmg.txt`

Germline variants whose gene appears in the gene list from the [American College of Medical Genetics (v2)](/datasources/acmg) will appear in this output. Columns included are:
- [Standardized features](#standardized-feature-columns)
- [Somatic and germline variants](#somatic-and-germline-variants)
- [Germline variants](#germline-variants)
- [Identifiers](#identifiers)

### Somatic cancers
Filename suffix: `.germline.cancer_related.txt`

Germline variants whose gene appears in appearing in [Molecular Oncology Almanac](/datasources/almanac), [Cancer Hotspots](/datasources/cancerhotspots), or [Cancer Gene Census](/datasources/cancergenecensus) will appear in this output. Columns included are:
- [Standardized features](#standardized-feature-columns)
- [Somatic and germline variants](#somatic-and-germline-variants)
- [Germline variants](#germline-variants)
- [Sorting somatic molecular features](#sorting-somatic-molecular-features)
- [Clinical and biological relevance](#clinical-and-biological-relevance)
- [Identifiers](#identifiers)

### Hereditary cancers
Filename suffix: `.germline.hereditary_cancers.txt`

Germline variants whose gene appears in a curated list of genes [related to hereditary cancer risk](/datasources/hereditary) will appear in this output. Columns included are:
- [Standardized features](#standardized-feature-columns)
- [Somatic and germline variants](#somatic-and-germline-variants)
- [Germline variants](#germline-variants)
- [Identifiers](#identifiers)

## Input metadata
Filename suffix: `.input-metadata.txt`

Table representation of the patient dictionary, which contains input details such as:
- `patient_id` (str) - the string label for the case 
- `reported_tumor_type` (str) - the tumor type string, as passed to the main function
- `stage` (str) - a free text description of the stage of disease
- `description` (str) - a free text description of the case patient
- `purity` (float) - tumor purity, as passed to the main function
- `ploidy` (float) - tumor ploidy, as passed to the main function
- `WGD` (boolean) - boolean for if whole-genome doubling is present, as passed to the main function
- `microsatellite_status` (str) - categorical label for microsatellite status: unknown, microsatellite stable, microsatellite instability low, microsatellite instability high

## Integrated summary
Filename suffix: `.integrated.summary.txt`

All genes catalogued in Molecular Oncology Almanac are listed along with a list of alterations that appear in somatic variants, germline variants, copy number alterations, and fusions per gene. This is written as a tab delimited file.

## Microsatellite Instability variants
Filename suffix: `.msi_variants.txt`

Molecular features associated with Microsatellite Instability are displayed in this tab delimited output. Specifically, molecular features involving a gene from the list below are included:
- _DOCK3_
- _ESRP1_
- _MLH1_ 
- _MSH2_
- _MSH6_
- _PMS2_
- _POLE_
- _POLE2_
- _PRDM2_

Columns present in this file are:
- [Standardized features](#standardized-feature-columns)
- [Somatic and germline variants](#somatic-and-germline-variants)
- [Datasource columns](#datasource-bins)
- [Identifiers](#identifiers)

## MOAlmanac execution json
Filename suffix: `.moalmanac-execution.json`

A JSON output that contains all runtime details for a single execution of moalmanac's main function. The dictionary contains the following keys:
- `config` - all settings present in `config.ini` for the process execution
- `execution_runtime` - start datetime of moalmanac main function as well as the elapsed seconds for execution
- `input_files` - paths to input files of case profile provided to moalmanac
- `input_datasources` - paths to input datasources as provided in `annotation-databases.ini`
- `input_metadata` - input metadata, passed to the main function with the `patient` dictionary
- `actionable` - dictionary representation of `.actionable.txt` output
- `germline_acmg` - dictionary representation of `.germline.acmg.txt` output
- `germline_cancer` - dictionary representation of `.germline.cancer_related.txt` output
- `germline_hereditary` - dictionary representation of `.germline.hereditary_cancers.txt` output
- `integrated` - dictionary representation of `.integrated.summary.txt` output
- `msi_variants` - dictioanry representation of `.msi_variants.txt` output
- `somatic_filtered` - dictionary representation of `.somatic.filtered.txt` output
- `somatic_scored` - dictionary representation of `.somatic.scored.txt` output
- `therapeutic_strategies` - dictionary representation of `.therapeutic_strategies.txt` output
- `tumor_mutational_burden` - dictionary representation of the `.mutational_burden.txt` output

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
- `percentile_tcga` (float, 0-100) - the percentile of `coding_mutational_burden_per_megabase` relative to TCGA from [Lawrence et al. 2013](/datasources/lawrence)
- `percentile_tcga_tissuetype` (float, 0-100) - the percentile of `coding_mutational_burden_per_megabase` relative to a matched TCGA tissue type, if the tumor type matched
- `high_burden?` (boolean, True or False) - boolean value for evaluated mutational burden

Molecular Oncology Almanac designates high mutational burden under two circumstances: 
- Mutations per Mb > 10 
- At least a mutational burden of 80th percentile of TCGA tumor type, if matched, or TCGA generally, if not matched.

## Preclinical efficacy
Filename suffix: `.preclinical_efficacy.txt`

Therapies listed in [actionable](#actionable) that have been evaluated on cancer cell lines through the Sanger Institute's GDSC are evaluated for efficacy in the presence and absence of the associated molecular feature. This is performed for relationships associated with therapeutic sensitivity. Columns include:
- `patient_id` (str) - the string associated with the given molecular profile (`--patient_id`)
- `feature_display` (str) - a condensed string which captures what the molecular feature is
- `tested_subfeature` (str) - Molecular Oncology Almanac will evalaute efficacy for the therapy on combinations of all feature information: the gene, feature type, alteration type, and alteration. For more information, see the publication's Methods section.
- `therapy_name` (str) - the name of the therapy associated with therapeutic sensitivity
- `n_mut` (int) - Number of cancer cell lines mutant for `tested_subfeature`
- `n_wt` (int) - Number of cancer cell lines wild type for `tested_subfeature`
- `n_mut_tested` (int) - Number of cancer cell lines mutant for `tested_subfeature` that were tested with `therapy_name`
- `n_wt_tested` (int) - Number of cancer cell lines wild type for `tested_subfeature` that were tested with `therapy_name`
- `mut_ic50_median` (float) - Median IC50 of cell lines mutant for `tested_subfeature` that were tested with `therapy_name`
- `mut_ic50_mean` (float) - Mean IC50 of cell lines mutant for `tested_subfeature` that were tested with `therapy_name`
- `mut_ic50_std` (float) - IC50 Standard deviation of cell lines mutant for `tested_subfeature` that were tested with `therapy_name`
- `wt_ic50_median` (float) - Median IC50 of cell lines wild type for `tested_subfeature` that were tested with `therapy_name`
- `wt_ic50_mean` (float) - Mean IC50 of cell lines wild type for `tested_subfeature` that were tested with `therapy_name`
- `wt_ic50_std` (float) - IC50 Standard deviation of cell lines wild type for `tested_subfeature` that were tested with `therapy_name`
- `pvalue_mww` (float) - P-value of Mann-Whitney-Wilcoxon test between mutant vs wild type IC50 values
- `statistic` (float) - Test statistic of Mann-Whitney-Wilcoxon test between mutant vs wild type IC50 values

An image is also generated for each relationship evaluated, named "`patient_id`.`feature_display`.`therapy_name`.png"

## Profile-to-cell line matchmaking
Filename suffix: `.matchmaker.txt`

Molecular Oncology Almanac performs profile-to-cell line matchmaking on the provided molecular profile, comparing it to a population of cancer cell lines using the genomic profile. For more information, see the publication's Methods section. This file contains three columns:
- `case` (str) - case profile being compared, the string associated with the given molecular profile (`--patient_id`)
- `comparison` (str) - the comparison profile being compared to
- `SNF: FDA & CGC` (float) - results from Similarity Network Fusion in ascending order. Rows higher in this list are more similar to the case profile.

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

Molecular features which are biologically relevant are listed without clinical association. Molecular features will appear here if the associated gene is catalogued in the Molecular Oncology Almanac but under a different feature type, variants are associated with microsatellite stability, and all present COSMIC v3.4 mutational signatures not associated with a clinical assertion are reported.

The last section of the report, comparison of molecular profile to cancer cell lines, displays results from Molecular Oncology Almanac's patient-to-cell line matchmaking module. **This will not appear in the report if `--disable_matchmaking` is passed as an argument**. The 5 most similar cancer cell lines to the provided profile are listed each listing the cell line name, sensitive therapies from GDSC, and clinically relevant features present. Users can click `[More details]` under each cell line's name for more details about a given cell line: aliases, sensitive therapies, clinically relevant molecular features, all somatic variants, copy number alterations, and fusions occuring in cancer gene census genes, and the 10 most sensitive therapies to the cancer cell line.

## Somatic molecular features
Somatic variants, copy number alterations, and fusions will either be evaluated or filtered by Molecular Oncology Almanac. Criteria for evaluation is as follows:
- Somatic variants: nonsynonymous variants
- Copy number alterations: Within specified [percentiles](/moalmanac/config.ini), if using total copy number
- Fusions: At least the [specified minimum spanning fragments](/moalmanac/config.ini)

### Somatic filtered
Filename suffix: `.somatic.filtered.txt`

Somatic variants, copy number alterations, and fusions not evaluated by Molecular Oncology Almanac will appear in this tab delimited output. Columns included are: 
- [Standardized features](#standardized-feature-columns)
- [Feature type specific columns](#feature-type-specific-columns)
- [Identifiers](#identifiers)

### Somatic scored
Filename suffix: `.somatic.scored.txt`

Somatic variants, copy number alterations, and fusions evaluated by Molecular Oncology Almanac will appear in this tab delimited output. Columns included are: 
- [Standardized features](#standardized-feature-columns)
- [Feature type specific columns](#feature-type-specific-columns)
- [Sorting somatic molecular features](#sorting-somatic-molecular-features)
- [Clinical and biological relevance](#clinical-and-biological-relevance)
- [Somatic variants overlap](#somatic-variants-overlap)
- [Identifiers](#identifiers)

## Validation sequencing image
Filename suffix: `.validation_overlap.png`

If somatic variants from both primary and validation sequencing are provided, the allelic fraction of variants appearing in primary sequencing from both sequencing sources will be plotted on a scatter plot. Variants will be colored and annotated on the scatter plot based on [configurable thresholds](/moalmanac/config.ini). See the publication's Methods section for more details.
