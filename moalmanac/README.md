# Molecular Oncology Almanac 
Molecular Oncology Almanac can be run by executing either `moalmanac.py` with [standard input formats](#standard-usage) or `simplified_input.py` with [simplified inputs](#simplified-input). Please follow the [installation instructions](../README.md#installation) before use.

## Standard usage
Molecular Oncology Almanac may be executed on any combination of input data but does require a patient_id to label output files. Additional settings can be set by modifying the [config.ini](config.ini) file and column names may be modified by editing the [colnames.ini](colnames.ini) file.

Required arguments:
```
    --patient_id            <string>    patient identifier
```

Optional arguments:
```
    --tumor_type            <string>    tumor ontology, default=Unknown
    --stage                 <string>    tumor stage, default=Unknown
    --snv_handle            <string>    handle for MAF file of somatic single nucleotide variants
    --indel_handle          <string>    handle for MAF file of somatic insertions and deletions
    --bases_covered_handle  <string>    handle for text file which contains the number of calcable somatic bases
    --called_cn_handle      <string>    handle for text file which contained genes and copy number calls, will be used over `--cnv_handle`
    --cnv_handle            <string>    handle for annotated seg file for somatic copy number
    --fusion_handle         <string>    handle for STAR fusion output, .final.abridged
    --germline_handle       <string>    handle for MAF file of germline single nucleotide variants and insertions and deletions
    --validation_handle     <string>    handle for MAF file of somatic single nucleotide variant called from validation sequencing
    --ms_status             <string>    microsatellite status as deemed by MSI sensor, MSI or MSS, default=Unknown
    --purity                <float>     tumor purity
    --ploidy                <float>     tumor ploidy
    --wgd                   <boolean>   specify the occurence of whole genome duplication
    --disable_matchmaking   <boolean>   remove patient-to-cell line matchmaking from report
    --description           <string>    description of patient
    --output-directory      <string>    specify location of produced outputs
```

Example:
```
python moalmanac.py \
    --patient_id "example" \
    --tumor_type "SKCM" \
    --stage "Metastatic" \
    --description "Example profile for interpretation with the Molecular Oncology Almanac" \
    --snv_handle "../example_data/example_patient.capture.somatic.snvs.maf" \
    --indel_handle "../example_data/example_patient.capture.somatic.indels.maf" \
    --bases_covered_handle "../example_data/example_patient.capture.somatic.coverage.txt" \
    --called_cn_handle "../example_data/example_patient.capture.somatic.called.cna.txt" \
    --fusion_handle "../example_data/example_patient.capture.somatic.seg.annotated" \
    --germline_handle "../example_data/example_patient.rna.star.fusions.txt" \
    --validation_handle "../example_data/example_patient.rna.somatic.snvs.maf" \
    --purity 0.85 \
    --ploidy 4.02 \
    --ms_status "msih"
    --wgd
```

These example inputs may also be processed by executing `run_example.py`. 

Please also view our additional documentation on the [descriptions of inputs](../docs/description-of-inputs.md) for more details about input file formats. 

## Simplified input
A simplified input of a single file for somatic variants, germline variants, called copy number alterations, and fusions may also be used for a minimal interpretation. This mode also allows for [MSI status](../docs/description-of-inputs.md#microsatellite-status) and [whole-genome doubling](../docs/description-of-inputs.md#whole-genome-doubling) to be considered. With this format, MOAlmanac will be unable to annotate with any datasources that rely on nucleotide position. 

As with the [standard usage](#standard-usage), additional settings can be set by modifying the [config.ini](config.ini) file and column names may be modified by editing the [colnames.ini](colnames.ini) file.

Input for simplified input is a tab delimited file with one genomic alteration per row based on MOAlmanac's [standardized feature columns](../docs/description-of-outputs.md#standardized-feature-columns). In short the following columns are expected,
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

This is also described in the [description of inputs](../docs/description-of-inputs.md#simplified-input) document within the [docs/](../docs/) folder of this repository.

Required arguments:
```
    --patient_id            <string>    patient identifier
```

Optional arguments:
```
    --tumor_type            <string>    tumor ontology, default=Unknown
    --stage                 <string>    tumor stage, default=Unknown
    --input                 <string>    handle for simplified tab delimited input table
    --ms_status             <string>    microsatellite status as deemed by MSI sensor, MSI or MSS, default=Unknown
    --purity                <float>     tumor purity
    --ploidy                <float>     tumor ploidy
    --wgd                   <boolean>   specify the occurence of whole genome duplication
    --description           <string>    description of patient
    --output-directory      <string>    specify location of produced outputs
```

Example:
```
python simplified_input.py \
    --patient_id "example" \
    --tumor_type "SKCM" \
    --stage "Metastatic" \
    --description "Example profile for interpretation with the Molecular Oncology Almanac" \
    --input "../example_data/example_patient.capture.somatic.snvs.maf" \
    --indel_handle "../example_data/example_patient.capture.somatic.indels.maf" \
    --bases_covered_handle "../example_data/example_patient.capture.somatic.coverage.txt" \
    --called_cn_handle "../example_data/example_patient.capture.somatic.called.cna.txt" \
    --fusion_handle "../example_data/example_patient.capture.somatic.seg.annotated" \
    --germline_handle "../example_data/example_patient.rna.star.fusions.txt" \
    --validation_handle "../example_data/example_patient.rna.somatic.snvs.maf" \
    --purity 0.85 \
    --ploidy 4.02 \
    --ms_status "msih"
    --wgd
```

## Configuration
MOAlmanac can be customized by modifying the [config.ini](config.ini) file and column names may be modified by editing the [colnames.ini](colnames.ini) file.

### config.ini
The configuration file [config.ini](config.ini) lets users change settings, thresholds, and input string values. The file contains the following sections,
- `function_toggle` allows users to enable or disable the [actionability report](../docs/description-of-outputs.md#report), [model_similarity](../docs/description-of-outputs.md#profile-to-cell-line-matchmaking), [mutational signature](../docs/description-of-outputs.md#mutational-signatures), and [preclinical efficacy](../docs/description-of-outputs.md#preclinical-efficacy) functions. 
- `versions` are string inputs to describe the MOAlmanac algorithm and database versions
- `exac` allows users to specify a threshold for [Allele frequency in ExAC](https://gnomad.broadinstitute.org/help/faf) to identify common variants
- `fusion` allows users to specify the minimum spanning fragments required for fusions
- `mutation` allows users to specify minimum values for coverage and allelic fraction for evaluation of somatic and germline variants
- `seg` allows users to specify thresholds for total copy number 
- `signatures` allows users to specify the minimum required contribution to consider mutational signatures
- `validation_sequencing` allows users to specify minimum power and allelic fraction to consider for variants from validation sequencing
- `feature_types` allows users to specify strings for considered feature types
- `databases` specifies file paths for databases used for annotation, found in the `moalmanac/databases/` folder
- `preclinical` specifies file paths for datasources used for preclinical functions, [model_similarity](../docs/description-of-outputs.md#profile-to-cell-line-matchmaking) and [preclinical efficacy](../docs/description-of-outputs.md#preclinical-efficacy)

### colnames.ini
The configuration file [colnames.ini](colnames.ini) lets users change strings associated with column names for input and output files. The file contains the following relevant sections,
- `input_data`, allows users to change column names for input data

Other sections in this configuration file are used internally to MOAlmanac for processing.

## Citation
If you find this tool or any code herein useful, please cite:  
> [Reardon, B., Moore, N.D., Moore, N.S., *et al*. Integrating molecular profiles into clinical frameworks through the Molecular Oncology Almanac to prospectively guide precision oncology. *Nat Cancer* (2021). https://doi.org/10.1038/s43018-021-00243-3](https://www.nature.com/articles/s43018-021-00243-3)

## Disclaimer - For research use only
DIAGNOSTIC AND CLINICAL USE PROHIBITED. DANA-FARBER CANCER INSTITUTE (DFCI) and THE BROAD INSTITUTE (Broad) MAKE NO REPRESENTATIONS OR WARRANTIES OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING, WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NONINFRINGEMENT OR VALIDITY OF ANY INTELLECTUAL PROPERTY RIGHTS OR CLAIMS, WHETHER ISSUED OR PENDING, AND THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER OR NOT DISCOVERABLE.

In no event shall DFCI or Broad or their Trustees, Directors, Officers, Employees, Students, Affiliates, Core Faculty, Associate Faculty and Contractors, be liable for incidental, punitive, consequential or special damages, including economic damages or injury to persons or property or lost profits, regardless of whether the party was advised, had other reason to know or in fact knew of the possibility of the foregoing, regardless of fault, and regardless of legal theory or basis. You may not download or use any portion of this program for any non-research use not expressly authorized by DFCI or Broad. You further agree that the program shall not be used as the basis of a commercial product and that the program shall not be rewritten or otherwise adapted to circumvent the need for obtaining permission for use of the program other than as specified herein.
