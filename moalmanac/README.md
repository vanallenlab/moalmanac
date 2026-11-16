# Molecular Oncology Almanac
Molecular Oncology Almanac can be run by executing `moalmanac.py` with [standard input formats](#standard-usage). Please follow the [installation instructions](../README.md#installation) before use.

## Standard usage
Molecular Oncology Almanac may be executed on any combination of input data but does require a patient_id to label output files. Additional settings can be set by modifying the [config.ini](#configini) file and column names may be modified by editing the [colnames.ini](#colnamesini) file. The [config.ini](config.ini) and [annotation-databases.ini](annotation-databases.ini) must be passed as arguments to moalmanac.py.

Required arguments:
```
    --patient_id            <string>    patient identifier
    --config                <string>    ini file that contains configuration details, config.ini
    --dbs                   <string>    ini file that contains database paths for annotation, annotation-databases.ini
```

Optional arguments:
```
    --tumor_type            <string>    tumor ontology, default=Unknown, MOAlmanac will attempt to map to both OncoTree codes and terms.
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
    --tmb                   <float>     tumor mutational burden
    --wgd                   <boolean>   specify the occurence of whole genome duplication
    --disable_matchmaking   <boolean>   remove patient-to-cell line matchmaking from report
    --description           <string>    description of patient
    --output-directory      <string>    specify location of produced outputs
    --preclinical-dbs       <string>    path to preclinical-databases.ini file
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
    --ms_status "msih" \
    --wgd \
    --config config.ini \
    --dbs annotation-databases.ini \
    --preclinical-dbs preclinical-databases.ini
```

These example inputs may also be processed by executing `run_example.py`.

Please also view our additional documentation on the [descriptions of inputs](../docs/description-of-inputs.md) for more details about input file formats.

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

### colnames.ini
The configuration file [colnames.ini](colnames.ini) lets users change strings associated with column names for input and output files. The file contains the following relevant sections,
- `input_data`, allows users to change column names for input data

Other sections in this configuration file are used internally to MOAlmanac for processing.

### annotation-databases.ini
The configuration file [annotation-databases.ini](annotation-databases.ini) lets users change paths to datasources being used to annotate genomic variants within the algorithm. This file contains a single `databases` section, with either a relative or absolute path being set to the `root` variable. By default, this points to the [datasources/](../datasources/) folder in the root directory of this repository.

### preclinical-databases.ini
Similar to `annotation-databases.ini`, the configuration file [preclinical-databases.ini](preclinical-databases.ini) lets users change paths to datasources being used to for preclinical comparison functions,  [model_similarity](../docs/description-of-outputs.md#profile-to-cell-line-matchmaking) and [preclinical efficacy](../docs/description-of-outputs.md#preclinical-efficacy).

## Citation
If you find this tool or any code herein useful, please cite:
> [Reardon, B., Moore, N.D., Moore, N.S., *et al*. Integrating molecular profiles into clinical frameworks through the Molecular Oncology Almanac to prospectively guide precision oncology. *Nat Cancer* (2021). https://doi.org/10.1038/s43018-021-00243-3](https://www.nature.com/articles/s43018-021-00243-3)

## Disclaimer

For purposes of this Disclaimer, the “Service” means the Molecular Oncology Almanac database, software, API, browsers, and related materials.

THIS DATABASE AND SOFTWARE HAVE NOT BEEN VALIDATED FOR DIAGNOSTIC OR CLINICAL USE, AND DANA-FARBER CANCER INSTITUTE (DFCI) MAKES NO WARRANTY OF FITNESS FOR DIAGNOSTIC OR CLINICAL USE. USERS ARE SOLELY RESPONSIBLE FOR ENSURING THAT THEIR ACCESS TO AND USE OF THE DATABASE, SOFTWARE, API, BROWSERS, AND RELATED MATERIALS COMPLIES WITH ALL APPLICABLE LAWS, REGULATIONS, INSTITUTIONAL POLICIES, AND ETHICAL REQUIREMENTS, INCLUDING ANY REQUIREMENTS APPLICABLE TO RESEARCH INVOLVING HUMAN SUBJECTS OR HUMAN-DERIVED DATA. DFCI MAKES NO REPRESENTATIONS OR WARRANTIES OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING, WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NONINFRINGEMENT OR VALIDITY OF ANY INTELLECTUAL PROPERTY RIGHTS OR CLAIMS, WHETHER ISSUED OR PENDING, AND THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER OR NOT DISCOVERABLE.

In no event shall DFCI or its trustees, directors, officers, employees, students, affiliates, core faculty, associate faculty and contractors be liable for incidental, punitive, consequential or special damages, including economic damages or injury to persons or property or lost profits, regardless of whether the party was advised, had other reason to know or in fact knew of the possibility of the foregoing, regardless of fault, and regardless of legal theory or basis. To the extent liability cannot be disclaimed, DFCI’s liability shall be limited to direct, actual damages resulting from DFCI’s failure to perform its obligations under the applicable terms.

The Molecular Oncology Almanac database is made available under the [Creative Commons Attribution 4.0 International License (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/), which permits use, including commercial use, provided appropriate credit is given. The software is licensed under the [GNU General Public License, Version 2 (GPL v2)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html). Certain third-party code or materials may be subject to separate license terms. For additional terms governing access to and use of the Service, please see the [Terms of Use](https://docs.moalmanac.org/about/terms-of-use/).
