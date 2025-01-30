# Annotated inputs
Molecular features formatted as described in the `formatted/` folder are to be annotated by [MOAlmanac](https://github.com/vanallenlab/moalmanac) before use. The inputs requirements are needed to run this step are,
- A formatted input for somatic variants, as described in the `formatted/` folder
- A formatted input for copy number alterations, as described in the `formatted/` folder
- A formatted input for fusions, as described in the `formatted/` folder

The `annotate-molecular-features.sh` script utilize these five requirements to annotate molecular features. These scripts have been modified from those in the [paper's Github repository](https://github.com/vanallenlab/moalmanac-paper/tree/main/analyses/preclinical), and furthermore subset to only include molecular features whose gene appears in either MOAlmanac, Cancer Hotspots, or Cancer Gene Census.

## Usage
Required arguments:
```
    INPUT_SOMATIC       <string>    Path to formatted somatic variants input
    INPUT_COPY_NUMBER   <string>    Path to formatted copy number alterations input
    INPUT_FUSION        <string>    Path to formatted fusions input
```

Optional arguments:
```
    OUTPUT_SOMATIC       <string>    Path to annotated somatic variants output, default: samples.variants.annotated.txt
    OUTPUT_COPY_NUMBER   <string>    Path to formatted copy number alterations output, default: samples.copy_numbers.annotated.txt
    OUTPUT_FUSION        <string>    Path to formatted fusions output, default: samples.fusions.annotated.txt
```

Example:
```bash
conda activate moalmanac
bash annotate-molecular-features.sh ../formatted/cell-lines.somatic-variants.txt \
                                    ../formatted/cell-lines.copy-numbers.txt \
                                    ../formatted/cell-lines.fusions.txt
```

This script will produce five output files, one each for somatic variants and copy number alterations and three for fusions. Output names for somatic variants, copy number alterations, and fusions can optionally be provided as arguments 4, 5, 6. 
```bash
conda activate moalmanac
bash annotate-molecular-features.sh ../formatted/samples.somatic_variants.txt \
                                    ../formatted/samples.copy_numbers.txt \
                                    ../formatted/samples.fusions.txt \
                                    example-samples.somatic-variants.txt \
                                    example-samples.copy-number-alterations.txt \
                                    example-samples.fusions.txt
```

Here, we provide detailed information on each output produced. 

[Return to top](#annotated-inputs)

## Annotated somatic variants
Somatic variants will be annotated with MOAlmanac, producing an output named by default as `cell-lines.somatic-variants.annotated.txt`. Details of the datasource columns (e.g. `cgc_bin`) can be found in the MOAlmanac documentation [description of outputs](https://github.com/vanallenlab/moalmanac/blob/main/docs/description-of-outputs.md#sorting-somatic-molecular-features).

### Fields
- `feature`, gene symbol associated with the variant
- `alteration_type`, variant consequence or classification resulting from the nucleotide change
- `alteration`, protein change associated with the variant using the one-letter amino-acid codes
- `sample_name`, string associated with each sample, molecular profile, or patient
- `feature_match_1`, 0 if the molecular feature's gene is not present in the database and 1 if it is
- `feature_match_2`, 0 if the molecular feature's gene is not catalogued under Somatic Variants and 1 if it is
- `feature_match_3`, 0 if the molecular feature's gene and alteration type is not catalogued under Somatic Variants and 1 if it is
- `feature_match_4`, 0 of the molecular feature's gene, alteration, and alteration is catalogued under Somatic Variants and 1 if it is
- `evidence`, [associated evidence](https://github.com/vanallenlab/moalmanac/blob/main/docs/description-of-outputs.md#associated-evidence-predictive-implication) with the molecular feature

[Return to top](#annotated-inputs)

## Annotated copy number alterations
Copy number alterations will be annotated with MOAlmanac, producing an output named by default as `cell-lines.copy-numbers.annotated.txt`. Details of the datasource columns (e.g. `cgc_bin`) can be found in the MOAlmanac documentation [description of outputs](https://github.com/vanallenlab/moalmanac/blob/main/docs/description-of-outputs.md#sorting-somatic-molecular-features).

### Fields
- `feature`, gene symbol associated with the copy number alteration
- `alteration_type`, direction of copy number alteration
- `alteration`, left blank
- `sample_name`, string associated with each sample, molecular profile, or patient
- `feature_match_1`, 0 if the molecular feature's gene is not present in the database and 1 if it is
- `feature_match_2`, 0 if the molecular feature's gene is not catalogued under Copy Number and 1 if it is
- `feature_match_3`, 0 if the molecular feature's gene and alteration type is not catalogued under Copy Number and 1 if it is
- `feature_match_4`, 0 for all rows
- `evidence`, [associated evidence](https://github.com/vanallenlab/moalmanac/blob/main/docs/description-of-outputs.md#associated-evidence-predictive-implication) with the molecular feature

[Return to top](#annotated-inputs)

## Annotated fusions, gene 1
Fusions are processed twice, once using provided feature partner values for each row and then with the two swapped. For each row, the "strongest" annotations for each row is also output. Annotating fusions based on the default configuration will produce an output named by default as `cell-lines.fusions.annotated.gene1.txt`. Details of the datasource columns (e.g. `cgc_bin`) can be found in the MOAlmanac documentation [description of outputs](https://github.com/vanallenlab/moalmanac/blob/main/docs/description-of-outputs.md#sorting-somatic-molecular-features).

### Fields
- `feature`, gene symbol associated with the fusion
- `partner`, gene symbol associated with the fusion's partner gene
- `sample_name`, string associated with each sample, molecular profile, or patient
- `feature_match_1`, 1 if feature is presented in MOAlmanac's gene list and 0 otherwise
- `feature_match_2`, 1 if feature is a gene associated with a Rearrangement within MOAlmanac and 0 otherwise
- `feature_match_3`, 1 if feature is a gene associated with a Rearrangement and Fusion within MOAlmanac and 0 otherwise
- `feature_match_4`, 1 if feature and the fusion partner are both catalogued as a fusion pair within MOAlmanac and 0 otherwise
- `feature_match`, valued 1-4 and reflects the sum of other feature_match columns, used for faster compute in matchmaking
- `evidence`, [associated evidence](https://github.com/vanallenlab/moalmanac/blob/main/docs/description-of-outputs.md#associated-evidence-predictive-implication) with the molecular feature

[Return to top](#annotated-inputs)

## Annotated fusions, gene 2
Fusions are processed twice, once using provided feature partner values for each row and then with the two swapped. This iteration annotates the fusions after swapping the `feature` and `partner` values. For each row, the "strongest" annotations for each row is also output. Annotating fusions based on the default configuration will produce an output named by default as `cell-lines.fusions.annotated.gene2.txt`. Details of the datasource columns (e.g. `cgc_bin`) can be found in the MOAlmanac documentation [description of outputs](https://github.com/vanallenlab/moalmanac/blob/main/docs/description-of-outputs.md#sorting-somatic-molecular-features).

### Fields
- `feature`, gene symbol associated with the fusion (`partner` in the gene1 output)
- `partner`, gene symbol associated with the fusion's partner gene (`feature` in the gene1 output)
- `sample_name`, string associated with each sample, molecular profile, or patient
- `feature_match_1`, 1 if feature is presented in MOAlmanac's gene list and 0 otherwise
- `feature_match_2`, 1 if feature is a gene associated with a Rearrangement within MOAlmanac and 0 otherwise
- `feature_match_3`, 1 if feature is a gene associated with a Rearrangement and Fusion within MOAlmanac and 0 otherwise
- `feature_match_4`, 1 if feature and the fusion partner are both catalogued as a fusion pair within MOAlmanac and 0 otherwise
- `feature_match`, valued 1-4 and reflects the sum of other feature_match columns, used for faster compute in matchmaking
- `evidence`, [associated evidence](https://github.com/vanallenlab/moalmanac/blob/main/docs/description-of-outputs.md#associated-evidence-predictive-implication) with the molecular feature

[Return to top](#annotated-inputs)

## Annotated fusions
Fusions are processed twice, once using provided feature partner values for each row and then with the two swapped. For each row, the "strongest" annotations for each row is produced, creating this output named by default as `cell-lines.fusions.annotated.txt`. Details of the datasource columns (e.g. `cgc_bin`) can be found in the MOAlmanac documentation [description of outputs](https://github.com/vanallenlab/moalmanac/blob/main/docs/description-of-outputs.md#sorting-somatic-molecular-features).

### Fields
- `feature`, gene symbol associated with the fusion 
- `partner`, gene symbol associated with the fusion's partner gene 
- `which_match`, which gene the match values are associated with 
- `sample_name`, string associated with each sample, molecular profile, or patient
- `feature_match_1`, 1 if feature is presented in MOAlmanac's gene list and 0 otherwise
- `feature_match_2`, 1 if feature is a gene associated with a Rearrangement within MOAlmanac and 0 otherwise
- `feature_match_3`, 1 if feature is a gene associated with a Rearrangement and Fusion within MOAlmanac and 0 otherwise
- `feature_match_4`, 1 if feature and the fusion partner are both catalogued as a fusion pair within MOAlmanac and 0 otherwise
- `feature_match`, valued 1-4 and reflects the sum of other feature_match columns, used for faster compute in matchmaking
- `evidence`, [associated evidence](https://github.com/vanallenlab/moalmanac/blob/main/docs/description-of-outputs.md#associated-evidence-predictive-implication) with the molecular feature

[Return to top](#annotated-inputs)
