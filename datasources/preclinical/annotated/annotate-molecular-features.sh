#!/bin/bash

INPUT_SOMATIC=${1-../formatted/cell-lines.somatic-variants.txt}
INPUT_COPY_NUMBER=${2-../formatted/cell-lines.copy-numbers.txt}
INPUT_FUSION=${3-../formatted/cell-lines.fusions.txt}
OUTPUT_SOMATIC=${4-cell-lines.somatic-variants.annotated.txt}
OUTPUT_COPY_NUMBER=${5-cell-lines.copy-numbers.annotated.txt}
OUTPUT_FUSION=${6-cell-lines.fusions.annotated.txt}
WD=$PWD
MOALMANAC_DIR=../../..

cp annotate-variants.py annotate-copy-numbers.py annotate-fusions.py "$MOALMANAC_DIR"
cp "$INPUT_SOMATIC" "$MOALMANAC_DIR"/cell-lines-somatic-variants.txt
cp "$INPUT_COPY_NUMBER" "$MOALMANAC_DIR"/cell-lines-copy-number-alterations.txt
cp "$INPUT_FUSION" "$MOALMANAC_DIR"/cell-lines-fusions.txt
cd "$MOALMANAC_DIR" || exit

python annotate-variants.py --input cell-lines-somatic-variants.txt --output "$OUTPUT_SOMATIC" --directory "$WD"
python annotate-copy-numbers.py --input cell-lines-copy-number-alterations.txt --output "$OUTPUT_COPY_NUMBER" --directory "$WD"
python annotate-fusions.py --input cell-lines-fusions.txt --output "$OUTPUT_FUSION" --directory "$WD"

rm annotate-variants.py annotate-copy-numbers.py annotate-fusions.py
rm cell-lines-somatic-variants.txt cell-lines-copy-number-alterations.txt cell-lines-fusions.txt
