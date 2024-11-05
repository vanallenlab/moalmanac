FROM vanallenlab/miniconda:3.11

WORKDIR /

RUN apt-get update && apt-get install -y

COPY requirements.txt /
RUN pip install --upgrade pip
RUN pip install -r requirements.txt

COPY example_data/ /example_data/
COPY example_output/ /example_output/

RUN mkdir /moalmanac/
RUN mkdir /moalmanac/datasources/
RUN mkdir /docs/

COPY moalmanac/test/ moalmanac/test/
COPY moalmanac/datasources/acmg/ /moalmanac/datasources/acmg/
COPY moalmanac/datasources/cancergenecensus/ /moalmanac/datasources/cancergenecensus/
COPY moalmanac/datasources/cancerhotspots/ /moalmanac/datasources/cancerhotspots/
COPY moalmanac/datasources/clinvar/ /moalmanac/datasources/clinvar/
COPY moalmanac/datasources/cosmic/ /moalmanac/datasources/cosmic/
COPY moalmanac/datasources/exac/ /moalmanac/datasources/exac/
COPY moalmanac/datasources/gsea_gene_sets/ /moalmanac/datasources/gsea_gene_sets/
COPY moalmanac/datasources/hereditary/ /moalmanac/datasources/hereditary/
COPY moalmanac/datasources/lawrence/ /moalmanac/datasources/lawrence/
COPY moalmanac/datasources/oncotree/ /moalmanac/datasources/oncotree/

COPY moalmanac/templates/ /moalmanac/templates/
COPY moalmanac/wrapper_deconstructsigs.sh moalmanac/run_deconstructsigs.R /moalmanac/
COPY moalmanac/*.py moalmanac/*.ini /moalmanac/

COPY moalmanac/datasources/moalmanac/ /moalmanac/datasources/moalmanac/

COPY moalmanac/datasources/preclinical/README.md /moalmanac/datasources/preclinical/README.md
COPY moalmanac/datasources/preclinical/generate-dictionary.ipynb /moalmanac/datasources/preclinical/generate-dictionary.ipynb
COPY moalmanac/datasources/preclinical/cell-lines.pkl /moalmanac/datasources/preclinical/cell-lines.pkl
COPY moalmanac/datasources/preclinical/formatted/almanac-gdsc-mappings.json /moalmanac/datasources/preclinical/formatted/almanac-gdsc-mappings.json
COPY moalmanac/datasources/preclinical/formatted/cell-lines.summary.txt /moalmanac/datasources/preclinical/formatted/cell-lines.summary.txt
COPY moalmanac/datasources/preclinical/annotated/cell-lines.somatic-variants.annotated.txt /moalmanac/datasources/preclinical/annotated/cell-lines.somatic-variants.annotated.txt
COPY moalmanac/datasources/preclinical/annotated/cell-lines.copy-numbers.annotated.txt  /moalmanac/datasources/preclinical/annotated/cell-lines.copy-numbers.annotated.txt
COPY moalmanac/datasources/preclinical/annotated/cell-lines.fusions.annotated.txt  /moalmanac/datasources/preclinical/annotated/cell-lines.fusions.annotated.txt
COPY moalmanac/datasources/preclinical/annotated/cell-lines.fusions.annotated.gene1.txt /moalmanac/datasources/preclinical/annotated/cell-lines.fusions.annotated.gene1.txt
COPY moalmanac/datasources/preclinical/annotated/cell-lines.fusions.annotated.gene2.txt /moalmanac/datasources/preclinical/annotated/cell-lines.fusions.annotated.gene2.txt
COPY moalmanac/datasources/preclinical/formatted/sanger.gdsc.txt /moalmanac/datasources/preclinical/formatted/sanger.gdsc.txt

COPY docs/* /docs/
COPY README.md /
COPY LICENSE /
COPY Dockerfile /
