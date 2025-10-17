FROM vanallenlab/miniconda:3.12

WORKDIR /

RUN apt-get update && apt-get install -y

COPY requirements.txt /
RUN pip install --upgrade pip
RUN pip install -r requirements.txt

COPY example_data/ /example_data/
COPY example_output/ /example_output/

RUN mkdir /moalmanac/
RUN mkdir /datasources
RUN mkdir /docs/

COPY moalmanac/test/ moalmanac/test/

COPY datasources/acmg/ /datasources/acmg/
COPY datasources/cancergenecensus/ /datasources/cancergenecensus/
COPY datasources/cancerhotspots/ /datasources/cancerhotspots/
COPY datasources/clinvar/ /datasources/clinvar/
COPY datasources/cosmic/CosmicMutantExport_empty.lite.txt /datasources/cosmic/
COPY datasources/cosmic/prepare_cosmic.py /datasources/cosmic/
COPY datasources/cosmic/README.md /datasources/cosmic/
COPY datasources/exac/ /datasources/exac/
COPY datasources/gsea_gene_sets/ /datasources/gsea_gene_sets/
COPY datasources/hereditary/ /datasources/hereditary/
COPY datasources/lawrence/ /datasources/lawrence/
COPY datasources/oncotree/ /datasources/oncotree/

COPY moalmanac/templates/ /moalmanac/templates/
COPY moalmanac/*.py moalmanac/*.ini /moalmanac/

COPY datasources/moalmanac/ /datasources/moalmanac/

COPY datasources/preclinical/README.md /datasources/preclinical/README.md
COPY datasources/preclinical/generate-dictionary.ipynb /datasources/preclinical/generate-dictionary.ipynb
COPY datasources/preclinical/cell-lines.pkl /datasources/preclinical/cell-lines.pkl
COPY datasources/preclinical/formatted/almanac-gdsc-mappings.json /datasources/preclinical/formatted/almanac-gdsc-mappings.json
COPY datasources/preclinical/formatted/cell-lines.summary.txt /datasources/preclinical/formatted/cell-lines.summary.txt
COPY datasources/preclinical/annotated/cell-lines.somatic-variants.annotated.txt /datasources/preclinical/annotated/cell-lines.somatic-variants.annotated.txt
COPY datasources/preclinical/annotated/cell-lines.copy-numbers.annotated.txt  /datasources/preclinical/annotated/cell-lines.copy-numbers.annotated.txt
COPY datasources/preclinical/annotated/cell-lines.fusions.annotated.txt  /datasources/preclinical/annotated/cell-lines.fusions.annotated.txt
COPY datasources/preclinical/annotated/cell-lines.fusions.annotated.gene1.txt /datasources/preclinical/annotated/cell-lines.fusions.annotated.gene1.txt
COPY datasources/preclinical/annotated/cell-lines.fusions.annotated.gene2.txt /datasources/preclinical/annotated/cell-lines.fusions.annotated.gene2.txt
COPY datasources/preclinical/formatted/sanger.gdsc.txt /datasources/preclinical/formatted/sanger.gdsc.txt

COPY docs/* /docs/
COPY README.md /
COPY LICENSE /
COPY Dockerfile /
