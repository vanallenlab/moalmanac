FROM vanallenlab/almanac:base

WORKDIR /

COPY requirements.txt /
RUN pip install -r requirements.txt

COPY example_data/ /example_data/
COPY example_output/ /example_output/

RUN mkdir /moalmanac/
RUN mkdir /moalmanac/datasources/
RUN mkdir /docs/

COPY moalmanac/test/ moalmanac/test/
COPY moalmanac/datasources/acmg/ /moalmanac/datasources/acmg/
COPY moalmanac/datasources/moalmanac/ /moalmanac/datasources/moalmanac/
COPY moalmanac/datasources/cancergenecensus/ /moalmanac/datasources/cancergenecensus/
COPY moalmanac/datasources/cancerhotspots/ /moalmanac/datasources/cancerhotspots/
COPY moalmanac/datasources/clinvar/ /moalmanac/datasources/clinvar/
COPY moalmanac/datasources/cosmic/ /moalmanac/datasources/cosmic/
COPY moalmanac/datasources/exac/ /moalmanac/datasources/exac/
COPY moalmanac/datasources/gsea_gene_sets/ /moalmanac/datasources/gsea_gene_sets/
COPY moalmanac/datasources/hereditary/ /moalmanac/datasources/hereditary/
COPY moalmanac/datasources/lawrence/ /moalmanac/datasources/lawrence/
COPY moalmanac/datasources/oncotree/ /moalmanac/datasources/oncotree/

COPY moalmanac/datasources/preclinical/almanac-gdsc-mappings.json /moalmanac/datasources/preclinical/almanac-gdsc-mappings.json
COPY moalmanac/datasources/preclinical/cell-lines.summary.txt /moalmanac/datasources/preclinical/cell-lines.summary.txt
COPY moalmanac/datasources/preclinical/ccle.variants.evaluated.txt /moalmanac/datasources/preclinical/ccle.variants.evaluated.txt
COPY moalmanac/datasources/preclinical/ccle.copy-numbers.evaluated.txt /moalmanac/datasources/preclinical/ccle.copy-numbers.evaluated.txt
COPY moalmanac/datasources/preclinical/sanger.fusions.evaluated.txt /moalmanac/datasources/preclinical/sanger.fusions.evaluated.txt
COPY moalmanac/datasources/preclinical/sanger.fusions.gene1.evaluated.txt /moalmanac/datasources/preclinical/sanger.fusions.gene1.evaluated.txt
COPY moalmanac/datasources/preclinical/sanger.fusions.gene2.evaluated.txt /moalmanac/datasources/preclinical/sanger.fusions.gene2.evaluated.txt
COPY moalmanac/datasources/preclinical/sanger.gdsc.txt /moalmanac/datasources/preclinical/sanger.gdsc.txt
COPY moalmanac/datasources/preclinical/cell-lines.pkl /moalmanac/datasources/preclinical/cell-lines.pkl

COPY moalmanac/templates/ /moalmanac/templates/
COPY moalmanac/wrapper_deconstructsigs.sh moalmanac/run_deconstructsigs.R /moalmanac/
COPY moalmanac/*.py moalmanac/*.ini /moalmanac/

COPY docs/* /docs/
COPY README.md /
COPY LICENSE /
COPY Dockerfile /
