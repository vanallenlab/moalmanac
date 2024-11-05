FROM vanallenlab/miniconda:3.11

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
COPY datasourcesacmg/ /datasourcesacmg/
COPY datasourcescancergenecensus/ /datasourcescancergenecensus/
COPY datasourcescancerhotspots/ /datasourcescancerhotspots/
COPY datasourcesclinvar/ /datasourcesclinvar/
COPY datasourcescosmic/ /datasourcescosmic/
COPY datasourcesexac/ /datasourcesexac/
COPY datasourcesgsea_gene_sets/ /datasourcesgsea_gene_sets/
COPY datasourceshereditary/ /datasourceshereditary/
COPY datasourceslawrence/ /datasourceslawrence/
COPY datasourcesoncotree/ /datasourcesoncotree/

COPY moalmanac/templates/ /moalmanac/templates/
COPY moalmanac/wrapper_deconstructsigs.sh moalmanac/run_deconstructsigs.R /moalmanac/
COPY moalmanac/*.py moalmanac/*.ini /moalmanac/

COPY datasourcesmoalmanac/ /datasourcesmoalmanac/

COPY datasourcespreclinical/README.md /datasourcespreclinical/README.md
COPY datasourcespreclinical/generate-dictionary.ipynb /datasourcespreclinical/generate-dictionary.ipynb
COPY datasourcespreclinical/cell-lines.pkl /datasourcespreclinical/cell-lines.pkl
COPY datasourcespreclinical/formatted/almanac-gdsc-mappings.json /datasourcespreclinical/formatted/almanac-gdsc-mappings.json
COPY datasourcespreclinical/formatted/cell-lines.summary.txt /datasourcespreclinical/formatted/cell-lines.summary.txt
COPY datasourcespreclinical/annotated/cell-lines.somatic-variants.annotated.txt /datasourcespreclinical/annotated/cell-lines.somatic-variants.annotated.txt
COPY datasourcespreclinical/annotated/cell-lines.copy-numbers.annotated.txt  /datasourcespreclinical/annotated/cell-lines.copy-numbers.annotated.txt
COPY datasourcespreclinical/annotated/cell-lines.fusions.annotated.txt  /datasourcespreclinical/annotated/cell-lines.fusions.annotated.txt
COPY datasourcespreclinical/annotated/cell-lines.fusions.annotated.gene1.txt /datasourcespreclinical/annotated/cell-lines.fusions.annotated.gene1.txt
COPY datasourcespreclinical/annotated/cell-lines.fusions.annotated.gene2.txt /datasourcespreclinical/annotated/cell-lines.fusions.annotated.gene2.txt
COPY datasourcespreclinical/formatted/sanger.gdsc.txt /datasourcespreclinical/formatted/sanger.gdsc.txt

COPY docs/* /docs/
COPY README.md /
COPY LICENSE /
COPY Dockerfile /
