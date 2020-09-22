FROM vanallenlab/miniconda:3.6

WORKDIR /

RUN echo "deb http://cran.rstudio.com/bin/linux/ubuntu xenial/" | tee -a /etc/apt/sources.list \
    && gpg --keyserver keyserver.ubuntu.com --recv-key E298A3A825C0D65DFD57CBB651716619E084DAB9 \
    && gpg -a --export E084DAB9 | apt-key add - \
    && apt-get update && apt-get install -y \
    aptitude \
    libcurl4-openssl-dev \
    libxml2-dev \
    r-base \
    r-base-dev \
    r-bioc-biocgenerics

RUN Rscript -e 'install.packages("RCurl", repos = "http://cran.rstudio.com/")' \
    && Rscript -e 'source("http://bioconductor.org/biocLite.R"); biocLite("BSgenome"); biocLite("BSgenome.Hsapiens.UCSC.hg19"); biocLite("GenomeInfoDb")' \
    && Rscript -e 'install.packages("reshape2", repos = "http://cran.rstudio.com/")' \
    && Rscript -e 'install.packages("deconstructSigs", repos = "http://cran.rstudio.com/")'

RUN mkdir /base-image
COPY Dockerfile /base-image
COPY README.md /base-image
