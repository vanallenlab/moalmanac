# Molecular Oncology Almanac

Molecular Oncology Almanac (MOAlmanac) is a clinical interpretation algorithm paired with an underlying knowledge base for precision oncology. The primary objective of MOAlmanac is to identify and associate molecular alterations with therapeutic sensitivity and resistance as well as disease prognosis. This is done for “first-order” genomic alterations -- individual events such as somatic variants, copy number alterations, fusions, and germline -- as well as “second-order” events -- those that are not associated with one single mutation, and may be descriptive of global processes in the tumor such as tumor mutational burden, microsatellite instability, mutational signatures, and whole-genome doubling. In addition to clinical insights, MOAlmanac will annotate and evaluate first-order events based on their presence in numerous other well established datasources as well as highlight connections between them. This method currently geared towards hg19/b37 reference files and whole-exome or targeted sequencing data.

This repository is specifically for the clinical interpretation method and there are several other resources within the Molecular Oncology Almanac ecosystem: 
- Molecular Oncology Almanac, the clinical interpretation algorithm for precision cancer medicine [[Github](https://github.com/vanallenlab/moalmanac)].
- [Molecular Oncology Almanac Browser](https://moalmanac.org), a website to browse our underlying knowledge base [[Github](https://github.com/vanallenlab/moalmanac-browser)].
- [Molecular Oncology Almanac Connector](https://chrome.google.com/webstore/detail/molecular-oncology-almana/jliaipolchffpaccagodphgjpfdpcbcm?hl=en), a Google Chrome extension to quickly suggest literature for cataloging [[Github](https://github.com/vanallenlab/moalmanac-extension)].
- Molecular Oncology Almanac Database, the content and release notes of our underlying knowledge base [[Github](https://github.com/vanallenlab/moalmanac-db)].
- [Molecular Oncology Almanac Portal](https://portal.moalmanac.org), a website to launch this method on the Broad Institute's Google Cloud platform called Terra [[Github](https://github.com/vanallenlab/moalmanac-portal)].

This method is also available on [Docker](https://hub.docker.com/repository/docker/vanallenlab/moalmanac) and [Terra](https://portal.firecloud.org/#methods/vanallenlab/moalmanac/). 

If you use this method, please cite our publication:
> [Reardon, B., Moore, N.D., Moore, N.S., *et al*. Integrating molecular profiles into clinical frameworks through the Molecular Oncology Almanac to prospectively guide precision oncology. *Nat Cancer* (2021). https://doi.org/10.1038/s43018-021-00243-3](https://www.nature.com/articles/s43018-021-00243-3)
