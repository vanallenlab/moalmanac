# Databases: Molecular Oncology Almanac
The Molecular Oncology Almanac uses its corresponding database to identify clinical relevant molecular features that may confer information regarding therapeutic sensitivity, resistance, or prognostic value. 

## About the Molecular Oncology Almanac
The Molecular Oncology Almanac attempts to capture the current body of knowledge on how genetic alterations affect clinical actionability. As the field of precision medicine grows, the quantity of research on how specific alterations affect therapeutic response and patient prognosis expands at an increasing rate. The Molecular Oncology Almanac seeks to curate this information from the literature, increasing the abilities of clinicians and researchers alike to rapidly assess the importance of an alteration. 

Several other services exist within the Molecular Oncology Almanac ecosystem. See [this repository's docs folder](/docs/) for more information.

## Usage: Formatting the database for use
This method uses a document-based format of the database, which is built using the [database repository](https://github.com/vanallenlab/moalmanac-db) and `create_almanac_db.py`. 

Arguments:
```
    --directory, -d     <string>    path to the moalmanac-db repository's content folder
    --version, -v       <string>    release name for the database content being used, should match the release from moalmanac-db
```

This should be run with this repository's virtual environment enabled. 
