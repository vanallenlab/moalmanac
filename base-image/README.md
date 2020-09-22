# Docker base image

The Molecular Oncology Almanac's docker also installs R in order to run [deconstructSigs](https://github.com/raerose01/deconstructSigs) as a subprocess. Here, we create a base image to be used so that we can rebuild the primary Docker without installing anything from ubuntu or for R.

## Installation
Build this docker with the tag base. 
```
docker build -t vanallenlab/almanac:base .
docker push vanallenlab/almanac:base
```
