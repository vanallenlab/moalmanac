#!/bin/bash

tag=$1

sed -i "" -e "s#colnames.ini#/moalmanac/colnames.ini#" moalmanac/config.py
sed -i "" -e "s#../datasources#/datasources#" moalmanac/annotation-databases.ini
sed -i "" -e "s#../datasources#/datasources#" moalmanac/preclinical-databases.ini
sed -i "" -e "s#exac.expanded.r1.empty.txt#exac.expanded.r1.txt#" moalmanac/annotation-databases.ini

docker buildx build --platform linux/amd64,linux/arm64/v8 -t vanallenlab/moalmanac:"${tag}" --push .
docker buildx build --platform linux/amd64,linux/arm64/v8 -t vanallenlab/moalmanac:latest --push .

sed -i "" -e "s#/moalmanac/colnames.ini#colnames.ini#" moalmanac/config.py
sed -i "" -e "s#/datasources#../datasources#" moalmanac/annotation-databases.ini
sed -i "" -e "s#/datasources#../datasources#" moalmanac/preclinical-databases.ini
sed -i "" -e "s#exac.expanded.r1.txt#exac.expanded.r1.empty.txt#" moalmanac/annotation-databases.ini
