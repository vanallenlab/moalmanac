#!/bin/bash

tag=$1

sed -i "" -e "s#colnames.ini#/moalmanac/colnames.ini#" moalmanac/config.py
sed -i "" -e "s#../datasources#/datasources#" moalmanac/annotation-databases.ini
sed -i "" -e "s#../datasources#/datasources#" moalmanac/preclinical-databases.ini

docker build -t vanallenlab/moalmanac:"${tag}" .
docker push vanallenlab/moalmanac:"${tag}"

docker build -t vanallenlab/moalmanac:latest .
docker push vanallenlab/moalmanac:latest

sed -i "" -e "s#/moalmanac/colnames.ini#colnames.ini#" moalmanac/config.py
sed -i "" -e "s#/datasources#../datasources#" moalmanac/annotation-databases.ini
sed -i "" -e "s#/datasources#../datasources#" moalmanac/preclinical-databases.ini
