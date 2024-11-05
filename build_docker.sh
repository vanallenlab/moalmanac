#!/bin/bash

tag=$1

sed -i "" -e "s#config.ini#/moalmanac/config.ini#" moalmanac/config.py
sed -i "" -e "s#colnames.ini#/moalmanac/colnames.ini#" moalmanac/config.py
sed -i "" -e "s#datasources/#/datasources/#" moalmanac/config.ini

docker build -t vanallenlab/moalmanac:"${tag}" .
docker push vanallenlab/moalmanac:"${tag}"

docker build -t vanallenlab/moalmanac:latest .
docker push vanallenlab/moalmanac:latest

sed -i "" -e "s#/moalmanac/config.ini#config.ini#" moalmanac/config.py
sed -i "" -e "s#/moalmanac/colnames.ini#colnames.ini#" moalmanac/config.py
sed -i "" -e "s#/datasources/#datasources/#" moalmanac/config.ini
