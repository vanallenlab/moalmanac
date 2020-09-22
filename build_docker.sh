#!/bin/bash

tag=$1

sed -i "" -e "s#config.ini#/moalmanac/config.ini#" moalmanac/config.py
sed -i "" -e "s#colnames.ini#/moalmanac/colnames.ini#" moalmanac/config.py
sed -i "" -e "s#datasources/#/moalmanac/datasources/#" moalmanac/config.ini
sed -i "" -e "s#wrapper_deconstructsigs.sh#/moalmanac/wrapper_deconstructsigs.sh#" moalmanac/features.py
sed -i "" -e "s#run_deconstructsigs.R#/moalmanac/run_deconstructsigs.R#" moalmanac/wrapper_deconstructsigs.sh

docker build -t vanallenlab/moalmanac:"${tag}" .

sed -i "" -e "s#/moalmanac/config.ini#config.ini#" moalmanac/config.py
sed -i "" -e "s#/moalmanac/colnames.ini#colnames.ini#" moalmanac/config.py
sed -i "" -e "s#/moalmanac/datasources/#datasources/#" moalmanac/config.ini
sed -i "" -e "s#/moalmanac/wrapper_deconstructsigs.sh#wrapper_deconstructsigs.sh#" moalmanac/features.py
sed -i "" -e "s#/moalmanac/run_deconstructsigs.R#run_deconstructsigs.R#" moalmanac/wrapper_deconstructsigs.sh

docker push vanallenlab/moalmanac:"${tag}"
