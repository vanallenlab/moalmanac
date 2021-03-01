#!/bin/bash

python -m unittest discover --pattern "*_tests.py" -v

rm "example"*".png"
