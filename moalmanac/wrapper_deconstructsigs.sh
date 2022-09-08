#!/bin/bash

patient_id=$1
snv_handle=$2
sample=$3
ref=$4
alt=$5
chr=$6
pos=$7
folder=$8

Rscript run_deconstructsigs.R ${patient_id} ${snv_handle} ${sample} ${ref} ${alt} ${chr} ${pos} ${folder}
