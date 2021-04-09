#!/bin/bash

args=("$@")

export SEED=${args[0]} 
export WEB=${args[1]} 
export LANDSCAPE=${args[2]}
export REMSPP=${args[3]}
export NUTSUPPLY=${args[4]}
export DIR=${args[5]}
export OUTPUT_DIR=${args[6]}
export INVASION_DIR=${args[7]}
export TS_DIR=${args[8]}
export TS_FILE=${args[9]}
export JOINT_BASH=${args[10]}
export ISOLATED_BASH=${args[11]}
export TIMESERIES_BASH=${args[12]}
export TEND_BASH=${args[13]}
export TEVAL_BASH=${args[14]}

cat << EOF
seed: $SEED
web: $WEB
landscape: $LANDSCAPE
nutrient supply: $NUTSUPPLY
joint scenario: $JOINT_BASH
isolated scenario: $ISOLATED_BASH
invasive species: $REMSPP
tend: $TEND_BASH
teval: $TEVAL_BASH
directory: $DIR
output directory: $OUTPUT_DIR
invasion input directory: $INVASION_DIR
joint scenario: $JOINT_BASH
isolated scenario: $ISOLATED_BASH
print timeseries: $TIMESERIES_BASH
timeseries directory: $TS_DIR
timeseries file: $TS_FILE
EOF

INVASION_BASH=0
export INVASION_BASH
./simulation $SEED $WEB $LANDSCAPE $REMSPP $NUTSUPPLY $INVASION_BASH $DIR $OUTPUT_DIR $INVASION_DIR $TS_DIR $TS_FILE $JOINT_BASH $ISOLATED_BASH $TIMESERIES_BASH $TEND_BASH $TEVAL_BASH 

Rscript make_invasion_input.R $OUTPUT_DIR $INVASION_DIR

INVASION_BASH=1
export INVASION_BASH
./simulation $SEED $WEB $LANDSCAPE $REMSPP $NUTSUPPLY $INVASION_BASH $DIR $OUTPUT_DIR $INVASION_DIR $TS_DIR $TS_FILE $JOINT_BASH $ISOLATED_BASH $TIMESERIES_BASH $TEND_BASH $TEVAL_BASH 
 


