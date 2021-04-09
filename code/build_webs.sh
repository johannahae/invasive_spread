#!/bin/bash

WEB_DIR="/homes/jh57masa/Desktop/invasive_spread/data/webs/"	

for WEB in 41 88 7 8 42 # 41 => web 1; 88 => web 2; 7 => web 3; 8 => web 4; 42 => web 5
do

for S in 21
do 

for SB in 5
do 

SC=$(($S-$SB))
SEED=$((54321+$WEB))

export SEED
export SB
export SC
export WEB
export WEB_DIR

./webs $WEB $SEED $SB $SC $WEB_DIR

done 
done 
done
