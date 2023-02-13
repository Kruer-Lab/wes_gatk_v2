#!/bin/bash

source settings_snv.sh

echo -e "Mother's fathers's and proband's directory: (example: bash main_trio.sh M_F186-001-U M_F186-002-U M_F186-003-A)"

#echo "$1"
#echo "$2"
#echo "$3"

source part1.sh $1 &
source part1.sh $2 &
source part1.sh $3 &
#source part1.sh $4 &


wait

source part2.sh $3
