#!/bin/bash

min_rsq=0
max_rsq=1
rsq_step=0.01

rm r_on_FAR_rsq.out

touch r_on_FAR_rsq.out

cur_rsq=$min_rsq
while (( $(echo "$cur_rsq <= $max_rsq" | bc -l) )); do
	echo "$cur_rsq"
	../grouping.py --wmax=10 --dmmin=100 --snmin=10 -s -r --rsqmmin=$cur_rsq ./split/* > ./r_on.out
	./calc_FAR.py r_on.out >> r_on_FAR_rsq.out
	cur_rsq=$(echo "$cur_rsq + $rsq_step" | bc)
done
