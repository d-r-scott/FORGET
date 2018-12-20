#!/bin/bash

min_snmin=10
max_snmin=12
snmin_step=0.1

rm simple_filter_FAR.out
rm r_off_FAR.out
rm r_on_FAR.out

touch simple_filter_FAR.out
touch r_off_FAR.out
touch r_on_FAR.out

cur_snmin=$min_snmin
while (( $(echo "$cur_snmin <= $max_snmin" | bc -l) )); do
	echo "$cur_snmin"
	../simple_filter.py --wmax=10 --dmmin=100 --snmin=$cur_snmin -s ./split/* > ./simple_filter.out
	../grouping.py --wmax=10 --dmmin=100 --snmin=12 -s ./split/* > ./r_off.out
	../grouping.py --wmax=10 --dmmin=100 --snmin=$cur_snmin -s -r --rsqmmin=0.5 ./split/* > ./r_on.out
	./calc_FAR.py simple_filter.out >> simple_filter_FAR.out
	./calc_FAR.py r_off.out >> r_off_FAR.out
	./calc_FAR.py r_on.out >> r_on_FAR.out
	cur_snmin=$(echo "$cur_snmin + $snmin_step" | bc)
done
