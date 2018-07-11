#!/bin/bash

index_start=0
if [ $# -ge 1 ]; then
	index_start=$1
fi

f="f_nn_output.txt"

i=0
skip_start=0
for a in `cat ${f}`
do
	i=`expr ${i} + 1`
	skip_start=`expr ${skip_start} + 1`
	if [ ${skip_start} -le ${index_start} ]; then
		continue
	fi

	#if [ ${i} -ge 30 ]; then
	if [ 1 ]; then
		echo -e "${a}\n" | nc -u -w0 192.168.11.2 12346
		echo -e "${a}\n"
		echo "${i}"
		#i=0
		./usleep 100
	fi
done
