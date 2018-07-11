#!/bin/bash

f="f_nn_output.txt"

if [ $1 -eq 0 ]; then
	data=`cat ./zero`
else
	data=`head -n $1 ${f} | tail -n 1`
fi
echo -e "${data}" | nc -u -w0 192.168.11.4 12347
echo -e "${data}"
