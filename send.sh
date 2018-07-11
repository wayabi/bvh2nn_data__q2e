#!/bin/sh

while [ 1 ]
do
	cat $1 | nc -u -w1 192.168.11.3 12346
	sleep 2
done
