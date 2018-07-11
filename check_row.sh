#!/bin/sh

fin="f_nn_input.txt"
fout="f_nn_output.txt"

cat ${fin} | head -n 1 | sed -e "s/,/\r/g" | wc
cat ${fin} | head -n 1
cat ${fout} | head -n 1 | sed -e "s/,/\r/g" | wc
cat ${fout} | head -n 1 
