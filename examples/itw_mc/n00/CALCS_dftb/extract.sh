#!/bin/bash
for file in *.vasp.log ; do
 e=$(grep 'Total Energy:' $file | awk '{print $5}' | tail -n1)
 echo "$e # $file"
done
