#!/bin/bash
for f in c*.vasp.log ; do x=$(grep 'Total Energy:' $f | tail -n1| awk '{print $5}'); echo "$x # $f" ; done > ENERGIES
