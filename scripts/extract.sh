#!/bin/bash
for f in *.vasp.gout ; do x=$(grep 'l e' $f | tail -n1| awk '{print $4}'); echo "$x # $f" ; done > ENERGIES
