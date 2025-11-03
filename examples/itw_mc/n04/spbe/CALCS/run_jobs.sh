#!/bin/bash
for i in *.vasp ; do bash vasp2gin.sh $i ; gulp < $i.gin > $i.gout ; done
