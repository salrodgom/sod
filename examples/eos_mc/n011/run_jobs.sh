#!/bin/bash
for i in *.vasp ; do 
 if [ ! -f $i.gout ] ; then
  echo "$i running"
  bash vasp2gin.sh $i ; gulp < $i.gin > $i.gout
 fi
done
