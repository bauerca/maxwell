#!/bin/bash

files=$(ls ./*.{h,hpp,cpp})
find=MxOrdinal
replace=MxIndex

for f in $files
do
  # find/replace
  cat $f | sed "s/$find/$replace/g" > $f.tmp

  # line deletion
  #cat $f | sed "/$find/d" > $f.tmp

  # overwrite original
  mv $f.tmp $f
done
