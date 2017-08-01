#!/usr/bin/env bash

#usage: ./parsePQR.sh <PQR file name>

if [ $# -ne 1 ]
then
    echo "Usage: ./parsePQR.sh <PQR file name>"
    exit -1
fi

grep 'ATOM\|HETATM' $1 | cut -c 30- > temp
wc -l < temp > $1.ext 
cat temp >> $1.ext
rm temp
