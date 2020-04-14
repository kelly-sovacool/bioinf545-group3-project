#!/bin/bash
infile=$1
temp=$(mktemp)
sed 's/[)('\'']//g' ${infile} > ${temp}
mv ${temp} ${infile}


