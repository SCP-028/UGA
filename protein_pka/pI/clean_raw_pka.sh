#!/bin/bash
# The header should be "residue | Position | pK | slope | 1000Chi2"
mkdir -p ./cleaned
for f in ./KAUST/*.pka; do
    awk '!/range$/ && !/^Residue/ {print FILENAME,$1,$3}' $f |  # dropped header
    awk '{$1=$1}1' OFS="," |  # change delimiter to comma
    sed 's/\.\/KAUST\/\(.\{4\}\).pka\(.*\)$/\1\2/g' > ./cleaned/${f:8}  # path to PDB ID
done
