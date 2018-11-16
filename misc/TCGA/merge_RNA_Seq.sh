#!/bin/bash
# Make sure to put this in the directory with all the downloaded directories from TCGA
ANNOTATION="/mnt/storage/yi/data/TCGA/annotation/fpkm_annot.csv"
OUTPUT="FPKM.tsv"
#########################################
# first gunzip all gz files
for f in ./*/*.gz; do
    gunzip "$f"
done

# get unique projects
awk -F, 'NR > 1 {print $3}' "$ANNOTATION" | sort | uniq > projects.out
while read -u 10 proj # for each project, -u 10 is for preventing reading from stdin in loop
do
    # get file names in project, variables used in awk have to be defined
    mapfile -t proj_files < <(awk -F, -v proj="$proj" '$3 == proj {print $1 "/" $2}' "$ANNOTATION" | sed 's/\.gz$//g')

    tmp=$(mktemp)
    # keep both columns in the first file
    cp "${proj_files[1]}" "$proj.$OUTPUT"
    echo -e "FILEID\t${proj_files[1]}" >> "$proj.$OUTPUT"
    # keep second column of other files
    for f in "${proj_files[@]:2}"; do
        cp "$f" "$tmp"
        echo -e "FILEID\t$f" >> "$tmp"
        paste "$proj.$OUTPUT" <(cut -d $'\t' -f2- "$tmp") > "$OUTPUT" && mv "$OUTPUT" "$proj.$OUTPUT"
    done
    tail -n 1 "$proj.$OUTPUT" > "$tmp" && head -n -1 "$proj.$OUTPUT" >> "$tmp" && mv "$tmp" "$proj.$OUTPUT"
    sed -i 's/\.htseq\.counts//g' "$proj.$OUTPUT"
done 10< projects.out && rm projects.out

