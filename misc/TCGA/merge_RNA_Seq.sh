#!/bin/bash
OUTPUT="counts.tsv"
# first gunzip all gz files
for f in ./*/*.gz; do
    gunzip "$f"
done

# merge 2nd column of counts files
tmp=$(mktemp)
counts_files=(*/*.counts)
## keep both columns of first file
cp "${counts_files[1]}" "$OUTPUT"
echo "FILEID\t${counts_files[1]}" >> "$OUTPUT"

## keep 2nd column of other files
for f in "${counts_files[@]:1}"; do
    cat "$f" > "$tmp"
    echo "FILEID\t$f" >> "$tmp"
    paste "$OUTPUT" <(cut -d $'\t' -f2- "$tmp") > "_$OUTPUT" && mv "_$OUTPUT" "$OUTPUT"
done
tail -n 1 "$OUTPUT" > "$tmp" && head -n -1 "$OUTPUT" >> "$tmp" && mv "$tmp" "$OUTPUT"
sed -i 's/\.htseq\.counts//g' "$OUTPUT"
