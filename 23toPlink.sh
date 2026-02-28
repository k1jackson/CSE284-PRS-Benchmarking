#!/bin/bash

input_dir="genomes"
output_dir="pink_data"

mkdir -p "$output_dir"

for file in "$input_dir"/*.txt; do
    id=$(basename "$file" .txt)
    echo $id
    plink1.9 --23file "$file" "$id" "$id" 0 -9 0 0 \
             --make-bed \
             --out "$output_dir/$id"
done

ls "$output_dir"/*.bed | sed 's/.bed//g' > "$output_dir".txt
plink1.9 --merge-list "$output_dir".txt --make-bed --out "$output_dir"

mkdir -p clean_"$output_dir"
while read filepath; do
    filename=$(basename "$filepath")
    plink1.9 --bfile "$filepath" --exclude "$output_dir"-merge.missnp --make-bed --out "clean_$output_dir/$filename"
done < "$output_dir".txt