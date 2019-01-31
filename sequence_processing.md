## Processing of 16S rRNA gene data
### Extracting barcodes
```
#QIIME v. 1.9.1

#validate mapping file
validate_mapping_file.py -m map.txt -o map_validate

#extract barcodes
extract_barcodes.py -f lane1_Undetermined.R1.fastq -r lane1_Undetermined.R2.fastq -l 12 -o bc_extacted


```


