## Processing of 16S rRNA gene data
### Extracting barcodes
```
#count reads in raw fastq files
grep -c '@' /Volumes/Untitled/BeeAmpliconsForPat/lane1_Undetermined.R1.fastq 
#303,222,273 reads

#QIIME v. 1.9.1
#module load qiime1

#validate mapping file
validate_mapping_file.py -m map.txt -o map_validate

#extract barcodes
extract_barcodes.py -f lane1_Undetermined.R1.fastq -l 12 -o bc_extracted -s "0:" -c barcode_in_label

#join paired ends
join_paired_ends.py -f lane1_Undetermined.R1.fastq -r lane1_Undetermined.R2.fastq -b bc_extracted/barcodes.fastq -o joined_reads

#count reads that could be joined
grep -c '@' /fastqjoin.join.fastq
#134,356,347

#demultiplex the reads
split_libraries_fastq.py -i joined_reads/fastqjoin.join.fastq -o demult_reads -b joined_reads/fastqjoin.join_barcodes.fastq --store_demultiplexed_fastq -m map.txt --barcode_type 12

#just R1
split_libraries_fastq.py -b bc_extracted/barcodes.fastq -o demult_reads --barcode_type 12 -m map.txt -i lane1_Undetermined.R1.fastq --store_demultiplexed_fastq

#split fastq file by sample (for submission to NCBI later, if needed)
split_sequence_file_on_sample_ids.py --file_type fastq -o fastq_split_by_sample -i demult_reads/
```


for i in *.fastq
do
grep -c "@" i > i_count.txt
done
