## Processing of 16S rRNA gene data
### Extracting barcodes
```
#count reads in raw fastq files
grep -c '@' /Volumes/Untitled/BeeAmpliconsForPat/lane1_Undetermined.R1.fastq 
#303,222,273 reads

#QIIME v. 1.9.1
#source activate qiime1

#validate mapping file
validate_mapping_file.py -m map.txt -o map_validate

#extract barcodes
extract_barcodes.py -f lane1_Undetermined.R1.fastq -l 12 -o bc_extracted -s "0:" -c barcode_in_label

#join paired ends
join_paired_ends.py -f lane1_Undetermined.R1.fastq -r lane1_Undetermined.R2.fastq -b bc_extracted/barcodes.fastq -o joined_reads

#count reads that could be joined
grep -c '@' joined_reads/fastqjoin.join.fastq
#134,356,347 reads

#demultiplex the reads
split_libraries_fastq.py -i joined_reads/fastqjoin.join.fastq -o demult_reads_paired -b joined_reads/fastqjoin.join_barcodes.fastq --store_demultiplexed_fastq -m map.txt --barcode_type 12 --rev_comp_mapping_barcodes

#just R1
split_libraries_fastq.py -b bc_extracted/barcodes.fastq -o demult_reads --barcode_type 12 -m map.txt -i lane1_Undetermined.R1.fastq --store_demultiplexed_fastq --rev_comp_mapping_barcodes

#just R2
split_libraries_fastq.py -b bc_extracted/barcodes.fastq -o demult_reads_R2 --barcode_type 12 -m map.txt -i lane1_Undetermined.R1.fastq --store_demultiplexed_fastq --rev_comp_mapping_barcodes

#count reads in dataset
grep -c '@' /Volumes/Untitled/BeeAmpliconsForPat/demult_reads/seqs.fastq 
#7652267

#split fastq file by sample (for submission to NCBI later, if needed)
split_sequence_file_on_sample_ids.py --file_type fastq -o fastq_split_by_sample_R1 -i demult_reads/seqs.fastq

split_sequence_file_on_sample_ids.py --file_type fastq -o fastq_split_by_sample_R2 -i demult_reads_R2/seqs.fastq
```

## USEARCH analysis
```
#get latest SILVA
wget https://www.arb-silva.de/fileadmin/silva_databases/qiime/Silva_132_release.zipwget https://www.arb-silva.de/fileadmin/silva_databases/qiime/Silva_132_release.zip

#check sample names
./usearch32 -fastx_get_sample_names demult_reads/seqs.fastq -output name_check.txt
#00:00 26Mb    100.0% 53,854 samples found

#fix sample names
sed 's/_/\./' demult_reads/seqs.fastq > demult_reads/seqs_fixed.fastq

#check fixed sample names
./usearch32 -fastx_get_sample_names demult_reads/seqs_fixed.fastq -output sample_name2.txt
#100.0% 78 samples found
```

```
#dereplicate sequences
./usearch32 -fastx_uniques test.fastq -fastqout test_derep.fq -sizeout
00:01 44Mb    100.0% Reading test.fastq
00:01 26Mb    100.0% DF                
00:01 27Mb   53854 seqs, 6268 uniques, 4759 singletons (75.9%)
00:01 27Mb   Min size 1, median 1, max 9705, avg 8.59
00:01 27Mb    100.0% Writing test_derep.fq

#remove singletons
./usearch32 -sortbysize test_derep.fq -fastqout test_nosingle.fq -minsize 2

#cluster against Silva 1.3.2
./usearch32 -usearch_global test_nosingle.fq -id 0.97 -strand plus -uc ref_seqs.uc -dbmatched closed_seqs.fna -notmatchedfq failed_seqs.fq -db SILVA_132_QIIME_release/rep_set/rep_set_16S_only/97/silva_132_97_16S.fna
#00:35 1.2Gb   100.0% Searching, 82.2% matched

#de novo OTU pick (UCHIIME chimera check)
./usearch32 -sortbysize failed_seqs.fq -fastqout failed_sort.fq

./usearch32 -cluster_otus failed_sort.fq -minsize 2 -otus denovo_otus.fasta -relabel OTU_dn_ -uparseout denovo_out.up
#00:00 5.0Mb   100.0% 31 OTUs, 41 chimeras

#catenate data
cat denovo_otus.fasta closed_seqs.fna > closed_denovo_seqs.fna

#map reads back to original fastq file
./usearch32 -usearch_global test.fastq -db closed_denovo_seqs.fna -strand plus -id 0.97 -uc otu_map.uc -otutabout OTU_table.txt
# 100.0% Searching, 98.0% matched
```

## QIIME analyses
```
#assign taxonomy
assign_taxonomy.py -i closed_denovo_seqs.fna -o taxonomy -r /Volumes/Untitled/BeeAmpliconsForPat/SILVA_132_QIIME_release/rep_set/rep_set_16S_only/97/silva_132_97_16S.fna -t /Volumes/Untitled/BeeAmpliconsForPat/SILVA_132_QIIME_release/taxonomy/16S_only/97/consensus_taxonomy_7_levels.txt
```
