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
split_libraries_fastq.py -b bc_extracted/barcodes.fastq -o demult_reads_R2 --barcode_type 12 -m map.txt -i lane1_Undetermined.R2.fastq --store_demultiplexed_fastq --rev_comp_mapping_barcodes

#count reads in dataset
grep -c '@' /Volumes/Untitled/BeeAmpliconsForPat/demult_reads/seqs.fastq 
#7,652,267

#split fastq file by sample (for submission to NCBI later, if needed)
split_sequence_file_on_sample_ids.py --file_type fastq -o fastq_split_by_sample_R1 -i demult_reads/seqs.fastq

split_sequence_file_on_sample_ids.py --file_type fastq -o fastq_split_by_sample_R2 -i demult_reads_R2/seqs.fastq
```

## USEARCH analysis
```
#get latest SILVA
#wget https://www.arb-silva.de/fileadmin/silva_databases/qiime/Silva_132_release.zipwget https://www.arb-silva.de/fileadmin/silva_databases/qiime/Silva_132_release.zip
```

## Dereplicate sequences
```
mkdir mergedfastq
./usearch64 -fastx_uniques seqs_fixed.fastq -fastqout mergedfastq/uniques_combined_merged.fastq -sizeout
```

## Remove Singeltons
```
./usearch64 -sortbysize mergedfastq/uniques_combined_merged.fastq -fastqout mergedfastq/nosigs_uniques_combined_merged.fastq -minsize 2
```

## Precluster Sequences
```
./usearch64 -cluster_fast mergedfastq/nosigs_uniques_combined_merged.fastq -centroids_fastq mergedfastq/denoised_nosigs_uniques_combined_merged.fastq -id 0.9 -maxdiffs 5 -abskew 10 -sizein -sizeout -sort size
```

## Reference-based OTU picking against the Silva v. 1.28
```
#pull Silva and extract it
# wget https://www.arb-silva.de/fileadmin/silva_databases/qiime/Silva_128_release.tgz

./usearch64 -usearch_global mergedfastq/denoised_nosigs_uniques_combined_merged.fastq -id 0.97 -db ./Silva_128_release/SILVA_128_QIIME_release/rep_set/rep_set_16S_only/97/97_otus_16S.fasta  -strand plus -uc mergedfastq/ref_seqs.uc -dbmatched mergedfastq/closed_reference.fasta -notmatchedfq mergedfastq/failed_closed.fq
```

## Sort by size and then de novo OTU picking on sequences that failed to hit GreenGenes
```
./usearch64 -sortbysize mergedfastq/failed_closed.fq -fastaout mergedfastq/sorted_failed_closed.fq

./usearch64 -cluster_otus mergedfastq/sorted_failed_closed.fq -minsize 2 -otus mergedfastq/denovo_otus.fasta -relabel OTU_dn_ -uparseout mergedfastq/denovo_out.up
```

## Combine the rep sets between de novo and reference-based OTU picking
```
cat mergedfastq/closed_reference.fasta mergedfastq/denovo_otus.fasta > mergedfastq/full_rep_set.fna
```

## Map rep_set back to pre-dereplicated sequences and make OTU tables
```
./usearch64 -usearch_global seqs_fixed.fastq -db mergedfastq/full_rep_set.fna  -strand plus -id 0.97 -uc mergedfastq/OTU_map.uc -otutabout OTU_table.txt
```

## Assign taxonomy with SILVA 
```
source activate qiime1
assign_taxonomy.py -i mergedfastq/full_rep_set.fna -o taxonomy -r Silva_128_release/SILVA_128_QIIME_release/rep_set/rep_set_16S_only/97/97_otus_16S.fasta -t 'Silva_128_release/SILVA_128_QIIME_release/taxonomy/16S_only/97/consensus_taxonomy_all_levels.txt'

#add tax to OTU table
biom convert -i OTU_table.txt -o OTU_table.biom --table-type='OTU table' --to-hdf5

biom add-metadata -i OTU_table.biom -o otu_table_tax.biom --observation-metadata-fp=taxonomy/full_rep_set_tax_assignments.txt --sc-separated=taxonomy --observation-header=OTUID,taxonomy

#summarize OTU table
biom summarize-table -i otu_table_tax.biom -o otu_table_summmary.txt
```
