#!/bin/bash

## demoise and merging ##
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs demux.qza \
  --p-trunc-len-f 260 \
  --p-trunc-len-r 220 \
  --p-trim-left-f 9 \
  --p-trim-left-r 9 \
  --p-n-threads 8 \
  --o-representative-sequences dada2-222-silva-138-rep-seqs.qza \
  --o-table dada2-222-silva-138-table.qza \
  --o-denoising-stats dada2-222-silva-138-stats.qza

## Silva 138 Reference database ##
qiime feature-classifier classify-sklearn \
  --i-classifier /mnt/shared/scratch/pyau/ref/silva-138-99-nb-classifier.qza \
  --i-reads dada2-222-silva-138-rep-seqs.qza \
  --p-n-jobs 1 \
  --o-classification dada2-222-silva-138-taxonomy.qza

qiime metadata tabulate \
  --m-input-file dada2-222-silva-138-taxonomy.qza \
  --o-visualization dada2-222-silva-138-taxonomy.qzv

## Taxonomy-based filtering of tables and sequences ## 
qiime taxa filter-table \
  --i-table dada2-222-silva-138-table.qza \
  --i-taxonomy dada2-222-silva-138-taxonomy.qza \
  --p-include p__ \
  --p-exclude mitochondria,chloroplast \
  --o-filtered-table dada2-222-silva-138-table-with-phyla-no-mitochondria-no-chloroplast.qza
  
qiime taxa filter-seqs \
  --i-sequences dada2-222-silva-138-rep-seqs.qza \
  --i-taxonomy dada2-222-silva-138-taxonomy.qza \
  --p-include p__ \
  --p-exclude mitochondria,chloroplast \
  --o-filtered-sequences dada2-222-silva-138-rep-seqs-with-phyla-no-mitochondria-no-chloroplast.qza
