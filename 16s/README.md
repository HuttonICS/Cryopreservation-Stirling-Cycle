# Cryopreservation Project 16s

This project uses QIIME2 with DADA2 and the Silva version 138 reference database to analyze 16s amplicon sequencing data. 

## Directory Structure

The main components of this project are:

- `qiime2/`: This directory contains the QIIME2 structure files generated based on the process from `16s-QIIME2-DADA2-Silva.md`. It includes the following files:
    - `430_327_213_rep-seqs.qza`
    - `430_327_213_rep-seqs-with-phyla-no-mitochondria-no-chloroplast.qza`
    - `430_327_213_stats.qza`
    - `430_327_213_table.qza`
    - `430_327_213_table-with-phyla-no-mitochondria-no-chloroplast.qza`
    - `430_327_213_taxonomy.qza`
    - `16s-meta-data.txt`

The numbering here means:
    - 430 nt long (maximum) for the sequence length
    - 237 nt - trimmed position for the forward sequences 
    - 213 nt - trimmed position for the reverse sequences

- `figaro/`: This directory contains the outcomes from the figaro prediction.

- `16s_figures/figure-markdown_strict/`: This directory contains figures for `16s-QIIME2-DADA2-Silva.md` visualisations.

- `16s.Rmd`: This is an R Markdown file that contains the code for the downstream analysis.

- `16s-QIIME2-DADA2-Silva.md`: This markdown file details the process of using QIIME2 with DADA2 and the UNITE v.9.0 reference database.

- `16s-R.md`: This is the markdown output from `16s.Rmd`, which can be viewed directly on GitHub.