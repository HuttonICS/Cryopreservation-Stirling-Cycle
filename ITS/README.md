# Cryopreservation project ITS

This project uses QIIME2 with DADA2 and the UNITE v.9.0 reference database to analyze ITS data. 

## Directory Structure

The main components of this project are:

- `qiime2/`: This directory contains the QIIME2 structure files generated based on the process from `ITS-QIIME2-DADA2-UNITE.md`. It includes the following files:
    - `taxonomy-its.qza`
    - `table-its.qza`
    - `stats-its.qza`
    - `rep-seqs-its.qza`
    - `rep-seqs-its-with-phyla-no-mitochondria-no-chloroplast.qza`
    - `table-its-with-phyla-no-mitochondria-no-chloroplast.qza`
    - `meta-data-ITS.txt`

- `figaro/`: This directory contains the outcomes from the figaro prediction.

- `ITS_figures/figure-markdown_strict/`: This directory contains figures for `16s-QIIME2-DADA2-Silva.md` visualisations.

- `ITS.Rmd`: This is an R Markdown file that contains the code for the analysis.

- `ITS-QIIME2-DADA2-UNITE.md`: This markdown file details the process of using QIIME2 with DADA2 and the UNITE v.9.0 reference database.

- `ITS-R.md`: This is the markdown output from `ITS.Rmd`, which can be viewed directly on GitHub.