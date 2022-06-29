# Cryopreservation of a soil microbiome using a Stirling Cycle approach -a genomic assessment

### Introduction
Collections of soil samples are kept by various institutions in either a refrigerated or occasionally frozen state, but conditions are not optimised to ensure the integrity of soil microbiome. In this study, we describe cryopreservation with a controlled rate cooler and estimate the genomic content of an exemplar soil sample before and after cryopreservation. Two methods of cryopreservation were applied and compared with control aliquots of soil. An optimised cryopreservation of soil samples is essential for the development of microbiome research in order to retain stable, functionally intact microbiomes.

### Table of contents
#### ITS Sequencing data
00-[(TXT) - Metadata](https://github.com/paytonyau/Cryopreservation-Stirling-Cycle/blob/main/ITS-meta-data.txt)

01-[(Shell Script) Training the Qiime2 classifier with-UNITE reference sequences](https://github.com/paytonyau/Cryopreservation-Stirling-Cycle/blob/main/ITS-UNITE-Rep-Seq.sh)

02-[(Shell Script) Download, QC, packing and trim the raw ITS sequencing data](https://github.com/paytonyau/Cryopreservation-Stirling-Cycle/blob/main/ITS-pre-processing.sh)

03-[(Document-zip) QC in multiqc file format](https://github.com/paytonyau/Cryopreservation-Stirling-Cycle/blob/main/ITS-multiqc.zip)

04-[(Document-qzv) Sequence quality - before cutadapt](https://github.com/paytonyau/Cryopreservation-Stirling-Cycle/blob/main/ITS-demux.qzv)

05-[(Document-qzv) Sequence quality - after cutadapt](https://github.com/paytonyau/Cryopreservation-Stirling-Cycle/blob/main/ITS-trimmed-demux.qzv)

06-[(Shell Script) ITS sequencing process using Qiime2 - DADA2](https://github.com/paytonyau/Cryopreservation-Stirling-Cycle/blob/main/ITS-DADA2.sh)

#### 16S Sequencing data
00 - [(TXT) - Metadata](https://github.com/paytonyau/Cryopreservation-Stirling-Cycle/blob/main/16s-meta-data.txt)

01-[(Shell Script) Download, QC and trim the raw 16s sequencing data](https://github.com/paytonyau/Cryopreservation-Stirling-Cycle/blob/main/16S-pre-processing.sh)

02-[(Document-zip) QC in multiqc file format](https://github.com/paytonyau/Cryopreservation-Stirling-Cycle/blob/main/16s-multiqc.zip)

03-[(Document-qzv) Sequence quality](https://github.com/paytonyau/Cryopreservation-Stirling-Cycle/blob/main/16s-demux.qzv)

04-[(Shell Script) 16S sequencing process using Qiime2 - DADA2](https://github.com/paytonyau/Cryopreservation-Stirling-Cycle/blob/main/16s-DADA2.sh)


####  Information
The preprint is available at https://agrirxiv.org/search-details/?pan=20210277652

Raw Sequencing data deposited on https://www.ncbi.nlm.nih.gov/bioproject/PRJEB46478/

The information in the .qzv files can be visualised via: https://view.qiime2.org/

