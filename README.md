
# Cryopreservation of a soil microbiome using a Stirling Cycle approach -a genomic assessment

### Introduction
Various institutions store collections of soil samples in either a refrigerated or sometimes frozen state. However, these conditions are not optimised to maintain the integrity of the soil microbiome. In this study, we explore cryopreservation using a controlled rate cooler and assess the genomic content of a representative soil sample before and after cryopreservation. We applied and compared two cryopreservation methods with control soil aliquots. Optimising the cryopreservation of soil samples is crucial for advancing microbiome research, as it helps preserve stable and functionally intact microbiomes.

The soil used for this project:
<img src="https://github.com/paytonyau/Cryopreservation-Stirling-Cycle/blob/main/soil.jpg" alt="soil texture" width="384" height="222"/>

![soil texture](https://github.com/paytonyau/Cryopreservation-Stirling-Cycle/blob/main/soil.jpg =100x40)

The study workflow:
```mermaid
graph LR
A[Soil] 
    A -->B[Control - no cryopreservation]
    A -->C[Plunge - Rapid cooler]
    A -->D[Rate - Stirling cycle cooler]
    B -->E[Recovery*]
    C -->E[Recovery*]
    D -->E[Recovery*]
    E -->F[DNA extraction]
    F -->G[Sequencing]
```
*The condition used for the recovery was 37°C for 5 minutes


### Table of contents
#### 16S Sequencing data
00-[(TXT) - Metadata](https://github.com/paytonyau/Cryopreservation-Stirling-Cycle/blob/main/16s/16s-meta-data.txt)

01-[(Markdown) Steps in proessing the ITS data](https://github.com/paytonyau/Cryopreservation-Stirling-Cycle/blob/main/16s-DADA2.md)

02-[(Rmarkdown) R data analysis](https://github.com/paytonyau/Cryopreservation-Stirling-Cycle/blob/main/16s/16s.Rmd)

02-[(Rmarkdown - pdf) R data analysis](https://github.com/paytonyau/Cryopreservation-Stirling-Cycle/blob/main/16s/16s.pdf)

03-[(Document-zip) QC in multiqc file format](https://github.com/paytonyau/Cryopreservation-Stirling-Cycle/blob/main/16s/multiqc/multiqc.zip)

04-[(Document-qiime2) qiime2 outcomes](https://github.com/paytonyau/Cryopreservation-Stirling-Cycle/blob/main/16s/qiime2)

05-[(Document-figaro) qiime2 outcomes](https://github.com/paytonyau/Cryopreservation-Stirling-Cycle/blob/main/16s/figaro)

#### ITS Sequencing data
00-[(TXT) - Metadata](https://github.com/paytonyau/Cryopreservation-Stirling-Cycle/blob/main/ITS/meta-data-ITS.txt)

01-[(Markdown) Steps in proessing the ITS data](https://github.com/paytonyau/Cryopreservation-Stirling-Cycle/blob/main/ITS-UNITE.md)

02-[(Rmarkdown) R data analysis](https://github.com/paytonyau/Cryopreservation-Stirling-Cycle/blob/main/ITS/ITS.Rmd)

02-[(Rmarkdown - pdf) R data analysis](https://github.com/paytonyau/Cryopreservation-Stirling-Cycle/blob/main/ITS/ITS.pdf)

03-[(Document-zip) QC in multiqc file format](https://github.com/paytonyau/Cryopreservation-Stirling-Cycle/blob/main/ITS/multiqc/multiqc.zip)

04-[(Document-qiime2) qiime2 outcomes](https://github.com/paytonyau/Cryopreservation-Stirling-Cycle/blob/main/ITS/qiime2)

05-[(Document-figaro) qiime2 outcomes](https://github.com/paytonyau/Cryopreservation-Stirling-Cycle/blob/main/ITS/figaro)

####  Information
The preprint is available at https://agrirxiv.org/search-details/?pan=20210277652

Raw Sequencing data deposited on https://www.ncbi.nlm.nih.gov/bioproject/PRJEB46478/

