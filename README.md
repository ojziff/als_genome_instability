## Meta-analysis of the ALS spectrum uncovers genome instability

Oliver Ziff 2022

This page contains scripts to analyse the data and reproduce the figures from the manuscript.

<img src="https://ojziff.github.io/ipsn_als_meta/figures/ipsn_meta_pipeline.png" height="400">

#### [Manuscript link](https://www.medrxiv.org/)

#### [Launch Shiny app](https://oliverziff.shinyapps.io/ipsn_als_meta/)


## Workbooks

[RNA-seq QC and motor neuron identities](https://ojziff.github.io/ipsn_als_meta/html/qc_identities.html) 

[ALS iPSN differential expression](https://ojziff.github.io/ipsn_als_meta/html/pan_als_expression.html) 

[Compare ALS genetic backgrounds](https://ojziff.github.io/ipsn_als_meta/html/compare_als_subgroups.html) 

[ALS postmortem differential expression](https://ojziff.github.io/ipsn_als_meta/html/postmortem_spinal_cord.html) 

[ALS iPSN differential splicing](https://ojziff.github.io/ipsn_als_meta/html/pan_als_splicing.html)

[Variants and gene fusions](https://ojziff.github.io/ipsn_als_meta/html/variants_fusions.html)

Rmarkdown scripts for each of these workbooks can be accessed [here](https://github.com/ojziff/ipsn_als_meta/blob/main/scripts/)


## Data availability

Raw RNA sequencing fastq files were accessed for each of the iPSN datasets in this table: 

| Reference             | Accession #           | Mutation              | ALS n | Control n | Library type | Paper URL                                    |
|-----------------------|-----------------------|-----------------------|-------|-----------|--------------|----------------------------------------------|
| Sareen et al, 2013    | GSE52202              | C9orf72               |     4 |         4 | polyA        | https://www.ncbi.nlm.nih.gov/pubmed/24154603 |
| Kiskinis et al, 2014  | GSE54409              | SOD1                  |     2 |         3 | polyA        | https://www.ncbi.nlm.nih.gov/pubmed/24704492 |
| Kapeli et al, 2016    | GSE77702              | FUS                   |     3 |         2 | polyA        | https://www.ncbi.nlm.nih.gov/pubmed/27378374 |
| Wang et al, 2017      | GSE95089              | SOD1                  |     2 |         2 | polyA        | https://pubmed.ncbi.nlm.nih.gov/28401346/    |
| De Santis et al, 2017 | GSE94888              | FUS                   |     3 |         3 | Ribo-zero    | https://www.ncbi.nlm.nih.gov/pubmed/28988989 |
| Bhinge et al, 2017    | PRJNA361408           | SOD1                  |     2 |         2 | Ribo-zero    | https://pubmed.ncbi.nlm.nih.gov/28366453/    |
| Luisier et al, 2018   | GSE98290              | VCP                   |     3 |         3 | polyA        | https://pubmed.ncbi.nlm.nih.gov/29789581/    |
| Abo-Rady et al, 2020  | GSE143743             | C9orf72               |     3 |         3 | polyA        | https://pubmed.ncbi.nlm.nih.gov/32084385/    |
| Dafinca et al. 2020   | GSE139144             | C9orf72               |     4 |         8 | polyA        | https://pubmed.ncbi.nlm.nih.gov/32330447/    |
|                       | GSE147544             | TARDBP                |     6 |         4 |              |                                              |
| Catanese et al., 2021 | GSE168831             | C9orf72               |     6 |         6 | polyA        | https://pubmed.ncbi.nlm.nih.gov/34125498/    |
|                       |                       | FUS                   |     6 |           |              |                                              |
| Smith et al, 2021     | PRJEB47567            | TARDBP                |     3 |         2 | polyA        | https://pubmed.ncbi.nlm.nih.gov/34660586/    |
| Hawkins et al, 2022   | GSE203168             | FUS                   |     2 |         2 | Ribo-zero    | https://pubmed.ncbi.nlm.nih.gov/35750046/                                     |
| Sommer et al, 2022    | GSE201407             | C9orf72               |     6 |         6 | polyA        | in press                                     |
| NeuroLINCS, 2022      | phs001231.v2.p1       | sporadic              |     8 |        14 | Ribo-zero    | https://pubmed.ncbi.nlm.nih.gov/34746695/    |
|                       |                       | SOD1                  |     6 |           |              |                                              |
|                       |                       | C9orf72               |    16 |           |              |                                              |
| Answer ALS, 2022      | AnswerALS data portal | sporadic              |   200 |        42 | Ribo-zero    | https://pubmed.ncbi.nlm.nih.gov/35115730/    |
|                       |                       | C9orf72               |    21 |           |              |                                              |
|                       |                       | SOD1                  |     8 |           |              |                                              |
|                       |                       | 6 other ALS mutations |     9 |           |              |                                              |

Raw fastq files were downloaded with [nf-core/fetchngs](https://nf-co.re/fetchngs) v1.7 and processed with [nf-core/rnaseq](https://nf-co.re/rnaseq) v3.8.1 utilising alignment with STAR and read quantification with salmon. 

Differential gene expression was performed using DESeq2 contrasting ALS versus control, adjusting for dataset and gender batch effects. 

Differential splicing was analysed using only polyA samples with MAJIQ v2.4 contrasting ALS versus control, adjusting for dataset and gender batch effects. 

Variant analysis was performed on Answer ALS samples only using [nf-core/rnavar](https://nf-co.re/rnavar) v1.0.0, which is based on GATK v4.2.6 short variant discovery workflow. Adjustment for read coverage was performed using a spline.

RNA fusion analysis was performed on paired-end datasets using [nf-core/rnafusion](https://nf-co.re/rnafusion) v2.0.0, utilising the STAR-Fusion workflow. Adjustment was made for dataset and read coverage (with a spline).

Schematics were created using [BioRender](https://biorender.com/).


