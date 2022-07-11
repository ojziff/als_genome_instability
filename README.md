## Meta-analysis of the ALS spectrum uncovers genome instability

Oliver Ziff 2022

### [Manuscript](https://www.medrxiv.org/)

### [Launch Shiny app](https://oliverziff.shinyapps.io/ipsn_als_meta/)

This page contains scripts to analyse the data and reproduce the figures from the meta-analysis of ALS iPSC-derived motor neurons.

<img src="https://github.com/ojziff/ipsn_als_meta/blob/main/figures/ipsn_meta_pipeline.png" height="400">

## [Rmarkdown workbooks](https://github.com/ojziff/ipsn_als_meta/blob/main/ipsc_mn_als_meta.Rmd)

Knitted HTMLs can be browsed with the links below.

### [RNA-seq QC and motor neuron identities](https://oliverziff.github.io/ipsn_als_meta/html/qc_identities.html) 

<img src="https://github.com/ojziff/ipsn_als_meta/blob/main/figures/ipsn_identities_pca_heatmap.png" height="700">

### [pan-ALS vs control differential expression](https://oliverziff.github.io/ipsn_als_meta/html/ipsn_als_expression.html) 

![Figure 2](https://github.com/ojziff/ipsn_als_meta/blob/main/figure_3.png)

### [Compare ALS genetic backgrounds](https://oliverziff.github.io/ipsn_als_meta/html/ipsn_genetic_backgrounds.html) 

![Figure 3](https://github.com/ojziff/ipsn_als_meta/blob/main/figure_3.png)

### [ALS postmortem spinal cord](https://oliverziff.github.io/ipsn_als_meta/html/postmortem_als.html) 

![Extended Data Figure 11](https://github.com/ojziff/ipsn_als_meta/blob/main/postmortem_als.png)

### [pan-ALS vs control differential splicing](https://oliverziff.github.io/ipsn_als_meta/html/ipsn_als_splicing.html)

![Figure 4](https://github.com/ojziff/ipsn_als_meta/blob/main/figure_4.png)

### [Variants and gene fusions]()

![Figure 5](https://github.com/ojziff/ipsn_als_meta/blob/main/figure_5.png)

## Data availability

iPSN RNA sequencing datasets used in this study:

| Reference             | Accession #           | Mutation              | ALS n | Control n | Library type | Paper URL                                    |
|-----------------------|-----------------------|-----------------------|-------|-----------|--------------|----------------------------------------------|
| Sareen et al, 2013    | GSE52202              | C9orf72               |     4 |         4 | polyA        | (https://www.ncbi.nlm.nih.gov/pubmed/24154603) |
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
| Hawkins et al, 2022   | GSE203168             | FUS                   |     2 |         2 | Ribo-zero    | in press                                     |
| Sommer et al, 2022    | GSE201407             | C9orf72               |     6 |         6 | polyA        | in press                                     |
| NeuroLINCS, 2022      | phs001231.v2.p1       | sporadic              |     8 |        14 | Ribo-zero    | https://pubmed.ncbi.nlm.nih.gov/34746695/    |
|                       |                       | SOD1                  |     6 |           |              |                                              |
|                       |                       | C9orf72               |    16 |           |              |                                              |
| Answer ALS, 2022      | AnswerALS data portal | sporadic              |   200 |        42 | Ribo-zero    | https://pubmed.ncbi.nlm.nih.gov/35115730/    |
|                       |                       | C9orf72               |    21 |           |              |                                              |
|                       |                       | SOD1                  |     8 |           |              |                                              |
|                       |                       | 6 other ALS mutations |     9 |           |              |                                              |

For each RNAseq dataset we process samples with nf-core/rnaseq v3.8.1 utilising alignment with STAR and read quantification with salmon. Differential gene expression was performed use DESeq2 and splicing with MAJIQ2, as per the rmarkdown scripts. Schematics were created using Biorender.com.


