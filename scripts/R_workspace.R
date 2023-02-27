# source(here(camp_path,"scripts/functions/R_workspace.R")) # run this script within R environment

# load packages ----------------------
# CRAN
library("tidyverse") # loads ggplot2, tibble, readr, tidyr, purr, dplyr, stringr, forcats
library("lubridate")
library("readxl")
library("vroom")
library("usethis")
library("rmarkdown")
library("data.table")
library("gsubfn")
library("devtools")
library("cowplot")
library("janitor")
library("snakecase")
library("magick")
library("ggplotify")
library("RColorBrewer")
library("scales")
library("ggpubr")
library("ggsci")
library("pheatmap")
library("circlize")
library("ggrepel")
library("colorRamps")
library("VennDiagram")
library("patchwork")
library("ggdendro")
library("dendextend")
library("ggforce")
library("gplots")
library("crayon")
library("ggvenn")
library("ggVennDiagram")
library("tidygraph")
library("ggraph")
library("ggupset")
library("raster")
library("ggrastr")
library("ggprism")
library("reshape2")
library("kableExtra")
library("jtools")

# Bioconductor
library("rhdf5")
library("DESeq2")
library("tximeta")
library("tximport")
library("splines")
library("BiocParallel")
library("rtracklayer")
library("ComplexHeatmap")
library("grid")
library("GO.db")
library("geneplotter")
library("genefilter") 
library("Gviz")
library("trackViewer")
library("GeneOverlap")
library("biovizBase")
library("EnhancedVolcano")
library("clusterProfiler")
library("org.Hs.eg.db")
library("limma") 
library("AnnotationDbi")
library("Glimma")
library("fgsea")
library("sva")
library("vsn")
library("RUVSeq")
library("gprofiler2")
library("rrvgo")
library("topGO") 
library("Biobase")
library("SummarizedExperiment")
library("rstatix")
library("GenomicFeatures")
library("BSgenome.Hsapiens.UCSC.hg38")
library("BSgenome.Mmusculus.UCSC.mm10")
library("BSgenome.Rnorvegicus.UCSC.rn6")
library("GenomicScores") # GC content
library("phastCons100way.UCSC.hg38")
library("DEP")
library("pcr")
library("plyranges")
library("seqinr")
library("GenomicRanges")
library("decoupleR")
library("progeny")
library("dorothea")
library("biomaRt")
library("gdata")
library("piano")
library("parallel")
library("snowfall")
library("GSEABase")
library("systemPipeR")
library("DEGreport")
library("GeneStructureTools")

# github
library("ClusterBurden") #install_github("adamwaring/ClusterBurden")
library("notNMD") # devtools::install_github('betsig/notNMD')
library("irlib") # devtools::install_github("jsha129/irlib")
library("GencoDymo") # remotes::install_github("monahton/GencoDymo") devtools::install_github("helixcn/seqRFLP")
library("VarCon") # remotes::install_github("caggtaagtat/VarCon")  https://bioconductor.org/packages/devel/bioc/manuals/VarCon/man/VarCon.pdf

options(stringsAsFactors = F)


# Gene annotations  ----------------------
# Tx2ens and gene2ens annotation
# mus_musculus.gtf <- import(here(camp_path,"/genomes/ensembl/Mus_musculus.GRCm38.101.gtf")) %>% as_tibble
# saveRDS(mus_musculus.gtf, here(camp_path,"/genomes/ensembl/Mus_musculus.GRCm38.101.gtf.RDS"))
# mus_musculus.gtf <- readRDS(here(camp_path,"/genomes/ensembl/Mus_musculus.GRCm38.101.gtf.RDS"))
# mus.musculus.gene2ens <- mus_musculus.gtf %>% dplyr::select(gene_id, gene_name) %>% unique()
# saveRDS(mus.musculus.gene2ens, here(camp_path,"/genomes/ensembl/Mus_musculus.GRCm38.101.gene2ens.RDS"))
mus.musculus.gene2ens <- readRDS(here(camp_path,"/genomes/ensembl/Mus_musculus.GRCm38.101.gene2ens.RDS")) %>% mutate(gene_name.mouse = gene_name, gene_name = toupper(gene_name))
# mus.musculus.tx2gene <- mus_musculus.gtf %>% dplyr::select(transcript_id, gene_id) %>% unique() %>% drop_na()
# saveRDS(mus.musculus.tx2gene, here(camp_path,"/genomes/ensembl/Mus_musculus.GRCm38.101.tx2gene.RDS"))
mus.musculus.tx2gene <- readRDS(here(camp_path,"/genomes/ensembl/Mus_musculus.GRCm38.101.tx2gene.RDS"))

# homo_sapiens.gtf <- import(here(camp_path,"/genomes/ensembl/Homo_sapiens.GRCh38.99.gtf")) %>% as_tibble
# transcript.gtf <- import(here(camp_path,"/genomes/ensembl/Homo_sapiens.GRCh38.99.gtf")) %>% as_tibble %>% distinct(transcript_id, .keep_all = TRUE) %>% mutate(coord = paste0("chr", seqnames, ":", start, "-", end))
# phast <- phastCons100way.UCSC.hg38
gtf <- readRDS(here(camp_path,"/genomes/ensembl/Homo_sapiens.GRCh38.99.gtf.coords.RDS"))
# gene_coordinates <- gtf %>% dplyr::mutate(gene_id = ensemblID, gene_coords = paste(.$seqnames, ":", .$start, "-", .$end, sep = "")) %>% dplyr::select(gene_coords, gene_id)
# gene2ens <- homo_sapiens.gtf %>% dplyr::select(gene_id, gene_name)
gene2ens <- readRDS(here(camp_path,"/genomes/ensembl/Homo_sapiens.GRCh38.99.gene2ens.RDS")) %>% unique
tx2gene <- readRDS(here(camp_path,"/genomes/ensembl/Homo_sapiens.GRCh38.99.chr_patch_hapl_scaff.tx2gene.rds")) %>% rename("transcript_id" = TXNAME, "gene_id" = GENEID) %>% unique # data.frame with 2 columns: ENST ID & ENSG ID
# homo_sapiens.gtf %>% dplyr::select(transcript_id, gene_name) %>% drop_na %>% unique %>% saveRDS(here(camp_path,"/genomes/ensembl/Homo_sapiens.GRCh38.99.enst2gene.rds"))
# enst2gene <- readRDS(here(camp_path,"/genomes/ensembl/Homo_sapiens.GRCh38.99.enst2gene.rds")) # data.frame with 2 columns: ENST ID & ENSG ID
# ensembl_transcript_biotype <- readRDS(here(camp_path,"/genomes/ensembl/Homo_sapiens.GRCh38.99.transcript_biotype.rds"))
# ensembl_gene_biotype <- readRDS(here(camp_path,"/genomes/ensembl/Homo_sapiens.GRCh38.99.gene_biotype.rds"))


# GO terms ----------------------
# gmtPathways() # reads in as vector read.gmt() # reads in as dataframe
gprofiler_database <- gmtPathways(here(camp_path,"/genomes/gene-ontology/gprofiler_full_hsapiens.name.edit.gmt")) # IDs removed to see term names
# https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp GO pathways
# curated.gene.sets.msigdb <- read.gmt(here(camp_path,"/genomes/gene-ontology/c2.all.v7.1.symbols.gmt")) # includes KEGG, REACTOME, PID, BIOCARTA, CGP
kegg.pathways.msigdb <- gmtPathways(here(camp_path,"/genomes/gene-ontology/c2.cp.kegg.v7.1.symbols.gmt")) #%>% mutate(term = gsub("KEGG|_", " ", term))
# reactome.pathways.msigdb <- read.gmt(here(camp_path,"/genomes/gene-ontology/c2.cp.reactome.v7.1.symbols.gmt")) %>% mutate(term = gsub("REACTOME|_", " ", term))
go.pathways.msigdb <- gmtPathways(here(camp_path,"/genomes/gene-ontology/c5.all.v7.0.symbols.gmt"))# %>% mutate(term = gsub("GO|_", " ", term))
# go.bp.pathways.msigdb<- gmtPathways(here(camp_path,"/genomes/gene-ontology/c5.bp.v7.1.symbols.gmt"))# %>% mutate(term = gsub("GO|_", " ", term))
# go.mf.pathways.msigdb<- gmtPathways(here(camp_path,"/genomes/gene-ontology/c5.mf.v7.1.symbols.gmt"))# %>% mutate(term = gsub("GO|_", " ", term))
# go.cc.pathways.msigdb <- gmtPathways(here(camp_path,"/genomes/gene-ontology/c5.cc.v7.1.symbols.gmt"))# %>% mutate(term = gsub("GO|_", " ", term))
hallmarks.msigdb <- read.gmt(here(camp_path,"/genomes/gene-ontology/h.all.v7.4.symbols.gmt")) %>% mutate(term = gsub("HALLMARK|_", " ", term))
# immune.pathways.msigdb <- read.gmt(here(camp_path,"/genomes/gene-ontology/c7.all.v7.1.symbols.gmt"))
# wikipathways.msigdb <- gmtPathways(here(camp_path,"/genomes/gene-ontology/c2.cp.wikipathways.v7.5.1.symbols.gmt"))# %>% mutate(term = gsub("GO|_", " ", term))
# hpo.msigdb <- gmtPathways(here(camp_path,"/genomes/gene-ontology/c5.hpo.v7.5.1.symbols.gmt"))# %>% mutate(term = gsub("GO|_", " ", term))
RNA_BINDING_PROTEINS <- read_excel(here(camp_path,"/genomes/gene-ontology/rbp_consensus_gerstberger_2014.xls"), sheet = "RBP table") %>% dplyr::pull(`gene name`)

# Cell type markers ----------------------

glia.genes <- c("S100B", "SOX9")
astrocyte.genes <- c("GFAP","AQP4","PLA2G7","SLC39A12","MLC1","DIO2","SLC14A1","ALDH1L1","ALDOC","TTPA","ACSBG1","CHRDL1","SLC4A4","SLC1A2","SLC25A18","SLC1A3","F3","PPP1R3G","FZD2","MERTK","EZR","EVA1A","GJB6","HAPLN1","RFX4","PAPSS2","SLC15A2","PPP1R3C","TLR3","ACOT11","ATP1A2","BMPR1B","PRODH","GLI3","TMEM47","SLC9A3R1","CTH","NTSR2","SLC7A10","VCAM1","FGFR3","CCDC80","ENTPD2","CYBRD1","KCNE5","FAM20A","TNC","TLCD1","S1PR1","CBS","PBXIP1","GRIN2C","ADHFE1","AGT","GLDC","SLC7A2","GJA1","PDK4","EGFR","SOX9","CLDN10","PLCD4","ID4","FMO1","EMP2","LONRF3","HTRA1","MGST1","THRSP")
panreactive.astrocyte.genes <- c("LCN2", "STEAP4", "S1PR3", "TIMP1", "HSPB1", "CXCL10", "CD44", "OSMR", "CP", "VIM", "GFAP", "SERPINA3", "ASPG") # "SERPINA5",
A1.astrocyte.genes <- c("C3", "HLA-A", "HLA-B", "HLA-C","HLA-E","HLA-F", "MICA", "H2-T23","H2-D1", #, "HCP5", "HLA-H", "HLA-G", "HLA-K", "HLA-L", "AL645929.2",
                        "GBP2", "AMIGO2", "SERPING1","GGTA1P","GGTA1", #"GLT6D1", "A3GALT2", "IRGC",
                        "FBLN5", "UGT1A1", "FKBP5", "PSMB8", "SRGN","IIGP1")#, "MX1")
A2.astrocyte.genes <- c("S100A10", "EMP1", "CLCF1", "TGM1", "PTX3", "SPHK1", "CD109", "PTGS2", "SLC10A6", "TM4SF1", "B3GNT5", "CD14")
astrocyte.subtype.genes <- list(panreactive.astrocyte.genes,A1.astrocyte.genes,A2.astrocyte.genes)
astrocyte.reactive.genes <- c(panreactive.astrocyte.genes, A1.astrocyte.genes, A2.astrocyte.genes) %>% unique
astrocyte.reactivity.markers = data.frame("gene_name" = astrocyte.reactive.genes) %>% 
  mutate(group = case_when(gene_name %in% A1.astrocyte.genes ~ "A1",
                           gene_name %in% A2.astrocyte.genes ~ "A2",
                           gene_name %in% panreactive.astrocyte.genes ~ "Pan-reactive"))
# astrocyte.reactivity.markers.mouse = data.frame("gene_name.mouse" = astrocyte.reactive.genes) %>% 
#   mutate(group = case_when(gene_name %in% A1.astrocyte.genes ~ "A1",
#                            gene_name %in% A2.astrocyte.genes ~ "A2",
#                            gene_name %in% panreactive.astrocyte.genes ~ "Pan-reactive"))

cytoskeleton = c("GFAP", "NES", "SYNM", "VIM")
metabolism = c("ALDOC", "FABP7", "MAOB", "TSPO")
chaperone = c("CRYAB", "HSPB1")
secreted = c("C3", "CHI3L1", "LCN2", "SERPINA3", "THBS1", "THBS2", "MT2A", "MT1E") #, "MT1F",  "MT1X", "MT1G", "MT1M", "MT3", "MT1H", "MT1A")
signalling = c("NFATC4","NFAT5", "NTRK2", "IL17R", "S100B", "SOX9", "STAT3", "IRAS1", "IRAS2") # "NFATC1", "NFATC2","NFATC3",
channels = c("SLC1A3", "SLC1A2", "KCNJ10")
escartin.markers = tibble(gene_name = c(cytoskeleton, metabolism, chaperone, secreted, signalling, channels), group = case_when(gene_name %in% cytoskeleton ~ "cytoskeleton", gene_name %in% metabolism ~ "metabolism", gene_name %in% chaperone ~ "chaperone", gene_name %in% secreted ~ "secreted", gene_name %in% signalling ~ "signalling", gene_name %in% channels ~ "channels"))

all.astrocyte.markers <- c("GRN", "PSEN1", "LRP1", "APP", "VIM", "IFNGR1", "GRN", "LAMC3", "PSEN1", "MT3", "TREM2", "CDK6", "IFNG", "LRP1", "PLP1", "IL1B", "NR1D1", "ADORA2A", "SMO", "LDLR", "GFAP", "KRAS", "TSPAN2", "IL6", "TLR4", "APP", "S100A8", "EIF2B5", "TTBK1", "EGFR", "ROR2", "FPR2", "LAMB2", "C1QA", "BACE2", "POU3F2", "DRD1", "ROR1", "MAPT", "NF1", "C5AR1", "DLL1", "AGER", "MIR181C", "MIR181B2", "MIR181B1", "TNF", "CNTF", "MIR142", "KRAS", "NF1", "CCR2", "GPR183", "WNT1", "FZD1", "CTNNB1", "HEXB", "CCL2", "CCR2", "APCDD1", "MMP14", "GPR183", "SCRIB", "CCL3", "IFNGR1", "GRN", "PSEN1", "TREM2", "IFNG", "LRP1", "IL1B", "NR1D1", "ADORA2A", "SMO", "LDLR", "IL6", "APP", "TTBK1", "EGFR", "FPR2", "C1QA", "BACE2", "MAPT", "C5AR1", "AGER", "MIR181C", "MIR181B2", "MIR181B1", "TNF", "CNTF", "MIR142", "SOX8", "PAX6", "VIM", "IFNGR1", "GRN", "LAMC3", "PSEN1", "MT3", "TREM2", "ABL1", "MAPK1", "MAPK3", "MAG", "CDK6", "PRPF19", "SOX6", "IFNG", "NR2E1", "HES1", "ID2", "EPHA4", "LRP1", "PLP1", "TTC21B", "SOX9", "IL1B", "NKX2-2", "BMP2", "NR1D1", "ADORA2A", "LIF", "SMO", "LDLR", "GFAP", "KRAS", "TSPAN2", "IL6ST", "SERPINE2", "IL6", "BIN1", "TLR4", "GCM1", "NTRK3", "MBD1", "APP", "S100A8", "EIF2B5", "TTBK1", "EGFR", "NOTCH1", "S100B", "TAL1", "SHH", "STAT3", "MAP2K1", "ROR2", "GPR37L1", "FPR2", "LAMB2", "ID4", "C1QA", "DAB1", "CLCF1", "PTPN11", "F2", "BACE2", "NOG", "CNTN2", "POU3F2", "DRD1", "ROR1", "MAPT", "NF1", "C5AR1", "HES5", "DLL1", "AGER", "MIR181C", "MIR181B2", "MIR181B1", "TNF", "CNTF", "MIR142", "MAG", "PRPF19", "NR2E1", "HES1", "ID2", "EPHA4", "BMP2", "NR1D1", "LIF", "LDLR", "IL6ST", "SERPINE2", "IL6", "BIN1", "NTRK3", "MBD1", "TTBK1", "NOTCH1", "GPR37L1", "ID4", "DAB1", "CLCF1", "F2", "NOG", "CNTN2", "NF1", "HES5", "MIR181C", "MIR181B2", "MIR181B1", "MIR142", "MAG", "PRPF19", "HES1", "ID2", "BMP2", "LIF", "IL6ST", "SERPINE2", "BIN1", "TTBK1", "NOTCH1", "CLCF1", "MIR142", "NR2E1", "NR1D1", "LDLR", "NTRK3", "MBD1", "GPR37L1", "ID4", "DAB1", "F2", "NOG", "NF1", "HES5", "MIR181C", "MIR181B2", "MIR181B1", "SOX8", "SOX9", "GCM1", "TAL1", "NR1D1", "LDLR", "IL6", "TTBK1", "MIR181C", "MIR181B2", "MIR181B1", "MIR142", "NR1D1", "LDLR", "MIR181C", "MIR181B2", "MIR181B1", "TTBK1", "MIR142", "CCR2", "GPR183", "KCNK2", "MT3", "EZR", "MLC1", "SLC1A2", "ATP1B2", "GFAP", "SYT4", "EIF2S1", "APP", "SLC7A11", "PINK1", "GRM2", "GJB2", "AQP4", "KCNJ10", "SLC17A8", "GRM3", "DMD", "ADGRG1", "MT3", "MLC1", "ATP1B2", "GFAP", "EIF2S1", "AQP4", "SLC17A8", "ADGRG1", "MLH1", "MSH2", "TSC2", "IFNG", "MSH3", "MSH6", "PMS2", "POT1", "APC", "IDH1", "BRCA2", "TP53", "ERBB2", "CDKN2A", "AIFM1", "TSC1", "IDH2", "NF2", "NF1")
all.astrocyte.markers <- c(all.astrocyte.markers,astrocyte.reactive.genes)
microglia.genes <- c("ITGAM", "CX3CR1", "CCL3", "CSFIR", "CCL4", "P2RY12", "C1QB", "PLEK", "GPR183")
myeloid.markers = c("TMEM119", "AIF1", "C1QA", "CSF1R", "ITGAM", "PTPRC", "CEMIP2", "RUNX1", "CD83", "TNF", "CCL3", "CCL2", "IL1A", "TLR2") #  "CX3CR1",
all.microglia.go <- c(go.pathways.msigdb$GO_MICROGLIA_DIFFERENTIATION, go.pathways.msigdb$GO_MICROGLIAL_CELL_PROLIFERATION, go.pathways.msigdb$GO_MICROGLIAL_CELL_MIGRATION, microglia.genes, myeloid.markers)
M1.microglia.genes <- c("IL6","IL12A","IL12B","IL17A","IL17B","IL17C", "IL17F","IL18","IL23A","MMP9","MMP12","CCL5","CCL20","CXCL1","CXCL9","CXCL10","IFNG",go.pathways.msigdb$GO_TUMOR_NECROSIS_FACTOR_MEDIATED_SIGNALING_PATHWAY,"TNFAIP1", "CD14","FCGR3A","FCGR3B","FCGR2A","FCGR2B","FCGR2C","CD40","CD86","HLA-DPA1",go.pathways.msigdb$GO_MHC_CLASS_II_PROTEIN_COMPLEX)
M2.microglia.genes <- c("CCL2","CCL17","CCL22","CCL24","IL4","IL10","IL13","CHIA","CHIT1","CHI3L1","RETN","RETNLA","RETNLB","BDNF","CSF1","FGFR1","FGFR3","FGFR2","FGFR4","FGFRL1","GDNF","IGF1",go.pathways.msigdb$GO_TRANSFORMING_GROWTH_FACTOR_BETA_ACTIVATION,go.pathways.msigdb$GO_TRANSFORMING_GROWTH_FACTOR_BETA_ACTIVATED_RECEPTOR_ACTIVITY,"NGF","NT","GRN","CD163","MRC1")

astrocyte.markers <- unique(c("S100B", "SOX9", "FGFR3", "GFAP", "AQP4", "VIM", "EDNRB","NFIX","PLA2G7","SLC39A12","MLC1","AU067665","DIO2","SLC14A1","ALDH1L1","CYP4F3","ALDOC","TTPA","ACSBG1","CHRDL1","GM266","SLC4A4","SLC1A2","SLC25A18","SLC1A3","F3","PPP1R3G","CYP4F12","LOC388419","FZD2","2900019G14RIK","MERTK","EZR","EVA1A","GJB6","HAPLN1","RFX4","PAPSS2","SLC15A2","PPP1R3C","TLR3","ACOT11","ATP1A2","BMPR1B","DEPDC2","PRODH","GLI3","TMEM47","SLC9A3R1","CTH","NTSR2","SLC7A10","VCAM1","FGFR3","CCDC80","ENTPD2","CYBRD1","KCNE5","FAM20A","KIAA1161","EG328479","TNC","TLCD1","S1PR1","CBS","PBXIP1","GRIN2C","A730056I06RIK","ADHFE1","AGT","GLDC","SLC7A2","BC055107","GJA1","PDK4","EGFR","SOX9","CLDN10","PLCD4","ID4","FMO1","EMP2","LONRF3","HTRA1","MGST1", "THRSP"))
motor.neuron.markers <- unique(c("OLIG2", "NEUROG2", "ISL1", "ISL2", "CHAT", "FOXP1", "LHX3", "HOXA5", "POU3F1", "TSHZ1", "PHOX2B", "OTX2", "MAFB", "NEFH", "TUBB3", "MAP2", "NCAM1", "NEUROD6","GLRA2","CCN3","PRDM8","SLA","CRH","GABRG2","HTR2C","HS3ST2","MAL2","9130024F11RIK","NECAB1","STMN2","GABRA5","NTS","A130090K04RIK","GABRA1","SATB2","GPR88","SYT1","4930544G21RIK","GDA","MYT1L","SLC17A6","A930034L06RIK","A930009L07RIK","NPAS4","CALB1","SLC12A5","EPHA7","VIP","MEF2C","SSTR2","TENM2","PCSK2","SNAP25","SCG2","PGM2L1","PLCXD3","DLX1","VSNL1","SYT4","NRG3","SCN2A","ICAM5","KCNF1","AI838057","CCK","VGF","TMEM130","CAMK4","ASPH","C030017B01RIK","SLC6A7","OLFR1344","ICA1L","MYO5B","NELL1","NEFM","NEFL","TTR","6330527O06RIK","CDH8","SV2B","GAP43","TRHDE","D11BWG0517E","CAMK2B","PENK","RGS4","LPL","CACNA1B","KCNC2","TTC9","L1CAM","9130430L19RIK","CLSTN2","NAPB","CYB561","NKX6-1","NKX6-2","CXADR","GLI2","NKX2-2","FOXA2","CHAT","MNX1","ISL1","FOXP1","LHX3","ETV4","LMO4","NOTCH1"))
ac.mn.gene.markers <- c(astrocyte.markers, motor.neuron.markers)

HLA_genes = gtf %>% filter(grepl("^HLA",gene_name)) %>% distinct(gene_name) %>% pull(gene_name)
nuclear_export.markers = unique(c(go.pathways.msigdb$GO_NUCLEAR_EXPORT, go.pathways.msigdb$GO_NUCLEAR_PORE, go.pathways.msigdb$GO_NUCLEAR_PORE_ORGANIZATION, go.pathways.msigdb$GO_TRANSCRIPTION_EXPORT_COMPLEX, go.pathways.msigdb$GO_NUCLEOCYTOPLASMIC_CARRIER_ACTIVITY, go.pathways.msigdb$GO_RNA_EXPORT_FROM_NUCLEUS))
nucleocytoplasmic_genes <- c(RNA_BINDING_PROTEINS, 
                             go.pathways.msigdb$GO_NUCLEAR_PORE, go.pathways.msigdb$GO_TRANSCRIPTION_EXPORT_COMPLEX, go.pathways.msigdb$GO_NUCLEAR_TRANSCRIBED_MRNA_CATABOLIC_PROCESS_NONSENSE_MEDIATED_DECAY, go.pathways.msigdb$GO_RNA_EXPORT_FROM_NUCLEUS, go.pathways.msigdb$GO_NUCLEOCYTOPLASMIC_CARRIER_ACTIVITY, go.pathways.msigdb$GO_NEGATIVE_REGULATION_OF_NUCLEOCYTOPLASMIC_TRANSPORT)
# nuclear & cyto markers #https://www.thermofisher.com/uk/en/home/life-science/antibodies/primary-antibodies/organelle-marker-antibodies/nuclear-marker-antibodies.html
#https://www.sciencedirect.com/science/article/pii/S2211124715013510?via%3Dihub#fig1
bahar_halpern.nuclear.markers = c("LMNA", "MALAT1", "NEAT1", "LMNA", "LMNB1","LMNB2", "LBR", "ASH2L", "PAF49","CENPA", "DDIT3", "EIF6", "ESR1", "HIF1A", "NUP98","MYC", "NOP2","SIRT1", "SRSF2", "KDM1A", "MLXIPL", "GCK",  
                                  "NLRP6", "GCGR", "POLR1G", "H2AX", "H2AZ1", "H4C1", "H4F3", "H3F3B", "HIST1H2BK", "HIST1H2BM", "HIST1H4A", "HIST2H2AA3", "HIST2H3A", "HIST1H2AC", "HIST1H1E", "HIST1H1C", "HIST4H4", "HIST1H2APS3", "HIST1H2APS5", "HIST2H2BD", "HIST2H2BF", "HIST2H2BK", "HIST2H2AK", "HIST1H2PS3")
bahar_halpern.cytoplasmic.markers = c("GAPDH", "KRT19", "TUBB", "HIF1A",  "VIM",  "NOS2","INS1","INS2", "ALDOA", "ACTB", "GLUL", "SLC2A2", "ASS1", "ALB","APOE", "PCK1")
bahar_halpern.nuclear_cytoplasmic.markers = c(bahar_halpern.nuclear.markers, bahar_halpern.cytoplasmic.markers)


ALS_RBPs = c("TARDBP", "SFPQ", "FUS", "ELAVL3", "HNRNPA1", "HNRNPA2B1", "MATR3", "ATXN2", "TAF15", "TIA1", "EWSR1", "ANG", "SETX", "ELP3") # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6909825/
ALS_GWAS = c("C9orf72", "UNC13A", "SOD1", "SCFD1", "MOBP", "RPSA", "KIF5A", "CFAP410","GPX3", "TNIP1", "SLC9A8", "SPATA2", "TBK1", "ERGIC1", "NEK1", "COG3", "PTPRN2", "KANK1","CAMTA1","IDE","LGALSL","NEFH")  # "HLA", https://www.nature.com/articles/s41588-021-00973-1/tables/1 C21orf2 == CFAP410
ALS_EWAS = c("DHCR24","FAM167B","LCK","VWA1","IL6R","TTC7A","RXFP4","SGMS2","HCLS1","GOLIM4","SGMS2","SLC7A11","MSMO1","PIK3R1","SLC26A8","GABBR1","FKBP5","SPIDR","C6orf223","SMURF1","RCL1","DNM1","NOLC1","EMG1","KRT2","LPCAT3","C1S","KRT2","RASA3","CILP","CCDC102A","IGF1R","PLEKHF1","ZFPM1","SIRT2","VMP1","NOL4L","ABCG1","TTC38") # https://www.science.org/doi/10.1126/scitranslmed.abj0264
ALS_modifiers = c("APOE", "VEGFA", "SMN2", "CHGB", "PPARGC1A","SLC52A2", "EPHA4")
sporadic_ALS_genes = c("NDUFS4","AC106707.1","ZC3H7B","AC023095.1","CCDC59","TXNP1","INPP5F","TNRC18","TOP2A","THRAP3","TRPM3","ATP10A","FAM184B","AC096747.1","NCS1","AC007690.1","AL033528.3","RN7SL33P","COX5A","AL161629.1","SLF1","LIPH","RPL5P16") #https://www.frontiersin.org/articles/10.3389/fgene.2022.851496/full
leigh_brown_ALS_FTD_genes = c("AR","PRNP","MAPT","ATN1","DNAJC5","DNAJC7","GLT8D1","DNMT1","APP","PSEN2","JPH3","ITM2B","TBP","PFN1","HTT","EPM2A","ATXN1","NOP56","NHLRC1", "UNC13A","UNC13B", "STMN2")
ALS_genes = unique(c(kegg.pathways.msigdb$KEGG_AMYOTROPHIC_LATERAL_SCLEROSIS_ALS, ALS_RBPs, ALS_GWAS, ALS_EWAS, ALS_modifiers, leigh_brown_ALS_FTD_genes, HLA_genes,"OPTN", "ALS2", "VAPB", "CHMP2B","FIG4", "UBQLN2","SQSTM1", "SIGMAR1", "VCP", "DCTN1", "PFN1", "CHCHD10", "TUBA4A", "ANXA11", "DAO", "ERBB4", "SPG11", "SPTLC1", "CCNF", "STMN2","IL18RAP")) # https://www.ncbi.nlm.nih.gov/books/NBK1116/?term=Amyotrophic+lateral+sclerosis https://www.ncbi.nlm.nih.gov/books/NBK1450/ https://www.sciencedirect.com/science/article/pii/S2352396421005892?via%3Dihub 
TDP43_spliced_genes = c("AKT3", "ANK3", "ANKRD20A5P", "APLP2", "ARHGEF11", "ATP8A2", "C2CD5", "CYFIP2 FAM193A", "KALRN", "CACNA1E", "CYP2C8", "FKBP15", "KCNQ2", "DLGAP4", "HP1BP3", "KIF3A", "EIF4G3", "ICA1", "KMT2C", "EPB41L1", "IMMT", "ERC2", "IQCB1", "CADPS", "CAMK2B", "CEP290", "CNOT6", "LOC100630923", "MADD", "SEPT7P2", "SETD5", "SH3KBP1", "MARK3", "PLEKHA1", "PTPRT", "MNAT1", "PLEKHA5", "RAP1GAP", "NCAM1", "POLDIP3", "RAPGEF4", "NFIX", "PPP6R3", "RAPGEF6", "SHISA9", "RNASET2", "STMN2", "SCN1A", "PACRGL", "PRUNE2", "PDE2A", "PTPRD", "STXBP1", "STXBP5L", "UNC79", "SYNE1", "UQCRC2", "SYT7", "VWA8", "WHSC1L1", "ZNF138", "ZNF292", "TRAPPC12", "TRRAP", "UNC13A")


focal.adhesion.go <- c("LAP3", "CD99", "LASP1", "ABCB4", "ITGA3", "ITGA2B", "RALA", "CD9", "MRC2", "TSPAN9", "PLAUR", "EHD3", "CAPN1", "FHL1", "VIM", "CD44", "ARHGAP31", "VCL", "TNC", "CTNNA1", "HSPA5", "LIMA1", "BCAR1", "CYBA", "SYNE2", "GDI2", "PPP1R12A", "NCKAP1", "RPL18", "CNN2", "SLC9A3R2", "TLE2", "RHOA", "FGFR3", "PABPC1", "CDC42", "MAP4K4", "RPL31", "ACTN1", "LIMS2", "PVR", "FERMT2", "CLASP1", "HACD3", "ACTB", "REXO2", "MCAM", "USP33", "APBB1IP", "ACTN2", "ITGA8", "FAP", "TNS1", "SENP1", "DNM2", "EPB41L2", "RAB21", "PTPRC", "ITGB5", "RPS5", "FAT1", "RAB10", "CD59", "CPNE3", "CTTN", "NOX4", "TRIP6", "ADD1", "CASS4", "ASAP3", "RPL6", "RPLP0", "PXN", "SLC9A1", "ICAM1", "ITGA6", "SNAP23", "EZR", "SORBS1", "JAK2", "NRP1", "MISP", "MAPK1", "TRIOBP", "PACSIN2", "RPL3", "MYH9", "ZFYVE21", "PROCR", "FERMT1", "HCK", "MAPRE1", "CD99L2", "ARHGEF7", "FLT1", "MAPK3", "CORO2B", "RPS16", "RPS19", "ITGB8", "CAV2", "CAV1", "HSPB1", "LIMK1", "ENG", "PDLIM1", "GIT1", "RPL19", "PFN1", "SLC6A4", "YWHAE", "LAMTOR3", "HSPA8", "LPXN", "CBL", "CD81", "RPS13", "PPFIBP1", "CORO1C", "TNS2", "TRPV4", "ARPC3", "NEDD9", "OPRM1", "WASF1", "PTK7", "HSPA9", "PDGFRB", "PRKAR2A", "ACTR3", "EPB41L5", "ITGB6", "ITGA4", "RPS15", "IL1RL1", "FHL2", "RND3", "RPL22", "RHOU", "ARHGEF2", "DOCK7", "CD46", "CNN3", "GNA13", "TEK", "SORBS3", "PTK2B", "CAT", "RPL5", "PLAU", "ADGRE5", "LRP1", "SDC4", "STX16", "RPS10", "AHNAK", "EFNB2", "RPL23", "VASP", "FLRT3", "RRAS", "AP006333.1", "AIF1L", "MAP2K2", "PTPN12", "CDC42EP1", "RAC2", "FLNC", "ARHGAP22", "PALLD", "AJUBA", "STARD8", "CNN1", "NECTIN2", "ACTN4", "ARPC1B", "PAK4", "AKAP12", "CAP1", "RPL27", "PPFIA1", "TNS4", "ITGB4", "FLOT2", "PTPRA", "DCTN4", "MPRIP", "KRAS", "RRAS2", "NUMB", "YWHAQ", "ANXA1", "TES", "AVIL", "FLNB", "LMO7", "LCP1", "TNS3", "RAC1", "ARPC5L", "TLN1", "HMGA1", "FLOT1", "MDC1", "ATAT1", "SDCBP", "RDX", "KIF23", "ITGA11", "RPLP1", "ADAM10", "BCAR3", "ACTR2", "ITGAV", "ARHGAP24", "SCARB2", "PARVG", "GIT2", "ITGB7", "IQGAP1", "TGFB1I1", "ARMC5", "CDH13", "RPS2", "PDPK1", "CLTC", "GRB7", "RPS11", "RPL13A", "EPHA2", "HSPG2", "RPS8", "MTF2", "DCAF6", "PRUNE1", "PIP5K1A", "S100A7", "ARF1", "RHOB", "FBLN7", "LIMD1", "PHLDB2", "LPP", "RPS3A", "ARHGAP26", "G3BP1", "GNA12", "EGFR", "SH3KBP1", "CASK", "MSN", "ZNF185", "RPL7", "GSN", "RPL7A", "RSU1", "CAPN5", "PAK1", "RPS3", "HYOU1", "ITGB1", "DIXDC1", "TWF1", "ADAM17", "DST", "ARL14EP", "TADA1", "GJA1", "DAB2", "THY1", "PGM5", "ENAH", "SORBS2", "RPL30", "MMP14", "FZD1", "SHROOM4", "CSRP1", "ACTC1", "ZYX", "ITGB2", "RPL8", "NPHS1", "ITGA5", "JAK1", "FBLIM1", "NEXN", "ARPC5", "DDR2", "NCSTN", "CAPN2", "XIRP2", "ARPC2", "NFASC", "CLASP2", "RPL9", "ANXA5", "ITGA2", "RPS14", "DLC1", "SLC4A2", "YWHAZ", "HNRNPK", "ARF6", "ILK", "HSP90B1", "B2M", "PPIB", "YWHAB", "MAPRE2", "PDIA3", "TPM4", "SRP68", "CTNNB1", "FAM107A", "IRF2", "XIRP1", "ADAM9", "SNTB2", "TM4SF20", "MAP2K1", "PTK2", "LIMS1", "ALCAM", "YWHAG", "PDCD6IP", "CDH2", "RPS9", "RPS7", "TLN2", "KLF11", "SNTB1", "BSG", "GNB2", "SYNPO2", "CORO1B", "CFL1", "RPL38", "DAG1", "PEAK1", "CSPG4", "JUP", "RPL4", "CSRP2", "YES1", "RHOG", "RPLP2", "CD151", "FLII", "PLEC", "GAK", "CALR", "FZD2", "NPM1", "ADGRB1", "FES", "CAV3", "RPS17", "FHL3", "ACTG1", "UBOX5", "FLRT2", "LMLN", "P4HB", "ATP6V0C", "PIP5K1C", "PPP1CC", "CHP1", "SPRY4", "NHS", "PEAK3", "FOCAD", "PARVB", "PPIA", "EVL", "AFAP1", "MME", "PDLIM7", "FLNA", "ANXA6", "IGF2R", "PCBP2", "SRC", "SVIL", "DPP4", "PARVA", "RPL37A", "RPL12", "MPZL1", "RPS4X", "ITGBL1", "RPL10A", "L1CAM", "TGM2", "LAYN", "HSPA1B", "HSPA1A", "ARL2", "PPP1CB", "RPS29", "ITGA1", "TSPAN4", "RPS18", "ALKBH6", "PI4KA", "SCARF2", "ACTN3", "LIMS4", "LIMS3", "ITGB3", "AC068234.1", "CYFIP1", "PRAG1", "MARCKS")
collagen.formation.go <- c("FOXC1","COL11A1","COL5A3","ADAMTS2","TGFB2","MMP11","CHADL","COMP","PLOD3","AEBP1","TGFBR1","COL1A1","COL12A1","LOX","LOXL3","FMOD","P4HA1","COLGALT1","PXDN","COL5A1","LOXL2","CYP1B1","EMILIN1","LOXL4","ADAMTS14","COL2A1","LUM","RB1","P3H4","DPT","SFRP2","SERPINH1","VIPAS39","ADAMTS3","ACAN","DDR2","COL1A2","ATP7A","GREM1","SERPINF2","TNXB","COL3A1","FOXC2","ANXA2","COL14A1","OPTC","NF1","COL11A2","COL5A2","SCX","MIR29B1","MIR29B2 MMP25","MRC2","FAP","COL19A1","ADAMTS2","MMP2","MMP11","MMP9","CST3","MMP15","VSIR","CTSD","MMP8","MMP19","PEPD","MMP24","CTSL","MMP7","MMP20","MMP27","MMP13","ADAMTS14","FURIN","CTSK","ADAM15","MMP3","ITGB1","MMP21","MMP16","ADAMTS3","MMP14","CTSS","CTSB","MMP10","MMP26","KLK6","PHYKPL","TMPRSS6","MMP23B","PRTN3","MMP1","COL13A1","MMP17","COL15A1","MMP12","MMP28","PRSS2","MMP25","MRC2","VIM","TRAM2","FAP","COL19A1","ADAMTS2","MMP2","P2RX7","P3H2","SUCO","MMP11","HIF1A","MMP9","CST3","RGCC","MMP15","SMPD3","TGFB1","PLOD3","ENG","VSIR","COL1A1","P3H3","TNS2","PPARD","PDGFRB","ERRFI1","RAP1A","P3H1","CTSD","MMP8","MYB","ARG1","CCN2","TGFB3","GOT1","MMP19","PEPD","AMELX","BMP4","MMP24","ID1","COL5A1","PPARG","CTSL","IL6","MMP7","MMP20","MMP27","MMP13","EMILIN1","ADAMTS14","FURIN","ARRB2","CBX8","P3H4","RCN3","CTSK","ADAM15","SERPINH1","MMP3","ITGB1","VIPAS39","MMP21","MMP16","ADAMTS3","MMP14","CREB3L1","CYGB","WNT4","CTSS","NPPC","IHH","UCN","ITGA2","COL1A2","CTSB","LARP6","SERPINB7","MFAP4","MMP10","MMP26","SERPINF2","KLK6","TNXB","PHYKPL","CIITA","F2","F2R","VPS33B","TMPRSS6","MMP23B","PRTN3","HDAC2","MMP1","COL13A1","MMP17","COL15A1","MIR149","MIR218-1","MIR218-2","SCX","MMP12","MMP28","PRSS2","MIR145","MIR92A1","MIR29B1","MIR29A","MIR21","MIR29B2","MIR92A2","COL4A4","COL1A1","COL4A2","ITGA11","UBASH3B","DDR2","ITGA2","SYK","COL4A3","OSCAR","COL4A1","COL4A5","COL4A6","DDR1","TLL1","COL9A2","COL23A1","LAMA3","LAMC2","COL11A1","COL17A1","P4HA2","COL5A3","COL4A4","COL19A1","PLOD1","COL16A1","ADAMTS2","P3H2","ITGA6","COL9A3","TLL2","MMP9","COL20A1","PCOLCE","PLOD3","COL1A1","P3H3","COL12A1","COL9A1","LOX","COL7A1","LOXL3","P3H1","P4HA1","COL10A1","COL21A1","LOXL1","COLGALT1","PXDN","COL5A1","ITGB4","LOXL2","COL4A2","CTSL","CTSV","MMP7","MMP20","MMP13","LOXL4","ADAMTS14","COL2A1","COL6A1","COL6A2","COL8A1","SERPINH1","P4HA3","MMP3","PLOD2","ADAMTS3","COL26A1","CTSS","COL6A3","PCOLCE2","COL1A2","CTSB","PPIB","BMP1","COL3A1","COL4A3","COL22A1","CRTAP","COL24A1","COL8A2","COL6A5","CD151","PLEC","COL18A1","P4HB","COL4A1","COL14A1","COL4A5","COL25A1","COL27A1","LAMB3","COL13A1","COL4A6","COLGALT2","COL11A2","COL5A2","COL15A1","COL6A6","COL28A1")
extracellular.matrix.go <- c("CFLAR","ST7","ITGAL","ITGA3","ITGA2B","NOX1","MMP25","DCN","CAPN1","PHLDB1","CD44","TNFRSF1B","IBSP","TIMP2","TLL1","VCAN","CDH1","TNC","GPM6B","COL9A2","ELN","LAMC3","COL23A1","LAMA3","FOXC1","LAMC2","COL11A1","COL17A1","TNFRSF1A","BCL3","CLASP1","NTN4","FSCN1","ICAM3","NFKB2","FBLN1","ITGA8","FAP","COL5A3","CPB2","COL4A4","COL19A1","ITGB5","ITGAE","COL16A1","B4GALT1","ERO1B","ADAMTS2","MMP2","NID2","KIF9","ICAM1","LAMB4","LAMB1","ITGA6","CCDC80","CMA1","COL9A3","TGFB2","TMEM38B","TLL2","ABL1","MADCAM1","MMP11","PDGFB","CHADL","CTSG","MMP9","FERMT1","CST3","LAMA1","TIMP1","RGCC","MMP15","HAS3","SMPD3","CRISPLD2","FOXF1","ERCC2","TGFB1","ICAM4","ICAM5","HAS1","COMP","HPN","ITGB8","CAV2","CAV1","DNAJB6","SERPINE1","PLOD3","AEBP1","MEGF9","TGFBR1","ECM2","ENG","SPOCK2","KAZALD1","SH3PXD2A","ICAM2","COL1A1","VTN","VWF","MYF5","COL12A1","ADTRP","COL9A1","NR2E1","SMOC2","LAMA4","LOX","SPARC","COL7A1","MPV17","ITGB6","ITGA4","LOXL3","FN1","TNR","QSOX1","EXOC8","NID1","MFAP2","MMP8","TTR","CCN2","SPP1","TGFBI","FMOD","PLG","SERAC1","RECK","P4HA1","PRDX4","MMP19","LRP1","COL10A1","MATN4","SOX9","TCF15","MMP24","CAPNS1","KDR","LOXL1","NCAN","COLGALT1","PXDN","COL5A1","LAMA5","POMT1","GFAP","RAMP2","MATN3","ITGB4","BCAN","HAPLN2","POSTN","MYH11","SPINK5","LOXL2","PDGFRA","COL4A2","ETS1","CTSL","ADAM19","ITGA7","AGT","LAMC1","LCP1","IL6","CTSV","FOXF2","FLOT1","SULF1","MMP7","MMP20","MMP27","MMP13","THBS1","ITGA11","ADAM10","CYP1B1","EMILIN1","LOXL4","ADAMTS14","ITGAV","FGF2","SLC39A8","FBN2","COL2A1","LUM","ITGB7","RB1","FBLN5","FURIN","ITGAX","GFOD2","P3H4","COL6A1","COL6A2","APP","HSPG2","CCN1","ITGA10","DPT","ADAMTSL4","CTSK","ADAM15","ITGA9","COL8A1","PHLDB2","SFRP2","HAPLN1","SCUBE3","CSGALNACT1","NOTCH1","ADAM12","GAS2","HSD17B12","SERPINH1","MMP3","ITGB1","VIPAS39","ADAM8","DSPP","DMP1","ABI3BP","WNT3A","MMP21","JAM2","ADAMTS5","MMP16","ADAMTS3","ITGAD","MMP14","MYO1E","CREB3L1","ACAN","F11R","ADAMTS4","SCUBE1","CARMIL2","ITGB2","MPZL3","FGFR4","NPHS1","ITGA5","PDPN","MATN1","VCAM1","DDR2","OLFML2B","CAPN2","CTSS","COL6A3","ELF3","IHH","FBLN2","CLASP2","ADAMTS9","PTX3","APBB2","MELTF","ITGA2","EGFLAM","KLKB1","COL1A2","TNFRSF11B","ATP7A","HTRA1","JAM3","SPINT1","FBN1","WDR72","MFAP4","MMP10","GREM1","SMAD3","MMP26","SPINT2","SERPINF2","KLK4","KLK2","KLK5","TNXB","BMP1","COL3A1","NPNT","CTRB1","CTRB2","COL4A3","KLK7","PTK2","COL22A1","ANTXR1","ITGAM","FSHR","HAS2","COL24A1","RXFP1","FGG","FGA","FGB","COL8A2","ANGPTL7","LAMB2","TPSAB1","BSG","EFEMP2","HPSE2","ADAMTS20","NDNF","DAG1","SH3PXD2B","A2M","FOXC2","RIC8A","VWA1","WASHC1","OTOL1","BGN","ANXA2","COL18A1","GAS6","WT1","FLRT2","OLFML2A","TMPRSS6","COL4A1","THSD4","COL14A1","COL4A5","AGRN","OPTC","MMP23B","NOXO1","SULF2","LAMA2","MMP1","NF1","COL27A1","CD47","LAMB3","PDGFA","COL13A1","ELANE","COL4A6","MFAP5","DPP4","ADAMTSL2","ERO1A","MMP17","EGFL6","COL11A2","COL5A2","COL15A1","DDR1","PRSS1","VIT","SERPINB5","COLQ","ITGA1","COL28A1","ATXN1L","TNF","MARCOL","CAPNS2","ITGB3","SCX","PECAM1","MMP12","MMP28","MIR98","PRSS2","MIR29B1","MIR29B2")
glutamate.uptake.go <- unique(c(go.pathways.msigdb$GO_GLUTAMATE_BINDING, go.pathways.msigdb$GO_GLUTAMATE_BIOSYNTHETIC_PROCESS, go.pathways.msigdb$GO_GLUTAMATE_GATED_CALCIUM_ION_CHANNEL_ACTIVITY, go.pathways.msigdb$GO_GLUTAMATE_METABOLIC_PROCESS, gprofiler_database$`glutamatergic synapse`, gprofiler_database$`glutamate reuptake`))
phagocytosis.go <- unique(c(go.pathways.msigdb$GO_REGULATION_OF_PHAGOCYTOSIS, go.pathways.msigdb$GO_PHAGOCYTOSIS))
all.reactivity.go <- c(focal.adhesion.go, collagen.formation.go, extracellular.matrix.go, go.pathways.msigdb$GO_RESPONSE_TO_WOUNDING,
                       RNA_BINDING_PROTEINS,
                       go.pathways.msigdb$GO_IMMUNE_EFFECTOR_PROCESS, go.pathways.msigdb$GO_CYTOKINE_MEDIATED_SIGNALING_PATHWAY, go.pathways.msigdb$GO_RESPONSE_TO_CYTOKINE, 
                       go.pathways.msigdb$GO_RESPONSE_TO_INTERFERON_ALPHA, go.pathways.msigdb$GO_RESPONSE_TO_INTERFERON_BETA, go.pathways.msigdb$GO_RESPONSE_TO_INTERFERON_GAMMA,
                       go.pathways.msigdb$GO_LIPOPOLYSACCHARIDE_MEDIATED_SIGNALING_PATHWAY, go.pathways.msigdb$GO_LIPOPOLYSACCHARIDE_IMMUNE_RECEPTOR_ACTIVITY,
                       gprofiler_database$`Cellular responses to stress`, all.astrocyte.markers, 
                       go.pathways.msigdb$GO_PHAGOCYTOSIS, glutamate.uptake.go, go.pathways.msigdb$GO_REGULATION_OF_PHAGOCYTOSIS,
                       go.pathways.msigdb$GO_EXCITATORY_SYNAPSE,go.pathways.msigdb$GO_NEUROFILAMENT, go.pathways.msigdb$GO_NEUROFILAMENT_CYTOSKELETON_ORGANIZATION,
                       go.pathways.msigdb$GO_NEUROINFLAMMATORY_RESPONSE, go.pathways.msigdb$GO_GLIAL_CELL_ACTIVATION, 
                       go.pathways.msigdb$GO_NEURON_APOPTOTIC_PROCESS, go.pathways.msigdb$GO_NEURON_CELL_CELL_ADHESION, go.pathways.msigdb$GO_NEURON_DEATH,
                       go.pathways.msigdb$GO_RESPONSE_TO_OXIDATIVE_STRESS,go.pathways.msigdb$GO_RESPONSE_TO_REACTIVE_OXYGEN_SPECIES, 
                       go.pathways.msigdb$GO_REACTIVE_OXYGEN_SPECIES_METABOLIC_PROCESS, go.pathways.msigdb$GO_NEURON_DEATH_IN_RESPONSE_TO_OXIDATIVE_STRESS, go.pathways.msigdb$GO_RESPONSE_TO_REACTIVE_OXYGEN_SPECIES,  
                       go.pathways.msigdb$GO_INTERMEDIATE_FILAMENT_CYTOSKELETON,go.pathways.msigdb$GO_INTERMEDIATE_FILAMENT_BASED_PROCESS, go.pathways.msigdb$GO_ACTIN_FILAMENT, go.pathways.msigdb$GO_CYTOSKELETON_ORGANIZATION,
                       go.pathways.msigdb$GO_MICROGLIA_DIFFERENTIATION, go.pathways.msigdb$GO_MICROGLIAL_CELL_PROLIFERATION, go.pathways.msigdb$GO_MICROGLIAL_CELL_MIGRATION, microglia.genes, myeloid.markers)

all.inflammation.go <- c(focal.adhesion.go, collagen.formation.go, extracellular.matrix.go, go.pathways.msigdb$GO_RESPONSE_TO_WOUNDING,
                         go.pathways.msigdb$GO_IMMUNE_EFFECTOR_PROCESS, go.pathways.msigdb$GO_CYTOKINE_MEDIATED_SIGNALING_PATHWAY, go.pathways.msigdb$GO_RESPONSE_TO_CYTOKINE, 
                         go.pathways.msigdb$GO_RESPONSE_TO_INTERFERON_ALPHA, go.pathways.msigdb$GO_RESPONSE_TO_INTERFERON_BETA, go.pathways.msigdb$GO_RESPONSE_TO_INTERFERON_GAMMA,
                         go.pathways.msigdb$GO_LIPOPOLYSACCHARIDE_MEDIATED_SIGNALING_PATHWAY, go.pathways.msigdb$GO_LIPOPOLYSACCHARIDE_IMMUNE_RECEPTOR_ACTIVITY,
                         gprofiler_database$`Cellular responses to stress`, 
                         go.pathways.msigdb$GO_PHAGOCYTOSIS, go.pathways.msigdb$GO_REGULATION_OF_PHAGOCYTOSIS,
                         go.pathways.msigdb$GO_NEUROINFLAMMATORY_RESPONSE, go.pathways.msigdb$GO_GLIAL_CELL_ACTIVATION, 
                         go.pathways.msigdb$GO_NEURON_APOPTOTIC_PROCESS, go.pathways.msigdb$GO_NEURON_CELL_CELL_ADHESION, go.pathways.msigdb$GO_NEURON_DEATH,
                         go.pathways.msigdb$GO_RESPONSE_TO_OXIDATIVE_STRESS,go.pathways.msigdb$GO_RESPONSE_TO_REACTIVE_OXYGEN_SPECIES, 
                         go.pathways.msigdb$GO_REACTIVE_OXYGEN_SPECIES_METABOLIC_PROCESS, go.pathways.msigdb$GO_NEURON_DEATH_IN_RESPONSE_TO_OXIDATIVE_STRESS, go.pathways.msigdb$GO_RESPONSE_TO_REACTIVE_OXYGEN_SPECIES,  
                         go.pathways.msigdb$GO_INTERMEDIATE_FILAMENT_CYTOSKELETON,go.pathways.msigdb$GO_INTERMEDIATE_FILAMENT_BASED_PROCESS, go.pathways.msigdb$GO_ACTIN_FILAMENT, go.pathways.msigdb$GO_CYTOSKELETON_ORGANIZATION,
                         go.pathways.msigdb$GO_MICROGLIA_DIFFERENTIATION, go.pathways.msigdb$GO_MICROGLIAL_CELL_PROLIFERATION, go.pathways.msigdb$GO_MICROGLIAL_CELL_MIGRATION)


astrocyte.reactivity.splicing.labels = c("VIM", "COL1A1", "MYL6", "DDX5", "FN1", "FLNA", "FUS", "SFPQ","OSMR", "NONO", "HLA-A", "HLA-B", "HLA-C", "STAT3", "RND3", "SLC26AA6","PIMREG", "RRP8", "RECK", "LOXL3","UBN1", "COL7A1","MX1",
                                         "TTBK1", "TSC1", "MBD1", "EMP1", "PRPF4", "PIF1", "AARS2", "PLSCR1", "NARPT", "CSRP1", "IRF7", "ID4", "AGER", "GBP2", "ANKZF1", "AGRN", "GTBP2", "TTBK1", "CCS", "TGFB1I1","COL1A1A", "TNC", 
                                         "RPL10", "NUP199", "LAMB2", "IDH1", "TP53", "LRP1")
astrocyte.reactivity.labels = c("VIM", "COL1A1", "MYL6", "FN1", "FLNA", "OSMR", "HLA-A", "HLA-B", "HLA-C", "STAT3", "RND3", "SLC26AA6","PIMREG", "RRP8", "RECK", "LOXL3","UBN1", "COL7A1","MX1",
                                "TTBK1", "TSC1", "MBD1", "EMP1", "PRPF4", "PIF1", "AARS2", "PLSCR1", "NARPT", "CSRP1", "IRF7", "ID4", "AGER", "GBP2", "ANKZF1", "AGRN", "GTBP2", "TTBK1", "CCS", "TGFB1I1","COL1A1A", "TNC", 
                                "RPL10", "NUP199", "LAMB2", "IDH1", "TP53", "LRP1",
                                "TINF2", "ILK", "MCAM", "FN1", "FBLN5", "CAPN2", "TGFB1I1", "CDKN1A", "ADAM19", "SERPINE1", "ITGA6", "ADAMTSL4", "ITGA5", "PDLIM7", "FLOT1", "DIXDC1", "BAG3", "PVR")  # overlap increased exp & decreased IR
astrocyte.reactivity.labels = c(astrocyte.reactivity.labels, astrocyte.reactive.genes)

astrocyte.reactivity.labels.ac_nuc <- c("VIM", "COL1A1", "MYL6", "FN1", "FLNA", "OSMR", "HLA-A", "HLA-B", "HLA-C", "STAT3", "RND3", "SLC26AA6",
                                        "PIMREG", "RRP8", "RECK", "MX1", "TTBK1", "MBD1", "EMP1", "PRPF4", "PLSCR1", "NARPT", "CSRP1", "IRF7", "ANKZF1", "GTBP2",   
                                        "CCS", "TGFB1I1", "COL1A1A", "RPL10", "NUP199", "LAMB2", "IDH1", "LRP1", "TINF2", "ILK" , "MCAM", "CAPN2", "ADAM19", "SERPINE1", "ITGA6", "PDLIM7", "FLOT1",  
                                        "BAG3", "LCN2", "STEAP4", "S1PR3", "TIMP1", "HSPB1", "CXCL10", "CD44", "CP", "GFAP", "SERPINA3", "ASPG", "C3", "HLA-F", "MICA", "H2-T23", "H2-D1", "AMIGO2", "SERPING1", "GGTA1P", "GGTA1", "UGT1A1",  
                                        "FKBP5", "PSMB8", "SRGN", "IIGP1", "S100A10", "CLCF1", "TGM1", "PTX3", "SPHK1", "CD109", "PTGS2", "SLC10A6", "TM4SF1", "B3GNT5", "CD14")  

astrocyte.reactivity.labels.ac_cyt <- c("VIM", "COL1A1", "FN1", "FLNA", "OSMR", "HLA-A", "HLA-B", "HLA-C", "STAT3", "RND3", "SLC26AA6","PIMREG", "RRP8", "RECK", "LOXL3", "COL7A1", "MX1", "EMP1",
                                        "NARPT", "CSRP1", "ID4", "AGER", "GBP2", "ANKZF1", "AGRN", "GTBP2",   "CCS", "TGFB1I1", "COL1A1A", "TNC", "NUP199", "LAMB2", "IDH1", "TP53", "ILK" ,    
                                        "MCAM", "FBLN5", "CAPN2", "CDKN1A", "ADAM19", "SERPINE1", "ADAMTSL4", "ITGA5", "PDLIM7", "FLOT1", "DIXDC1",  "BAG3", "LCN2", "STEAP4", "S1PR3", "TIMP1", "HSPB1", "CXCL10", "CD44", "CP", "GFAP", "SERPINA3",
                                        "ASPG", "C3", "HLA-E", "MICA", "H2-T23", "H2-D1", "AMIGO2", "SERPING1", "GGTA1P", "GGTA1", "UGT1A1",  "FKBP5", "SRGN", "IIGP1", "S100A10", "CLCF1", "TGM1", "PTX3", "SPHK1", "CD109", "PTGS2", "SLC10A6", "TM4SF1", "B3GNT5", "CD14")  

astrocyte.reactivity.labels.ac_nucIR_cytDGE <- c("VIM", "COL1A1", "FN1", "HLA-A", "HLA-B", "HLA-C", "SLC26AA6","PIMREG", "RRP8", "RECK", "MX1", "TTBK1", "EMP1", "PRPF4",   
                                                 "NARPT", "CSRP1", "IRF7", "GBP2", "ANKZF1", "GTBP2","CCS", "TGFB1I1", "COL1A1A", "RPL10", "NUP199", "TP53", "TINF2", "ILK",    
                                                 "MCAM", "CAPN2", "ADAM19", "SERPINE1", "ITGA6", "ADAMTSL4", "FLOT1","BAG3", "PVR", "LCN2", "STEAP4", "S1PR3", "TIMP1", "HSPB1", "CXCL10", "CD44", "CP", "GFAP", "SERPINA3",
                                                 "ASPG", "C3", "HLA-E", "HLA-F", "MICA", "H2-T23", "H2-D1", "AMIGO2", "SERPING1", "GGTA1P", "GGTA1", "UGT1A1",  "FKBP5", "PSMB8", "SRGN", "IIGP1", "S100A10", "CLCF1", "TGM1", "PTX3", "SPHK1", "CD109", "PTGS2", "SLC10A6", "TM4SF1", "B3GNT5", "CD14")  

astrocyte.reactivity.labels.ac_c9orf72 <- c("VIM", "COL1A1", "FN1", "FLNA", "OSMR", "HLA-A", "HLA-B", "HLA-C", "STAT3", "RND3", "SLC26AA6","PIMREG", "RRP8", "RECK", "UBN1", "MX1", "TTBK1", "TSC1", "MBD1", "EMP1", "PRPF4",   
                                            "PIF1", "AARS2", "PLSCR1", "NARPT", "CSRP1", "IRF7", "ID4", "AGER", "GBP2", "ANKZF1", "AGRN", "GTBP2","CCS", "TGFB1I1", "COL1A1A", "TNC", "RPL10", "NUP199", "LAMB2", "IDH1", "TP53", "LRP1", "TINF2", "ILK" ,    
                                            "MCAM", "FBLN5", "CAPN2", "CDKN1A", "ADAM19", "SERPINE1", "ITGA6", "ADAMTSL4", "ITGA5", "PDLIM7", "FLOT1", "DIXDC1", "BAG3", "PVR", "LCN2", "STEAP4", "S1PR3", "TIMP1", "HSPB1", "CXCL10", "CD44", "CP", "GFAP", "SERPINA3",
                                            "ASPG", "C3", "HLA-E", "HLA-F", "MICA", "H2-T23", "H2-D1", "AMIGO2", "SERPING1", "GGTA1P", "GGTA1", "UGT1A1",  "FKBP5", "PSMB8", "SRGN", "IIGP1", "S100A10", "CLCF1", "TGM1", "SPHK1", "CD109", "PTGS2", "SLC10A6", 
                                            "TM4SF1", "B3GNT5", "CD14")  

astrocyte.reactivity.labels.ac_sod1 <-  c("COL1A1", "MYL6", "FN1", "FLNA", "OSMR", "HLA-A", "HLA-B", "HLA-C", "STAT3", "RND3", "SLC26AA6",
                                          "PIMREG", "RRP8", "RECK", "LOXL3", "UBN1", "COL7A1", "MX1", "TTBK1", "TSC1", "MBD1", "EMP1", "PRPF4",   
                                          "PIF1", "AARS2", "PLSCR1", "NARPT", "CSRP1", "IRF7", "ID4", "AGER", "GBP2", "ANKZF1", "AGRN", "GTBP2",   
                                          "CCS", "TGFB1I1", "COL1A1A", "TNC", "RPL10", "NUP199", "LAMB2", "IDH1", "TP53", "LRP1", "TINF2", "ILK" ,    
                                          "FBLN5", "CAPN2", "CDKN1A", "ADAM19", "SERPINE1", "ITGA6", "ADAMTSL4", "ITGA5", "PDLIM7", "DIXDC1",  
                                          "BAG3", "PVR", "LCN2", "STEAP4", "S1PR3", "TIMP1", "HSPB1", "CXCL10", "CD44", "CP", "GFAP", "SERPINA3",
                                          "ASPG", "C3", "HLA-E", "HLA-F","H2-T23", "H2-D1", "AMIGO2", "SERPING1", "GGTA1P", "GGTA1", "UGT1A1",  
                                          "FKBP5", "PSMB8", "SRGN", "IIGP1", "S100A10", "CLCF1", "TGM1", "PTX3", "SPHK1", "CD109", "PTGS2", "SLC10A6", 
                                          "TM4SF1", "B3GNT5", "CD14")  
astrocyte.reactivity.labels.ac_a1 <- c("FN1", "ILK", "OSMR", "TINF2", "IRF7", "MCAM", "SOD2", "PSMB7", "PSMB10", "DNAJC2", "KMT2E", "COL27A1", "ERF")

astrocyte.reactivity.labels.ac_tdp43 <- c("FLOT1", "MCAM", "MBD1", "ANKZF1", "AARS2", "PRRT1", "RAB34", "DUSP1", "SPG11", "ABCD1", "CFP", "TRIP6", "ADAMTS20", "PTPN6", "ITGA10", "ADAMTS9", "BCAN", "NSMF", "MAP7", "PLOD3")

astrocyte.reactivity.labels.ac_protein <- c("FN1", "ILK", "OSMR", "TINF2", "IRF7", "MCAM", "SOD2", "PSMB7", "PSMB10", "DNAJC2", "KMT2E", "COL27A1", "ERF", "ATP6V1F", "RFTN1", "IDH1", "HNRNPU", "VCP", "ACTB", "COL5A2", "HSPB1", "PDLIM2", "PURA", "SVIL", "RAP1A", "TGFB1I1", "HLA-E", "ITGA3", "RELA", "TUBA1C")
RNA.capping <-c("POLR2J","MNAT1","POLR2B","NCBP3","POLR2E","POLR2F","RNMT","POLR2C","ERCC2","POLR2I","GTF2H1","GTF2H3","RNGTT","NCBP2","GTF2F1","CDK7","CCNH","NCBP1","CMTR1","TGS1","POLR2D","GTF2H2","POLR2K","ERCC3","POLR2H","POLR2G","RAMAC","POLR2L","CMTR2","POLR2A","GTF2F2","SUPT5H","GTF2H4","RAMACL","GTF2H5")
exon.junction.complex <- c("UPF1","THRAP3","SMG6","TDRD3","PNN","R3HCC1","CASC3","MAGOHB","UPF3B","SRSF1","EIF4A3","SAP18","UPF2","MAGOH","R3HCC1L","UPF3A","PYM1","RNPS1","RBM8A","SRRM1","PRPF8","EIF4A3","MLN51","MAGOH","Y14","CASC3","EIF4A3","MAGOH","RBM8A","PYM1","NXF1","DDX39B")
RNA.polyadenylation <- c("PAF1","ZC3H3","PABPC1","CPSF1","SNRPA","YTHDC1","PAPOLA","ZC3H14","PABPC1L","CSTF2","RNF40","CPEB3","MTPAP","CPSF6","TENT4A","PAPOLG","CPSF3","TENT4B","SYMPK","CCNT1","GRSF1","CDC73","WDR33","CDK9","PNPT1","APP","FIP1L1","TUT1","CPSF7","RNF20","SSU72","TENT2","VIRMA","PCF11","CPSF2","LEO1","NUDT21","AHCYL1","CLP1","CSTF3","HSF1","SUPT5H","CTR9","NELFE","SCAF8","CPEB1","PAPOLB","SSU72P8","AP002990.1","SSU72P4","SSU72P5","SSU72P2","SSU72P7","SSU72P3","NCBP2","NCBP1","","CPSF4","NUDT21","CSTF2T","CPSF4")
# Transport of Mature Transcript to Cytoplasm	NUP160	TPR	THOC3	ZC3H11A	NDC1	U2AF2	NUP133	CPSF1	NUP37	THOC1	NUP50	AAAS	NUP188	POLDIP3	THOC5	SRSF5	RAE1	LUZP4	NUP93	CASC3	NUP88	MAGOHB	NUP107	SRSF9	SRSF3	NUP155	NCBP2	SRSF7	SRSF4	SRSF11	CPSF3	GLE1	NUP43	FYTTD1	DDX39A	SRSF6	NUP153	UPF3B	NUP85	THOC2	SYMPK	NUP214	THOC6	NUP210	NXT1	SRRM1	NUP42	SRSF1	WDR33	NCBP1	NUP54	DHX38	EIF4A3	FIP1L1	EIF4E	RANBP2	NUP205	SEC13	U2AF1	CHTOP	CPSF4	U2AF1L4	SRSF2	NXF1	MAGOH	NUP35	THOC7	SLBP	SLU7	CPSF2	CDC40	POM121	DDX39B	SARNP	RNPS1	NUP62	RBM8A	NXF2	NXF2B	POM121C
minor.spliceosome.go <- c("RNU4ATAC","RNU6ATAC","PHF5A","SF3B1","SF3B2","SF3B3","SF3B4","SF3B5","SF3B6","SNRPB","SNRPD1","SNRPD2","SNRPD3","SNRPE","SNRPF","SNRPG","SNRPN","DHX15","PDCD7","RNPC3","RNU11","RNU12","SNRNP25","SNRNP35","SNRNP48","YBX1","ZCRB1","ZMAT5","ZRSR2","CD2BP2","DDX23","EFTUD2","PRPF6","PRPF8","RNU5A-1","SNRNP200","SNRNP40","TXNL4A")
minor.splicing.go <- unique(c(gprofiler_database$`mRNA Splicing - Minor Pathway`, minor.spliceosome.go))
RNA.capping <-c("POLR2J","MNAT1","POLR2B","NCBP3","POLR2E","POLR2F","RNMT","POLR2C","ERCC2","POLR2I","GTF2H1","GTF2H3","RNGTT","NCBP2","GTF2F1","CDK7","CCNH","NCBP1","CMTR1","TGS1","POLR2D","GTF2H2","POLR2K","ERCC3","POLR2H","POLR2G","RAMAC","POLR2L","CMTR2","POLR2A","GTF2F2","SUPT5H","GTF2H4","RAMACL","GTF2H5")
exon.junction.complex <- c("UPF1","THRAP3","SMG6","TDRD3","PNN","R3HCC1","CASC3","MAGOHB","UPF3B","SRSF1","EIF4A3","SAP18","UPF2","MAGOH","R3HCC1L","UPF3A","PYM1","RNPS1","RBM8A","SRRM1","PRPF8","EIF4A3","MLN51","MAGOH","Y14","CASC3","EIF4A3","MAGOH","RBM8A","PYM1","NXF1","DDX39B")
RNA.polyadenylation <- c("PAF1","ZC3H3","PABPC1","CPSF1","SNRPA","YTHDC1","PAPOLA","ZC3H14","PABPC1L","CSTF2","RNF40","CPEB3","MTPAP","CPSF6","TENT4A","PAPOLG","CPSF3","TENT4B","SYMPK","CCNT1","GRSF1","CDC73","WDR33","CDK9","PNPT1","APP","FIP1L1","TUT1","CPSF7","RNF20","SSU72","TENT2","VIRMA","PCF11","CPSF2","LEO1","NUDT21","AHCYL1","CLP1","CSTF3","HSF1","SUPT5H","CTR9","NELFE","SCAF8","CPEB1","PAPOLB","SSU72P8","AP002990.1","SSU72P4","SSU72P5","SSU72P2","SSU72P7","SSU72P3","NCBP2","NCBP1","","CPSF4","NUDT21","CSTF2T","CPSF4")
exosome.complex <- c("MTREX","EXOSC7","EXOSC5","DIS3","KHSRP","GTPBP1","EXOSC3","EXOSC8","EXOSC9","EXOSC2","WDR74","ZFC3H1","MPHOSPH6","PNPT1","NVL","DIS3L2","CARHSP1","SUPV3L1","DIS3L","EXOSC1","EXOSC10","EXOSC4","C1D","EXOSC6")
# response.to.typeIinterferon <- c("SP100","IFI35","IP6K2","UBE2K","MAVS","OAS1","MUL1","HSP90AB1","TTLL12","SAMHD1","CACTIN","TYK2","CDC37","OAS3","OAS2","PTPN6","WNT5A","STAT1","IRF6","IFIT3","IFIT2","IFNA6","IFNA8","EGR1","TRIM6","ZBP1","IRF1","IRF3","IFI6","IRF5","BST2","SHFL","XAF1","RSAD2","OASL","RNASEL","IFNA21","IRF4","NLRC5","IRF8","IFITM3","IFNAR1","IFNA5","IFNA16","MX1","IFNAR2","ADAR","JAK1","GBP2","DCST1","ABCE1","METTL3","IFI27","FADD","IRF2","TRIM56","STAT2","IFNB1","ISG20","MYD88","PTPN2","PTPN11","SETD2","SHMT2","MX2","TBK1","IRAK1","USP18","IFITM2","IRF7","YTHDF3","IFIT1","IFITM1","IFNA10","ISG15","IFNA2","PTPN1","IFNA1","YTHDF2","CNOT7","PSMB8","HLA-C","HLA-E","HLA-G","HLA-F","HLA-A","TREX1","IRF9","IFNA7","IFNA14","IFNA13","HLA-B","IFNA17","IFNA4","LSM14A","MMP12","IKBKE","MIR21")	
# response.to.typeIinterferon <- c("SP100","IFI35","IP6K2","UBE2K","MAVS","OAS1","MUL1","HSP90AB1","TTLL12","SAMHD1","CACTIN","TYK2","CDC37","OAS3","OAS2","PTPN6","WNT5A","STAT1","IRF6","IFIT3","IFIT2","IFNA6","IFNA8","EGR1","TRIM6","ZBP1","IRF1","IRF3","IFI6","IRF5","BST2","SHFL","XAF1","RSAD2","OASL","RNASEL","IFNA21","IRF4","NLRC5","IRF8","IFITM3","IFNAR1","IFNA5","IFNA16","MX1","IFNAR2","ADAR","JAK1","GBP2","DCST1","ABCE1","METTL3","IFI27","FADD","IRF2","TRIM56","STAT2","IFNB1","ISG20","MYD88","PTPN2","PTPN11","SETD2","SHMT2","MX2","TBK1","IRAK1","USP18","IFITM2","IRF7","YTHDF3","IFIT1","IFITM1","IFNA10","ISG15","IFNA2","PTPN1","IFNA1","YTHDF2","CNOT7","PSMB8","HLA-C","HLA-E","HLA-G","HLA-F","HLA-A","TREX1","IRF9","IFNA7","IFNA14","IFNA13","HLA-B","IFNA17","IFNA4","LSM14A","MMP12","IKBKE","MIR21")	
# p53.signalling.pathway = unique(c(gprofiler_database$`TP53 Network`, go.pathways.msigdb$GO_P53_BINDING, go.pathways.msigdb$GO_POSITIVE_REGULATION_OF_SIGNAL_TRANSDUCTION_BY_P53_CLASS_MEDIATOR, go.pathways.msigdb$GO_POSITIVE_REGULATION_OF_DNA_DAMAGE_RESPONSE_SIGNAL_TRANSDUCTION_BY_P53_CLASS_MEDIATOR, go.pathways.msigdb$GO_POSITIVE_REGULATION_OF_INTRINSIC_APOPTOTIC_SIGNALING_PATHWAY_BY_P53_CLASS_MEDIATOR))
p53.signalling.pathway = unique(c(go.pathways.msigdb$GO_SIGNAL_TRANSDUCTION_BY_P53_CLASS_MEDIATOR, kegg.pathways.msigdb$KEGG_P53_SIGNALING_PATHWAY))
p53.DDR.repair_gene_set = unique(c(p53.signalling.pathway, gprofiler_database$`DNA Damage Response`,gprofiler_database$`DNA repair`))


# RNAInter database
# sbatch -N 1 -c 4 --mem=10G -t 1:00:00 --wrap="wget 'http://www.rna-society.org/rnainter/download/RNA-Protein.zip'" --mail-type=ALL,ARRAY_TASKS --mail-user=oliver.ziff@crick.ac.uk --job-name=RNAinter
# sbatch -N 1 -c 4 --mem=10G -t 2:00:00 --wrap="wget 'http://www.rna-society.org/rnainter/download/Full%20Version.zip'" --mail-type=ALL,ARRAY_TASKS --mail-user=oliver.ziff@crick.ac.uk --job-name=RNAinterFULL
# sbatch -N 1 -c 4 --mem=10G -t 1:00:00 --wrap="unzip 'Full Version.zip'" --mail-type=ALL,ARRAY_TASKS --mail-user=oliver.ziff@crick.ac.uk --job-name=RNAinter
# sbatch -N 1 -c 4 --mem=10G -t 1:00:00 --wrap="unzip 'RNA-Protein.zip'" --mail-type=ALL,ARRAY_TASKS --mail-user=oliver.ziff@crick.ac.uk --job-name=RNAinter
# RNAInter_interaction_full <- read.table("/camp/lab/luscomben/home/users/ziffo/genomes/rbp-databases/rnainter/RNAInter_interaction_full.txt", header=FALSE, sep = "\t")
# RNAInter_RNA_Protein.txt <- read.table("/camp/lab/luscomben/home/users/ziffo/genomes/rbp-databases/rnainter/RNA-Protein.txt", header=FALSE, sep = "\t") %>% as_tibble() %>% select(gene_name = V2, gene_id = V3, biotype = V4, gene.organism = V5, protein_name = V6, protein_id = V7, protein_type = V8, protein.organism = V9, score = V10)
# RNAInter_RNA_Protein.human <- RNAInter_RNA_Protein.txt %>% filter(gene.organism == "Homo sapiens", protein.organism == "Homo sapiens")
# saveRDS(RNAInter_RNA_Protein.human, "/camp/lab/luscomben/home/users/ziffo/genomes/rbp-databases/rnainter/rna_protein.human.rds")
# RNAInter_RNA_Protein.human <- readRDS("/camp/lab/luscomben/home/users/ziffo/genomes/rbp-databases/rnainter/rna_protein.human.rds")
# Sys.setenv("http_proxy" = "http://my.proxy.org:9999") # options(RCurlOptions = list(proxy="uscache.kcc.com:80",proxyuserpwd="------:-------"))
# httr::set_config(httr::config(ssl_verifypeer = FALSE)) # httr::set_config(httr::config(ssl_cipher_list = "DEFAULT@SECLEVEL=1"))
# ensembl99 <- useEnsembl(biomart = 'genes', dataset = 'hsapiens_gene_ensembl', version = 99) # Select biomaRt database (genes) and dataset (hsapiens)
# ens2entrez <- getBM(attributes =  c("ensembl_gene_id",  "entrezgene_id"),values = TRUE,mart = ensembl99,useCache = FALSE) %>% select(ensemblID = ensembl_gene_id, entrezID = entrezgene_id) # # NCBI entrez gene ID to Ensembl Gene ID
# gtf <- gtf %>% left_join(ens2entrez)
# RNAInter_RNA_Protein <- RNAInter_RNA_Protein.human %>%  filter(!gene_name %in% gtf$gene_name | gene_id %in% gtf$ensemblID | gene_id %in% gtf$entrezID) # 1,990,466 interactions can be matched based on gene identifiers (many missing e.g. ENSG00000124835)
# RNAInter_RNA_RBP <- RNAInter_RNA_Protein.human %>%  filter(protein_name %in% RNA_BINDING_PROTEINS) %>% # 22,501,074 interactions, 1,471 RBPs
# # filter(protein_type == "RBP")  # some important RBPs are annotated as TFs eg FUS and SFPQ so cant rely on protein_type == 192 RBPs and 36,280 genes
# filter(!gene_name %in% gtf$gene_name | gene_id %in% gtf$ensemblID | gene_id %in% gtf$entrezID) # 1,990,466 interactions can be matched based on gene identifiers (many missing e.g. ENSG00000124835)
# saveRDS(RNAInter_RNA_RBP, "/camp/lab/luscomben/home/users/ziffo/genomes/rbp-databases/rnainter/rna_rbp.human.rds")
# RNAInter_RNA_RBP <- readRDS("/camp/lab/luscomben/home/users/ziffo/genomes/rbp-databases/rnainter/rna_rbp.human.rds")


# Generic Functions ------------------
make_matrix <- function(df,rownames = NULL){
  my_matrix <-  as.matrix(df)
  if(!is.null(rownames))
    rownames(my_matrix) = rownames
  my_matrix
}

tibble_to_matrix <- function(tbl, ..., row_names=NULL){
  cols <- rlang::enquos(...)
  mat <- as.matrix(dplyr::select(tbl, !!! cols))
  if (!is.null(row_names)){
    if (length(row_names) == 1 & row_names[1] %in% colnames(tbl)){
      row_names <- tbl[[row_names]]
    }
    rownames(mat) <- row_names
  }
  return(mat)
}

scale2 <- function(x, na.rm = FALSE) (x - mean(x, na.rm = na.rm)) / sd(x, na.rm)

my_summarise <- function(data, var, ...) {
  data %>% 
    group_by(...) %>%
    summarise(
      "mean_{{var}}" := mean({{ var }}, na.rm = TRUE), "median_{{var}}" := median({{ var }}, na.rm = TRUE), 
      "sd_{{var}}" := sd({{ var }}, na.rm = TRUE), "Q1_{{var}}" := quantile({{ var }}, probs = 0.25), "Q3_{{var}}" := quantile({{ var }}, probs = 0.75),
      "min_{{var}}" := min({{ var }}, na.rm = TRUE), "max_{{var}}" := max({{ var }}, na.rm = TRUE), "sum_{{var}}" := sum({{ var }}, na.rm = TRUE),
      "n_{{var}}" := n(),
      "se_{{var}}" := plotrix::std.error({{ var }})
    )
}


left_join_NA <- function(x, y, ...) {
  left_join(x = x, y = y, by = ...) %>% 
    mutate_each(funs(replace(., which(is.na(.)), 0)))
}



# transpose_tibble <- function(tbl, col_names, id_col = "columns"){
#   col_names <- rlang::enquo(col_names)
#   tibble_to_matrix(tbl, -!!col_names, row_names = dplyr::pull(tbl, !!col_names)) %>%
#     t() %>% dplyr::as_tibble(rownames = id_col) %>% return()
# }

transpose_tibble <- function(tbl,old_rownames,new_rownames){
  tbl %>% column_to_rownames(old_rownames) %>% t() %>% as_tibble(rownames = new_rownames) %>% return()
}
# transpose_tibble(vsd_counts.phag_glut, old_rownames = "gene_name", new_rownames = "group")



ens_to_gene <- function(ens) {
  gene <- replace(ens, ens %in% as.character(gene2ens$gene_id), as.character(gene2ens$gene_name[gene2ens$gene_id %in% ens]))
  # gene <- gene2ens %>% dplyr::filter(gene_id == ens) %>% dplyr::pull(gene_name) %>% as.character()
  return(gene)
}
gene_to_ens <- function(gene) {
  ens <- gene2ens %>% dplyr::filter(gene_name == gene) %>% dplyr::pull(gene_id) %>% as.character()
  print(ens)
}
ens_to_gene_list <- function(ens) {
  gene <- list()
  for (i in seq_along(ens)) {
    gene[[i]] <- replace(ens[[i]], ens[[i]] %in% as.character(gene2ens$gene_id), as.character(gene2ens$gene_name[gene2ens$gene_id %in% ens[[i]]]))
  }
  return(gene)
}

# Sys.setenv("http_proxy" = "http://my.proxy.org:9999") # Configure biomaRt for CAMP. NB this will stop gprofiler2 working
# Sys.unsetenv("http_proxy") # unset proxy to fix gprofiler2
# httr::set_config(httr::config(ssl_verifypeer = FALSE)) # httr::set_config(httr::config(ssl_cipher_list = "DEFAULT@SECLEVEL=1"))

# human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
# mouse2human = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = mus_musculus.gtf$gene_name, mart = mouse, attributes = c("hgnc_symbol"), martL = human, uniqueRows=T) # useCache = FALSE solves no such table: metadata error
# mouse2human = readRDS(here(camp_path,"/genomes/ensembl/mouse2human.rds"))

# Convert mouse to human gene names https://www.r-bloggers.com/2016/10/converting-mouse-to-human-gene-names-with-biomart-package/
# convertMouseGeneList <- function(x){
#   require("biomaRt")
#   human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
#   mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
#   genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = mus_musculus.gtf$gene_name, mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
#   humanx <- unique(genesV2[, 2])
#   # Print the first 6 genes found to the screen
#   print(head(humanx))
#   return(humanx)
# }

# genbank to ensembl
# genbank = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# genbank2ensembl = getLDS(filters = "refseq_mrna", values = mus_musculus.gtf$gene_name, mart = mouse, attributes = c("hgnc_symbol"), martL = human, uniqueRows=T)
# searchAttributes(mart = human, pattern = "refseq")




theme_oz <- function (xaxis_arrow = TRUE) { 
  if(xaxis_arrow == TRUE) {
  theme_bw(base_size=8) %+replace% # theme_bw(base_size=8, base_family="Helvetica") %+replace% 
    theme(
      panel.grid = element_blank(),
      strip.background = element_blank(),
      #panel.border = element_blank(),
      axis.line = element_line(arrow = arrow(length = unit(2,"mm"))),
      # axis.line.y = element_line(arrow = arrow(length = unit(2,"mm"))),
      #text = element_text(color = "black"), 
      strip.text = element_text(color = "black"),
      axis.text = element_text(colour = "black"),
      axis.ticks = element_line(colour = "black"),
      panel.background  = element_blank(),
      plot.background = element_rect(fill="white", colour=NA), 
      legend.background = element_rect(fill="transparent", colour=NA),
      legend.key = element_rect(fill="transparent", colour=NA), 
      strip.text.x = element_text(face = "bold", margin = margin(t = 2,r = 0,b = 2,l=0)))
  } else if (xaxis_arrow == FALSE) {
    theme_bw(base_size=8) %+replace% # theme_bw(base_size=8, base_family="Helvetica") %+replace% 
      theme(
        panel.grid = element_blank(),
        strip.background = element_blank(),
        #panel.border = element_blank(),
        # axis.line = element_line(arrow = arrow(length = unit(2,"mm"))),
        axis.line.y = element_line(arrow = arrow(length = unit(2,"mm"))),
        #text = element_text(color = "black"), 
        strip.text = element_text(color = "black"),
        axis.text = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        panel.background  = element_blank(),
        plot.background = element_rect(fill="white", colour=NA), 
        legend.background = element_rect(fill="transparent", colour=NA),
        legend.key = element_rect(fill="transparent", colour=NA), 
        strip.text.x = element_text(face = "bold", margin = margin(t = 2,r = 0,b = 2,l=0)))
  }
}

blank_theme <- function () { theme_minimal() + theme(axis.title.x = element_blank(), axis.title.y = element_blank(), panel.border = element_blank(), panel.grid=element_blank(), axis.ticks = element_blank(), plot.title=element_text(size=14, face="bold")) + theme(axis.text.x=element_blank()) }


# Gene Expression Functions ------------------
DESeq.analysis <- function(metadata, design = ~1, contrast = NULL, combine.contrasts = FALSE, LRT_full = NULL, LRT_reduced = NULL, join_multiqc = TRUE, species = "human", run_deseq = TRUE, 
                           unique_names = "sample", salmon_names = "sample", irfinder_names = "sample", # each sample may have a different filename for salmon and irfinder. For the dataframe, sample rownames need to be unique (unique_names)
                           gene.level = TRUE, transcript.level = FALSE, vsd_tx = FALSE,run_irfinder = FALSE, irfinder_design = NULL, irfinder_contrast_variable = "condition", irfinder_contrast_name = NULL,
                           plot.path.prefix = NULL, plot.figures = FALSE){
  
  # check unique_names are unique
  metadata$duplicated.names = metadata %>% select({{ unique_names }}) %>% duplicated() 
  if (TRUE %in% metadata$duplicated.names){
    duplicated.names = metadata %>% filter(duplicated.names == TRUE) %>% pull({{ unique_names }}) # print duplicated names
    print(duplicated.names)
    stop("unique_names variable has these duplicated names in")
  }
  
  # join multiqc_general_stats.txt
  if(join_multiqc == FALSE){
    print("join_multiqc == FALSE. Skipping join.")
  } else if (FALSE %in% file.exists(paste0(metadata$database_dir,"/nfcore/multiqc/star_salmon/multiqc_data/multiqc_general_stats.txt"))){
    print("multiqc general stats not found for all samples. Skipping join.")
  } else {    
    print("joining multiqc general stats to metadata")
    multiqc_stat_paths = metadata %>% distinct(database_dir) %>% mutate(multiqc_stat_paths = paste0(database_dir,"/nfcore/multiqc/star_salmon/multiqc_data/multiqc_general_stats.txt")) %>% pull(multiqc_stat_paths)
    multiqc_stats.tsv = multiqc_stat_paths %>% map(read_tsv, show_col_types = FALSE)
    names(multiqc_stats.tsv) =  metadata %>% distinct(database_dir,.keep_all=TRUE) %>% pull(database_dir)
    multiqc_data <- map_dfr(multiqc_stats.tsv, bind_rows, .id = "database_dir") %>% clean_names() %>% drop_na(samtools_mqc_generalstats_samtools_flagstat_total)
    metadata = metadata %>% left_join(multiqc_data) # need unique sample names - only works if sample & dataset are unique. If multiple dataset names within a single multiqc_general_stats.txt then will only join the first dataset name
  }
  
  if(run_deseq == TRUE){
    # salmon files
    metadata.salmon_files = metadata %>% mutate(file_salmon = file.path(database_dir, "nfcore/star_salmon", metadata[[salmon_names]], "quant.sf")) %>% pull(file_salmon)
    names(metadata.salmon_files) = metadata %>% pull({{unique_names}}) # sample filename in nfcore outdir
    rownames(metadata) <- metadata %>% pull({{unique_names}})  # unique sample name
    if (FALSE %in% file.exists(metadata.salmon_files)){
      missing.bams = metadata.salmon_files[file.exists(metadata.salmon_files) == FALSE] # print any missing files to console, these will cause tximport to fail.
      print(missing.bams)
      stop("these count files do not exist")
    }
    
    if(species == "human"){
      cat(blue("using human gene names\n"))
      transcript_gene=tx2gene
      gene_ensembl=gene2ens
    } else if(species == "mouse"){
      cat(blue("using mouse gene names\n"))
      transcript_gene=mus.musculus.tx2gene
      gene_ensembl=mus.musculus.gene2ens
    }
    
    if(gene.level == TRUE){
      if(is.null(LRT_full) == TRUE){
        cat(red(paste("Running DESeq2 with Wald test: Design =",design,". Results Contrast =", contrast, ".\n")))
        dds <- tximport(metadata.salmon_files, type="salmon", tx2gene=transcript_gene) %>% DESeqDataSetFromTximport(colData = metadata, design = design) %>% DESeq(test = "Wald")
      } else if (!is.null(LRT_full)){
        cat(red(paste("Running DESeq2 with Likelihood ratio test: Full model =", LRT_full, ". Reduced model =", LRT_reduced, ". Results Contrast =", contrast, ".\n")))
        dds <- tximport(metadata.salmon_files, type="salmon", tx2gene=transcript_gene) %>% DESeqDataSetFromTximport(colData = metadata, design = LRT_full) %>% DESeq(test = "LRT", reduced = LRT_reduced)
        res <- DESeq2::results(dds) %>% as_tibble(rownames = "gene_id") %>% left_join(gene_ensembl) %>% arrange(-abs(stat)) 
      }
      vsd <- vst(dds, blind=TRUE)
      # vsd.counts <- as_tibble(assay(vsd), rownames = "gene_id") %>% left_join(gene_ensembl)
      print(resultsNames(dds))
      
      if(length(contrast) == 1 | combine.contrasts == TRUE) {
        if(!is.null(contrast) & combine.contrasts == FALSE) { 
          res <- DESeq2::results(dds, name = contrast) %>% as_tibble(rownames = "gene_id") %>% left_join(gene_ensembl) %>% arrange(-abs(stat)) 
          # lfcshrink <- DESeq2::lfcShrink(dds, coef = contrast, type = "ashr") %>% as_tibble(rownames = "gene_id") %>% left_join(gene_ensembl) %>% arrange(pvalue) 
        } else if (!is.null(contrast) & combine.contrasts == TRUE) {
          res <- DESeq2::results(dds, contrast = contrast) %>% as_tibble(rownames = "gene_id") %>% left_join(gene_ensembl) %>% arrange(-abs(stat)) 
          # lfcshrink <- DESeq2::lfcShrink(dds, coef = contrast, type = "ashr") %>% as_tibble(rownames = "gene_id") %>% left_join(gene_ensembl) %>% arrange(pvalue) 
        }
        if(!is.null(contrast)){ cat(blue(paste("Significant events: padj < 0.05 = ", nrow(filter(res, padj < 0.05)), ", pvalue < 0.05 = ", nrow(filter(res, pvalue < 0.05)),"\n"))) }
      } else if(length(contrast) > 1 & combine.contrasts == FALSE){
        cat(green(paste(length(contrast), "contrasts specified")))
        res <- list() # lfcshrink <- list()
        for(i in seq_along(contrast)) {
          res[[i]] <- DESeq2::results(dds, name = contrast[[i]]) %>% as_tibble(rownames = "gene_id") %>% left_join(gene_ensembl) %>% arrange(-abs(stat)) 
          names(res)[i] <- contrast[[i]]
          # lfcshrink[[i]] <- DESeq2::lfcShrink(dds, coef = contrast[[i]], type = "ashr") %>% as_tibble(rownames = "gene_id") %>% left_join(gene_ensembl) %>% arrange(pvalue) 
          # names(lfcshrink)[i] <- contrast[[i]]
        }
      }
    }
    if(transcript.level == TRUE){
      cat(blue("Running transcript level analysis\n"))
      dds.tx <- tximport(metadata.salmon_files, type="salmon", txOut = TRUE) %>% DESeqDataSetFromTximport(colData = metadata, design = design) %>% DESeq()
      if(length(contrast) == 1 | combine.contrasts == TRUE){
        if(!is.null(contrast) & combine.contrasts == FALSE){ res.tx <- DESeq2::results(dds.tx, name = contrast) %>% as_tibble(rownames = "transcript_id") %>% left_join(transcript_gene, by="transcript_id") %>% left_join(gene_ensembl) %>% arrange(-abs(stat)) 
        } else if (!is.null(contrast) & combine.contrasts == TRUE) {res.tx <- DESeq2::results(dds.tx, contrast = contrast) %>% as_tibble(rownames = "transcript_id") %>% left_join(transcript_gene, by="transcript_id") %>% left_join(gene_ensembl) %>% arrange(-abs(stat)) }
      } else if(length(contrast) > 1 & combine.contrasts == FALSE){
        res.tx <- list()
        for(i in seq_along(contrast)) {
          res.tx[[i]] <- DESeq2::results(dds.tx, name = contrast[[i]]) %>% as_tibble(rownames = "transcript_id") %>% left_join(tx2gene, by="transcript_id") %>% arrange(-abs(stat)) 
          names(res.tx)[i] <- contrast[[i]]
        }
      }
      if(vsd_tx == TRUE){ vsd.tx <- vst(dds.tx, blind=TRUE) }
    }
  }
  if(run_irfinder == TRUE){
    cat(blue("Gene Expression Complete\nRunning IRFinder...\n"))
    irfinder = IRFinder.analysis(metadata = metadata, sample.names = irfinder_names, contrast_variable = irfinder_contrast_variable, design = irfinder_design, contrast = irfinder_contrast_name, ge.res = res, animal = species)
  }
  if(!is.null(plot.path.prefix)){
   write_csv(res, file = paste0(metadata$database_dir[1],"/expression/deseq2/results/",plot.path.prefix,"_deseq_results.csv"), na = "")
    if(plot.figures == TRUE){
      cat(blue("Gene Expression Complete\nPlotting figures...\n"))
      # PCA
      pcaData <- plotPCA(vsd, intgroup=c("sample", "database_dir", "condition"), returnData=TRUE) %>% as_tibble()
      pca <- ggplot(pcaData, aes(PC1, PC2, colour=condition)) + geom_point(size=1) +
        theme_oz() + theme(legend.text=element_text(size=8), legend.title=element_text(size=7), legend.position = "top", legend.spacing = unit(0.01, 'cm')) + 
        scale_colour_manual(values = c("firebrick2", "dodgerblue2", "forestgreen", "gold2")) + 
        xlab(paste0("PC1: ",round(100 * attr(pcaData, "percentVar"))[1],"% variance")) +  ylab(paste0("PC2: ",round(100 * attr(pcaData, "percentVar"))[2],"% variance")) + #coord_fixed() +
        ggforce::geom_mark_ellipse(aes(label = condition), show.legend = FALSE) + geom_text_repel(aes(label=sample), size = 2.4, max.overlaps = 50) +
        guides(color=guide_legend(title="", keywidth=0.1, keyheight=0.15, default.unit="inch", reverse = TRUE, nrow = 1)) + ggtitle("")
      ggsave(pca, filename = paste0(metadata$database_dir[1],"/expression/deseq2/figures/",plot.path.prefix,"_pca.pdf"), height = 6, width = 6)
      
      # Volcano
      volcano.als.labels <- res %>% filter(padj < 0.05, gene_name %in% c(ALS_genes, RNA_BINDING_PROTEINS, p53.signalling.pathway, gprofiler_database$`DNA Damage Response`)) %>% top_n(20, wt = -pvalue) %>% distinct(gene_name) %>% pull(gene_name)
      volcano.labels <- res %>% filter(!grepl("^AC|^AL", gene_name)) %>% top_n(20, wt = -pvalue) %>% pull(gene_name)
      volcano <- volcano_plot(res, numerator = levels(metadata$condition)[2], denomenator = levels(metadata$condition)[1], title = plot.path.prefix, dpi = 300, gene_labels = c(volcano.als.labels, volcano.labels), distinct_labels = TRUE)
      ggsave(volcano, filename = paste0(metadata$database_dir[1],"/expression/deseq2/figures/",plot.path.prefix,"_volcano.pdf"), height = 6, width = 6)
      
      # GO terms
      if(nrow(filter(res, padj < 0.05))>2) {
        cat(blue("GO terms...\n"))
        if(species =="human"){ gost <- res %>% filter(padj < 0.05) %>% group_split(stat > 0) %>% map("gene_id") %>% gost(organism = "hsapiens", sources = c("GO", "KEGG", "REAC", "CORUM")) } else if (species == "mouse"){ gost <- res %>% filter(padj < 0.05) %>% group_split(stat > 0) %>% map("gene_id") %>% gost(organism = "mmusculus", sources = c("GO", "KEGG", "REAC", "CORUM")) }
          if(!is.null(gost)) {
          gost_filt <- gost$result %>% arrange(-log10(p_value)) %>% mutate(query = case_when(query == "query_1" ~ levels(metadata$condition)[1], query == "query_2" ~ levels(metadata$condition)[2]))
          # gost_query1 <- gost_filt %>% filter(query == levels(metadata$condition)[1]) %>% top_n(10, wt = -p_value)
          # gost_query2 <- gost_filt %>% filter(query == levels(metadata$condition)[2]) %>% top_n(10, wt = -p_value)
          # if(!is.null(gost_query1) & !is.null(gost_query2)) {
          #   go_terms <- plot_grid(curated_profiles(gprofiles = gost_query2, title = plot.path.prefix) + xlab(""), curated_profiles(gprofiles = gost_query1, cols = "dodgerblue2"), nrow = 2, align = "v")
          #   ggsave(go_terms, filename = paste0(metadata$database_dir[1],"/expression/deseq2/figures/",plot.path.prefix,"_go_terms.pdf"), height = 8, width = 6)
          # }
          write_csv(gost_filt, file = paste0(metadata$database_dir[1],"/expression/deseq2/results/",plot.path.prefix,"_go_terms.csv"), na = "")
        }
      }
      
      # GSEA https://yulab-smu.top/biomedical-knowledge-mining-book/universal-api.html
      cat(blue("GSEA...\n"))
      library(msigdbr)
      msigdbr_go = msigdbr(species = "Homo sapiens", category = "C5") # C5 = GO pathways       # msigdbr_hs = msigdbr(species = "Homo sapiens")
      msigdbr_go.t2g = msigdbr_go %>% select(gs_name, ensembl_gene) %>% filter(!grepl("^HP_",gs_name))
      geneList <- res %>% drop_na(stat) %>% pull(stat)
      names(geneList) <- res %>% drop_na(stat) %>% pull(gene_id) %>% as.character()
      geneList = sort(geneList, decreasing = TRUE)
      gsea <- GSEA(geneList, TERM2GENE = msigdbr_go.t2g) %>% as_tibble() %>% arrange(-abs(NES)) 
      write_csv(as_tibble(gsea), file = paste0(metadata$database_dir[1],"/expression/deseq2/results/",plot.path.prefix,"_gsea.csv"), na = "")
      
      # Progeny
      cat(blue("decoupleR...\n"))
      if(species =="human"){ net_progeny <- get_progeny(organism = 'human', top = 100) } else if (species == "mouse"){ net_progeny <- get_progeny(organism = 'mouse', top = 100) }
      res.matrix <- res %>% select(gene_name, stat) %>% drop_na(stat) %>% distinct(gene_name, .keep_all = TRUE) %>% tibble_to_matrix(stat, row_names = "gene_name")
      progeny.contrast_acts <- run_wmean(mat=res.matrix, net=net_progeny, .source='source', .target='target', .mor='weight', times = 1000, minsize = 5) %>% filter(statistic == 'norm_wmean') %>% mutate(p.signif = case_when(p_value < 0.003 ~ "***", p_value < 0.009 ~ "***",p_value < 0.01 ~ "**",p_value < 0.05 ~ "*", TRUE ~ ""), group1 = source, group2 = source)
      progeny_pathways <- ggplot(progeny.contrast_acts, aes(x=reorder(source, score), y = score)) + geom_bar(aes(fill = score), stat = "identity") +
        scale_fill_gradient2(low = "darkblue", high = "indianred", mid = "whitesmoke", midpoint = 0) + 
        theme_oz(xaxis_arrow = FALSE) + theme(axis.text.x = element_text(angle = 90), legend.position = "none", plot.title = element_text(hjust = 0.5)) +
        labs(title = paste("Signaling Pathways:", plot.path.prefix), x = "", y = "Normalised enrichment") +
        stat_pvalue_manual(progeny.contrast_acts, label = "p.signif", y.position = max(progeny.contrast_acts$score), xmin = "source", xmax = NULL, size = 4, hide.ns = TRUE) 
      ggsave(progeny_pathways, filename = paste0(metadata$database_dir[1],"/expression/deseq2/figures/",plot.path.prefix,"_progeny_pathways.pdf"), height = 6, width = 6)
      
      # Dorothea
      if(species =="human"){ net_dorothea <- get_dorothea(organism='human', levels=c('A', 'B', 'C')) } else if (species == "mouse"){ net_dorothea <- get_dorothea(organism='mouse', levels=c('A', 'B', 'C')) }
      dorothea.contrast_acts <- run_wmean(mat=res.matrix, net=net_dorothea, .source='source', .target='target', .mor='mor', times = 1000, minsize = 5) %>% filter(statistic == 'norm_wmean') %>% 
        mutate(p.signif = case_when(p_value < 0.05 ~ "*", TRUE ~ ""), sig = case_when(p_value < 0.05 & abs(score) < 4 ~ "sig", p_value < 0.05 & abs(score) >= 4 ~ "sig_strong", p_value >= 0.05 ~ "non_sig"), direction = ifelse(score > 0, "up", "down"), class = paste(sig, direction))
      dorothea.contrast_acts.labels.up = dorothea.contrast_acts %>% top_n(score, n = 5) %>% pull(source)
      dorothea.contrast_acts.labels.down = dorothea.contrast_acts %>% top_n(-score, n = 5) %>% pull(source)
      dorothea_volcano = ggplot(dorothea.contrast_acts, aes(x = score, y = -log10(p_value))) + geom_point(aes(colour = class)) +
        scale_colour_manual(values = c("non_sig up" = "gray", "non_sig down" = "gray",  "sig up" = "#F36F6F", "sig_strong up" = "#EB4445", "sig down" = "#A6CEE3", "sig_strong down" = "#79B1D3")) +
        labs(y = expression(-log[10]~P~value), x = "Normalised enrichment", title = paste("Transcription Factors:",plot.path.prefix)) + guides(colour = "none") +
        theme_oz() +  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5), panel.border = element_blank(), axis.ticks = element_line(colour = "black") ) +
        geom_text_repel(fontface = "italic", data = filter(dorothea.contrast_acts, source %in% c(dorothea.contrast_acts.labels.up, dorothea.contrast_acts.labels.down)), aes(label = source), max.overlaps = Inf, min.segment.length = unit(0, "lines"), size = 2.5) +  
        geom_vline(xintercept = 0, linetype = 'dotted') + ylim(c(0,2.8))
      ggsave(dorothea_volcano, filename = paste0(metadata$database_dir[1],"/expression/deseq2/figures/",plot.path.prefix,"_dorothea_tf.pdf"), height = 6, width = 6)
      write_csv(dorothea.contrast_acts, file = paste0(metadata$database_dir[1],"/expression/deseq2/results/",plot.path.prefix,"_dorothea_results.csv"), na = "")
    }
  }
  
  cat(blue("Complete. Returning results...\n"))
  if(gene.level == TRUE & transcript.level == TRUE & run_irfinder == TRUE & (!is.null(contrast) | !is.null(LRT_full))){ return(list(dds=dds,vsd=vsd,res=res,res.tx=res.tx,irfinder=irfinder$irfinder,ir.ge=irfinder$ir.ge)) # lfcshrink=lfcshrink,,vsd.counts=vsd.counts
  } else if(gene.level == TRUE & transcript.level == TRUE & run_irfinder == FALSE & (!is.null(contrast) | !is.null(LRT_full))){ return(list(dds=dds,vsd=vsd,res=res,res.tx=res.tx)) #lfcshrink=lfcshrink,,vsd.counts=vsd.counts
  } else if(gene.level == TRUE & transcript.level == FALSE & run_irfinder == TRUE & (!is.null(contrast) | !is.null(LRT_full))){ return(list(dds=dds, vsd=vsd, res=res, irfinder=irfinder$irfinder, ir.ge=irfinder$ir.ge)) # lfcshrink=lfcshrink,,vsd.counts=vsd.counts
  } else if(gene.level == TRUE & transcript.level == FALSE & run_irfinder == FALSE & (!is.null(contrast) | !is.null(LRT_full))){ return(list(dds=dds, vsd=vsd, res=res)) # lfcshrink=lfcshrink,,vsd.counts=vsd.counts
  } else if(gene.level == TRUE & transcript.level == FALSE & run_irfinder == TRUE & is.null(contrast) == TRUE){ return(list(dds=dds,vsd=vsd,irfinder=irfinder$irfinder,ir.ge=irfinder$ir.ge)) #,vsd.counts=vsd.counts
  } else if(gene.level == TRUE & transcript.level == FALSE & run_irfinder == FALSE & is.null(contrast) == TRUE){ return(list(dds=dds,vsd=vsd)) #,vsd.counts=vsd.counts
  } else if(gene.level == FALSE & transcript.level == TRUE & vsd_tx == FALSE & run_irfinder == FALSE & (!is.null(contrast) | !is.null(LRT_full))){ return(list(res.tx=res.tx)) 
  } else if(gene.level == FALSE & transcript.level == TRUE & vsd_tx == TRUE & run_irfinder == FALSE & (!is.null(contrast) | !is.null(LRT_full))){ return(list(vsd.tx=vsd.tx,res.tx=res.tx)) 
  } else if(gene.level == FALSE & transcript.level == TRUE & vsd_tx == TRUE & run_irfinder == FALSE & is.null(contrast) == TRUE){ return(list(vsd.tx=vsd.tx)) }
}



gprofiler2rrvgo <- function(deseq.res, query1 = "Down", query2 = "Up", query3 = NULL, query4 = NULL, GO_cat = "BP", sig.filter = "padj", species = "human"){
  print(paste0("filtering significant genes with ", sig.filter))
  
  if(species == "human"){
    if(sig.filter == "padj"){gost.res <- deseq.res %>% filter(padj < 0.05) %>% group_split(log2FoldChange > 0) %>%  map("gene_id") %>% gost(organism = "hsapiens")
    } else if(sig.filter == "pvalue"){gost.res <- deseq.res %>% filter(pvalue < 0.05) %>% group_split(log2FoldChange > 0) %>%  map("gene_id") %>% gost(organism = "hsapiens")
    } else { gost.res <- deseq.res %>% filter(significant == "both") %>% group_by(direction) %>% distinct(gene_id, .keep_all = TRUE) %>% ungroup %>% group_split(direction) %>% map("gene_id") %>% gost(organism = "hsapiens") }
  } else if(species == "mouse") {
    if(sig.filter == "padj"){gost.res <- deseq.res %>% filter(padj < 0.05) %>% group_split(log2FoldChange > 0) %>%  map("gene_id") %>% gost(organism = "mmusculus")
    } else if(sig.filter == "pvalue"){gost.res <- deseq.res %>% filter(pvalue < 0.05) %>% group_split(log2FoldChange > 0) %>%  map("gene_id") %>% gost(organism = "mmusculus")
    } else { gost.res <- deseq.res %>% filter(significant == "both") %>% group_by(direction) %>% distinct(gene_id, .keep_all = TRUE) %>% ungroup %>% group_split(direction) %>% map("gene_id") %>% gost(organism = "mmusculus") }
  }
  if(is.null(query2) & is.null(query3) & is.null(query4)) {gost.res_filt <- gost.res$result %>% filter(str_detect(source, paste0("GO:",GO_cat))) %>% arrange(-log10(p_value)) %>% mutate(query = case_when(query == "query_1" ~ query1), term_name = as.factor(term_name))
  } else if(is.null(query3) & is.null(query4)) {gost.res_filt <- gost.res$result %>% filter(str_detect(source, paste0("GO:",GO_cat))) %>% arrange(-log10(p_value)) %>% mutate(query = case_when(query == "query_1" ~ query1, query == "query_2" ~query2), term_name = as.factor(term_name))
  } else if(is.null(query4)) {gost.res_filt <- gost.res$result %>% filter(str_detect(source, paste0("GO:",GO_cat))) %>% arrange(-log10(p_value)) %>% mutate(query = case_when(query == "query_1" ~ query1, query == "query_2" ~query2, query == "query_3" ~query3), term_name = as.factor(term_name))
  } else {gost.res_filt <- gost.res$result %>% filter(str_detect(source, paste0("GO:",GO_cat))) %>% arrange(-log10(p_value)) %>% mutate(query = case_when(query == "query_1" ~ query1, query == "query_2" ~query2, query == "query_3" ~query3, query == "query_4" ~query4), term_name = as.factor(term_name))} 
  
  gost.res_filt.query1 <- gost.res_filt %>% filter(query == query1)
  cat(blue(paste0("\n", query1, " gost\n")))
  paste('"',unique(gost.res_filt.query1$term_name),'",', sep = "", collapse = '\n') %>% cat()
  cat(blue("\nremove term redundancy with revigo rrvgo\n")) 
  simMatrix.query1 <-  calculateSimMatrix(gost.res_filt.query1$term_id, orgdb="org.Hs.eg.db",ont=GO_cat, method="Rel")
  scores.query1 <- setNames(-log10(gost.res_filt.query1$p_value), gost.res_filt.query1$term_id)
  reducedTerms.query1 <- reduceSimMatrix(simMatrix.query1, scores.query1, threshold=0.7, orgdb="org.Hs.eg.db")
  cat(blue("\nDownregulated rrvgo\n"))
  print(reducedTerms.query1)
  reducedTerms.query1.filt = reducedTerms.query1 %>% arrange(-score) %>% distinct(parentTerm, .keep_all = TRUE)
  if(!is.null(query2)){
    gost.res_filt.query2 <- gost.res_filt %>% filter(query == query2)
    cat(blue(paste0("\n", query2, " gost\n")))
    paste('"',unique(gost.res_filt.query2$term_name),'",', sep = "", collapse = '\n') %>% cat()
    simMatrix.query2 <-  calculateSimMatrix(gost.res_filt.query2$term_id, orgdb="org.Hs.eg.db",ont="BP", method="Rel")
    scores.query2 <- setNames(-log10(gost.res_filt.query2$p_value), gost.res_filt.query2$term_id)
    reducedTerms.query2 <- reduceSimMatrix(simMatrix.query2, scores.query2, threshold=0.7, orgdb="org.Hs.eg.db")
    cat(blue("\nUpregulated rrvgo\n"))
    print(reducedTerms.query2)
    reducedTerms.query2.filt = reducedTerms.query2 %>% arrange(-score) %>% distinct(parentTerm, .keep_all = TRUE)
  }
  if(!is.null(query3)){
    gost.res_filt.query3 <- gost.res_filt %>% filter(query == query3)
    cat(blue(paste0("\n", query3, " gost\n")))
    paste('"',unique(gost.res_filt.query3$term_name),'",', sep = "", collapse = '\n') %>% cat()
    simMatrix.query3 <-  calculateSimMatrix(gost.res_filt.query3$term_id, orgdb="org.Hs.eg.db",ont="BP", method="Rel")
    scores.query3 <- setNames(-log10(gost.res_filt.query3$p_value), gost.res_filt.query3$term_id)
    reducedTerms.query3 <- reduceSimMatrix(simMatrix.query3, scores.query3, threshold=0.7, orgdb="org.Hs.eg.db")
    cat(blue("\nUpregulated rrvgo\n"))
    print(reducedTerms.query3)
    reducedTerms.query3.filt = reducedTerms.query3 %>% arrange(-score) %>% distinct(parentTerm, .keep_all = TRUE)
  }
  if(!is.null(query4)){
    gost.res_filt.query4 <- gost.res_filt %>% filter(query == query4)
    cat(blue(paste0("\n", query4, " gost\n")))
    paste('"',unique(gost.res_filt.query4$term_name),'",', sep = "", collapse = '\n') %>% cat()
    simMatrix.query4 <-  calculateSimMatrix(gost.res_filt.query4$term_id, orgdb="org.Hs.eg.db",ont="BP", method="Rel")
    scores.query4 <- setNames(-log10(gost.res_filt.query4$p_value), gost.res_filt.query4$term_id)
    reducedTerms.query4 <- reduceSimMatrix(simMatrix.query4, scores.query4, threshold=0.7, orgdb="org.Hs.eg.db")
    cat(blue("\nUpregulated rrvgo\n"))
    print(reducedTerms.query4)
    reducedTerms.query4.filt = reducedTerms.query4 %>% arrange(-score) %>% distinct(parentTerm, .keep_all = TRUE)
  }
  if(is.null(query2) & is.null(query3) & is.null(query4)) {
    reducedTerms_combined <- mutate(reducedTerms.query1.filt, query = query1)
  } else if(is.null(query3) & is.null(query4)) {
    reducedTerms_combined <- bind_rows(mutate(reducedTerms.query1.filt, query = query1), mutate(reducedTerms.query2.filt, query = query2)) %>% mutate(query = factor(query, levels = c(query2, query1)))
  } else if (is.null(query4)) {
    reducedTerms_combined <- bind_rows(mutate(reducedTerms.query1.filt, query = query1), mutate(reducedTerms.query2.filt, query = query2), mutate(reducedTerms.query3.filt, query = query3)) %>% mutate(query = factor(query, levels = c(query3, query2, query1)))
  } else {
    reducedTerms_combined <- bind_rows(mutate(reducedTerms.query1.filt, query = query1), mutate(reducedTerms.query2.filt, query = query2), mutate(reducedTerms.query3.filt, query = query3), mutate(reducedTerms.query4.filt, query = query4)) %>% mutate(query = factor(query, levels = c(query4, query3, query2, query1)))
  }
  return(reducedTerms_combined)
}


rrvgo_plot <- function(gprofiler2rrvgo = reducedTerms_combined, n_terms = 10, remove_terms = NULL, colours = 4, cols =  c("firebrick2", "dodgerblue2", "forestgreen", "gold2"), label_angle = 0){
  gprofiles.filt <- gprofiler2rrvgo %>% arrange(score) %>% distinct(parentTerm, .keep_all = TRUE) %>% mutate(parentTerm = factor(parentTerm, levels = parentTerm)) %>% group_by(query) %>% top_n(n = n_terms, wt = score) %>% filter(!parentTerm %in% remove_terms)
  if(colours < 5){
    ggplot(gprofiles.filt, aes( x = parentTerm, y = score, fill = query) ) + 
      geom_col(aes(y = ( (score)/sum(score) * 100.0) )) + # makes more space to axis labels
      scale_fill_manual(values = cols) +
      facet_grid(query ~ ., space = "free", scales = "free") +
      coord_flip() + theme_classic() + guides(fill="none") +
      ylab(label = expression(log[10]~(P~value)) ) + xlab("") +
      theme(strip.text.y = element_text(size=8, face = "bold",angle = label_angle), axis.title = element_text(size = 9), legend.position = "none", plot.subtitle = element_text(size = 10)) + # legend.position = c(0.8, 0.2), 
      scale_y_continuous(expand = c(0, 0))
  } else{
    ggplot(gprofiles.filt, aes( x = parentTerm, y = score, fill = query) ) + 
      geom_col(aes(y = ( (score)/sum(score) * 100.0) )) + # makes more space to axis labels
      scale_fill_manual(values = get_palette("npg", 10)) +
      facet_grid(query ~ ., space = "free", scales = "free") +
      coord_flip() + theme_classic() + guides(fill="none") +
      ylab(label = expression(log[10]~(P~value)) ) + xlab("") +
      theme(strip.background = element_blank(), strip.text.y = element_text(size=8, face = "bold",angle = label_angle), axis.title = element_text(size = 9), legend.position = "none", plot.subtitle = element_text(size = 10)) + # legend.position = c(0.8, 0.2), 
      scale_y_continuous(expand = c(0, 0))
  }
}


curated_profiles <- function(gprofiles = IRevents_gprofiles_curated, colours = 4, cols =  c("firebrick2", "dodgerblue2", "forestgreen", "gold2"), label_angle = 90, title = NULL, subtitle = NULL, theme = "classic"){
  gprofiles <- gprofiles %>% distinct(term_name, .keep_all = TRUE) %>%  arrange( -log10(p_value) ) %>%  mutate(term_name = factor(term_name, levels = term_name))
  plot = ggplot(gprofiles, aes(x =  -log10(p_value), y = term_name, fill = query) ) + 
    geom_col(aes(x = ( (-log10(p_value))/sum(-log10(p_value)) * 100.0) )) + # makes more space to axis labels
    facet_grid(query ~ ., space = "free", scales = "free") + guides(fill="none") + 
    scale_x_continuous(expand = c(0, 0)) +  labs(x = expression(log[10]~FDR), y = "", title = title, subtitle = subtitle)
  if(theme == "oz"){ plot = plot + theme_oz() } else{ plot = plot + theme_classic() }
  plot = plot + theme(strip.text.y = element_text(size=8, face = "bold", angle = label_angle), axis.text = element_text(size = 8), axis.title = element_text(size = 9), legend.position = "none", 
        plot.title.position = "plot", plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5), axis.line.y = element_line(arrow = arrow(length = unit(0,"mm"))))
  # , plot.subtitle.position = "plot", plot.subtitle = element_text(size = 10)) + # legend.position = c(0.8, 0.2), 
  if(colours < 5){  plot = plot + scale_fill_manual(values = cols)  } else{   plot = plot + scale_fill_manual(values = get_palette("npg", 10))  }
  return(plot)
}


signed_z <- function(log2fc, pvalue){
  z <- ifelse(test = log2fc > 0,
              yes = qnorm(1 - (pvalue / 2)),
              no = qnorm(pvalue / 2) )
  return(z)
}


nest_batch <- function(inner, within, set_to=".") {
  n_instances <- apply(table(inner, within)!=0, 1, sum)
  if (any(n_instances>1)) {
    stop(paste(names(n_instances)[n_instances>1], collapse=", "), " appear in multiple parents")
  }
  levels_in <- levels(inner)
  # each inner level has in a unique group, find it
  corresponding_group <-within[match(levels_in, inner)]
  # for each group, find index of first inner level
  which_inner <- match(levels(within), corresponding_group)
  # and set it to the common value
  if (any(levels(inner)[-which_inner]==set_to )) { # we're going to duplicate an existing level so best stop
    stop("One of the batches is already called ", set_to, ". Please use a non-existing level")
  }
  levels(inner)[which_inner] <- set_to
  inner
}

recode_within <- function(inner, ...) {
  within <- do.call(interaction, alist(...))
  tab <- table(inner, within)!=0 # which batches are in which nest
  if (any(rowSums(tab)>1)) {
    stop("Some inner levels appear in multiple outer groups.")
  }
  factor(apply(tab, 2, cumsum)[cbind(as.character(inner),as.character(within))]) # cumsum to get incrementing index within group.
}


res_top <- function(r, i=1) {
  ind <- which(r$padj<0.05) # indexes of significant genes
  ind[order(abs(r$log2FoldChange)[ind])][i]
}


# Gene Expression Plots ------------------------------------------------------
ma_plot <- function(deseq_res, title = NULL, subtitle = NULL, gene_labels = NULL, distinct_labels = TRUE, arrow_labels = FALSE, numerator = NULL, denomenator = NULL, ymax = 10, xmax = 100000, xpos = 0.5, dpi = 300, significance_threshold = "padj"){
  deseq_res <- deseq_res %>% mutate(sig = case_when(!!sym(significance_threshold) < 0.05 & abs(log2FoldChange) < 1 ~ "sig", !!sym(significance_threshold) < 0.05 & abs(log2FoldChange) >= 1 ~ "sig_strong", !!sym(significance_threshold) >= 0.05 ~ "non_sig"), 
                                    direction = ifelse(log2FoldChange > 0, "up", "down"), class = paste(sig, direction), log2FoldChange = case_when(log2FoldChange > xmax ~ Inf, log2FoldChange < -xmax ~ -Inf, TRUE ~ log2FoldChange), pvalue = case_when(-log10(pvalue) > ymax ~ 10^-ymax, TRUE ~ pvalue))
  
  plot <- ggplot(data = deseq_res, aes(x = baseMean, y = log2FoldChange)) +  # NB if pvalue == NA then will not show
    ggrastr::rasterise(geom_point(aes(colour = class), size = 0.5), dpi = dpi) +
    scale_colour_manual(values = c("non_sig up" = "gray", "non_sig down" = "gray", 
                                   "sig up" = "#F36F6F",
                                   "sig_strong up" = "#EB4445",
                                   "sig down" = "#A6CEE3",#"#4F8FC4",
                                   "sig_strong down" = "#79B1D3")) +
    labs(x = expression(log[2]~base~mean~expression), y = paste0("log2 fold change (",numerator," vs ",denomenator,")"), title = title, subtitle = subtitle) + guides(colour = "none") +
    scale_y_continuous(expand = c(0,0), limits = c(-ymax,ymax)) +scale_x_continuous(limits = c(0,xmax)) +
    theme_oz() +  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5), panel.border = element_blank(), axis.ticks = element_line(colour = "black") ) +
    geom_hline(yintercept = 0, linetype = 'dotted')
  
  
  if(!is.null(gene_labels)){
    if(distinct_labels == TRUE){
      plot <- plot + geom_text_repel(fontface = "italic", data = distinct(arrange(filter(deseq_res, gene_name %in% gene_labels), wt = -abs(stat)), gene_name, .keep_all = TRUE), aes(x = baseMean, y = log2FoldChange, label = gene_name), max.overlaps = Inf, min.segment.length = unit(0, "lines"), size = 2.3)
    } else if(distinct_labels == FALSE){
      plot <- plot + geom_text_repel(fontface = "italic", data = arrange(filter(deseq_res, gene_name %in% gene_labels), wt = -abs(stat)), aes(x = baseMean, y = log2FoldChange, label = gene_name), max.overlaps = Inf, min.segment.length = unit(0, "lines"),size = 2.3)
    }
  }
  if(arrow_labels == TRUE){
    plot <- plot + 
      geom_segment(aes(y = ymax*0.25, yend = ymax*0.95, x= xmax*0.88, xend= xmax * 0.88), arrow=arrow(length=unit(0.3,"cm"))) + geom_segment(aes(y = -ymax*0.25, yend = -ymax*0.95, x= xmax*0.88, xend= xmax*0.88),arrow=arrow(length=unit(0.3,"cm"))) +
      annotate("text", y = ymax*0.6, x = xmax*0.93, label = paste0("Up in ",numerator), size = 3) + annotate("text", y = -ymax*0.6, x = xmax*0.93, label = paste0("Up in ",denomenator), size = 3)
  }
  return(plot)
}


volcano_plot <- function(deseq_res, title = NULL, subtitle = NULL, gene_labels = NULL, transcript_labels = NULL, distinct_labels = FALSE, ymax = 16.5, xmax = 3, xpos = 0.5, numerator = NULL, denomenator = NULL, arrow_labels = TRUE, 
                         dpi = 300, significance_threshold = "padj", facet_by = NULL, facet_nrow = NULL, facet_ncol = NULL, facet_labeller = "label_value"){
  deseq_res %<>% mutate(sig = case_when(!!sym(significance_threshold) < 0.05 & abs(log2FoldChange) < 1 ~ "sig", !!sym(significance_threshold) < 0.05 & abs(log2FoldChange) >= 1 ~ "sig_strong", !!sym(significance_threshold) >= 0.05 ~ "non_sig"), 
                        direction = ifelse(log2FoldChange > 0, "up", "down"), class = paste(sig, direction), log2FoldChange = case_when(log2FoldChange > xmax ~ Inf, log2FoldChange < -xmax ~ -Inf, TRUE ~ log2FoldChange), pvalue = case_when(-log10(pvalue) > ymax ~ 10^-ymax, TRUE ~ pvalue))
  # de_tally <- deseq_res %>% group_by(sig, direction, class) %>% count() %>% filter(sig != "non_sig") %>% drop_na() %>% mutate(position = ifelse(sig == "sig", xpos, xmax-1), position = ifelse( direction == "down", -1 * position, position), n = formatC(n, format="f", big.mark=",", digits=0)) #%>% 
  # mutate(pvalue = case_when(-log10(pvalue) < ymax ~ Inf, TRUE ~ pvalue)) %>% # include dots & labels for when Pvalue is out of coordinates of ymax
  plot <- ggplot(deseq_res, aes(x = log2FoldChange, y = -log10(pvalue))) +  #geom_point(aes(colour = class ), size = 0.5) +
    ggrastr::rasterise(geom_point(aes(colour = class), size = 0.5), dpi = dpi) +
    scale_colour_manual(values = c("non_sig up" = "gray", "non_sig down" = "gray", 
                                   "sig up" = "#F36F6F",
                                   "sig_strong up" = "#EB4445",
                                   "sig down" = "#A6CEE3",#"#4F8FC4",
                                   "sig_strong down" = "#79B1D3")) +
    labs(y = expression(-log[10]~P~value), x = paste0("log2 fold change (",numerator," vs ",denomenator,")"), title = title, subtitle = subtitle) +
    guides(colour = "none") +
    scale_y_continuous(expand = c(0,0), limits = c(0,ymax)) + scale_x_continuous(limits = c(-xmax,xmax)) +
    theme_oz() +  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5), panel.border = element_blank(), axis.ticks = element_line(colour = "black") )
    # geom_text(fontface = "bold", data = de_tally, aes(x = position, y = ymax - 0.5, label = n, colour = class), size = 2.5 ) +

  if(!is.null(gene_labels)){
    if(distinct_labels == TRUE){
      plot <- plot + geom_text_repel(fontface = "italic", data = distinct(arrange(filter(deseq_res, gene_name %in% gene_labels), wt = pvalue), gene_name, .keep_all = TRUE), aes(x = log2FoldChange, y = -log10(pvalue), label = gene_name), max.overlaps = Inf, min.segment.length = unit(0, "lines"), size = 2.3)
    } else if(distinct_labels == FALSE){
      plot <- plot + geom_text_repel(fontface = "italic", data = filter(deseq_res, gene_name %in% gene_labels), aes(x = log2FoldChange, y = -log10(pvalue), label = gene_name), max.overlaps = Inf, min.segment.length = unit(0, "lines"),size = 2.3)
    }
  } else if(!is.null(transcript_labels)){
    plot <- plot + geom_text_repel(fontface = "italic", data = filter(deseq_res, transcript_id %in% transcript_labels), aes(x = log2FoldChange, y = -log10(pvalue), label = gene_name), max.overlaps = 15, min.segment.length = unit(0, "lines"), size = 2.3)
  }
  
  if(arrow_labels == TRUE){
    plot <- plot + 
      geom_segment(aes(x = xmax*0.25, xend = xmax*0.95, y= ymax*0.88, yend= ymax * 0.88), arrow=arrow(length=unit(0.3,"cm"))) + geom_segment(aes(x = -xmax*0.25, xend = -xmax*0.95, y= ymax*0.88, yend= ymax*0.88),arrow=arrow(length=unit(0.3,"cm"))) +
      annotate("text", x = xmax*0.6, y = ymax*0.93, label = paste0("Up in ",numerator), size = 3) + annotate("text", x = -xmax*0.6, y = ymax*0.93, label = paste0("Up in ",denomenator), size = 3)
  }
  
  return(plot)
}


# volcano_plot_deseq <- 
#   function(data, labels_list = NULL){
#     if(is.null(labels_list)){labels_list = data %>% arrange(padj) %>% distinct(gene_name, .keep_all = TRUE) %>% top_n(-20, wt = pvalue) %>% pull(gene_name)}
#     labels_plot <- data  %>% arrange(padj) %>% distinct(gene_name, .keep_all = TRUE) %>% filter(gene_name %in% labels_list, padj < 0.05)
#     data2 <- data %>% mutate(direction = case_when(padj < 0.05 & log2FoldChange > 0 ~ "upregulated", padj < 0.05 & log2FoldChange < 0 ~ "downregulated", TRUE ~ "unchanged"))
#     ggplot(data = data2, aes(x = log2FoldChange, y =  -log10(padj), colour = direction)) +  geom_point(size = 0.5) + #xlim(1,5) + ylim(-3,3) +
#       scale_colour_manual( values = c(unchanged = "darkgray", upregulated = "firebrick2", downregulated = "dodgerblue2") ) +  theme_classic() +  theme(legend.title = element_blank(), legend.position = "none") + 
#       ylab(expression( -log[10]~FDR )) + 
#       geom_text_repel(data = labels_plot, aes(x = log2FoldChange, y = -log10(padj), label = gene_name), fontface = "italic", colour = "black", min.segment.length = 0, box.padding = 0.5, point.padding = 0.5, size=2.6, max.overlaps = Inf) +
#       geom_hline(yintercept = -log10(0.05), linetype = 2) + geom_vline(xintercept =0, linetype = 2)
#   }
# 
# # volcano plotting undajusted p value
# volcano_plot_deseq.pvalue <-  
#   function(data = mn_d25.res, labels_list = NULL){
#     if(is.null(labels_list)){labels_list = data %>% top_n(-20, wt = pvalue) %>% pull(gene_name)}
#     labels_plot <- data %>% filter(gene_name %in% labels_list, pvalue < 0.05) %>% arrange(pvalue) %>% distinct(gene_name, .keep_all = TRUE)
#     data2 <- data %>% mutate(direction = case_when(pvalue < 0.05 & log2FoldChange > 0 ~ "upregulated", pvalue < 0.05 & log2FoldChange < 0 ~ "downregulated", TRUE ~ "unchanged"))
#     ggplot(data = data2, aes(x = log2FoldChange, y =  -log10(pvalue), colour = direction)) +  geom_point(size = 0.5) + #xlim(1,5) + ylim(-3,3) +
#       scale_colour_manual( values = c(unchanged = "darkgray", upregulated = "firebrick2", downregulated = "dodgerblue2") ) +  
#       theme_classic() +  theme(legend.title = element_blank(), legend.position = "none", axis.text=element_text(size=8), axis.title=element_text(size=10)) + 
#       ylab(expression( -log[10]~P~value )) + 
#       geom_text_repel(data = labels_plot, aes(x = log2FoldChange, y = -log10(pvalue), label = gene_name), fontface = "italic", colour = "black", min.segment.length = 0, box.padding = 0.5, point.padding = 0.5, size=2.6, max.overlaps = Inf) +
#       geom_hline(yintercept = -log10(0.05), linetype = 2) + geom_vline(xintercept =0, linetype = 2)
#   }



violin_plot <- function(data, color_by = "condition", continuous = "log2FoldChange", cols = c("dodgerblue2","firebrick2","forestgreen","gold2","purple"), null_value = NULL, ylims = NULL, #c(-2,2), 
                        plot_stats = TRUE, stats_test = "t_test", one_sample = FALSE, stat_ypos = NULL, stat_label = "p.signif", adjust_p = TRUE, tip_length = 0.01, plot_violin = TRUE, violin_adjust = 1, 
                        facet_by = NULL, facet_by2 = NULL, facet_nrow = NULL, facet_ncol = NULL, facet_labeller = "label_value", facet_order = NULL, ylabel = NULL, xlabel = "", title = NULL, subtitle = NULL, arrow_labels = TRUE, arrow_label_x = 0.94, numerator = NULL, denomenator = NULL, dpi = 300){
  
  if(!is.null(facet_by2)){
    print(paste0("Faceting by ",facet_by," and ",facet_by2))
    data %<>% dplyr::select(group = all_of(color_by), lfc = all_of(continuous), facet = all_of(facet_by), facet2 = all_of(facet_by2), everything())
  } else if(!is.null(facet_by)){
    print(paste0("Faceting by ",facet_by))
    data %<>% dplyr::select(group = all_of(color_by), lfc = all_of(continuous), facet = all_of(facet_by), everything())
    if(!is.null(facet_order)){data %<>% mutate(facet = factor(facet, levels = facet_order))}
  } else { data <- data %>% dplyr::select(group = all_of(color_by), lfc = all_of(continuous), everything()) }
  
  plot <- ggplot(data, aes(x = group, y = lfc, colour = group))
  if(plot_violin == TRUE){ plot = plot + ggrastr::rasterise(geom_violin(adjust = violin_adjust), dpi = {{dpi}}) }
  plot = plot + ggrastr::rasterise(geom_sina(adjust = violin_adjust, size = 0.5), dpi = {{dpi}}) + 
    ggrastr::rasterise(geom_boxplot(aes(fill = group), colour = "black", width=0.2, alpha = 0.5, outlier.shape = NA), dpi = dpi) + 
    scale_colour_manual(values = {{ cols }}) + scale_fill_manual(values = {{ cols }}) +
    theme_oz(xaxis_arrow = FALSE) + theme(legend.position = "none", plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5), panel.border = element_blank(), axis.line.x = element_line(arrow = arrow(length = unit(0,"mm"))),
) +#guides(colour = "none") +
    labs(y = paste0("log2 fold change (",numerator," vs ",denomenator,")"), x = xlabel, title = title, subtitle = subtitle) 
  if(!is.null(ylims)){ plot <- plot + coord_cartesian(ylim = {{ ylims }}) }
  if(!is.null(null_value)){ plot <- plot + geom_hline(yintercept = {{ null_value }}, linetype = "dashed", size = 0.5) }
  if(!is.null(facet_by2)) { plot <- plot + facet_grid(facet2 ~ facet, scales = "free", labeller = facet_labeller) 
    } else if(!is.null(facet_by)) { plot <- plot + facet_wrap(~facet, scales = "free", nrow = {{ facet_nrow }}, ncol = {{ facet_ncol }}, labeller = facet_labeller) }
  if(!is.null(ylabel)){ plot <- plot + labs(y = ylabel) }
  if(arrow_labels == TRUE){ plot <- plot + 
    annotation_custom(grid::textGrob(label = paste0(numerator), x = unit({{ arrow_label_x }}, "npc"), y = unit(0.8, "npc"), gp = grid::gpar(fontsize = 8, fontface=1), rot = 90)) + 
    annotation_custom(grid::textGrob(label = paste0(denomenator), x = unit({{ arrow_label_x }}, "npc"), y = unit(0.2, "npc"), gp = grid::gpar(fontsize = 8, fontface=1), rot = 90)) + 
    annotation_custom(grid::linesGrob(arrow=arrow(type="open", ends="first", length=unit(3,"mm")), gp=gpar(col="black", lwd=1), y = c(0.95,0.65), x = unit(0.97, "npc"))) +
    annotation_custom(grid::linesGrob(arrow=arrow(type="open", ends="first", length=unit(3,"mm")), gp=gpar(col="black", lwd=1), y = c(0.05,0.35), x = unit(0.97, "npc"))) }
  
  if(plot_stats == TRUE){     
    
    if(!is.null(facet_by2)){
      if(one_sample == TRUE){ cohens_effect_size = data %>% group_by(facet, facet2, group) %>% cohens_d(lfc ~ 1, var.equal = TRUE) 
      if(stats_test == "t_test"){ print(paste("Using one sample t_test"))
        stat <- data %>% group_by(facet, facet2, group) %>% filter(lfc != Inf, lfc != -Inf) %>% t_test(lfc ~ 1, mu = {{null_value}}) %>% add_significance() %>% add_xy_position(x = "group")    
      } else if(stats_test == "wilcox_test"){ print(paste("Using one sample wilcox_test"))
        stat <- data %>% group_by(facet, facet2, group) %>% filter(lfc != Inf, lfc != -Inf) %>% wilcox_test(lfc ~ 1, mu = {{null_value}}) %>% add_significance() %>% add_xy_position(x = "group")
      }
      } else if (one_sample == FALSE){ cohens_effect_size = data %>% group_by(facet, facet2) %>% cohens_d(lfc ~ group, var.equal = TRUE) 
      if(stats_test == "t_test"){ print(paste("Using two sample t_test"))
        stat <- data %>% group_by(facet, facet2) %>% filter(lfc != Inf, lfc != -Inf) %>% t_test(lfc ~ group) %>% add_significance() %>% add_xy_position(x = "group")     
      } else if(stats_test == "wilcox_test"){ print(paste("Using two sample wilcox_test"))
        stat <- data %>% group_by(facet, facet2) %>% filter(lfc != Inf, lfc != -Inf) %>% wilcox_test(lfc ~ group) %>% add_significance() %>% add_xy_position(x = "group")
      }
      } 
    } else if(!is.null(facet_by)){
      if(one_sample == TRUE){ cohens_effect_size = data %>% group_by(facet, group) %>% cohens_d(lfc ~ 1, var.equal = TRUE) 
      if(stats_test == "t_test"){ print(paste("Using one sample t_test"))
        stat <- data %>% group_by(facet, group) %>% filter(lfc != Inf, lfc != -Inf) %>% t_test(lfc ~ 1, mu = {{null_value}}) %>% add_significance() %>% add_xy_position(x = "group")    
      } else if(stats_test == "wilcox_test"){ print(paste("Using one sample wilcox_test"))
        stat <- data %>% group_by(facet, group) %>% filter(lfc != Inf, lfc != -Inf) %>% wilcox_test(lfc ~ 1, mu = {{null_value}}) %>% add_significance() %>% add_xy_position(x = "group")
      }
      } else if (one_sample == FALSE){ cohens_effect_size = data %>% group_by(facet) %>% cohens_d(lfc ~ group, var.equal = TRUE) 
      if(stats_test == "t_test"){ print(paste("Using two sample t_test"))
        stat <- data %>% group_by(facet) %>% filter(lfc != Inf, lfc != -Inf) %>% t_test(lfc ~ group) %>% add_significance() %>% add_xy_position(x = "group")     
      } else if(stats_test == "wilcox_test"){ print(paste("Using two sample wilcox_test"))
        stat <- data %>% group_by(facet) %>% filter(lfc != Inf, lfc != -Inf) %>% wilcox_test(lfc ~ group) %>% add_significance() %>% add_xy_position(x = "group")
      }
      } 
    } else if (is.null(facet_by) == TRUE) { 
      if(one_sample == TRUE){ cohens_effect_size = data %>% group_by(group) %>% cohens_d(lfc ~ 1, var.equal = TRUE) 
        if(stats_test == "t_test"){ print(paste("Using one sample t_test"))
          stat <- data %>% group_by(group) %>% filter(lfc != Inf, lfc != -Inf) %>% t_test(lfc ~ 1, mu = {{null_value}}) %>% add_significance() #%>% add_xy_position()
      } else if(stats_test == "wilcox_test"){ print(paste("Using one sample wilcox_test"))
        stat <- data %>% group_by(group) %>% filter(lfc != Inf, lfc != -Inf) %>% wilcox_test(lfc ~ 1, mu = {{null_value}}) %>% add_significance() #%>% add_xy_position()
      }
      } else if (one_sample == FALSE){ cohens_effect_size = data %>% cohens_d(lfc ~ group, var.equal = TRUE) 
      if(stats_test == "t_test"){ print(paste("Using two sample t_test"))
        stat <- data %>% filter(lfc != Inf, lfc != -Inf) %>% t_test(lfc ~ group) %>% add_significance() #%>% add_xy_position()      
      } else if(stats_test == "wilcox_test"){ print(paste("Using two sample wilcox_test"))
        stat <- data %>% filter(lfc != Inf, lfc != -Inf) %>% wilcox_test(lfc ~ group) %>% add_significance() #%>% add_xy_position()
      }
      }
    }
    if(adjust_p == FALSE){ stat %<>% mutate( p.signif = case_when(p<0.001~"****", p<0.001~"***", p<0.01~"**", p<0.05~"*", TRUE ~ "ns")) %>% select(-contains("adj")) }
    if(one_sample == TRUE){ plot <- plot + stat_pvalue_manual(stat, label = {{stat_label}}, y.position = stat_ypos, xmin = "group", xmax = NULL, size = 5, hide.ns = TRUE)
    } else if (one_sample == FALSE){
      if(!is.null(stat_ypos)){ stat %<>% mutate(y.position = stat_ypos) } else { stat %<>% mutate(y.position = 0.95*ylims[2]) }
      plot <- plot + stat_pvalue_manual(stat, label = {{stat_label}}, tip.length = {{tip_length}}, size = 5, hide.ns = TRUE)
    }
    print(stat, n=Inf)
    print(cohens_effect_size, n=Inf)
  }
  return(plot)
}


# sinaplot.one.sample <- 
#   function(data = boxplot.lfc.merged, groups = "category", continuous = "log2FoldChange", flip = "yes", tip.length = 0.01, labels = "gene_name", colours = 3, stat.y.position = 1.00, label_number = 10, stats.test = "t_test", null_value = 0, width = 0.2){
#     data2 <- data %>% dplyr::select(group = groups, lfc = continuous, label = labels)
#     print(data2)
#     print(paste("Stats test used = ", stats.test))
#     if(stats.test == "t_test"){
#       stats_test <- data2 %>% dplyr::select(group, lfc) %>% group_by(group) %>% filter(lfc != Inf) %>% filter(lfc != -Inf) %>% t_test(lfc ~ 1, mu = null_value) %>% add_significance() %>% mutate(y.position = stat.y.position)
#       print(stats_test)
#     } else{
#       stats_test <- data2 %>% dplyr::select(group, lfc) %>% group_by(group) %>% wilcox_test(lfc ~ 1, mu = null_value) %>% add_significance() %>% mutate(y.position = stat.y.position)
#       print(stats_test)
#     }
#     if(flip == "yes"){
#       print(paste("Flipping coords"))
#       if(colours < 5){
#         ggplot(data2, aes(x = group, y= lfc)) + 
#           geom_violin(aes(colour = group)) + 
#           geom_sina(size = 0.5, aes(colour = group)) + 
#           geom_boxplot(data = data2, aes(fill = group), colour = "black", width=width, outlier.shape = NA) +
#           scale_colour_manual(values = c("firebrick2", "dodgerblue2", "forestgreen", "gold2")) + scale_fill_manual(values = c("firebrick2", "dodgerblue2", "forestgreen", "gold2")) +
#           theme_classic() + theme(legend.position = "none", axis.text =element_text(size=10), axis.title = element_text(size=10)) +
#           geom_hline(yintercept =null_value, linetype = 2) + coord_flip() +
#           geom_text_repel(data = top_n(data2, label_number, abs(lfc)), aes(x = group, y = lfc, label = label), fontface = "italic", colour = "black", min.segment.length = 0, box.padding = 0.5, point.padding = 0, size=2.7, max.overlaps = Inf) +
#           stat_pvalue_manual(stats_test, label = "p.signif", tip.length = tip.length, xmin = "group", xmax = NULL, hide.ns = TRUE, size = 3) #, position_dodge(0.8)) # add manually p values
#       } else{
#         ggplot(data2, aes(x = group, y= lfc, color = group)) + 
#           geom_violin(aes(colour = group)) + 
#           geom_sina(size = 0.5, aes(colour = group)) + 
#           geom_boxplot(data = data2, aes(fill = group), colour = "black", width=width, outlier.shape = NA) +
#           scale_colour_manual(values = get_palette("npg",colours)) + scale_fill_manual(values = get_palette("npg",colours)) +
#           theme_classic() + theme(legend.position = "none", axis.text =element_text(size=10), axis.title = element_text(size=10)) +
#           geom_hline(yintercept =null_value, linetype = 2) + coord_flip() +
#           geom_text_repel(data = top_n(data2, label_number, abs(lfc)), aes(x = group, y = lfc, label = label), fontface = "italic", colour = "black", min.segment.length = 0, box.padding = 0.5, point.padding = 0, size=2.7, max.overlaps = Inf) +
#           stat_pvalue_manual(stats_test, label = "p.signif", tip.length = tip.length, xmin = "group", xmax = NULL, hide.ns = TRUE, size = 3) #, position_dodge(0.8)) # add manually p values
#       }
#     } else{
#       print(paste("Not flipping coords"))
#       if(colours < 5){
#         ggplot(data2, aes(x = group, y= lfc, color = group)) + 
#           geom_violin() +  #aes(colour = group)
#           geom_sina(size = 0.5) + # , aes(colour = group)
#           geom_boxplot(data = data2, aes(fill = group), colour = "black", width=width, outlier.shape = NA) + 
#           scale_colour_manual(values = c("firebrick2", "dodgerblue2", "forestgreen", "gold2")) + scale_fill_manual(values = c("firebrick2", "dodgerblue2", "forestgreen", "gold2")) +
#           theme_classic() + theme(legend.position = "none", axis.text = element_text(size=10), axis.title = element_text(size=10)) +
#           geom_hline(yintercept = null_value, linetype = 2) + #coord_flip() +
#           geom_text_repel(data = top_n(data2, label_number, abs(lfc)), aes(x = group, y = lfc, label = label), fontface = "italic", colour = "black", min.segment.length = 0, box.padding = 0.5, point.padding = 0, size=2.7, max.overlaps = Inf) +
#           stat_pvalue_manual(stats_test, label = "p.signif", tip.length = tip.length, xmin = "group", xmax = NULL, hide.ns = TRUE, size = 3) #, position_dodge(0.8)) # add manually p values
#       } else{
#         ggplot(data2, aes(x = group, y= lfc, color = group)) + 
#           geom_violin(aes(colour = group)) + 
#           geom_sina(size = 0.5, aes(colour = group)) + 
#           geom_boxplot(data = data2, aes(fill = group), colour = "black", width=width, outlier.shape = NA) +
#           scale_colour_manual(values = get_palette("npg",colours)) + scale_fill_manual(values = get_palette("npg",colours)) +
#           theme_classic() + theme(legend.position = "none", axis.text = element_text(size=10), axis.title = element_text(size=10)) +
#           geom_hline(yintercept = null_value, linetype = 2) + #coord_flip() +
#           geom_text_repel(data = top_n(data2, label_number, abs(lfc)), aes(x = group, y = lfc, label = label), fontface = "italic", colour = "black", min.segment.length = 0, box.padding = 0.5, point.padding = 0, size=2.7, max.overlaps = Inf) +
#           stat_pvalue_manual(stats_test, label = "p.signif", tip.length = tip.length, xmin = "group", xmax = NULL, hide.ns = TRUE, size = 3) #, position_dodge(0.8)) # add manually p values
#       }
#     }
#   }


density_plot <- function(data, color_by = "condition", cols = c("dodgerblue2", "firebrick2", "forestgreen"), title = NULL, subtitle = NULL, ymax = NULL, xlims = c(-2,2), xpos = 0.5, numerator = NULL, denomenator = NULL, arrow_labels = TRUE, dpi = 300){
  mu <- data %>% group_by({{ color_by }}) %>% summarise(grp.mean=mean(log2FoldChange, na.rm = TRUE))
  plot <- ggplot(data, aes(x = log2FoldChange, color = {{ color_by }}, fill = {{ color_by }})) +
    ggrastr::rasterise(geom_density(alpha = 0.5), dpi = dpi) +
    scale_colour_manual(values = {{ cols }}) + scale_fill_manual(values = {{ cols }}) +
    geom_vline(data = mu, aes(xintercept=grp.mean, color={{ color_by }}), linetype="dashed") +
    geom_vline(xintercept =0, linetype = 2) +
    theme_oz() + theme(legend.position = "none", plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5), panel.border = element_blank(), axis.ticks = element_line(colour = "black") ) +#guides(colour = "none") +
    labs(y = "", x = paste0("log2 fold change (",numerator," vs ",denomenator,")"), title = title, subtitle = subtitle) +
    scale_y_continuous(expand = c(0,0), limits = c(0,ymax)) +
    scale_x_continuous(limits = xlims)
  if(arrow_labels == TRUE){
    plot <- plot + 
      annotation_custom(grid::textGrob(label = paste0("Up in ",numerator), x = unit(0.8, "npc"), y = unit(0.94, "npc"), gp = grid::gpar(fontface=1, fontsize = 8))) + 
      annotation_custom(grid::textGrob(label = paste0("Up in ",denomenator), x = unit(0.2, "npc"), y = unit(0.94, "npc"), gp = grid::gpar(fontface=1, fontsize = 8))) + 
      annotation_custom(grid::linesGrob(arrow=arrow(type="open", ends="first", length=unit(3,"mm")), gp=gpar(col="black", lwd=1), x = c(0.9,0.7), y = unit(0.85, "npc"))) +
      annotation_custom(grid::linesGrob(arrow=arrow(type="open", ends="first", length=unit(3,"mm")), gp=gpar(col="black", lwd=1), x = c(0.1,0.3), y = unit(0.85, "npc"))) 
  }
  return(plot)
}


upset_plot <- function(data, title = NULL, subtitle = NULL, ymax = NULL, overlap_by = "mutation", overlap_level = "gene_id", order_counts = "degree", order_group_list = NULL, plot_y_labels = TRUE,
                       split_direction = FALSE, split_up = "up", split_down = "down", significance_threshold = "padj", plot_numbers = TRUE, dot_colour = NULL, colour_by = NULL, label_colours = NULL){
  if(split_direction == TRUE) {
    data %<>% filter(!!sym(significance_threshold) < 0.05) %>% mutate(direction = case_when(stat > 0 ~ split_up, TRUE ~ split_down), group = paste(!!sym(overlap_by), direction)) %>% select(group, !!sym(overlap_level)) %>% group_by(!!sym(overlap_level)) %>% summarise(group = list(group))
  } else {data %<>% filter(!!sym(significance_threshold) < 0.05) %>% select({{overlap_by}}, !!sym(overlap_level)) %>% group_by(!!sym(overlap_level)) %>% summarise(group = list(!!sym(overlap_by)))}
  if(is.null(ymax)) { ymax = data }
  plot =  ggplot(data, aes(x = group)) + geom_bar() +
    geom_text(stat='count', aes(label=after_stat(count)), vjust=-0.2, hjust = 0.5, size = 3, angle = 0) +
    scale_y_continuous(breaks = NULL, lim = c(0, ymax), name = "") +
    theme_classic(base_size=8) +  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5), legend.title = element_blank(), legend.position = "none", axis.text=element_text(size=8), axis.title=element_text(size=10)) + labs(x = "", y = "") +
    theme_combmatrix(combmatrix.label.make_space = TRUE) + labs(title = title, subtitle = subtitle)
  if(order_counts == "degree"){plot = plot + scale_x_upset(order_by = "degree")} else {plot = plot + scale_x_upset(order_by = "freq")} # degree or freq
  if(!is.null(order_group_list)){plot = plot + scale_x_upset(sets = order_group_list)}
  if(plot_numbers == TRUE){plot = plot + geom_text(stat='count', aes(label=after_stat(count)), vjust=-0.2, hjust = 0.5, size = 3, angle = 0)}
  if(!is.null(dot_colour)){ plot = plot + theme_combmatrix(combmatrix.panel.point.color.fill = dot_colour) } #, combmatrix.panel.line.size = 0, combmatrix.label.make_space = FALSE)
  if(!is.null(label_colours)){ plot = plot + theme(axis.text.y = element_text(colour = label_colours))}
  # if(plot_y_labels == FALSE){ plot = plot + theme(axis.text.y = element_blank())}
  return(plot)
}




plot_de_correlation <- function(model_list1, model_list2, mutation1, mutation2, col = "stat", plot_corr_line = TRUE, plot_r_value = TRUE, plot_p_value = FALSE, 
                                gene_labels = NULL, transcript_labels = NULL, lsv_junction_labels = NULL, join_by = c("gene_name", "gene_id"), distinct_labels = TRUE, title = NULL, subtitle = NULL ){
  res <- dplyr::left_join(model_list1[[mutation1]], model_list2[[mutation2]],  by = join_by, suffix = c(".1", ".2")) 
  x_string <- paste0(col, ".1")
  y_string <- paste0(col, ".2")
  plot <- ggplot(res, aes_string(x = x_string, y = y_string )) + #geom_point(size = 1, alpha = 0.1) + 
    labs(x = gsub("_", " ", mutation1), y = gsub("_", " ", mutation2), title = title, subtitle = subtitle ) + 
    geom_hline(yintercept = 0, linetype = 3) + geom_vline(xintercept = 0, linetype = 3) + geom_abline(slope =1 , intercept = 0, linetype = 3) +    
    geom_bin2d(bins = 100) +
    scale_fill_continuous(type = "viridis") +
    guides(fill = "none") + theme_bw() + theme_oz() + theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
  if(plot_r_value == TRUE & plot_p_value == FALSE){ plot = plot + ggpubr::stat_cor(method = "pearson", aes(label = ..r.label..) )
  } else if(plot_r_value == TRUE & plot_p_value == TRUE){ plot = plot + ggpubr::stat_cor(method = "pearson") }
  if(plot_corr_line == TRUE){ plot = plot + geom_smooth(method='lm') }
  
  if(!is.null(gene_labels)){
    if(distinct_labels == TRUE){
      plot <- plot + geom_point(data = distinct(arrange(filter(res, gene_name %in% gene_labels), wt = paste0(col, ".1")), gene_name, .keep_all = TRUE), aes_string(x = x_string, y = y_string), colour = "red", size = 0.5) +
        geom_text_repel(fontface = "italic", data = distinct(arrange(filter(res, gene_name %in% gene_labels), wt = paste0(col, ".1")), gene_name, .keep_all = TRUE), aes_string(x = x_string, y = y_string, label = "gene_name"), colour = "red", max.overlaps = 15, min.segment.length = unit(0, "lines"), size = 2.3)
    } else if(distinct_labels == FALSE){
      plot <- plot + geom_point(data = filter(res, gene_name %in% gene_labels), aes_string(x = x_string, y = y_string), colour = "red", size = 0.5) +
        geom_text_repel(fontface = "italic", data = filter(res, gene_name %in% gene_labels), aes_string(x = x_string, y = y_string, label = "gene_name"), colour = "red", max.overlaps = 15, min.segment.length = unit(0, "lines"), size = 2.3)
    }
  } else if(!is.null(transcript_labels)){
    plot <- plot + geom_point(data = filter(res, transcript_id %in% transcript_labels), aes_string(x = x_string, y = y_string), colour = "red", size = 0.5) +
      geom_text_repel(fontface = "italic", data = filter(res, transcript_id %in% transcript_labels), aes_string(x = x_string, y = y_string, label = "gene_name"), colour = "red", max.overlaps = 15, min.segment.length = unit(0, "lines"), size = 2.3)
  } else if(!is.null(lsv_junction_labels)){
    plot <- plot + geom_point(data = filter(res, lsv_junction_id %in% lsv_junction_labels), aes_string(x = x_string, y = y_string), colour = "red", size = 0.5) +
      geom_text_repel(fontface = "italic", data = filter(res, lsv_junction_id %in% lsv_junction_labels), aes_string(x = x_string, y = y_string, label = "gene_name"), colour = "red", max.overlaps = 15, min.segment.length = unit(0, "lines"), size = 2.3)
  }
  return(plot)
}


get_lower_tri<-function(cormat){ # Get lower triangle of the correlation matrix
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
get_upper_tri <- function(cormat){ # Get upper triangle of the correlation matrix
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

plot_vsd_gene_correlation <- function(vsd_counts, gene, variable, metadata, log10_convert = FALSE, filter_zeros = TRUE){
  require(ggrepel)
  min_vst = min(vsd_counts[,2])
  vsd_counts.long = vsd_counts %>% filter(gene_name == {{gene}} | gene_id == {{gene}}) %>% pivot_longer(!c(gene_id, gene_name), names_to = "sample", values_to = "vsd") %>% mutate(vsd_zerod = vsd - min_vst) %>% left_join(select(metadata, sample, metric = all_of(variable)))
  if(filter_zeros == TRUE){vsd_counts.long = vsd_counts.long %>% filter(vsd > min_vst) }
  if(log10_convert == FALSE) {
    ggplot(vsd_counts.long, aes(x = metric, y = vsd_zerod)) + 
      labs(x = gsub("_", " ", variable), y = paste({{ gene }},"normalised expression") ) + 
      ggpubr::stat_cor(method = "pearson", aes(label = ..r.label..) ) +
      geom_smooth(method='lm') + 
      geom_point(size = 1) + # , alpha = 0.1
      geom_bin2d(bins = 100) +
      scale_fill_continuous(type = "viridis") +
      guides(fill = "none") +
      theme_bw() + theme_oz()
  } else {
    ggplot(vsd_counts.long, aes(x = metric, y = log10(vsd_zerod))) + 
      labs(x = gsub("_", " ", variable), y = paste({{ gene }},"log10 normalised expression") ) + 
      ggpubr::stat_cor(method = "pearson", aes(label = ..r.label..) ) +
      geom_smooth(method='lm') + 
      geom_point(size = 1) + # , alpha = 0.1
      geom_bin2d(bins = 100) +
      scale_fill_continuous(type = "viridis") +
      guides(fill = "none") +
      theme_bw() + theme_oz() + scale_y_log10()
  }
}


progeny_scores <- function(data, statistic = "stat",species="Human"){
  progeny_score = data %>% select(gene_name, all_of(statistic)) %>% drop_na(all_of(statistic)) %>% distinct(gene_name, .keep_all = TRUE) %>% column_to_rownames(var = "gene_name") %>% as.matrix() %>% 
    progeny(scale=TRUE, organism={{species}}, top = 100, perm = 10000, z_scores = TRUE) %>% t() %>% as_tibble(rownames = "Pathway") %>% mutate(Pathway = factor(Pathway))
  return(progeny_score)
}

# data(dorothea_hs, package = "dorothea")
# regulons <- dorothea_hs %>%  dplyr::filter(confidence %in% c("A", "B","C"))
dorothea_scores <- function(data, statistic = "stat",tf_n=25){
  dorothea_score = data %>% select(gene_name, all_of(statistic)) %>% drop_na(all_of(statistic)) %>% distinct(gene_name, .keep_all = TRUE) %>% column_to_rownames(var = "gene_name") %>% as.matrix() %>% 
    dorothea::run_viper(regulons, options =  list(minsize = 5, eset.filter = FALSE, cores = 4, verbose = TRUE, nes = TRUE)) %>% as_tibble(rownames = "gene_name") %>% dplyr::top_n({{tf_n}}, wt = abs(stat)) %>% mutate(gene_name = factor(gene_name))
  return(dorothea_score)
}


display_venn <- function(x, ...){
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
}


plotPCA.san <- function(object, intgroup = "condition", ntop = 500, PC_x = 1, PC_y = 2, returnData = FALSE) {
  rv <- rowVars(assay(object))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca <- prcomp(t(assay(object)[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(object)))) {  stop("the argument 'intgroup' should specify columns of colData(dds)") }
  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop = FALSE])
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = " : "))
  }
  else {
    colData(object)[[intgroup]]
  }
  d <- data.frame(PCX = pca$x[, PC_x], PCY = pca$x[, PC_y], group = group,  intgroup.df, name = colData(object)[,1])
  if (returnData) {
    attr(d, "percentVar") <- percentVar[PC_x:PC_y]
    return(d)
  }
  ggplot(data = d, aes_string(x = "PCX", y = "PCY", color = "group", label = "name")) + geom_point(size = 3) + 
    xlab(paste0("PC",PC_y,": ", round(percentVar[1] * 100), "% variance")) + ylab(paste0("PC",PC_y,": ", round(percentVar[2] * 100), "% variance")) + 
    coord_fixed() + 
    geom_text_repel(size=3) 
  
}

# Splicing Functions ----------------------
IRFinder.analysis <- function(metadata, sample.names = irfinder_names, design = ~ Condition + Condition:IRFinder, contrast_variable = "condition", contrast = NULL, ge.res = NULL, animal = "human"){
  experiment = metadata %>% mutate(file_paths.dir = file.path(database_dir,"splicing/irfinder",metadata[[sample.names]],"IRFinder-IR-dir.txt"), file_paths.nondir = file.path(database_dir,"splicing/irfinder",metadata[[sample.names]],"IRFinder-IR-nondir.txt")) %>% 
    rename(UniqueNames = {{sample.names}}) #SampleNames = {{sample.names}}, 
  if(FALSE %in% file.exists(experiment$file_paths.dir) == FALSE){
    cat(blue("\nImporting stranded IRFinder output: IRFinder-IR-dir.txt\n"))
    metaList = DESeqDataSetFromIRFinder(filePaths=experiment$file_paths.dir, designMatrix=experiment, designFormula = design)
  }else{
    cat(blue("\nImporting non-stranded IRFinder output: IRFinder-IR-nondir.txt\n"))
    metaList = DESeqDataSetFromIRFinder(filePaths=experiment$file_paths.nondir, designMatrix=experiment, designFormula = design)
  }
  cat(blue(paste0("Design: ~",design[2],"\n")))
  dds = DESeq(metaList$DESeq2Object)
  # vsd <- vst(dds, blind=TRUE) # useful for PCA plot of IR but often gives error less than 'nsub' rows with mean normalized count > 5. so use use varianceStabilizingTransformation directly
  # vsd.counts <- as_tibble(assay(vsd), rownames = "gene_id")  %>% left_join(gene2ens)
  print(resultsNames(dds))
  if(length(contrast_variable) == 1){
    denom = metadata %>% pull(contrast_variable) %>% fct_drop() %>% fct_unique() %>% first()
    numer = metadata %>% pull(contrast_variable) %>% fct_drop() %>% fct_unique() %>% last()
    cat(blue(paste0("IRFinder DESeq: Denomenator = ", denom, ", Numerator = ",numer,"\n")))
    res.denom = DESeq2::results(dds, name = paste0(contrast_variable,denom,".IRFinderIR")) # tests IR reads vs spliced reads in denomenator
    res.numer = DESeq2::results(dds, name = paste0(contrast_variable,numer,".IRFinderIR")) # get IR ratio in the numerator samples
    res = DESeq2::results(dds, contrast=list(paste0(contrast_variable,numer,".IRFinderIR"),  paste0(contrast_variable,denom,".IRFinderIR")))
    ir.res = cleanIRFinder(mutant_vs_ctrl = res, mutant = res.numer, ctrl = res.denom, species = animal)
  } else if(length(contrast_variable) == 2){
    cat(red(paste("Multi-factor design examining ", contrast_variable[1], "and ", contrast_variable[2],"\n")))
    res = DESeq2::results(dds, name = contrast)
    ir.res = cleanIRFinder.multidesign(IRFinder = res, species = animal)
  }
  if(is.null(ge.res) == TRUE){ 
    return(list(irfinder=ir.res)) # ir.vsd=vsd
  } else { 
    ir.ge.res = summariseIRjoinGE(cleanIRFinder = ir.res, GE = ge.res)
    return(list(irfinder=ir.res, ir.ge=ir.ge.res)) # ir.vsd=vsd, 
  }
}


# using output from IRFinder DESeq2 GLM with replicates https://github.com/williamritchie/IRFinder/wiki/Generalized-Linear-Model-with-Replicates
DESeqDataSetFromIRFinder = function(filePaths,designMatrix,designFormula){
  print("reading in IRFinder files with read_tsv")
  irfinder.tsv = filePaths %>% map(read_tsv, col_types = list(Chr="c",Name="c",Strand="c",Warnings="c",.default="n"))#, show_col_types = FALSE) # ! Can't combine `ExcludedBases` <double> and `ExcludedBases` <character>.
  names(irfinder.tsv) = designMatrix$UniqueNames # needs to be unique names e.g. dataset_sample
  irtab <- map_dfr(irfinder.tsv, bind_rows, .id = "UniqueNames") %>% 
    mutate(IntronDepth = round(IntronDepth), SpliceExact = round(SpliceExact), MaxSplice = round(pmax(SpliceLeft, SpliceRight)), irnames = paste0(Name,"/",Chr,":",Start,"-",End,":",Strand))
  # introndepth.duplicate = irtab %>% pivot_wider(irnames, names_from = UniqueNames, values_from = IntronDepth, names_glue = "intronDepth.{UniqueNames}", values_fn = length) %>% count(across(where(is.numeric)))
  # irtab %>% pivot_wider(irnames, names_from = UniqueNames, values_from = IntronDepth, names_glue = "intronDepth.{UniqueNames}", values_fill = list(IntronDepth = 0))
  IntronDepth = irtab %>% pivot_wider(irnames, names_from = UniqueNames, values_from = IntronDepth, names_glue = "intronDepth.{UniqueNames}") %>% column_to_rownames("irnames")
  SpliceExact = irtab %>% pivot_wider(irnames, names_from = UniqueNames, values_from = SpliceExact, names_glue = "totalSplice.{UniqueNames}") %>% column_to_rownames("irnames")
  MaxSplice = irtab %>% pivot_wider(irnames, names_from = UniqueNames, values_from = MaxSplice, names_glue = "maxSplice.{UniqueNames}") %>% column_to_rownames("irnames")
  group = bind_rows(mutate(designMatrix, IRFinder = "IR"), mutate(designMatrix, IRFinder = "Splice")) %>% mutate(IRFinder = factor(IRFinder, levels=c("Splice","IR")))
  counts.IRFinder = bind_cols(IntronDepth,MaxSplice) %>% drop_na()
  dd = DESeqDataSetFromMatrix(countData = counts.IRFinder, colData = group, design = designFormula)
  sizeFactors(dd)=rep(1,nrow(group))
  final=list(dd,IntronDepth,SpliceExact,MaxSplice)
  names(final)=c("DESeq2Object","IntronDepth","SpliceDepth","MaxSplice")
  return(final)
}


cleanIRFinder <- function(mutant_vs_ctrl = irfinder.res.ac_nuc_vcp_vs_ctrl, mutant = irfinder.res_ac_nuc.vcp, ctrl = irfinder.res_ac_nuc.ctrl, deseq_res = NULL, mass_spec = NULL, species = "human"){
  irfinder.mutant <- mutant %>% as_tibble(rownames = "Intron.GeneName.GeneID.Coords") %>% 
    mutate(IR_vs_Splice.mutant = 2^.$log2FoldChange, IRratio.mutant = IR_vs_Splice.mutant/(1+IR_vs_Splice.mutant), baseMean.mutant = baseMean) %>% dplyr::select(Intron.GeneName.GeneID.Coords, IRratio.mutant, baseMean.mutant) 
  irfinder.ctrl <- ctrl %>% as_tibble(rownames = "Intron.GeneName.GeneID.Coords") %>% 
    mutate(IR_vs_Splice.ctrl = 2^.$log2FoldChange, IRratio.ctrl = IR_vs_Splice.ctrl/(1+IR_vs_Splice.ctrl), baseMean.ctrl = baseMean) %>% dplyr::select(Intron.GeneName.GeneID.Coords, IRratio.ctrl, baseMean.ctrl) 
  IRFinder <- mutant_vs_ctrl %>% as_tibble(rownames = "Intron.GeneName.GeneID.Coords") %>% 
    left_join(irfinder.mutant, by = "Intron.GeneName.GeneID.Coords") %>% left_join(irfinder.ctrl, by = "Intron.GeneName.GeneID.Coords") %>% 
    mutate(gene_name = str_split_fixed(Intron.GeneName.GeneID.Coords, "/", 4)[,1], gene_id = str_split_fixed(Intron.GeneName.GeneID.Coords, "/", 4)[,2], 
           coords = str_split_fixed(Intron.GeneName.GeneID.Coords, "/",4)[,4], Chr = str_split_fixed(coords, ":", 3)[,1], Start_End = str_split_fixed(coords, ":", 3)[,2], Start = as.numeric(str_split_fixed(Start_End, "-", 2)[,1]), End = as.numeric(str_split_fixed(Start_End, "-", 2)[,2]), 
           Direction = str_split_fixed(coords, ":", 3)[,3], intron_type = str_split_fixed(Intron.GeneName.GeneID.Coords, "/", 4)[,3], intron_length = GenomicRanges::width(as.GenomicRange(Intron.GeneName.GeneID.Coords)),
           
           IRratio.diff = IRratio.mutant - IRratio.ctrl, IRdirection = case_when(IRratio.diff > 0 ~ "up", IRratio.diff < 0 ~ "down", TRUE ~ "unchanged"), 
           
           retained = case_when((IRratio.mutant >= 0.1 | IRratio.ctrl >= 0.1) ~ "retained", TRUE ~ "spliced"),
           # reliable = case_when(baseMean < 10 ~ "unreliable", TRUE ~ "reliable"), reliable_retained = case_when(reliable == "reliable" & retained == "retained" ~ "retained", TRUE ~ "spliced"),
           
           p0.05 = case_when(pvalue < 0.05 ~ "significant", TRUE ~ "none_significant"), p0.05 = factor(p0.05, levels = c("significant", "none_significant")),
           p0.05_IRdirection = case_when(IRdirection == "up" & pvalue < 0.05 ~ "up", IRdirection == "down" & pvalue < 0.05 ~ "down", TRUE ~ "none_significant"),
           # p0.05_reliable = case_when(pvalue < 0.05 & reliable == "reliable" ~ "significant", TRUE ~ "none_significant"),
           # p0.05_reliable_IRdirection = case_when(IRdirection == "up" & p0.05_reliable == "significant" ~ "IR up", IRdirection == "down" & p0.05_reliable  == "significant" ~ "IR down", TRUE ~ "none_significant"),
           
           padj0.05 = case_when(padj < 0.05 ~ "significant", TRUE ~ "none_significant"), padj0.05 = factor(padj0.05, levels = c("significant", "none_significant")),
           padj0.05_IRdirection = case_when(IRdirection == "up" & padj < 0.05 ~ "up", IRdirection == "down" & padj < 0.05 ~ "down", TRUE ~ "none_significant")) %>%
    # padj0.05_reliable = case_when(padj < 0.05 & reliable == "reliable" ~ "significant", TRUE ~ "none_significant"),
    # padj0.05_reliable_IRdirection = case_when(IRdirection == "up" & padj0.05_reliable == "significant" ~ "IR up", IRdirection == "down" & padj0.05_reliable  == "significant" ~ "IR down", TRUE ~ "none_significant")) %>% 
    select(-Start_End)
  print(paste(nrow(filter(IRFinder, pvalue < 0.05)), "significant IR events pvalue < 0.05"))
  print(table(IRFinder$p0.05_IRdirection))
  print(paste(nrow(filter(IRFinder, padj < 0.05)), "significant IR events padj < 0.05"))
  print(table(IRFinder$padj0.05_IRdirection))
  gr_introns <- as.GenomicRange(IRFinder$Intron.GeneName.GeneID.Coords)
  seqlevelsStyle(gr_introns) <- "UCSC"
  if(species == "human"){
    print("GC content based on human")
    IRFinder$gc_content <- GCcontent(Hsapiens, gr_introns)
  }else if(species == "mouse"){
    print("GC content based on mouse")
    IRFinder$gc_content <- GCcontent(Mmusculus, gr_introns)
  }else if(species == "rat"){
    print("GC content based on rat")
    IRFinder$gc_content <- GCcontent(Rnorvegicus, gr_introns)
  }
  if(species != "human"){
    IRFinder <- IRFinder %>% mutate(gene_name = toupper(gene_name))  # covert gene_names to upper case for matching to Deseq results & Human gene names
  }
  # IRFinder$phast_cons <- gscores(phast, gr_introns) # too slow to keep in function
  # if(!is.null(deseq_res)){
  #   deseq_res <- deseq_res %>% rename(ge.log2FoldChange = log2FoldChange, ge.padj = padj, ge.pvalue = pvalue, ge.baseMean = baseMean) %>% dplyr::select(-lfcSE, -stat)
  #   IRFinder <- IRFinder %>% left_join(deseq_res, by = c("gene_id", "gene_name")) %>% dplyr::mutate(ge.direction = case_when(ge.log2FoldChange > 0 ~ "up", TRUE ~ "down"),
  #                                          direction = case_when(ge.log2FoldChange > 0 & IRratio.diff < 0 ~ "IR down & expression up",
  #                                                                ge.log2FoldChange > 0 & IRratio.diff > 0 ~ "IR up & expression up",
  #                                                                ge.log2FoldChange < 0 & IRratio.diff < 0 ~ "IR down & expression down",
  #                                                                ge.log2FoldChange < 0 & IRratio.diff > 0 ~ "IR up & expression down"))
  # }
  # if(!is.null(mass_spec)){
  #   IRFinder <- IRFinder %>%  left_join(mass_spec, by = c("gene_name"="name"))
  #   IRFinder <- IRFinder %>% dplyr::mutate(massspec_direction = case_when(vcp_vs_ctrl_ratio > 0 ~ "up", vcp_vs_ctrl_ratio < 0 ~ "down"),
  #                                          massspec_ir_direction = case_when(vcp_vs_ctrl_ratio > 0 & IRratio.diff < 0 ~ "IR down & protein up", vcp_vs_ctrl_ratio > 0 & IRratio.diff > 0 ~ "IR up & protein up",
  #                                                                            vcp_vs_ctrl_ratio < 0 & IRratio.diff < 0 ~ "IR down & protein down", vcp_vs_ctrl_ratio < 0 & IRratio.diff > 0 ~ "IR up & protein down"),
  #                                          massspec_expression_ir_direction = case_when(direction == "IR down & expression up" & massspec_direction == "up" ~ "IR down & expression up & protein up",
  #                                                                                       direction == "IR down & expression up" & massspec_direction == "down" ~ "IR down & expression up & protein down", 
  #                                                                                       direction == "IR down & expression down" & massspec_direction == "up" ~ "IR down & expression up & protein down", 
  #                                                                                       direction == "IR down & expression down" & massspec_direction == "down" ~ "IR down & expression down & protein down",
  #                                                                                       direction == "IR up & expression down" & massspec_direction == "up" ~ "IR up & expression down & protein up",
  #                                                                                       direction == "IR up & expression up" & massspec_direction == "up" ~ "IR up & expression up & protein up", 
  #                                                                                       direction == "IR up & expression up" & massspec_direction == "down" ~ "IR up & expression up & protein down",
  #                                                                                       direction == "IR up & expression down" & massspec_direction == "down" ~ "IR up & expression down & protein up",
  #                                                                                       TRUE ~ "other"))
  # }
  return(IRFinder)
}



# normalise IR per gene and left_join to GE results
summariseIRjoinGE <- function(cleanIRFinder, GE){
  cleanIRFinder.summarised <- cleanIRFinder %>% group_by(gene_id, gene_name) %>% 
    summarise(all_intron_lengths = sum(intron_length), introns = dplyr::n(), 
              IR.log2FoldChange = sum( log2FoldChange * (intron_length / all_intron_lengths)), # normalise for intron length. # IRratio.diff.norm = sum( IRratio.diff * (intron_length / all_intron_lengths)), 
              IR.stat = sum( stat * (intron_length / all_intron_lengths)),
              IR.pvalue = min(pvalue, na.rm = TRUE), IR.padj = min(padj, na.rm = TRUE), 
              IR.baseMean = sum( baseMean * (intron_length / all_intron_lengths))) %>% ungroup %>% 
    mutate(IR.direction = case_when(IR.log2FoldChange > 0 ~ "IR up", IR.log2FoldChange < 0 ~ "IR down")) %>% 
    dplyr::select(gene_id, gene_name, IR.log2FoldChange, IR.stat, IR.padj, IR.pvalue, IR.baseMean, IR.direction) # , IRratio.diff.norm
  cleanIRFinder.summarised.GE <- GE %>% left_join(cleanIRFinder.summarised, by = c("gene_id", "gene_name")) %>% 
    mutate(GE.direction = case_when(log2FoldChange > 0 ~ "GE up", log2FoldChange < 0 ~ "GE down")) %>% 
    rename(GE.log2FoldChange = log2FoldChange, GE.stat = stat, GE.baseMean = baseMean, GE.padj = padj, GE.pvalue = pvalue) %>%
    mutate(significant = case_when(GE.pvalue < 0.05 & IR.pvalue < 0.05 ~ "both", GE.pvalue < 0.05 ~ "GE only", IR.pvalue < 0.05 ~ "IR only", TRUE ~ "None"), 
           direction = case_when(IR.direction == "IR up" & GE.direction == "GE down" ~ "IR up mRNA down", IR.direction == "IR up" & GE.direction == "GE up" ~ "IR up mRNA up",
                                 IR.direction == "IR down" & GE.direction == "GE down" ~ "IR down mRNA down", IR.direction == "IR down" & GE.direction == "GE up" ~ "IR down mRNA up")) %>%
    dplyr::select(gene_id, gene_name, GE.log2FoldChange, GE.stat, GE.padj, GE.pvalue, GE.baseMean, IR.log2FoldChange, IR.stat, IR.padj, IR.pvalue, IR.baseMean, direction, significant, IR.direction, GE.direction) # , IRratio.diff.norm
  print(paste("Normalising", nrow(cleanIRFinder), "IR events in", nrow(GE), "GE events"))
  return(cleanIRFinder.summarised.GE)
}




cleanIRFinder.multidesign <- function(IRFinder = irfinder.ac.nuc_vs_cyt.vcp_vs_ctrl, species = "human"){
  IRFinder <- IRFinder %>% as_tibble(rownames = "Intron.GeneName.GeneID.Coords") %>%
    mutate(gene_name = str_split_fixed(.$Intron.GeneName.GeneID.Coords, "/", 4)[,1], gene_id = str_split_fixed(.$Intron.GeneName.GeneID.Coords, "/", 4)[,2], 
           IRdirection = case_when(log2FoldChange > 0 ~ "up", log2FoldChange < 0 ~ "down", TRUE ~ "unchanged"), 
           reliable = case_when(baseMean < 10 ~ "unreliable", TRUE ~ "reliable"), 
           p0.05 = case_when(pvalue < 0.05 ~ "significant", TRUE ~ "none_significant"), p0.05 = factor(p0.05, levels = c("significant", "none_significant")),
           p0.05_IRdirection = case_when(IRdirection == "up" & pvalue < 0.05 ~ "up", IRdirection == "down" & pvalue < 0.05 ~ "down", TRUE ~ "none_significant"),
           p0.05_reliable = case_when(pvalue < 0.05 & reliable == "reliable" ~ "significant", TRUE ~ "none_significant"),
           p0.05_reliable_IRdirection = case_when(IRdirection == "up" & p0.05_reliable == "significant" ~ "up", IRdirection == "down" & p0.05_reliable  == "significant" ~ "down", TRUE ~ "none_significant"),
           padj0.05 = case_when(padj < 0.05 ~ "significant", TRUE ~ "none_significant"), padj0.05 = factor(padj0.05, levels = c("significant", "none_significant")),
           padj0.05_IRdirection = case_when(IRdirection == "up" & padj < 0.05 ~ "up", IRdirection == "down" & padj < 0.05 ~ "down", TRUE ~ "none_significant"),
           padj0.05_reliable = case_when(padj < 0.05 & reliable == "reliable" ~ "significant", TRUE ~ "none_significant"),
           padj0.05_reliable_IRdirection = case_when(IRdirection == "up" & padj0.05_reliable == "significant" ~ "up", IRdirection == "down" & padj0.05_reliable  == "significant" ~ "down", TRUE ~ "none_significant"),
           
           intron_type = str_split_fixed(.$Intron.GeneName.GeneID.Coords, "/", 4)[,3], intron_length = GenomicRanges::width(as.GenomicRange(Intron.GeneName.GeneID.Coords)), 
           coords = str_split_fixed(.$Intron.GeneName.GeneID.Coords, "/",4)[,4], Chr = str_split_fixed(.$Intron.GeneName.GeneID.Coords, "/|:", 6)[,4], Start = as.numeric(str_split_fixed(.$Intron.GeneName.GeneID.Coords, "/|:|-", 7)[,5]), 
           End = as.numeric(str_split_fixed(.$Intron.GeneName.GeneID.Coords, "/|:|-", 7)[,6]), Direction = str_split_fixed(.$Intron.GeneName.GeneID.Coords, ":", 3)[,3])
  print(paste(nrow(filter(IRFinder, pvalue < 0.05)), "significant IR events P < 0.05"))
  print(table(IRFinder$p0.05_IRdirection))
  print(paste(nrow(filter(IRFinder, padj < 0.05)), "significant IR events Padj < 0.05"))
  print(table(IRFinder$padj0.05_IRdirection))
  gr_introns <- as.GenomicRange(IRFinder$Intron.GeneName.GeneID.Coords)
  seqlevelsStyle(gr_introns) <- "UCSC"
  if(species == "human"){
    print("GC content based on human")
    IRFinder$gc_content <- GCcontent(Hsapiens, gr_introns)
  }else if(species == "mouse"){
    print("GC content based on mouse")
    IRFinder$gc_content <- GCcontent(Mmusculus, gr_introns)
  }else if(species == "rat"){
    print("GC content based on rat")
    IRFinder$gc_content <- GCcontent(Rnorvegicus, gr_introns)
  }
  if(species != "human"){
    IRFinder <- IRFinder %>% mutate(gene_name = toupper(gene_name))  # covert gene_names to upper case for matching to Deseq results & Human gene names
  }
  # IRFinder$phast_cons <- gscores(phast, gr_introns) # too slow to keep in function
  return(IRFinder)
}


# using output from Small Amount of Replicates via Audic and Claverie Test https://github.com/williamritchie/IRFinder/wiki/Small-Amounts-of-Replicates-via-Audic-and-Claverie-Test
cleanIRFinder.lowreplicates <- function(IRFinder = ac_who_vcp_vs_ctrl.IRFinder, deseq_res = NULL, mass_spec = NULL, species = "Hsapiens"){
  IRFinder <- IRFinder %>% mutate(gene_name = str_split_fixed(IRFinder$Intron.GeneName.GeneID, "/", 3)[,1], gene_id = str_split_fixed(IRFinder$Intron.GeneName.GeneID, "/", 3)[,2], intron_type = str_split_fixed(IRFinder$Intron.GeneName.GeneID, "/", 3)[,3],
                                  coords = paste(.$Chr, ":", .$Start, "-", .$End, sep = ""), Intron.GeneName.GeneID.Coords = paste(Intron.GeneName.GeneID, "/", coords, ":", Direction, sep = ""),
                                  intron_length = GenomicRanges::width(as.GenomicRange(Intron.GeneName.GeneID.Coords)),
                                  # reliable / stable intron expression
                                  A.IR.coverage = A.SplicesMax + A.IntronDepth,
                                  B.IR.coverage = B.SplicesMax + B.IntronDepth,
                                  IR.coverage = (A.IR.coverage + B.IR.coverage)/2,
                                  # reliable = case_when(A.IRok != "-" ~ "unreliable", B.IRok != "-" ~ "unreliable",  TRUE ~ "reliable"),
                                  reliable = case_when(A.IRok %in% c("LowCover","LowSplicing") ~ "unreliable", B.IRok %in% c("LowCover","LowSplicing") ~ "unreliable",  TRUE ~ "reliable"),
                                  # reliable = case_when(((A.SplicesMax + A.IntronDepth) > 100) & ((B.SplicesMax + B.IntronDepth) > 100) ~ "reliable", TRUE ~ "unreliable"), # based on https://github.com/lbroseus/SIRFindeR/blob/master/vignettes/SIRFindeR.pdf
                                  # reliable_lenient = case_when(((A.SplicesMax + A.IntronDepth) > 10) & ((B.SplicesMax + B.IntronDepth) > 10) ~ "reliable", TRUE ~ "unreliable"),
                                  retained = case_when((A.IRratio >= 0.1 | B.IRratio >= 0.1) ~ "retained", TRUE ~ "spliced"),
                                  reliable_retained = case_when(reliable == "reliable" & retained == "retained" ~ "retained", TRUE ~ "spliced"),
                                  # differential IR
                                  ## effect size
                                  IRratio.diff = A.IRratio - B.IRratio, # VCP - CTRL as per Dadi at IRFinder https://github.com/williamritchie/IRFinder/issues/99
                                  IR.lfc = log(IRFinder$A.IRratio/(1-IRFinder$A.IRratio)) -  log(IRFinder$B.IRratio/(1-IRFinder$B.IRratio)) %>% replace(., is.na(.), 0), # as per lucile broseus https://github.com/lbroseus/SIRFindeR/issues/1#issuecomment-668611143
                                  IR.l2fc = log2(IRFinder$A.IRratio / IRFinder$B.IRratio) %>% replace(., is.na(.), 0),
                                  IRdirection = case_when(IRratio.diff > 0 ~ "up", TRUE ~ "down"), # if positive then IRratio greater in VCP, negative then IRratio greater in CTRL
                                  ## significance
                                  p0.05 = case_when(p.diff < 0.05 ~ "significant", TRUE ~ "none_significant"), # classify significant IR events by p.diff < 0.05
                                  p0.05_IRdirection = case_when(IRdirection == "up" & p.diff < 0.05 ~ "IR up", IRdirection == "down" & p.diff < 0.05 ~ "IR down", TRUE ~ "none_significant"),
                                  p0.05_reliable = case_when(p.diff < 0.05 & reliable == "reliable" ~ "significant", TRUE ~ "none_significant"),
                                  p0.05_reliable_IRdirection = case_when(IRdirection == "up" & p0.05_reliable == "significant" ~ "IR up",
                                                                         IRdirection == "down" & p0.05_reliable  == "significant" ~ "IR down", TRUE ~ "none_significant"))
  # # adjust for multiple comparisons
  # padj = p.adjust(p.diff),
  # padj0.05 = case_when(padj < 0.05 ~ "significant", TRUE ~ "none_significant"), # classify significant IR events by padj < 0.05
  # padj0.05_IRdirection = case_when(IRdirection == "up" & padj < 0.05 ~ "IR up", IRdirection == "down" & padj < 0.05 ~ "IR down", TRUE ~ "none_significant"),
  # padj0.05_reliable = case_when(padj < 0.05 & reliable == "reliable" ~ "significant", TRUE ~ "none_significant"),
  # padj0.05_reliable_IRdirection = case_when(IRdirection == "up" & padj0.05_reliable == "significant" ~ "IR up",
  #                                           IRdirection == "down" & padj0.05_reliable  == "significant" ~ "IR down",TRUE ~ "none_significant"))
  gr_introns <- as.GenomicRange(IRFinder$Intron.GeneName.GeneID.Coords)
  seqlevelsStyle(gr_introns) <- "UCSC"
  if(species == "Hsapiens"){
    print("GC content based on Hsapiens")
    IRFinder$gc_content <- GCcontent(Hsapiens, gr_introns)
  }else if(species == "Mmusculus"){
    print("GC content based on Mmusculus")
    IRFinder$gc_content <- GCcontent(Mmusculus, gr_introns)
  }else if(species == "Rnorvegicus"){
    print("GC content based on Rnorvegicus")
    IRFinder$gc_content <- GCcontent(Rnorvegicus, gr_introns)
  }
  if(species != "Hsapiens"){
    IRFinder <- IRFinder %>% mutate(gene_name = toupper(gene_name))  # covert gene_names to upper case for matching to Deseq results & Human gene names
  }
  # IRFinder$phast_cons <- gscores(phast, gr_introns) # too slow to keep in function
  if(!is.null(deseq_res)){
    IRFinder <- IRFinder %>% left_join(deseq_res, by = c("gene_id", "gene_name")) # left join deseq2 results
    IRFinder <- IRFinder %>% dplyr::mutate(DGE.direction = case_when(log2FoldChange > 0 ~ "up", TRUE ~ "down"),
                                           direction = case_when(log2FoldChange > 0 & IRratio.diff < 0 ~ "IR down & expression up",
                                                                 log2FoldChange > 0 & IRratio.diff > 0 ~ "IR up & expression up",
                                                                 log2FoldChange < 0 & IRratio.diff < 0 ~ "IR down & expression down",
                                                                 log2FoldChange < 0 & IRratio.diff > 0 ~ "IR up & expression down"))
  }
  if(!is.null(mass_spec)){
    IRFinder <- IRFinder %>%  left_join(mass_spec, by = c("gene_name"="name"))
    IRFinder <- IRFinder %>% dplyr::mutate(massspec_direction = case_when(vcp_vs_ctrl_ratio > 0 ~ "up", vcp_vs_ctrl_ratio < 0 ~ "down"),
                                           massspec_ir_direction = case_when(vcp_vs_ctrl_ratio > 0 & IRratio.diff < 0 ~ "IR down & protein up", vcp_vs_ctrl_ratio > 0 & IRratio.diff > 0 ~ "IR up & protein up",
                                                                             vcp_vs_ctrl_ratio < 0 & IRratio.diff < 0 ~ "IR down & protein down", vcp_vs_ctrl_ratio < 0 & IRratio.diff > 0 ~ "IR up & protein down"),
                                           massspec_expression_ir_direction = case_when(direction == "IR down & expression up" & massspec_direction == "up" ~ "IR down & expression up & protein up",
                                                                                        direction == "IR down & expression up" & massspec_direction == "down" ~ "IR down & expression up & protein down", 
                                                                                        direction == "IR down & expression down" & massspec_direction == "up" ~ "IR down & expression up & protein down", 
                                                                                        direction == "IR down & expression down" & massspec_direction == "down" ~ "IR down & expression down & protein down",
                                                                                        direction == "IR up & expression down" & massspec_direction == "up" ~ "IR up & expression down & protein up",
                                                                                        direction == "IR up & expression up" & massspec_direction == "up" ~ "IR up & expression up & protein up", 
                                                                                        direction == "IR up & expression up" & massspec_direction == "down" ~ "IR up & expression up & protein down",
                                                                                        direction == "IR up & expression down" & massspec_direction == "down" ~ "IR up & expression down & protein up",
                                                                                        TRUE ~ "other"))
  }
  return(IRFinder)
}


# clean the pooled replicate quantified IRFinder output (not differential IR). Needed to calculate absolute IR per condition.
cleanIRFinderQuant <- function(IRFinder = ac_who_ctrl.IRFinder, species = "Hsapiens"){
  IRFinder <- IRFinder %>% mutate(gene_name = str_split_fixed(IRFinder$Name, "/", 3)[,1], gene_id = str_split_fixed(IRFinder$Name, "/", 3)[,2], intron_type = str_split_fixed(IRFinder$Name, "/", 3)[,3],
                                  coords = paste(.$Chr, ":", .$Start, "-", .$End, sep = ""), Name.Coords = paste(Name, "/", coords, ":", Strand, sep = ""),
                                  intron_length = GenomicRanges::width(as.GenomicRange(Name.Coords)),
                                  retained = case_when(IRratio >= 0.1 ~ "retained", TRUE ~ "spliced"),
                                  ExonDepth = case_when(SpliceRight > SpliceLeft ~ SpliceRight, TRUE ~ SpliceLeft), # aka SplicesMax
                                  # IR.coverage = ExonDepth + IntronDepth,
                                  # reliable = case_when(Warnings == "-" ~ "reliable", TRUE ~ "unreliable"), # case_when((ExonDepth + IntronDepth) > 100 ~ "reliable", TRUE ~ "unreliable"),
                                  reliable = case_when(Warnings %in% c("LowCover","LowSplicing") ~ "unreliable", TRUE ~ "reliable"), # case_when((ExonDepth + IntronDepth) > 100 ~ "reliable", TRUE ~ "unreliable"),
                                  # reliable_lenient = case_when((ExonDepth + IntronDepth) > 10 ~ "reliable", TRUE ~ "unreliable"),
                                  reliable_retained = case_when(reliable == "reliable" & retained == "retained" ~ "reliable-retained", TRUE ~ "unreliable-or-spliced"))
  # reliable_lenient_retained = case_when(reliable_lenient == "reliable" & retained == "retained" ~ "reliable-retained", TRUE ~ "unreliable-or-spliced"))
  gr_introns <- as.GenomicRange(IRFinder$Name.Coords)
  seqlevelsStyle(gr_introns) <- "UCSC"
  if(species == "Hsapiens"){
    print("GC content based on Hsapiens")
    IRFinder$gc_content <- GCcontent(Hsapiens, gr_introns)
  }else if(species == "Mmusculus"){
    print("GC content based on Mmusculus")
    IRFinder$gc_content <- GCcontent(Mmusculus, gr_introns)
  }else if(species == "Rnorvegicus"){
    print("GC content based on Rnorvegicus")
    IRFinder$gc_content <- GCcontent(Rnorvegicus, gr_introns)
  }
  if(species != "Hsapiens"){
    IRFinder <- IRFinder %>% mutate(gene_name = toupper(gene_name))  # covert gene_names to upper case for matching to Deseq results & Human gene names
  }
  return(IRFinder)
}

# Compare the mutant and control pooled replicate quantified IRFinder output (not differential IR)
cleanIRFinderQuantCompared <- function(mutant = ac_who_vcp.IR, control = ac_who_ctrl.IR, differentialIR= NULL, deseq_res = NULL, mass_spec = NULL){
  IRFinder.control <- control %>% dplyr::select(gene_id, gene_name, coords, Chr, Start, End, Name, Strand, intron_length, gc_content,
                                                "B.IRratio" = IRratio, "B.retained" = retained, "B.reliable" = reliable, "B.reliable_retained" = reliable_retained, "B.Coverage" = Coverage, "B.IntronDepth" = IntronDepth, 
                                                "B.SpliceLeft" = SpliceLeft, "B.SpliceRight" = SpliceRight, "B.SpliceExact" = SpliceExact, "B.Warnings" = Warnings, "B.intron_type" = intron_type)
  IRFinder.mutant <- mutant %>% dplyr::select(gene_id, gene_name, coords, Chr, Start, End, Name, Strand, intron_length, gc_content,
                                              "A.IRratio" = IRratio, "A.retained" = retained, "A.reliable" = reliable, "A.reliable_retained" = reliable_retained, "A.Coverage" = Coverage, "A.IntronDepth" = IntronDepth, 
                                              "A.SpliceLeft" = SpliceLeft, "A.SpliceRight" = SpliceRight, "A.SpliceExact" = SpliceExact, "A.Warnings" = Warnings, "A.intron_type" = intron_type)
  IRFinder <- IRFinder.control %>% full_join(IRFinder.mutant, by = c("gene_id", "gene_name", "coords", "Chr", "Start", "End", "Name", "Strand", "intron_length", "gc_content"))
  IRFinder <- IRFinder %>% mutate(
    # reliable event
    A.IR.coverage = max(A.SpliceLeft, A.SpliceRight) + A.IntronDepth,
    B.IR.coverage = max(B.SpliceLeft, B.SpliceRight) + B.IntronDepth,
    IR.coverage = (A.IR.coverage + B.IR.coverage)/2,
    reliable = case_when(A.Warnings %in% c("LowCover","LowSplicing") ~ "unreliable", B.Warnings %in% c("LowCover","LowSplicing") ~ "unreliable",  TRUE ~ "reliable"),
    retained = case_when((A.IRratio >= 0.1 | B.IRratio >= 0.1) ~ "retained", TRUE ~ "spliced"),
    reliable_retained = case_when(reliable == "reliable" & retained == "retained" ~ "retained", TRUE ~ "spliced"),
    # compare IR between conditions
    IRratio.diff = A.IRratio - B.IRratio, # Mutant - CTRL as per Dadi at IRFinder https://github.com/williamritchie/IRFinder/issues/99
    IR.lfc = log(IRFinder$A.IRratio/(1-IRFinder$A.IRratio)) -  log(IRFinder$B.IRratio/(1-IRFinder$B.IRratio)) %>% replace(., is.na(.), 0), # as per lucile broseus https://github.com/lbroseus/SIRFindeR/issues/1#issuecomment-668611143
    IR.l2fc = log2(IRFinder$A.IRratio / IRFinder$B.IRratio) %>% replace(., is.na(.), 0),
    IR.direction = case_when(IRratio.diff > 0 ~ "up", TRUE ~ "down")) # if positive then IRratio greater in VCP, negative then IRratio greater in CTRL
  if(!is.null(differentialIR)){
    IRFinder.differential <- differentialIR %>% dplyr::select(Chr, Start, End, "Name" = Intron.GeneName.GeneID, "Strand" = Direction, p.diff)
    IRFinder <- IRFinder %>% left_join(IRFinder.differential, by = c("Chr", "Start", "End", "Name", "Strand")) # left join Differential IR results
    IRFinder <- IRFinder %>% dplyr::mutate(p0.05 = case_when(p.diff < 0.05 ~ "significant", TRUE ~ "none_significant"), # classify significant IR events by p.diff < 0.05
                                           # p0.05_IRr0.1 = case_when(p.diff < 0.05 & (A.IRratio >= 0.1 | B.IRratio >= 0.1) ~ "significant", TRUE ~ "none_significant"), # significant IR events p.diff < 0.05 and IRratio 0.1 in either sample
                                           # p0.05_dIRr0.1 = case_when(p.diff < 0.05 & (abs(IRratio.diff) >= 0.1) ~ "significant", TRUE ~ "none_significant"), # significant IR events p.diff < 0.05 and delta IRratio of 0.1 between samples
                                           p0.05_IRdirection = case_when(IR.direction == "up" & p.diff < 0.05 ~ "IR up",
                                                                         IR.direction == "down" & p.diff < 0.05 ~ "IR down",
                                                                         TRUE ~ "none_significant"),
                                           # p0.05_dIRr0.1_IRdirection = case_when(IRdirection == "up" & p.diff < 0.05 & (abs(IRratio.diff) >= 0.1) ~ "IR up",
                                           # IRdirection == "down" & p.diff < 0.05 & (abs(IRratio.diff) >= 0.1) ~ "IR down",
                                           # TRUE ~ "none_significant"),
                                           p0.05_reliable = case_when(p.diff < 0.05 & reliable == "reliable" ~ "significant", TRUE ~ "none_significant"),
                                           # p0.05_reliable_lenient = case_when(p.diff < 0.05 & reliable_lenient == "reliable" ~ "significant", TRUE ~ "none_significant"),
                                           p0.05_reliable_IR.direction = case_when(IR.direction == "up" & p0.05_reliable == "significant" ~ "IR up",
                                                                                   IR.direction == "down" & p0.05_reliable  == "significant" ~ "IR down",
                                                                                   TRUE ~ "none_significant"))
  }
  if(!is.null(deseq_res)){
    IRFinder <- IRFinder %>% left_join(deseq_res, by = c("gene_id", "gene_name")) # left join deseq2 results
    IRFinder <- IRFinder %>% dplyr::mutate(DGE.direction = case_when(log2FoldChange > 0 ~ "up", TRUE ~ "down"),
                                           direction = case_when(log2FoldChange > 0 & IRratio.diff < 0 ~ "IR down & expression up",
                                                                 log2FoldChange > 0 & IRratio.diff > 0 ~ "IR up & expression up",
                                                                 log2FoldChange < 0 & IRratio.diff < 0 ~ "IR down & expression down",
                                                                 log2FoldChange < 0 & IRratio.diff > 0 ~ "IR up & expression down"))
  }
  if(!is.null(mass_spec)){
    IRFinder <- IRFinder %>%  left_join(mass_spec, by = c("gene_name"="name"))
    IRFinder <- IRFinder %>% dplyr::mutate(massspec_direction = case_when(vcp_vs_ctrl_ratio > 0 ~ "up", TRUE ~ "down"),
                                           massspec_ir_direction = case_when(vcp_vs_ctrl_ratio >= 0 & IRratio.diff <= 0 ~ "IR down & protein up", vcp_vs_ctrl_ratio > 0 & IRratio.diff > 0 ~ "IR up & protein up",
                                                                             vcp_vs_ctrl_ratio < 0 & IRratio.diff < 0 ~ "IR down & protein down", vcp_vs_ctrl_ratio < 0 & IRratio.diff > 0 ~ "IR up & protein down"),
                                           massspec_expression_ir_direction = case_when(direction == "IR down & expression up" & massspec_direction == "up" ~ "IR down & expression up & protein up",
                                                                                        direction == "IR down & expression up" & massspec_direction == "down" ~ "IR down & expression up & protein down", 
                                                                                        direction == "IR down & expression down" & massspec_direction == "up" ~ "IR down & expression up & protein down", 
                                                                                        direction == "IR down & expression down" & massspec_direction == "down" ~ "IR down & expression down & protein down",
                                                                                        direction == "IR up & expression down" & massspec_direction == "up" ~ "IR up & expression down & protein up",
                                                                                        direction == "IR up & expression up" & massspec_direction == "up" ~ "IR up & expression up & protein up", 
                                                                                        direction == "IR up & expression up" & massspec_direction == "down" ~ "IR up & expression up & protein down",
                                                                                        direction == "IR up & expression down" & massspec_direction == "down" ~ "IR up & expression down & protein up",
                                                                                        TRUE ~ "other"))
  }
  return(IRFinder)
}




# IRFinder.analysis <- function(metadata = treated.who.metadata, sample.names = "sample", condition = "stimulated", ctrl = "A0", mut = "A1", batch = "NA", file.var = "Sample limsid", ge.res = treated.who.a1_vs_a0$res, 
#                               irfinder.dir = "/camp/lab/luscomben/home/shared/projects/patani-collab/astrocyte-stimulated-bulk-rnaseq/splicing/irfinder"){
#   if(batch == "NA"){
#     experiment = metadata %>% select(SampleNames = sample.names, Condition = condition, file.name = file.var) %>% remove_rownames() %>% mutate(Condition = factor(Condition,levels=c(ctrl, mut)))
#     file_paths = experiment %>% mutate(file_paths = case_when(file.exists(file.path(irfinder.dir,file.name,"IRFinder-IR-dir.txt")) ~ file.path(irfinder.dir,file.name,"IRFinder-IR-dir.txt"), TRUE ~ file.path(irfinder.dir,file.name,"IRFinder-IR-nondir.txt"))) %>% pull(file_paths)
#     print("Importing stranded IRFinder output: IRFinder-IR-dir.txt or IRFinder-IR-nondir.txt")
#       metaList = DESeqDataSetFromIRFinder(filePaths=file_paths, designMatrix=experiment, designFormula =~ Condition + Condition:IRFinder)
#   } else{
#     print(paste0("Accounting for batch: ", batch))
#     experiment = metadata %>% select(SampleNames = sample.names, Condition = condition, Batch = batch, file.name = file.var) %>% remove_rownames() %>% mutate(Condition = factor(Condition,levels=c(ctrl, mut)))
#     file_paths = experiment %>% mutate(file_paths = case_when(file.exists(file.path(irfinder.dir,file.name,"IRFinder-IR-dir.txt")) ~ file.path(irfinder.dir,file.name,"IRFinder-IR-dir.txt"), TRUE ~ file.path(irfinder.dir,file.name,"IRFinder-IR-nondir.txt"))) %>% pull(file_paths)
#     print("Importing stranded IRFinder output: IRFinder-IR-dir.txt or IRFinder-IR-nondir.txt")
#     metaList = DESeqDataSetFromIRFinder(filePaths=file_paths, designMatrix=experiment, designFormula = ~ Batch + Condition + Condition:IRFinder)
#   }
#   dds = DESeq(metaList$DESeq2Object)
#   res.ctrl = DESeq2::results(dds, name = paste0("Condition",ctrl,".IRFinderIR")) # tests IR reads vs spliced reads in ctrl
#   res.mut = DESeq2::results(dds, name = paste0("Condition",mut,".IRFinderIR")) # get IR ratio in the mut samples
#   res = DESeq2::results(dds, contrast=list(paste0("Condition",mut,".IRFinderIR"),  paste0("Condition",ctrl,".IRFinderIR")))
#   ir.res = cleanIRFinder(mutant_vs_ctrl = res, mutant = res.mut, ctrl = res.ctrl)
#   ir.ge.res = summariseIRjoinGE(cleanIRFinder = ir.res, GE = ge.res)
#   return(list(irfinder=ir.res,ir.ge=ir.ge.res))
# }



mergeIRFinderQuantDESeq2 <- function(cleanIRFinderQuantCompared = NULL, deseq_res = NULL, mass_spec = NULL){
  cleanIRFinderQuantCompared.GE <- cleanIRFinderQuantCompared  %>% filter(!A.Warnings %in% c("LowCover","LowSplicing"), !B.Warnings %in% c("LowCover","LowSplicing")) %>%
    group_by(gene_name) %>% dplyr::summarise(all_intron_lengths = sum(intron_length), introns = dplyr::n(), 
                                             IRratio.diff.norm = sum( IRratio.diff * (intron_length / all_intron_lengths)), 
                                             IRratio.diff.mean = mean(IRratio.diff, na.rm = TRUE),
                                             IRratio.diff.max = IRratio.diff[which.max(abs(IRratio.diff))],
                                             # IRl2fc.norm = sum( IR.l2fc * (intron_length / all_intron_lengths)),
                                             # IRl2fc.mean = mean(IR.l2fc, na.rm = TRUE),
                                             # IRl2fc.max = IR.l2fc[which.max(abs(IR.l2fc))],
                                             IR.pval = min(p.diff, na.rm = TRUE), IR.p0.05_reliable = max(p0.05_reliable, na.rm = TRUE)) %>% ungroup %>% mutate(IR.pval = case_when(IR.pval == Inf ~ NaN, TRUE ~ IR.pval))
  if(!is.null(deseq_res)){
    deseq_res <- deseq_res %>% group_by(gene_name) %>% slice_max(baseMean, n = 1, with_ties = FALSE)
    cleanIRFinderQuantCompared.GE <- cleanIRFinderQuantCompared.GE %>% 
      full_join(filter(dplyr::select(deseq_res, gene_name, "mRNA.lfc" = log2FoldChange, "mRNA.FDR" = padj, "mRNA.baseMean" = baseMean)))
  }
  if(!is.null(mass_spec)){
    cleanIRFinderQuantCompared.GE <- cleanIRFinderQuantCompared.GE %>% 
      full_join(dplyr::select(mass_spec, "gene_name" = name, protein.lfc = vcp_vs_ctrl_ratio, "protein.pval" = vcp_vs_ctrl_p.val), by = "gene_name")
  }
  return(cleanIRFinderQuantCompared.GE)
}


read_majiq_deltapsi <- function(deltapsi.tsv, grp1 = "als", grp2 = "ctrl"){
  if(!file.exists(deltapsi.tsv)){ stop("MAJIQ deltapsi.tsv file not found") }
  read_tsv({{deltapsi.tsv}}) %>% 
    separate_rows(c("mean_dpsi_per_lsv_junction", "probability_changing","probability_non_changing",paste0({{grp1}},"_mean_psi"),paste0({{grp2}},"_mean_psi"),"junctions_coords"), sep = ";", convert = TRUE) %>% 
    mutate(gene_id = gsub("gene:","",gene_id), lsv_id = gsub("gene:","",lsv_id), padj = p.adjust(probability_changing)) %>% arrange(-abs(mean_dpsi_per_lsv_junction)) %>% left_join(gene2ens, by = "gene_id") %>% select(gene_name, everything())
}


read_voila_het <- function(voila.tsv, grp1 = "als", grp2 = "ctrl"){
  if(!file.exists(voila.tsv)){ stop("MAJIQ voila.tsv file not found") }
  read_tsv({{voila.tsv}}, skip = 15) %>% clean_names() %>%
    separate_rows(c("tnom", "ttest","wilcoxon","tnom_quantile","ttest_quantile","wilcoxon_quantile","tnom_score" ,"changing","nonchanging", paste0({{grp1}},"_median_psi"), paste0({{grp2}},"_median_psi"), paste0({{grp1}},"_percentile25_psi"), 
                    paste0({{grp2}},"_percentile25_psi"), paste0({{grp1}},"_percentile75_psi"), paste0({{grp2}},"_percentile75_psi"), "de_novo_junctions","junctions_coords"), sep = ";", convert = TRUE) %>%  #,"exons_coords"
    mutate(gene_id = gsub("gene:","",gene_id), lsv_id = gsub("gene:","",lsv_id),
           # tnom_padj = p.adjust(tnom), ttest_padj = p.adjust(ttest), wilcoxon_padj = p.adjust(wilcoxon), 
           deltapsi = !!sym(paste0({{grp1}},"_median_psi")) - !!sym(paste0({{grp2}},"_median_psi")), changing_dpsi0.2 = as.logical(changing), changing_dpsi0.1 = case_when(tnom < 0.05 & abs(deltapsi) > 0.1 ~ TRUE, TRUE ~ FALSE)) %>% # calculate dPSI = grp1_median_psi - grp2_median_psi
    group_by(lsv_id) %>% mutate(junction_id = row_number()) %>% ungroup %>% mutate(junction_id = paste(lsv_id,junction_id,sep="_"), lsv_junction_id = paste0(lsv_id,":j:",junctions_coords)) %>% # add junction_id unique variable for each lsv based on junction_coords
    select(gene_name, gene_id, lsv_id, lsv_junction_id, junction_id, deltapsi, tnom, ttest, wilcoxon, ends_with("_median_psi"), changing_dpsi0.2, changing_dpsi0.1, everything()) # ends_with("_padj"), 
}

read_voila_deltapsi <- function(voila.tsv, grp1 = "als", grp2 = "ctrl"){
  if(!file.exists(voila.tsv)){ stop("MAJIQ voila.tsv file not found") }
  read_tsv({{voila.tsv}}, skip = 10, col_types = list(num_exons="d")) %>% 
    separate_rows(c("mean_dpsi_per_lsv_junction", "probability_changing","probability_non_changing",paste0({{grp1}},"_mean_psi"),paste0({{grp2}},"_mean_psi"),"de_novo_junctions","junctions_coords"), sep = ";", convert = TRUE) %>% 
    mutate(gene_id = gsub("gene:","",gene_id), lsv_id = gsub("gene:","",lsv_id), log10_test_stat = -log10(1-probability_changing), 
           changing_dpsi0.1 = case_when(probability_changing > 0.9 & abs(mean_dpsi_per_lsv_junction) > 0.1 ~ TRUE, TRUE ~ FALSE), changing_dpsi0.2 = case_when(probability_changing > 0.9 & abs(mean_dpsi_per_lsv_junction) > 0.2 ~ TRUE, TRUE ~ FALSE)) %>% 
    group_by(lsv_id) %>% mutate(junction_id = row_number()) %>% ungroup %>% mutate(junction_id = paste(lsv_id,junction_id,sep="_"), lsv_junction_id = paste0(lsv_id,":j:",junctions_coords)) %>% # add junction_id unique variable for each lsv
    arrange(-abs(mean_dpsi_per_lsv_junction)) %>% select(gene_name, gene_id, lsv_id, lsv_junction_id, junction_id, deltapsi = mean_dpsi_per_lsv_junction, everything())
}

sep <- function(...) {
  dots <- list(...)
  n <- stringr::str_count(dots[[1]][[dots[[2]]]], "\\d+")
  separate_(..., into = sprintf("%s_%d", dots[[2]], 1:n), sep = ";", convert = TRUE)
}

read_majiq_voila.wide <- function(voila.tsv, grp1 = "als", grp2 = "ctrl"){
  if(!file.exists(voila.tsv)){ stop("MAJIQ voila.tsv file not found") }
  read_tsv({{voila.tsv}}, skip = 10) %>% 
    Reduce(f = sep, x = (c("mean_dpsi_per_lsv_junction", "probability_changing","probability_non_changing",paste0({{grp1}},"_mean_psi"),paste0({{grp2}},"_mean_psi"),"de_novo_junctions","junctions_coords"))) %>% 
    mutate(gene_id = gsub("gene:","",gene_id), lsv_id = gsub("gene:","",lsv_id), padj = p.adjust(probability_non_changing_1)) %>% arrange(-abs(mean_dpsi_per_lsv_junction_1))}




read_voila_modulizer <- function(moduilzer_path, skip_rows = 7, piechart = TRUE){
  alternative_intron.tsv = read_tsv(file.path(moduilzer_path,"alternative_intron.tsv"), skip = {{ skip_rows }}, col_types = list(module_id="c",gene_id="c",gene_name="c",seqid="c",strand="c",lsv_id="c",event_id="c",reference_exon_coord="c",spliced_with="c",spliced_with_coord="c",junction_name="c",junction_coord="c",module_event_combination="c", other_junctions="c", other_exons="c", exons_skipped_coords="c", exon1_coord="c", exon2_coordinate="c", complex="l", denovo="l",event_non_changing="l",event_changing="l",junction_changing="l",intron="l",.default="n"))  %>% clean_names() %>% mutate(modulizer = "alternative_intron", modulizer_simplified = "intron_retention")
  
  alt3and5prime.tsv = read_tsv(file.path(moduilzer_path,"alt3and5prime.tsv"), skip = {{ skip_rows }}, col_types = list(module_id="c",gene_id="c",gene_name="c",seqid="c",strand="c",lsv_id="c",event_id="c",reference_exon_coord="c",spliced_with="c",spliced_with_coord="c",junction_name="c",junction_coord="c",module_event_combination="c", other_junctions="c", other_exons="c", exons_skipped_coords="c", exon1_coord="c", exon2_coordinate="c", complex="l", denovo="l",event_non_changing="l",event_changing="l",junction_changing="l",intron="l",.default="n")) %>% clean_names() %>% mutate(modulizer = "alt3and5prime", modulizer_simplified = "alt3and5prime")
  alt5prime.tsv = read_tsv(file.path(moduilzer_path,"alt5prime.tsv"), skip = {{ skip_rows }}, col_types = list(module_id="c",gene_id="c",gene_name="c",seqid="c",strand="c",lsv_id="c",event_id="c",reference_exon_coord="c",spliced_with="c",spliced_with_coord="c",junction_name="c",junction_coord="c",module_event_combination="c", other_junctions="c", other_exons="c", exons_skipped_coords="c", exon1_coord="c", exon2_coordinate="c", complex="l", denovo="l",event_non_changing="l",event_changing="l",junction_changing="l",intron="l",.default="n"))  %>% clean_names() %>% mutate(modulizer = "alt5prime", modulizer_simplified = "alt5prime")
  p_alt5prime.tsv = read_tsv(file.path(moduilzer_path,"p_alt5prime.tsv"), skip = {{ skip_rows }}, col_types = list(module_id="c",gene_id="c",gene_name="c",seqid="c",strand="c",lsv_id="c",event_id="c",reference_exon_coord="c",spliced_with="c",spliced_with_coord="c",junction_name="c",junction_coord="c",module_event_combination="c", other_junctions="c", other_exons="c", exons_skipped_coords="c", exon1_coord="c", exon2_coordinate="c", complex="l", denovo="l",event_non_changing="l",event_changing="l",junction_changing="l",intron="l",.default="n"))  %>% clean_names() %>% mutate(modulizer = "p_alt5prime", modulizer_simplified = "alt5prime")
  alt3prime.tsv = read_tsv(file.path(moduilzer_path,"alt3prime.tsv"), skip = {{ skip_rows }}, col_types = list(module_id="c",gene_id="c",gene_name="c",seqid="c",strand="c",lsv_id="c",event_id="c",reference_exon_coord="c",spliced_with="c",spliced_with_coord="c",junction_name="c",junction_coord="c",module_event_combination="c", other_junctions="c", other_exons="c", exons_skipped_coords="c", exon1_coord="c", exon2_coordinate="c", complex="l", denovo="l",event_non_changing="l",event_changing="l",junction_changing="l",intron="l",.default="n"))  %>% clean_names() %>% mutate(modulizer = "alt3prime", modulizer_simplified = "alt3prime")
  p_alt3prime.tsv = read_tsv(file.path(moduilzer_path,"p_alt3prime.tsv"), skip = {{ skip_rows }}, col_types = list(module_id="c",gene_id="c",gene_name="c",seqid="c",strand="c",lsv_id="c",event_id="c",reference_exon_coord="c",spliced_with="c",spliced_with_coord="c",junction_name="c",junction_coord="c",module_event_combination="c", other_junctions="c", other_exons="c", exons_skipped_coords="c", exon1_coord="c", exon2_coordinate="c", complex="l", denovo="l",event_non_changing="l",event_changing="l",junction_changing="l",intron="l",.default="n"))  %>% clean_names() %>% mutate(modulizer = "p_alt3prime", modulizer_simplified = "alt3prime")
  
  alternate_first_exon.tsv = read_tsv(file.path(moduilzer_path,"alternate_first_exon.tsv"), skip = {{ skip_rows }}, col_types = list(module_id="c",gene_id="c",gene_name="c",seqid="c",strand="c",lsv_id="c",event_id="c",reference_exon_coord="c",spliced_with="c",spliced_with_coord="c",junction_name="c",junction_coord="c",module_event_combination="c", other_junctions="c", other_exons="c", exons_skipped_coords="c", exon1_coord="c", exon2_coordinate="c", complex="l", denovo="l",event_non_changing="l",event_changing="l",junction_changing="l",intron="l",.default="n"))  %>% clean_names() %>% mutate(modulizer = "alternate_first_exon", modulizer_simplified = "alternate_first_exon")
  p_alternate_first_exon.tsv = read_tsv(file.path(moduilzer_path,"p_alternate_first_exon.tsv"), skip = {{ skip_rows }}, col_types = list(module_id="c",gene_id="c",gene_name="c",seqid="c",strand="c",lsv_id="c",event_id="c",reference_exon_coord="c",spliced_with="c",spliced_with_coord="c",junction_name="c",junction_coord="c",module_event_combination="c", other_junctions="c", other_exons="c", exons_skipped_coords="c", exon1_coord="c", exon2_coordinate="c", complex="l", denovo="l",event_non_changing="l",event_changing="l",junction_changing="l",intron="l",.default="n"))  %>% clean_names() %>% mutate(modulizer = "p_alternate_first_exon", modulizer_simplified = "alternate_first_exon")
  alternate_last_exon.tsv = read_tsv(file.path(moduilzer_path,"alternate_last_exon.tsv"), skip = {{ skip_rows }}, col_types = list(module_id="c",gene_id="c",gene_name="c",seqid="c",strand="c",lsv_id="c",event_id="c",reference_exon_coord="c",spliced_with="c",spliced_with_coord="c",junction_name="c",junction_coord="c",module_event_combination="c", other_junctions="c", other_exons="c", exons_skipped_coords="c", exon1_coord="c", exon2_coordinate="c", complex="l", denovo="l",event_non_changing="l",event_changing="l",junction_changing="l",intron="l",.default="n"))  %>% clean_names() %>% mutate(modulizer = "alternate_last_exon", modulizer_simplified = "alternate_last_exon")
  p_alternate_last_exon.tsv =read_tsv(file.path(moduilzer_path,"p_alternate_last_exon.tsv"), skip = {{ skip_rows }}, col_types = list(module_id="c",gene_id="c",gene_name="c",seqid="c",strand="c",lsv_id="c",event_id="c",reference_exon_coord="c",spliced_with="c",spliced_with_coord="c",junction_name="c",junction_coord="c",module_event_combination="c", other_junctions="c", other_exons="c", exons_skipped_coords="c", exon1_coord="c", exon2_coordinate="c", complex="l", denovo="l",event_non_changing="l",event_changing="l",junction_changing="l",intron="l",.default="n"))  %>% clean_names() %>% mutate(modulizer = "p_alternate_last_exon", modulizer_simplified = "alternate_last_exon")
  
  tandem_cassette.tsv = read_tsv(file.path(moduilzer_path,"tandem_cassette.tsv"), skip = {{ skip_rows }}, col_types = list(module_id="c",gene_id="c",gene_name="c",seqid="c",strand="c",lsv_id="c",event_id="c",reference_exon_coord="c",spliced_with="c",spliced_with_coord="c",junction_name="c",junction_coord="c",module_event_combination="c", other_junctions="c", other_exons="c", exons_skipped_coords="c", exon1_coord="c", exon2_coordinate="c", complex="l", denovo="l",event_non_changing="l",event_changing="l",junction_changing="l",intron="l",.default="n"))  %>% clean_names() %>% mutate(modulizer = "tandem_cassette", modulizer_simplified = "skipped_exon")
  cassette.tsv = read_tsv(file.path(moduilzer_path,"cassette.tsv"), skip = {{ skip_rows }}, col_types = list(module_id="c",gene_id="c",gene_name="c",seqid="c",strand="c",lsv_id="c",event_id="c",reference_exon_coord="c",spliced_with="c",spliced_with_coord="c",junction_name="c",junction_coord="c",module_event_combination="c", other_junctions="c", other_exons="c", exons_skipped_coords="c", exon1_coord="c", exon2_coordinate="c", complex="l", denovo="l",event_non_changing="l",event_changing="l",junction_changing="l",intron="l",.default="n"))  %>% clean_names() %>% mutate(modulizer = "cassette_exon", modulizer_simplified = "skipped_exon")
  mutually_exclusive.tsv = read_tsv(file.path(moduilzer_path,"mutually_exclusive.tsv"), skip = {{ skip_rows }}, col_types = list(module_id="c",gene_id="c",gene_name="c",seqid="c",strand="c",lsv_id="c",event_id="c",reference_exon_coord="c",spliced_with="c",spliced_with_coord="c",junction_name="c",junction_coord="c",module_event_combination="c", other_junctions="c", other_exons="c", exons_skipped_coords="c", exon1_coord="c", exon2_coordinate="c", complex="l", denovo="l",event_non_changing="l",event_changing="l",junction_changing="l",intron="l",.default="n"))  %>% clean_names() %>% mutate(modulizer = "mutually_exclusive_exon", modulizer_simplified = "mutually_exclusive_exon")
  multi_exon_spanning.tsv = read_tsv(file.path(moduilzer_path,"multi_exon_spanning.tsv"), skip = {{ skip_rows }}, col_types = list(module_id="c",gene_id="c",gene_name="c",seqid="c",strand="c",lsv_id="c",event_id="c",reference_exon_coord="c",spliced_with="c",spliced_with_coord="c",junction_name="c",junction_coord="c",module_event_combination="c", other_junctions="c", other_exons="c", exons_skipped_coords="c", exon1_coord="c", exon2_coordinate="c", complex="l", denovo="l",event_non_changing="l",event_changing="l",junction_changing="l",intron="l",.default="n"))  %>% clean_names() %>% mutate(modulizer = "multi_exon_spanning", modulizer_simplified = "multi_exon_spanning")
  
  orphan_junction.tsv = read_tsv(file.path(moduilzer_path,"orphan_junction.tsv"), skip = {{ skip_rows }}, col_types = list(module_id="c",gene_id="c",gene_name="c",seqid="c",strand="c",lsv_id="c",event_id="c",reference_exon_coord="c",spliced_with="c",spliced_with_coord="c",junction_name="c",junction_coord="c",module_event_combination="c", other_junctions="c", other_exons="c", exons_skipped_coords="c", exon1_coord="c", exon2_coordinate="c", complex="l", denovo="l",event_non_changing="l",event_changing="l",junction_changing="l",intron="l",.default="n"))  %>% clean_names() %>% mutate(modulizer = "orphan_junction", modulizer_simplified = "orphan_junction")
  other.tsv = read_tsv(file.path(moduilzer_path,"other.tsv"), skip = {{ skip_rows }}, col_types = list(module_id="c",gene_id="c",gene_name="c",seqid="c",strand="c",lsv_id="c",event_id="c",reference_exon_coord="c",spliced_with="c",spliced_with_coord="c",junction_name="c",junction_coord="c",module_event_combination="c", other_junctions="c", other_exons="c", exons_skipped_coords="c", exon1_coord="c", exon2_coordinate="c", complex="l", denovo="l",event_non_changing="l",event_changing="l",junction_changing="l",intron="l",.default="n"))  %>% clean_names() %>% mutate(modulizer = "other", modulizer_simplified = "other")
  modulizer_bind = bind_rows(alt3and5prime.tsv, alt5prime.tsv, alternate_last_exon.tsv, cassette.tsv, mutually_exclusive.tsv, other.tsv, p_alt5prime.tsv, p_alternate_last_exon.tsv, alt3prime.tsv, alternate_first_exon.tsv, alternative_intron.tsv,  multi_exon_spanning.tsv, 
                             orphan_junction.tsv, p_alt3prime.tsv, p_alternate_first_exon.tsv, tandem_cassette.tsv) %>% #heatmap.tsv, summary.tsv, junctions.tsv
    select(modulizer, everything()) 
  modulizer_bind.summary = modulizer_bind %>% count(modulizer, modulizer_simplified) %>% mutate(modulizer = fct_reorder(modulizer, n, .desc = TRUE)) %>% arrange(desc(modulizer)) %>% mutate(pct = n/sum(n) *100, pct.label = case_when(pct > 5 ~ as.character(round(pct,1)), TRUE ~ ""), ypos = cumsum(pct) - 0.5*pct) #%>% print
  modulizer_bind.summary %>% print
  modulizer_simplfied.summary = modulizer_bind %>% count(modulizer_simplified) %>% mutate(modulizer_simplified = fct_reorder(modulizer_simplified, n, .desc = TRUE)) %>% arrange(desc(modulizer_simplified)) %>% mutate(pct = n/sum(n) *100, pct.label = case_when(pct > 5 ~ as.character(round(pct,1)), TRUE ~ ""), ypos = cumsum(pct) - 0.5*pct) #%>% print
  modulizer_simplfied.summary %>% print
  if(piechart == TRUE){  # Modulizer pie chart
    # ggplot(modulizer_bind.summary, aes(x = "", y = pct, fill = fct_reorder(modulizer, n, .desc = TRUE))) + geom_bar(width = 1, stat = "identity", color="white") + coord_polar("y", start=0) + theme_void() + theme(legend.title=element_blank()) + geom_text(aes(y = ypos, label = pct.label), color = "white", size=6)
    ggplot(modulizer_simplfied.summary, aes(x = "", y = pct, fill = fct_reorder(modulizer_simplified, n, .desc = TRUE))) + geom_bar(width = 1, stat = "identity", color="white") + coord_polar("y", start=0) + theme_void() + theme(legend.title=element_blank()) + geom_text(aes(y = ypos, label = pct.label), color = "white", size=6)
  }
  return(modulizer_bind)
}
# summary.tsv = read_tsv(file.path(moduilzer_path,"summary.tsv"), skip = {{ skip_rows }}, col_types = list(module_id="c",gene_id="c",gene_name="c",seqid="c",strand="c",lsv_id="c",event_id="c",reference_exon_coord="c",spliced_with="c",spliced_with_coord="c",junction_name="c",junction_coord="c",module_event_combination="c", other_junctions="c", other_exons="c", exons_skipped_coords="c", exon1_coord="c", exon2_coordinate="c", complex="l", denovo="l",event_non_changing="l",event_changing="l",junction_changing="l",intron="l",.default="n")) %>% clean_names() %>% mutate(modulizer = "summary")
# heatmap.tsv = read_tsv(file.path(moduilzer_path,"heatmap.tsv"), skip = {{ skip_rows }}, col_types = list(module_id="c",gene_id="c",gene_name="c",seqid="c",strand="c",lsv_id="c",event_id="c",reference_exon_coord="c",spliced_with="c",spliced_with_coord="c",junction_name="c",junction_coord="c",module_event_combination="c", other_junctions="c", other_exons="c", exons_skipped_coords="c", exon1_coord="c", exon2_coordinate="c", complex="l", denovo="l",event_non_changing="l",event_changing="l",junction_changing="l",intron="l",.default="n"))  %>% clean_names() %>% mutate(modulizer = "heatmap")
# junctions.tsv = read_tsv(file.path(moduilzer_path,"junctions.tsv"), skip = {{ skip_rows }}, col_types = list(module_id="c",gene_id="c",gene_name="c",seqid="c",strand="c",lsv_id="c",event_id="c",reference_exon_coord="c",spliced_with="c",spliced_with_coord="c",junction_name="c",junction_coord="c",module_event_combination="c", other_junctions="c", other_exons="c", exons_skipped_coords="c", exon1_coord="c", exon2_coordinate="c", complex="l", denovo="l",event_non_changing="l",event_changing="l",junction_changing="l",intron="l",.default="n"))  %>% clean_names() %>% mutate(modulizer = "junctions")




majiq_upset_plot <- function(data, title = NULL, subtitle = NULL, ymax = 1000, overlap_by = "mutation", overlap_level = "lsv_junction_id", order_counts = "degree", order_group_list = NULL, 
                             split_direction = FALSE, split_up = "up", split_down = "down", significance_threshold = "probability_changing", deltapsi_threshold = 0.2, plot_numbers = TRUE, dot_colour = NULL, colour_by = NULL, label_colours = NULL){
  if(significance_threshold == "probability_changing"){
    if(split_direction == TRUE) {
      data %<>% filter(!!sym(significance_threshold) > 0.9, abs(deltapsi) > {{deltapsi_threshold}}) %>% mutate(direction = case_when(deltapsi > 0 ~ split_up, TRUE ~ split_down), group = paste(!!sym(overlap_by), direction)) %>% select(group, !!sym(overlap_level)) %>% group_by(!!sym(overlap_level)) %>% summarise(group = list(group))
    } else {data %<>% filter(!!sym(significance_threshold) > 0.9, abs(deltapsi) > {{deltapsi_threshold}}) %>% select({{overlap_by}}, !!sym(overlap_level)) %>% group_by(!!sym(overlap_level)) %>% summarise(group = list(!!sym(overlap_by)))}
  } else {
    if(split_direction == TRUE) {
      data %<>% filter(!!sym(significance_threshold) < 0.05, abs(deltapsi) > {{deltapsi_threshold}}) %>% mutate(direction = case_when(deltapsi > 0 ~ split_up, TRUE ~ split_down), group = paste(!!sym(overlap_by), direction)) %>% select(group, !!sym(overlap_level)) %>% group_by(!!sym(overlap_level)) %>% summarise(group = list(group))
    } else {data %<>% filter(!!sym(significance_threshold) < 0.05, abs(deltapsi) > {{deltapsi_threshold}}) %>% select({{overlap_by}}, !!sym(overlap_level)) %>% group_by(!!sym(overlap_level)) %>% summarise(group = list(!!sym(overlap_by)))}
    
  }
  plot =  ggplot(data, aes(x = group)) + geom_bar() +
    geom_text(stat='count', aes(label=after_stat(count)), vjust=-0.2, hjust = 0.5, size = 3, angle = 0) +
    scale_y_continuous(breaks = NULL, lim = c(0, ymax), name = "") +
    theme_classic(base_size=8) +  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5), legend.title = element_blank(), legend.position = "none", axis.text=element_text(size=8), axis.title=element_text(size=10)) + labs(x = "", y = "") +
    theme_combmatrix(combmatrix.label.make_space = TRUE) + labs(title = title, subtitle = subtitle)
  if(order_counts == "degree"){plot = plot + scale_x_upset(order_by = "degree")} else {plot = plot + scale_x_upset(order_by = "freq")} # degree or freq
  if(!is.null(order_group_list)){plot = plot + scale_x_upset(sets = order_group_list)}
  if(plot_numbers == TRUE){plot = plot + geom_text(stat='count', aes(label=after_stat(count)), vjust=-0.2, hjust = 0.5, size = 3, angle = 0)}
  if(!is.null(dot_colour)){ plot = plot + theme_combmatrix(combmatrix.panel.point.color.fill = dot_colour) } #, combmatrix.panel.line.size = 0, combmatrix.label.make_space = FALSE)
  if(!is.null(label_colours)){ plot = plot + theme(axis.text.y = element_text(colour = label_colours))}
  return(plot)
}

cleanVASTTOOLS <- function(diff = NULL, comp = untx_vcp_vs_ctrl_compare_all_events){
  if(!is.null(diff)){
    txt = comp %>% full_join(diff, by = c("GENE", "EVENT")) %>% clean_names %>%
      mutate(type = case_when(grepl("HsaEX", event) ~ "exon skipping",grepl("HsaINT", event) ~ "intron retention", grepl("HsaALTA", event) ~ "alternate 3'", grepl("HsaALTD", event) ~ "alternate 5'"),
             compare.direction = case_when(d_psi > 0 ~ "up", d_psi < 0 ~ "down"), 
             d_psi0.1 = case_when(abs(d_psi) > 0.1 ~ "significant", TRUE ~ "none_significant"), 
             d_psi0.1_direction = case_when(d_psi > 0.1 ~ "up", d_psi < -0.1 ~ "down", TRUE ~ "none_significant"),
             diff.direction = case_when(e_d_psi > 0 ~ "up", e_d_psi < 0 ~ "down"), # if positive then IRratio greater in VCP, negative then IRratio greater in CTRL
             MV0.1 = case_when(mv_d_psi_at_0_95 > 0.1 ~ "significant", TRUE ~ "none_significant"), 
             MV0.1_dPSI0.1 = case_when(mv_d_psi_at_0_95 > 0.1 & abs(e_d_psi) >= 0.1 ~ "significant", TRUE ~ "none_significant"), # significant IR events p.diff < 0.05 and delta IRratio of 0.1 between samples
             MV0.1_direction = case_when(e_d_psi > 0 & mv_d_psi_at_0_95 > 0.1 ~ "up", e_d_psi < 0 & mv_d_psi_at_0_95 > 0.1 ~ "down", TRUE ~ "none_significant"),
             MV0.1_dPSI0.1_direction = case_when(e_d_psi > 0 & mv_d_psi_at_0_95 > 0.1 & (abs(e_d_psi) >= 0.1) ~ "up", e_d_psi < 0 & mv_d_psi_at_0_95 > 0.1 & (abs(e_d_psi) >= 0.1) ~ "down", TRUE ~ "none_significant"),
             chr = gsub("chr", "", str_split_fixed(coord, ":", 2)[,1]), start = as.numeric(str_split_fixed(coord, "(:|-)", 3)[,2]), end = as.numeric(str_split_fixed(coord, "-", 2)[,2])) %>% 
      select(gene, event, chr, start, end, coord, length, full_co, everything())
  } else {
    txt = comp %>% clean_names %>%
      mutate(type = case_when(grepl("HsaEX", event) ~ "exon skipping",grepl("HsaINT", event) ~ "intron retention", grepl("HsaALTA", event) ~ "alternate 3'", grepl("HsaALTD", event) ~ "alternate 5'"),
             compare.direction = case_when(d_psi > 0 ~ "up", d_psi < 0 ~ "down"), 
             d_psi0.1 = case_when(abs(d_psi) > 0.1 ~ "significant", TRUE ~ "none_significant"), 
             d_psi0.1_direction = case_when(d_psi > 0.1 ~ "up", d_psi < -0.1 ~ "down", TRUE ~ "none_significant"),
             chr = gsub("chr", "", str_split_fixed(coord, ":", 2)[,1]), start = as.numeric(str_split_fixed(coord, "(:|-)", 3)[,2]), end = as.numeric(str_split_fixed(coord, "-", 2)[,2])) %>% 
      select(gene, event, chr, start, end, coord, length, full_co, everything())
  }
  return(txt)
}


# labels = c("VIM", "COL1A1", "MYL6", "DDX5", "FN1", "FLNA", "HNRNPK", "FUS", "HNRNPL", "SFPQ", "EMP1", "HNRNPH1", "OSMR", "NONO", "HLA-A", "HLA-B", "HLA-C", "STAT3", "ANXA2", "SPARC", "HNRNPU", "HNRNPA2B1", "B2M", "PPIA", "LRP1","SRRM2", "SRSF5", "TNC", "SLC44A2", "CXCL2")


# Splicing Plots --------------------------------------

histogram_plot_irfinder <- 
  function(data = irfinder.ac_nuc_vcp_vs_ctrl, reliabilityThr = 10, legend = TRUE){
    data <- data %>% mutate(reliable.threshold = case_when(baseMean < reliabilityThr ~ "unreliable", TRUE ~ "reliable"))
    if(legend == TRUE){
      ggplot(filter(data, abs(IRratio.diff) > 0), aes(x = IRratio.diff, fill = reliable.threshold)) + 
        geom_histogram(position="dodge",binwidth=0.03) +  theme_bw() + 
        theme(strip.background = element_rect(colour="black", fill="white"), panel.grid = element_blank(), legend.position = "top", legend.box = "horizontal", legend.margin=margin(0.2,0.2,0.2,0.2), legend.box.margin=margin(-5,-5,-5,-5)) + labs(fill = expression(Coverage)) + 
        geom_vline(xintercept = 0.1, linetype = 3, colour = "darkgrey")  + geom_vline(xintercept = -0.1, linetype = 3, colour = "darkgrey")  + geom_vline(xintercept = 0, linetype = 1, colour = "darkgrey")  + 
        geom_hline(yintercept = 0, linetype = 1, colour = "black") +
        scale_fill_manual( values = c("dodgerblue2", "firebrick2") ) + scale_colour_manual( values = c("dodgerblue2", "firebrick2") ) +
        coord_cartesian(xlim=c(-0.5,0.5)) + xlab(label = expression(Delta~IR~ratio) ) + ylab(label = expression(IR~event~count) )
    } else{
      ggplot(filter(data, abs(IRratio.diff) > 0), aes(x = IRratio.diff, fill = reliable.threshold)) + 
        geom_histogram(position="dodge",binwidth=0.03) +  theme_bw() + 
        theme(strip.background = element_rect(colour="black", fill="white"), panel.grid = element_blank(), legend.position = "none") + 
        geom_vline(xintercept = 0.1, linetype = 3, colour = "darkgrey")  + geom_vline(xintercept = -0.1, linetype = 3, colour = "darkgrey")  + geom_hline(yintercept = 0, linetype = 1, colour = "black") +
        scale_fill_manual( values = c("dodgerblue2", "firebrick2") ) + scale_colour_manual( values = c("dodgerblue2", "firebrick2") )+
        coord_cartesian(xlim=c(-0.5,0.5)) + xlab(label = expression(Delta~IR~ratio) ) + ylab(label = expression(IR~event~count) )
    }
  }


histogram_plot_irfinder.log2FoldChange <- 
  function(data = irfinder.ac_nuc_vcp_vs_ctrl, reliabilityThr = 10, legend = TRUE, y.position = 100000){
    data <- data %>% mutate(reliable.threshold = case_when(baseMean < reliabilityThr ~ "unreliable", TRUE ~ "reliable"))
    if(legend == TRUE){
      ggplot(filter(data, abs(log2FoldChange) > 0), aes(x = log2FoldChange, fill = reliable.threshold)) + 
        geom_histogram(position="dodge",binwidth=0.1) +  theme_bw() + 
        theme(strip.background = element_rect(colour="black", fill="white"), panel.grid = element_blank(), legend.position = "top", legend.box = "horizontal", legend.margin=margin(0.2,0.2,0.2,0.2), legend.box.margin=margin(-5,-5,-5,-5)) + labs(fill = expression(Coverage)) + 
        geom_vline(xintercept = 0.1, linetype = 3, colour = "darkgrey")  + geom_vline(xintercept = -0.1, linetype = 3, colour = "darkgrey")  + geom_vline(xintercept = 0, linetype = 1, colour = "darkgrey")  + 
        geom_hline(yintercept = 0, linetype = 1, colour = "black") +
        scale_fill_manual( values = c("dodgerblue2", "firebrick2") ) + scale_colour_manual( values = c("dodgerblue2", "firebrick2") ) +
        coord_cartesian(xlim=c(-3,3)) + xlab(label = expression(LFC~IR) ) + ylab(label = expression(IR~event~count) ) #+
      # geom_segment(aes(x = 1, xend = 6, y= y.position, yend= y.position), arrow=arrow(length=unit(0.3,"cm")), color = "black") + geom_segment(aes(x = 1, xend = -6, y= y.position, yend= y.position),arrow=arrow(length=unit(0.3,"cm")), color = "black") +
      # annotate("text", x = 5, y = y.position + 10000, label = "IR Up", colour = "black") +   annotate("text", x = -5, y =  y.position + 10000, label = "IR Down", colour = "black")
    } else{
      ggplot(filter(data, abs(log2FoldChange) > 0), aes(x = log2FoldChange, fill = reliable.threshold)) + 
        geom_histogram(position="dodge",binwidth=0.1) +  theme_bw() + 
        theme(strip.background = element_rect(colour="black", fill="white"), panel.grid = element_blank(), legend.position = "none") + 
        geom_vline(xintercept = 0.1, linetype = 3, colour = "darkgrey")  + geom_vline(xintercept = -0.1, linetype = 3, colour = "darkgrey")  + geom_hline(yintercept = 0, linetype = 1, colour = "black") +
        scale_fill_manual( values = c("dodgerblue2", "firebrick2") ) + scale_colour_manual( values = c("dodgerblue2", "firebrick2") )+
        coord_cartesian(xlim=c(-3,3)) + xlab(label = expression(LFC~IR) ) + ylab(label = expression(IR~event~count) ) +
        geom_segment(aes(x = 1, xend = 6, y= y.position, yend= y.position), arrow=arrow(length=unit(0.3,"cm")), color = "black") + geom_segment(aes(x = 1, xend = -6, y= y.position, yend= y.position),arrow=arrow(length=unit(0.3,"cm")), color = "black") +
        annotate("text", x = 5, y = y.position + 10000, label = "IR Up", colour = "black") +   annotate("text", x = -5, y =  y.position + 10000, label = "IR Down", colour = "black")
    }
  }


histogram_plot_irfinder.lowreplicates <- 
  function(ctrl = ac_who_ctrl.IR, mutant = ac_who_vcp.IR, ctrl.name = "ctrl", mutant.name = "VCP", legend = TRUE){
    ctrl.filt <- ctrl %>% dplyr::select(Chr, Start, End, Name, Strand, IRratio, gene_name, reliable) %>% mutate(condition = ctrl.name)
    mutant.filt <- mutant %>% dplyr::select(Chr, Start, End, Name, Strand, IRratio, gene_name, reliable) %>% mutate(condition = mutant.name)
    ctrl_and_mutant.filt <- ctrl.filt %>% bind_rows(mutant.filt) %>% mutate(condition = factor(condition, levels = c(ctrl.name, mutant.name)))
    if(legend == TRUE){
      ggplot(filter(ctrl_and_mutant.filt, IRratio > 0), aes(x = IRratio, fill = reliable)) + 
        geom_histogram(position="dodge",binwidth=0.03) +  theme_bw() + 
        theme(strip.background = element_rect(colour="black", fill="white"), panel.grid = element_blank(), legend.position = "top", legend.box = "horizontal", legend.margin=margin(0.2,0.2,0.2,0.2), legend.box.margin=margin(-5,-5,-5,-5)) + labs(fill = expression(Coverage)) + 
        geom_vline(xintercept = 0.1, linetype = 3, colour = "darkgrey")  + geom_hline(yintercept = 0, linetype = 1, colour = "black") +
        facet_grid(. ~ condition) + scale_fill_manual( values = c("dodgerblue2", "firebrick2") ) + scale_colour_manual( values = c("dodgerblue2", "firebrick2") )
    } else{
      ggplot(filter(ctrl_and_mutant.filt, IRratio > 0), aes(x = IRratio, fill = reliable)) + 
        geom_histogram(position="dodge",binwidth=0.03) +  theme_bw() + 
        theme(strip.background = element_rect(colour="black", fill="white"), panel.grid = element_blank(), legend.position = "none") + 
        geom_vline(xintercept = 0.1, linetype = 3, colour = "darkgrey")  + geom_hline(yintercept = 0, linetype = 1, colour = "black") +
        facet_grid(. ~ condition)  + scale_fill_manual( values = c("dodgerblue2", "firebrick2") ) + scale_colour_manual( values = c("dodgerblue2", "firebrick2") )
    }
  }



IRratio.diff_histogram_plot_irfinder.lowreplicates <- 
  function(data = ac_who_vcp_vs_ctrl, reliabilityThr = 100){
    if(reliabilityThr == 100){
      ggplot(data, aes(x = IRratio.diff, fill = reliable)) + 
        geom_histogram(colour="black", position="dodge",binwidth=0.03) +  theme_bw() + 
        theme(strip.background = element_rect(colour="black", fill="white"), panel.grid = element_blank(), legend.position = "top") + labs(fill = expression(Coverage)) + 
        geom_vline(xintercept = 0.1, linetype = 2, colour = "darkred") + geom_vline(xintercept = -0.1, linetype = 2, colour = "darkred")  + geom_vline(xintercept = 0, linetype = 1, colour = "black") +
        scale_fill_manual( values = c("dodgerblue", "firebrick2") )
    } else{
      ggplot(data, aes(x = IRratio.diff, fill = reliable_lenient)) + 
        geom_histogram(colour="black", position="dodge",binwidth=0.01) +  theme_bw() + 
        theme(strip.background = element_rect(colour="black", fill="white"), panel.grid = element_blank(), legend.position = "top") + labs(fill = expression(Coverage)) + 
        geom_vline(xintercept = 0.1, linetype = 2, colour = "darkred") + geom_vline(xintercept = -0.1, linetype = 2, colour = "darkred")  + geom_vline(xintercept = 0, linetype = 1, colour = "black")  +
        scale_fill_manual( values = c("dodgerblue", "firebrick2") )
    }
  }


reliabile_plot_irfinder.lowreplicates <- 
  function(data = ac_who_ctrl, reliabilityThr = 100){
    ggplot(filter(data, IRratio > 0), aes( x = IRratio, y = (ExonDepth + IntronDepth),  label = gene_name)) + 
      theme(strip.background = element_rect(colour="black", fill="white"), panel.grid = element_blank()) + theme_bw() +
      geom_point(colour = "purple", alpha = 0.5, size = 0.5) +  theme_bw() +  xlab( expression(IR~ratio)  ) +  ylab( expression(log[2]~Intron~Exon~abundance) ) + scale_y_continuous(trans='log2') +
      geom_hline(yintercept = log2(reliabilityThr), linetype = 3, colour = "darkred") + geom_vline(xintercept = 0.1, linetype = 3, colour = "darkred")
  }


ma_plot_irfinder <- 
  function(data = mn.nuc_vs_cyt.vcp_vs_ctrl$irfinder, labels_list = mn.nuc_vs_cyt.vcp_vs_ctrl$irfinder.ma.labels, padj = TRUE){
    data <- data %>% mutate(p0.05_reliable_IRdirection = factor(p0.05_reliable_IRdirection, levels = c("up", "down", "none_significant")))
    if(padj == TRUE){
      labels_maplot <- data %>% filter(gene_name %in% labels_list, padj0.05_reliable == "significant") %>% arrange(-(abs(log2FoldChange))) %>% distinct(gene_name, .keep_all = TRUE)
      ggplot(data = data, aes(x = log10(baseMean), y = log2FoldChange, colour = padj0.05_reliable_IRdirection)) +  # NB if padj== NA then will not show
        geom_point(size = 0.5) +
        scale_colour_manual(values = c("up"="firebrick2", "down"="dodgerblue2", "none_significant"="darkgray") ) +  
        theme_bw() +  theme(strip.background = element_rect(colour="black", fill="white"), panel.grid = element_blank(), legend.position = "none") +
        xlab(expression( Log[10]~IR~event~coverage )) + ylab(expression(log[2]~Fold~Change)) + 
        geom_text_repel(data = labels_maplot, aes(x = log10(baseMean), y = log2FoldChange, label = gene_name), fontface = "italic", colour = "black", min.segment.length = 0, box.padding = 0.5, point.padding = 0, size=2.7, max.overlaps = Inf) +
        geom_hline(yintercept = 0, linetype = 1) + geom_vline(xintercept = log10(10), linetype = 3, colour = "darkgrey")
    }
    else{
      labels_maplot <- data %>% filter(gene_name %in% labels_list, p0.05_reliable == "significant") %>% arrange(-(abs(log2FoldChange))) %>% distinct(gene_name, .keep_all = TRUE)
      ggplot(data = data, aes(x = log10(baseMean), y = log2FoldChange, colour = p0.05_reliable_IRdirection)) +  # NB if pvalue == NA then will not show
        geom_point(size = 0.5) +
        scale_colour_manual(values = c("up"="firebrick2", "down"="dodgerblue2", "none_significant"="darkgray") ) +  
        theme_bw() +  theme(strip.background = element_rect(colour="black", fill="white"), panel.grid = element_blank(), legend.position = "none") +
        xlab(expression( Log[10]~IR~event~coverage )) + ylab(expression(log[2]~Fold~Change)) + 
        geom_text_repel(data = labels_maplot, aes(x = log10(baseMean), y = log2FoldChange, label = gene_name), fontface = "italic", colour = "black", min.segment.length = 0, box.padding = 0.5, point.padding = 0, size=2.7, max.overlaps = Inf) +
        geom_hline(yintercept = 0, linetype = 1) + geom_vline(xintercept = log10(10), linetype = 3, colour = "darkgrey")
    }
  }

ma_plot_irfinder.IRratio.diff <- 
  function(data = irfinder.ac_nuc_vcp_vs_ctrl, labels_list = irfinder.ac_nuc_vcp_vs_ctrl.ma.labels){
    labels_maplot <- data %>% filter(gene_name %in% labels_list) %>% filter(p0.05_reliable == "significant") %>% arrange(-(abs(IRratio.diff))) %>% distinct(gene_name, .keep_all = TRUE)
    ggplot(data = data, aes(x = log10(baseMean), y = IRratio.diff, colour = p0.05_reliable_IRdirection)) +  # NB if pvalue == NA then will not show
      geom_point(size = 0.5) +
      scale_colour_manual( values = c("firebrick2", "dodgerblue2", "darkgray") ) +  theme_bw() +  
      theme(panel.grid = element_blank(), legend.title = element_blank(), legend.position = "none") +
      xlab(expression( Log[10]~IR~event~coverage )) + ylab(expression(Delta~IR~ratio)) + 
      geom_text_repel(data = labels_maplot, aes(x = log10(baseMean), y = IRratio.diff, label = gene_name), fontface = "italic", colour = "black", min.segment.length = 0, box.padding = 0.5, point.padding = 0, size=2.7, max.overlaps = Inf) +
      geom_hline(yintercept = 0, linetype = 1) + geom_vline(xintercept = log10(10), linetype = 3, colour = "darkgrey")
  }


ma_plot_irfinder.lowreplicates <- 
  function(data = ac_ctrl_cyt_vs_nuc, xlabs = "Log[10] exon + intron coverage", labels_list = astrocyte.reactivity.labels){
    ggmaplot <- data %>% mutate(baseMean = log10(IR.coverage))
    labels_maplot <- data %>% mutate(baseMean = log10(IR.coverage)) %>% filter(p0.05_reliable == "significant") %>% filter(gene_name %in% labels_list) %>% arrange( -(abs(IRratio.diff)) ) %>% distinct(gene_name, .keep_all = TRUE)
    ggplot(data = ggmaplot, aes(x = baseMean, y = IRratio.diff, colour = p0.05_reliable_IRdirection)) +  geom_point(size = 0.5) + #xlim(1,5) + ylim(-3,3) +
      scale_colour_manual( values = c("firebrick2", "dodgerblue2", "darkgray") ) +  theme_bw() +  
      theme(panel.grid = element_blank(), axis.title.y = element_text(size = 10), legend.position = "none") +
      xlab(expression( Log[10]~coverage )) + 
      geom_text_repel(data = labels_maplot, aes(x = baseMean, y = IRratio.diff, label = gene_name), fontface = "italic", colour = "black", min.segment.length = 0, box.padding = 0.5, point.padding = 0,
                      size=2.5, max.overlaps = Inf) +
      geom_hline(yintercept = 0, linetype = 2) + geom_vline(xintercept = 1, linetype = 2)
  }


volcano_plot_irfinder <- 
  function(data = ac_ctrl_cyt_vs_nuc, xlabs = "Log[10] exon + intron coverage", labels_list = astrocyte.reactivity.labels, reliability = "reliable"){
    if(reliability == "reliable"){
      labels_plot <- data %>% filter(p0.05_reliable == "significant") %>% filter(gene_name %in% labels_list) %>% arrange( -(abs(IRratio.diff)) ) %>% distinct(gene_name, .keep_all = TRUE)
      ggplot(data = data, aes(x = IRratio.diff, y =  -log10(p.diff), colour = p0.05_reliable_IRdirection)) +  geom_point(size = 0.5) + #xlim(1,5) + ylim(-3,3) +
        scale_colour_manual( values = c("firebrick2", "dodgerblue", "darkgray") ) +  theme_bw() +  theme(panel.grid = element_blank(), legend.title = element_blank(), legend.position = "none") +
        ylab(expression( -log[10]~FDR )) + #xlab( expression( IR~log[2]~fold~change~CTRL~cytoplasmic:nuclear)) +
        geom_text_repel(data = labels_plot, aes(x = IRratio.diff, y = -log10(p.diff), label = gene_name), fontface = "italic", colour = "black", min.segment.length = 0, box.padding = 0.5, point.padding = 0, max.overlaps = Inf) +
        geom_hline(yintercept =1.30, linetype = 2) + geom_vline(xintercept =0, linetype = 2)
    } else if(reliability == "reliable_lenient"){
      labels_plot <- data %>% filter(p0.05_reliable_lenient == "significant") %>% filter(gene_name %in% labels_list) %>% arrange( -(abs(IRratio.diff)) ) %>% distinct(gene_name, .keep_all = TRUE)
      ggplot(data = data, aes(x = IRratio.diff, y =  -log10(p.diff), colour = p0.05_reliable_lenient_IRdirection)) +  geom_point(size = 0.5) + #xlim(1,5) + ylim(-3,3) +
        scale_colour_manual( values = c("firebrick2", "dodgerblue", "darkgray") ) +  theme_bw() +  theme(panel.grid = element_blank(), legend.title = element_blank(), legend.position = "none") +
        ylab(expression( -log[10]~FDR )) + #xlab( expression( IR~log[2]~fold~change~CTRL~cytoplasmic:nuclear)) +
        geom_text_repel(data = labels_plot, aes(x = IRratio.diff, y = -log10(p.diff), label = gene_name), fontface = "italic", colour = "black", min.segment.length = 0, box.padding = 0.5, point.padding = 0, max.overlaps = Inf) +
        geom_hline(yintercept =1.30, linetype = 2) + geom_vline(xintercept =0, linetype = 2)
    } else{
      labels_plot <- data %>% filter(p0.05 == "significant") %>% filter(gene_name %in% labels_list) %>% arrange( -(abs(IRratio.diff)) ) %>% distinct(gene_name, .keep_all = TRUE)
      ggplot(data = data, aes(x = IRratio.diff, y =  -log10(p.diff), colour = p0.05_reliable_IRdirection)) +  geom_point(size = 0.5) + #xlim(1,5) + ylim(-3,3) +
        scale_colour_manual( values = c("firebrick2", "dodgerblue", "darkgray") ) +  theme_bw() +  theme(panel.grid = element_blank(), legend.title = element_blank(), legend.position = "none") +
        ylab(expression( -log[10]~FDR )) + #xlab( expression( IR~log[2]~fold~change~CTRL~cytoplasmic:nuclear)) +
        geom_text_repel(data = labels_plot, aes(x = IRratio.diff, y = -log10(p.diff), label = gene_name), fontface = "italic", colour = "black", min.segment.length = 0, box.padding = 0.5, point.padding = 0, max.overlaps = Inf) +
        geom_hline(yintercept =1.30, linetype = 2) + geom_vline(xintercept =0, linetype = 2)
    }
  }

volcano_plot_irfinder_p0.05 <- 
  function(data = ac_ctrl_cyt_vs_nuc, xlabs = "Log[10] exon + intron coverage", labels_list = labels){
    labels_plot <- data %>% filter(p0.05 == "significant") %>% filter(gene_name %in% labels_list) %>% arrange( -(abs(IRratio.diff)) ) %>% distinct(gene_name, .keep_all = TRUE)
    ggplot(data = data, aes(x = IRratio.diff, y =  -log10(p.diff), colour = p0.05_IRdirection)) +  geom_point(size = 0.5) + #xlim(1,5) + ylim(-3,3) +
      scale_colour_manual( values = c("firebrick2", "dodgerblue", "darkgray") ) +  theme_bw() +  theme(panel.grid = element_blank(), legend.title = element_blank(), legend.position = "none") +
      ylab(expression( -log[10]~FDR )) + #xlab( expression( IR~log[2]~fold~change~CTRL~cytoplasmic:nuclear)) +
      geom_text_repel(data = labels_plot, aes(x = IRratio.diff, y = -log10(p.diff), label = gene_name), fontface = "italic", colour = "black", min.segment.length = 0, box.padding = 0.5, point.padding = 0, max.overlaps = Inf) +
      geom_hline(yintercept =1.30, linetype = 2) + geom_vline(xintercept =0, linetype = 2)
  }


# volcano_plot_irfinder(data = ac_ctrl_cyt_vs_nuc, labels_list = labels) + xlab( expression( IR~log[2]~fold~change~CTRL~cytoplasmic:nuclear )) + #xlim(1,5) + ylim(-3,3)


scatter_plot_irfinder <-
  function(data = ctrl_vcp_scatter, labels_list = labels, x = IR.lfc_ctrl, y = IR.lfc_vcp){
    # labels_plot <- data %>% filter(gene_name %in% labels) %>% arrange( -(abs(x)) ) %>% distinct(gene_name, .keep_all = TRUE)
    ggplot(data, aes( x = x, y = y)) + geom_point(colour = "purple", alpha = 0.5, size = 0.5) +  theme_bw() +  
      # geom_text_repel(data = labels_plot, aes(x = x, y = y, label = gene_name), fontface = "italic", colour = "black", min.segment.length = 0, box.padding = 0.5, point.padding = 0.5, max.overlaps = Inf) +
      theme(strip.background = element_rect(colour="black", fill="white"), panel.grid = element_blank()) +
      geom_hline(yintercept = 0, linetype = 3, colour = "darkgrey") + geom_vline(xintercept = 0, linetype = 3, colour = "darkgrey") + 
      geom_smooth(data = data, aes(x = x, y = y), method = "lm", se = FALSE, show.legend = TRUE, colour = "black") # + geom_abline(slope = -0.11, intercept = 0, linetype = 1)
  }



majiq_deltapsi_volcano_plot <- function(majiq.deltapsi, title = NULL, subtitle = NULL, junction_labels = NULL, gene_labels = NULL, distinct_labels = TRUE, effect_size = "deltapsi", significance_value = "probability_changing", ymax = 5, xmax = 1, xpos = 0.5, 
                                        numerator = NULL, denomenator = NULL, arrow_labels = TRUE, dpi = 300){
  
  majiq.deltapsi %<>% #rename(effect_size = all_of(effect_size), significance_value = all_of(significance_value), significance_threshold = all_of(significance_threshold)) %>%
    mutate(sig = case_when(!!sym(significance_value) > 0.9 & abs(!!sym(effect_size)) >= 0.4 ~ "sig_strong", !!sym(significance_value) > 0.9 & abs(!!sym(effect_size)) > 0.2 ~ "sig",  TRUE ~ "non_sig"), direction = ifelse(!!sym(effect_size) > 0, "up", "down"), class = paste(sig, direction), 
           effect_size = case_when(!!sym(effect_size) > xmax ~ Inf, !!sym(effect_size) < -xmax ~ -Inf, TRUE ~ !!sym(effect_size)), log10_test_stat = ifelse(is.infinite(log10_test_stat), ymax, log10_test_stat)) #case_when(-log10(1-!!sym(significance_value)) > ymax ~ 10^-ymax, TRUE ~ -log10(1-!!sym(significance_value))))
  print(count(majiq.deltapsi, class))
  
  plot <- ggplot(majiq.deltapsi, aes(x = effect_size, y = log10_test_stat)) +  #geom_point(aes(colour = class ), size = 0.5) +
    ggrastr::rasterise(geom_point(aes(colour = class), size = 0.5), dpi = dpi) +
    scale_colour_manual(values = c("non_sig up" = "gray", "non_sig down" = "gray",  "sig up" = "#F36F6F", "sig_strong up" = "#EB4445",  "sig down" = "#A6CEE3","sig_strong down" = "#79B1D3")) + #"#4F8FC4",  
    labs(y = expression(-log[10]~Test~Statistic), x = paste0("Delta PSI (",numerator," - ",denomenator,")"), title = title, subtitle = subtitle) +
    guides(colour = "none") +  scale_y_continuous(expand = c(0,0), limits = c(0,ymax)) + scale_x_continuous(limits = c(-xmax,xmax)) +
    theme_oz() +  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5), panel.border = element_blank(), axis.ticks = element_line(colour = "black") )
  
  if(!is.null(junction_labels)){
    plot <- plot + geom_text_repel(fontface = "italic", data = filter(majiq.deltapsi, lsv_junction_id %in% junction_labels, log10_test_stat > 1), aes(x = effect_size, y = log10_test_stat, label = gene_name), max.overlaps = Inf, min.segment.length = unit(0, "lines"),size = 2.3)
  }
  if(!is.null(gene_labels)){
    if(distinct_labels == TRUE){
      plot <- plot + geom_text_repel(fontface = "italic", data = distinct(arrange(filter(majiq.deltapsi, gene_name %in% gene_labels), wt = -log10_test_stat), gene_name, .keep_all = TRUE), aes(x = effect_size, y = log10_test_stat, label = gene_name), max.overlaps = Inf, min.segment.length = unit(0, "lines"), size = 2.3)
    } else if(distinct_labels == FALSE){
      plot <- plot + geom_text_repel(fontface = "italic", data = arrange(filter(majiq.deltapsi, gene_name %in% gene_labels), wt = -log10_test_stat), aes(x = effect_size, y = log10_test_stat, label = gene_name), max.overlaps = Inf, min.segment.length = unit(0, "lines"),size = 2.3)
    }
  }
  if(arrow_labels == TRUE){
    plot <- plot + 
      geom_segment(aes(x = xmax*0.3, xend = xmax*0.9, y= ymax*0.85, yend= ymax * 0.85), arrow=arrow(length=unit(0.3,"cm"))) + geom_segment(aes(x = -xmax*0.3, xend = -xmax*0.9, y= ymax*0.85, yend= ymax*0.85),arrow=arrow(length=unit(0.3,"cm"))) +
      annotate("text", x = xmax*0.6, y = ymax*0.90, label = paste0("Up in ",numerator), size = 3) + annotate("text", x = -xmax*0.6, y = ymax*0.90, label = paste0("Up in ",denomenator), size = 3)
  }
  return(plot)
}

majiq_het_volcano_plot <- function(voila.het, title = NULL, subtitle = NULL, junction_labels = NULL, gene_labels = NULL, distinct_labels = TRUE, effect_size = "deltapsi", significance_value = "tnom", significance_threshold = "tnom", ymax = 16.5, xmax = 1, xpos = 0.5, 
                                   numerator = NULL, denomenator = NULL, arrow_labels = TRUE, dpi = 300){
  voila.het <- voila.het %>% # rename(effect_size = all_of(effect_size), significance_value = all_of(significance_value), significance_threshold = all_of(significance_threshold)) %>%
    mutate(sig = case_when(!!sym(significance_threshold) < 0.05 & abs(!!sym(effect_size)) >= 0.2 ~ "sig_strong", !!sym(significance_threshold) < 0.05 & abs(!!sym(effect_size)) >= 0.1 ~ "sig", TRUE ~ "non_sig"),
           direction = ifelse(!!sym(effect_size) > 0, "up", "down"), class = paste(sig, direction), 
           effect_size = case_when(!!sym(effect_size) > xmax ~ Inf, effect_size < -xmax ~ -Inf, TRUE ~ !!sym(effect_size)), significance_value = case_when(-log10(!!sym(significance_value)) > ymax ~ 10^-ymax, TRUE ~ !!sym(significance_value)))
  print(count(voila.het, sig))
  # de_tally <- voila.het %>% group_by(sig, direction, class) %>% count() %>% filter(sig != "non_sig") %>% drop_na() %>%
  #   mutate(position = ifelse(sig == "sig", xpos, xmax-1), position = ifelse( direction == "down", -1 * position, position), n = formatC(n, format="f", big.mark=",", digits=0)) #%>% 
  # mutate(significance_value = case_when(-log10(significance_value) < ymax ~ Inf, TRUE ~ significance_value)) %>% # include dots & labels for when significance_value is out of coordinates of ymax
  
  plot <- ggplot(voila.het, aes(x = effect_size, y = -log10(significance_value))) +  #geom_point(aes(colour = class ), size = 0.5) +
    ggrastr::rasterise(geom_point(aes(colour = class), size = 0.5), dpi = dpi) +
    scale_colour_manual(values = c("non_sig up" = "gray", "non_sig down" = "gray", 
                                   "sig up" = "#F36F6F",
                                   "sig_strong up" = "#EB4445",
                                   "sig down" = "#A6CEE3",#"#4F8FC4",
                                   "sig_strong down" = "#79B1D3")) +
    labs(y = expression(-log[10]~Test~Statistic), x = paste0("Delta PSI (",numerator," - ",denomenator,")"), title = title, subtitle = subtitle) +
    guides(colour = "none") +
    scale_y_continuous(expand = c(0,0), limits = c(0,ymax)) +
    theme_oz() +  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5), panel.border = element_blank(), axis.ticks = element_line(colour = "black") ) +
    # geom_text(fontface = "bold", data = de_tally, aes(x = position, y = ymax - 0.5, label = n, colour = class), size = 2.5 ) +
    scale_x_continuous(limits = c(-xmax,xmax))
  
  if(!is.null(junction_labels)){
    plot <- plot + geom_text_repel(fontface = "italic", data = filter(voila.het, lsv_junction_id %in% junction_labels, significance_value < 0.05), aes(x = effect_size, y = -log10(significance_value), label = gene_name), max.overlaps = Inf, min.segment.length = unit(0, "lines"),size = 2.3)
  }
  
  if(!is.null(gene_labels)){
    if(distinct_labels == TRUE){
      plot <- plot + geom_text_repel(fontface = "italic", data = distinct(arrange(filter(voila.het, gene_name %in% gene_labels, significance_value < 0.05), wt = -abs(effect_size)), gene_name, .keep_all = TRUE), aes(x = effect_size, y = -log10(significance_value), label = gene_name), max.overlaps = Inf, min.segment.length = unit(0, "lines"), size = 2.3)#
    } else if(distinct_labels == FALSE){
      plot <- plot + geom_text_repel(fontface = "italic", data = filter(voila.het, gene_name %in% gene_labels, significance_value < 0.05), aes(x = effect_size, y = -log10(significance_value), label = gene_name), max.overlaps = Inf, min.segment.length = unit(0, "lines"),size = 2.3)
    }
  }
  if(arrow_labels == TRUE){
    plot <- plot + 
      geom_segment(aes(x = xmax*0.25, xend = xmax*0.95, y= ymax*0.88, yend= ymax * 0.88), arrow=arrow(length=unit(0.3,"cm"))) + geom_segment(aes(x = -xmax*0.25, xend = -xmax*0.95, y= ymax*0.88, yend= ymax*0.88),arrow=arrow(length=unit(0.3,"cm"))) +
      annotate("text", x = xmax*0.6, y = ymax*0.93, label = paste0("Up in ",numerator), size = 3) + annotate("text", x = -xmax*0.6, y = ymax*0.93, label = paste0("Up in ",denomenator), size = 3)
  }
  return(plot)
}



dpsi_majiq_scatter <- function(model_list1, model_list2, mutation1, mutation2, col = "deltapsi", plot_corr_line = TRUE, gene_labels = NULL){
  require(ggrepel)
  res <- dplyr::left_join(model_list1[[mutation1]], model_list2[[mutation2]],  by = c("gene_name", "junction_id") , suffix = c(".1", ".2") ) 
  x_string <- paste0(col, ".1")
  y_string <- paste0(col, ".2")
  plot <- res %>%
    ggplot(aes_string(x = x_string, y = y_string )) + 
    labs(x = gsub("_", " ", mutation1), y = gsub("_", " ", mutation2) ) + 
    theme_bw() + 
    ggpubr::stat_cor(method = "pearson", aes(label = ..r.label..) ) +
    geom_hline(yintercept = 0, linetype = 3) + geom_vline(xintercept = 0, linetype = 3) +
    geom_abline(slope =1 , intercept = 0, linetype = 3) +    
    #geom_point(size = 1, alpha = 0.1) + 
    geom_bin2d(bins = 100) +
    scale_fill_continuous(type = "viridis") +
    guides(fill = "none") +
    #xlim(-10,10) + ylim(-8,8) +
    theme_oz()
  
  if(plot_corr_line == TRUE){ plot = plot + geom_smooth(method='lm') }
  
  if(!is.null(gene_labels)){
    plot <- plot +
      geom_point(data = filter(res, gene_name %in% gene_labels), aes_string(x = x_string, y = y_string), colour = "red") +
      geom_text_repel(data = filter(res, gene_name %in% gene_labels), aes_string(x = x_string, y = y_string, label = "gene_name"), colour = "red",size = 2.3)
  }
  return(plot)
}


# CLIP xlinks score
# irfinder.introns <- ipsc_mn_als_datasets$irfinder %>% distinct(Intron.GeneName.GeneID.Coords, intron_length) # 249,624 unique introns
# irfinder.introns.gr.NCBI <- as.GenomicRange(irfinder.introns$Intron.GeneName.GeneID.Coords)
# seqlevelsStyle(irfinder.introns.gr.NCBI) <- "NCBI"
# 
# # HepG2 eCLIP
# clip_HepG2 <- list.files("/camp/lab/luscomben/home/users/ziffo/genomes/rbp-databases/clip-bedgraphs/HepG2/")
# clip_HepG2 <- clip_HepG2[!grepl("3nt_peaks",clip_HepG2)] # only assess raw crosslink files without any peak calling for enrichment
# clip_HepG2 <- gsub(".merged.bed.gz", "", clip_HepG2)
# clip_HepG2_score <- list()
# for (i in seq_along(clip_HepG2)) { ## sum the iclip xlink scores for each RBP in HepG2 in a for loop
#   print(clip_HepG2[[i]])
#   bed <- import(paste("/camp/lab/luscomben/home/users/ziffo/genomes/rbp-databases/clip-bedgraphs/HepG2/", clip_HepG2[[i]], ".merged.bed.gz",sep=""))
#   seqlevelsStyle(bed) <- "NCBI" # remove "chr"
#   clip_HepG2_score[[i]] <-  GenomicRanges::countOverlaps(query = irfinder.introns.gr.NCBI, subject=bed)
#   names(clip_HepG2_score)[i] <- clip_HepG2[[i]]
# }
# saveRDS(clip_HepG2_score, "/camp/lab/luscomben/home/users/ziffo/genomes/rbp-databases/clip-bedgraphs/intron.clip_HepG2_score.rds")
# # clip_HepG2_score <- readRDS("/camp/lab/luscomben/home/users/ziffo/genomes/rbp-databases/clip-bedgraphs/intron.clip_HepG2_score.rds")
# 
# # K562 eCLIP
# clip_K562 <- list.files("/camp/lab/luscomben/home/users/ziffo/genomes/rbp-databases/clip-bedgraphs/K562/") # 90 RBPs
# clip_K562 <- clip_K562[!grepl("3nt_peaks",clip_K562)] # only assess raw crosslink files without any peak calling for enrichment
# clip_K562 <- gsub(".merged.bed.gz", "", clip_K562)
# clip_K562 <- clip_K562[!clip_K562 %in% clip_HepG2] # keep only 42 / 90  RBPs without scores in HepG2 cells
# clip_K562_score <- list()
# for (i in seq_along(clip_K562)) {
#   print(clip_K562[[i]])
#   bed <- import(paste("/camp/lab/luscomben/home/users/ziffo/genomes/rbp-databases/clip-bedgraphs/K562/", clip_K562[[i]], ".merged.bed.gz", sep=""))
#   seqlevelsStyle(bed) <- "NCBI" # remove "chr"
#   clip_K562_score[[i]] <-  GenomicRanges::countOverlaps(query = irfinder.introns.gr.NCBI, subject=bed)
#   names(clip_K562_score)[i] <- clip_K562[[i]]
# }
# saveRDS(clip_K562_score, "/camp/lab/luscomben/home/users/ziffo/genomes/rbp-databases/clip-bedgraphs/intron.clip_K562_score.rds")
# # clip_K562_score <- readRDS("/camp/lab/luscomben/home/users/ziffo/genomes/rbp-databases/clip-bedgraphs/intron.clip_K562_score.rds")
# 
# # Hela-hg19 iCLIP
# gr_introns <- as.GenomicRange(ipsc_mn_als_datasets$irfinder$Intron.GeneName.GeneID.Coords)
# seqlevelsStyle(gr_introns) <- "UCSC"
# library(rtracklayer) # https://genviz.org/module-01-intro/0001/06/02/liftoverTools/
# ch = import.chain("/camp/lab/luscomben/home/users/ziffo/genomes/ucsc/hg38ToHg19.over.chain") # import chain for liftOver http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
# gr_introns.UCSC.hg19 = liftOver(gr_introns, ch)
# 
# clip_Hela <- list.files("/camp/lab/luscomben/home/users/ziffo/genomes/rbp-databases/clip-bedgraphs/Hela_hg19/xlinks/") # 21 RBPs
# clip_Hela <- gsub(".bed", "", clip_Hela)
# clip_Hela <- clip_Hela[!clip_Hela %in% c(clip_HepG2, clip_K562)] # keep only 12 / 21  RBPs without scores in HepG2 or K562 cells
# clip_Hela_score <- list()
# for (i in seq_along(clip_Hela)) {
#   print(clip_Hela[[i]])
#   bed <- import(paste0("/camp/lab/luscomben/home/users/ziffo/genomes/rbp-databases/clip-bedgraphs/Hela_hg19/xlinks/", clip_Hela[[i]], ".bed"))
#   clip_Hela_score[[i]] <-  GenomicRanges::countOverlaps(query = gr_introns.UCSC.hg19, subject=bed)
#   names(clip_Hela_score)[i] <- clip_Hela[[i]]
# }
# saveRDS(clip_Hela_score, "/camp/lab/luscomben/home/users/ziffo/genomes/rbp-databases/clip-bedgraphs/intron.clip_Hela_score.rds")
# # clip_Hela_score <- readRDS("/camp/lab/luscomben/home/users/ziffo/genomes/rbp-databases/clip-bedgraphs/intron.clip_Hela_score.rds")
# 
# clip_HepG2_score <- readRDS("/camp/lab/luscomben/home/users/ziffo/genomes/rbp-databases/clip-bedgraphs/intron.clip_HepG2_score.rds")
# clip_K562_score <- readRDS("/camp/lab/luscomben/home/users/ziffo/genomes/rbp-databases/clip-bedgraphs/intron.clip_K562_score.rds")
# clip_Hela_score <- readRDS("/camp/lab/luscomben/home/users/ziffo/genomes/rbp-databases/clip-bedgraphs/intron.clip_Hela_score.rds")
# 
# iclip_score.HepG2_K562_Hela.normalised <-  ipsc_mn_als_datasets$irfinder %>% select(Intron.GeneName.GeneID.Coords) %>% # add gene_id to join DESeq2
#   bind_cols(bind_rows(clip_HepG2_score)) %>% # collapse iclip score list into a single df. Columns = RBPs. Rows = genes
#   bind_cols(bind_rows(clip_K562_score)) %>% #select(all_of(c(clip_HepG2, clip_K562)))
#   bind_cols(bind_rows(clip_Hela_score)) %>%
#   mutate(across(c(clip_HepG2, clip_K562, clip_Hela), function(x) (x / ipsc_mn_als_datasets$irfinder$intron_length))) %>% # normalise all  scores for intron length 
#   mutate(mean_iclip_score = rowMeans(select(., all_of(c(clip_HepG2, clip_K562, clip_Hela))), na.rm = TRUE), # take mean of RBP score for each intron
#          sum_iclip_score = rowSums(select(., all_of(c(clip_HepG2, clip_K562, clip_Hela))), na.rm = TRUE)) # take sum of RBP score for each intron - give same pattern of results as rowSums
# saveRDS(iclip_score.HepG2_K562_Hela.normalised, "/camp/lab/luscomben/home/users/ziffo/genomes/rbp-databases/clip-bedgraphs/iclip_score.HepG2_K562_Hela.intron.normalised.rds")


# Coverage Plots ---------------
# ggbio http://www.sthda.com/english/wiki/ggbio-visualize-genomic-data   https://www.bioconductor.org/packages/release/bioc/vignettes/ggbio/inst/doc/ggbio.pdf
# Gviz is better # https://blog.liang2.tw/posts/2016/01/plot-seq-depth-gviz/ http://www.sthda.com/english/wiki/gviz-visualize-genomic-data # Dave Tang blog: https://davetang.org/muse/2013/10/03/using-gviz/
options(ucscChromosomeNames=FALSE) # crucial to ensure AlignmentsTrack plot works
# gtrack <- GenomeAxisTrack() # genome axis track
# itrack <- IdeogramTrack(genome = "hg38", chromosome = 1) # add chromosome ideogram # Error in readHTMLTable(url)[[1L]] : subscript out of bounds
# Transcript annotations
### TURN ON ONLY WHEN NEEDED
# txdb <- makeTxDbFromGFF(here(camp_path,"/genomes/ensembl/Homo_sapiens.GRCh38.99.chr_patch_hapl_scaff.gtf"), format="gtf") # make txdb from GTF - turn on only when needed
# seqlevels(txdb, pruning.mode="coarse") <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y") # removes all the other chromsome annotations & scaffolds which are also in the BAM

getGeneIDsFromTxDb <- function(gr, txdb){
  stopifnot(is(gr, "GRanges"))
  stopifnot(length(gr)>0)
  stopifnot(is(txdb, "TxDb"))
  if(length(gr)>1){
    warning("The length of gr is greater than 1. Only first genomic location will be used.")
    gr <- gr[1]
  }
  genes <- genes(txdb, columns="gene_id")
  genes <- subsetByOverlaps(genes, gr)
  return(genes$gene_id)
}

# Track viewer https://www.bioconductor.org/packages/release/bioc/vignettes/trackViewer/inst/doc/trackViewer.html#introduction
track_viewer_nuc_cyt <- function(level = "gene", gene_name = NULL, transcript_id = NULL, cond1.nuc.bam = ctrl_nuclear.bam, cond1.cyt.bam = ctrl_cytoplasm.bam, cond1.label = "CTRL", cond2.cyt.bam = NULL, cond2.nuc.bam = NULL, cond2.label = "MUT", 
                                 cond3.cyt.bam = NULL, cond3.nuc.bam = NULL, cond3.label = "MUT2", cond4.cyt.bam = NULL, cond4.nuc.bam = NULL, cond4.label = "MUT3",  cond5.cyt.bam = NULL, cond5.nuc.bam = NULL, cond5.label = "MUT4", 
                                 clip_track = NULL, plot_transcript_number = 8, project_path, limit.max.cond = NULL, limit.min.cond = NULL, ylab.size = 0.8){
  print(paste("plotting at ",level, " level..."))
  if(level == "gene") { # gene level
    gtf.gene <- gtf %>% filter(gene_name == {{ gene_name }}) %>% group_by(ensemblID, gene_name, strand, seqnames) %>% summarise(start = min(start), end = max(end), width = max(width))
    print(gtf.gene)
    if (nrow(gtf.gene) == 0){
      stop("ensembl ID not found in GTF")
    }
    gene_name <-  gtf.gene$gene_name
    thechr <-  gtf.gene$seqnames  # without chr before number
    st <- gtf.gene$start
    en <- gtf.gene$end
    strand <- gtf.gene$strand
  } else if (level == "transcript"){ # transcript level
    gtf.tx <- gtf %>% filter(transcript == transcript_id) %>% group_by(transcript, gene_name, strand, seqnames) %>% summarise(start = min(start), end = max(end), width = max(width))
    print(gtf.tx)
    if (nrow(gtf.tx) == 0){
      stop("transcript_id not found in GTF")
    }
    gene_name <-  gtf.tx$gene_name
    thechr <-  gtf.tx$seqnames  # without chr before number
    st <- gtf.tx$start
    en <- gtf.tx$end
    strand <- gtf.tx$strand
  }
  gr <- GRanges(thechr, IRanges(st - 100, en + 100), strand = strand) # Granges for retained intron with padding either side
  print(paste("Building Gene model..."))
  gene <- geneModelFromTxdb(txdb, org.Hs.eg.db, gr = gr, strand = as.character(strand))   # Build Gene model
  
  print(paste("importing ",cond1.label," BAM & calculating coverage..."))
  cond1 <- importBam(cond1.nuc.bam, cond1.cyt.bam, ranges=gr, pairs = TRUE) # import BAMs
  cond1$dat <- coverageGR(cond1$dat)    # calculate coverage
  n.conditions=1
  if(!is.null(cond2.cyt.bam) & !is.null(cond2.nuc.bam)) { # plot cond2
    print(paste("importing ",cond2.label," BAM & calculating coverage..."))
    cond2 <- importBam(cond2.nuc.bam, cond2.cyt.bam, ranges=gr, pairs = TRUE)
    cond2$dat <- coverageGR(cond2$dat)
    n.conditions=2
  }
  if(!is.null(cond3.cyt.bam) & !is.null(cond3.nuc.bam)) { # plot cond3
    print(paste("importing ",cond3.label," BAM & calculating coverage..."))
    cond3 <- importBam(cond3.nuc.bam, cond3.cyt.bam, ranges=gr, pairs = TRUE)
    cond3$dat <- coverageGR(cond3$dat)
    n.conditions=3
  }
  if(!is.null(cond4.cyt.bam) & !is.null(cond4.nuc.bam)) { # plot cond4
    print(paste("importing ",cond4.label," BAM & calculating coverage..."))
    cond4 <- importBam(cond4.nuc.bam, cond4.cyt.bam, ranges=gr, pairs = TRUE)
    cond4$dat <- coverageGR(cond4$dat)
    n.conditions=4
  }
  if(!is.null(cond5.cyt.bam) & !is.null(cond5.nuc.bam)) { # plot cond5
    print(paste("importing ",cond5.label," BAM & calculating coverage..."))
    cond5 <- importBam(cond5.nuc.bam, cond5.cyt.bam, ranges=gr, pairs = TRUE)
    cond5$dat <- coverageGR(cond5$dat)
    n.conditions=5
  }
  if(!is.null(clip_track)) { # plot CLIP peaks
    if(file.exists(paste0(here(camp_path,"/genomes/rbp-databases/clip-bedgraphs/HepG2/"), clip_track, ".merged.bed.gz"))){ print(paste0(clip_track, " clip bed found in HepG2"))
      CLIP <- import(paste0(here(camp_path,"/genomes/rbp-databases/clip-bedgraphs/HepG2/"), clip_track, ".merged.bed.gz"), "BED")
    }else if(file.exists(paste0(here(camp_path,"/genomes/rbp-databases/clip-bedgraphs/K562/"), clip_track, ".merged.bed.gz"))){ print(paste0(clip_track, " clip bed found in K562"))
      CLIP <- import(paste0(here(camp_path,"/genomes/rbp-databases/clip-bedgraphs/K562/"), clip_track, ".merged.bed.gz"), "BED")
    } else if(file.exists(paste0(here(camp_path,"/genomes/rbp-databases/clip-bedgraphs/Hela_hg19/xlinks/"), clip_track, ".bed"))){ print(paste0(clip_track, " clip bed found in Hela_hg19. liftOver to hg38."))
      CLIP <- import(paste0(here(camp_path,"/genomes/rbp-databases/clip-bedgraphs/Hela_hg19/xlinks/"), clip_track, ".bed"), "BED") # hg19
      # hg38toHg19 for Hela
      require(rtracklayer) # https://genviz.org/module-01-intro/0001/06/02/liftoverTools/
      ch = import.chain(here(camp_path,"genomes/rbp-databases/clip-bedgraphs/hg19ToHg38.over.chain")) # import chain for liftOver http://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz
      seqlevelsStyle(CLIP) <- "UCSC"
      CLIP = liftOver(CLIP, ch)
      CLIP = unlist(CLIP)
    }
    seqlevelsStyle(CLIP) <- "NCBI"
    CLIP <- new("track", dat=CLIP, type="data", format="BED")
    n.conditions=n.conditions+1
  }
  
  print(paste("Setting track arguments..."))
  if(!is.null(clip_track) & !is.null(cond5.cyt.bam) & !is.null(cond5.nuc.bam)) { optSty <- optimizeStyle(trackList(cond1, cond2, cond3, cond4, cond5, CLIP, gene)) 
  }else if(!is.null(cond5.cyt.bam) & !is.null(cond5.nuc.bam)) { optSty <- optimizeStyle(trackList(cond1, cond2, cond3, cond4, cond5, gene)) 
  }else if(!is.null(cond4.cyt.bam) & !is.null(cond4.nuc.bam)) { optSty <- optimizeStyle(trackList(cond1, cond2, cond3, cond4, gene)) 
  }else if(!is.null(cond3.cyt.bam) & !is.null(cond3.nuc.bam)) { optSty <- optimizeStyle(trackList(cond1, cond2, cond3, gene)) 
  }else if(!is.null(cond2.cyt.bam) & !is.null(cond2.nuc.bam)) { optSty <- optimizeStyle(trackList(cond1, cond2, gene)) 
  }else if(!is.null(cond1.cyt.bam) & !is.null(cond1.nuc.bam)) { optSty <- optimizeStyle(trackList(cond1, gene)) }
  trackList <- optSty$tracks
  viewerStyle <- optSty$style  # viewerStyle <- trackViewerStyle()
  setTrackViewerStyleParam(viewerStyle, "xaxis", FALSE)
  setTrackViewerStyleParam(viewerStyle, "margin", c(.01, .10, -0.1, .01))
  
  for(i in 1:n.conditions) {  # optimise coverage plots
    setTrackXscaleParam(trackList[[i]], "gp", list(cex=0.6))
    setTrackStyleParam(trackList[[i]], "height", .80/n.conditions) # height of each condition
    setTrackStyleParam(trackList[[i]], "color", c("firebrick2", "dodgerblue2"))
    setTrackStyleParam(trackList[[i]], "ylabgp", list(cex=ylab.size))
    if (!is.null(limit.max.cond)){  setTrackStyleParam(trackList[[i]], "ylim", c(limit.min.cond, limit.max.cond)) }  ## Adjust the limit of y-axis
  }
  
  if(level == "gene") { # gene level: plot multiple transcript tracks
    if (length(trackList) > (plot_transcript_number + n.conditions) ){
      print(paste("There are ", length(trackList) - n.conditions, " ensembl transcripts. Plotting only the first ",plot_transcript_number,". Removing these:")) 
      for(i in (plot_transcript_number + n.conditions):length(trackList)) {
        print(trackList[[i]]$name) 
      }
      idx = c((plot_transcript_number + n.conditions):length(trackList))
      trackList <- trackList[-idx] # only plot first n transcript tracks - remove remaining transcript tracks
    }
    for(i in (n.conditions+1):length(trackList)) {  # modify the transcript tracks
      setTrackStyleParam(trackList[[i]], "height", .17/plot_transcript_number) # height of transcript tracks
      setTrackStyleParam(trackList[[i]], "ylabgp", list(cex=0)) # hide transcript names
      setTrackStyleParam(trackList[[i]], "color", "black")
    }
  } else if (level == "transcript"){ # transcript level: plot single transcript track
    if(!is.null(clip_track) & !is.null(cond5.cyt.bam) & !is.null(cond5.nuc.bam)) { trackList <- trackList[c("cond1", "cond2", "cond3", "cond4", "cond5", "CLIP", transcript_id)]  # keep only the relevant transcript_id track
    }else if(!is.null(cond5.cyt.bam) & !is.null(cond5.nuc.bam)) { trackList <- trackList[c("cond1", "cond2", "cond3", "cond4", "cond5", transcript_id)]  # keep only the relevant transcript_id track
    }else if(!is.null(cond4.cyt.bam) & !is.null(cond4.nuc.bam)) { trackList <- trackList[c("cond1", "cond2", "cond3", "cond4", transcript_id)] 
    }else if(!is.null(cond3.cyt.bam) & !is.null(cond3.nuc.bam)) { trackList <- trackList[c("cond1", "cond2", "cond3", transcript_id)] 
    }else if(!is.null(cond2.cyt.bam) & !is.null(cond2.nuc.bam)) { trackList <- trackList[c("cond1", "cond2", transcript_id)] 
    }else if(!is.null(cond1.cyt.bam) & !is.null(cond1.nuc.bam)) { trackList <- trackList[c("cond1", transcript_id)] }
    setTrackStyleParam(trackList[[(n.conditions+1)]], "height", .13) # height of single transcript track
    setTrackStyleParam(trackList[[(n.conditions+1)]], "ylabgp", list(cex=0)) # hide transcript name
    setTrackStyleParam(trackList[[(n.conditions+1)]], "color", "black")
  }
  if(!is.null(clip_track)){# modify the CLIP track
    setTrackStyleParam(trackList[[n.conditions]], "height", .30/n.conditions) # height of transcript tracks
    setTrackStyleParam(trackList[[n.conditions]], "color", "forestgreen")
    setTrackStyleParam(trackList[[n.conditions]], "ylim", c(0, 1))## set limit of y-axis
    setTrackStyleParam(trackList[[n.conditions]], "ylabgp", list(cex=0.6)) # CLIP name size
  }
  
  if(!is.null(cond1.cyt.bam) & !is.null(cond1.nuc.bam)) { names(trackList)[1] = cond1.label } # rename conditions
  if(!is.null(cond2.cyt.bam) & !is.null(cond2.nuc.bam)) { names(trackList)[2] = cond2.label }
  if(!is.null(cond3.cyt.bam) & !is.null(cond3.nuc.bam)) { names(trackList)[3] = cond3.label }
  if(!is.null(cond4.cyt.bam) & !is.null(cond4.nuc.bam)) { names(trackList)[4] = cond4.label }
  if(!is.null(cond5.cyt.bam) & !is.null(cond5.nuc.bam)) { names(trackList)[5] = cond5.label }
  if(!is.null(clip_track)) { names(trackList)[n.conditions] = paste0(clip_track,"\nCLIP") }
  
  if(level == "gene") { # gene level file names
    print(paste0("Saving figure to: ", project_path, "/expression/coverage/", gene_name,".", n.conditions, "_conditions.nuc_cyt_coverage_plot.png"))
    png(file = paste0(project_path, "/expression/coverage/", gene_name,".", n.conditions, "_conditions.nuc_cyt_coverage_plot.png"), height = 3.5, width = 6.5, units ="in", res = 300, type="cairo")
  } else if (level == "transcript"){ # transcript level file names
    print(paste0("Saving figure to: ", project_path, "/expression/coverage/", gene_name,".",transcript_id,".", n.conditions, "_conditions.nuc_cyt_coverage_plot.png"))
    png(file = paste0(project_path, "/expression/coverage/", gene_name,".",transcript_id,".", n.conditions, "_conditions.nuc_cyt_coverage_plot.png"), height = 3.3, width = 6.5, units ="in", res = 300, type="cairo")
  }
  vp <- viewTracks(trackList, gr=gr, viewerStyle=viewerStyle)
  grid.text("Nuclear reads", x=.85, y=.85, just="bottom", gp=gpar(fontsize=6, col = "firebrick2"))
  grid.text("Cytoplasm reads", x=0.85, y=.01, just="bottom", gp=gpar(fontsize=6, col = "dodgerblue2"))
  grid.text(gene_name, x=.06, y=.92, just="bottom", gp=gpar(fontsize=10))
  dev.off()
}



# homo_sapiens.gtf <- import(here(camp_path,"/genomes/ensembl/Homo_sapiens.GRCh38.99.gtf")) %>% as_tibble
# saveRDS(homo_sapiens.gtf, here(camp_path, "genomes/ensembl/Homo_sapiens.GRCh38.99.rds"))
# homo_sapiens.gtf = readRDS(here(camp_path, "genomes/ensembl/Homo_sapiens.GRCh38.99.rds"))
# gtf.3utr <- homo_sapiens.gtf %>% filter(type == "three_prime_utr", transcript_biotype == "protein_coding")  %>% group_by(gene_id) %>% top_n(end, n = 2) %>% select(transcript_id, gene_id, gene_name, seqnames, start, end, width, strand, ccds_id) %>% unique()
# gtf.exon <- homo_sapiens.gtf %>% filter(type == "exon", transcript_biotype == "protein_coding") %>% group_by(gene_id) %>% top_n(end, n = 2) %>% select(transcript_id, gene_id, gene_name, seqnames, start, end, width, strand, ccds_id) %>% unique()
# gtf.apa <- homo_sapiens.gtf %>% filter(type %in% c("three_prime_utr", "exon"), transcript_biotype == "protein_coding") %>% group_by(gene_id) %>% top_n(end, n = 2) %>% select(transcript_id, gene_id, gene_name, seqnames, start, end, width, strand, ccds_id) %>% unique()
# saveRDS(gtf.apa, here(camp_path,"/genomes/ensembl/Homo_sapiens.GRCh38.99.gtf.apa.coords.rds"))
# gtf.apa <- readRDS(here(camp_path,"/genomes/ensembl/Homo_sapiens.GRCh38.99.gtf.apa.coords.rds"))

# track_viewer_nuc_cyt_apa <- function(gene.name = "ZSWIM7", mut.nuc.bam = mn_nuc_vcp.bam, ctrl.nuc.bam = mn_nuc_ctrl.bam, mut.cyt.bam = mn_cyt_vcp.bam, ctrl.cyt.bam = mn_cyt_ctrl.bam,
#                                      plot_transcript_number = 6, print2screen = TRUE, save = TRUE, limit.max.ctrl = "NA", limit.min.ctrl = "NA", limit.max.vcp = "NA", limit.min.vcp = "NA"){
#   gtf.gene <- gtf.apa %>% filter(gene_name == gene.name) %>% group_by(gene_id, gene_name, strand, seqnames) %>% summarise(start = min(start), end = max(end), width = max(width))
#   print(gtf.gene)
#   # cat(blue(gtf.gene))
#   if (nrow(gtf.gene) == 0){
#     stop("ensembl ID not found in GTF")
#   }
#   gene_name <-  gtf.gene$gene_name
#   thechr <-  gtf.gene$seqnames  # without chr before number
#   st <- gtf.gene$start
#   en <- gtf.gene$end
#   strand <- gtf.gene$strand
#   utr_length = gtf.gene$width
#   # xpad <- utr_length[1]*2.0
#   gr <- GRanges(thechr, IRanges(st - utr_length[1]*0.1, en + utr_length[1]*0.1), strand = strand)
#   # gr <- GRanges(thechr, IRanges(st - 50, en + 400), strand = strand) # Granges for retained intron with padding either side
#   
#   CTRL <- importBam(ctrl.nuc.bam, ctrl.cyt.bam, ranges=gr, pairs = TRUE) # import BAMs
#   VCP <- importBam(mut.nuc.bam, mut.cyt.bam, ranges=gr, pairs = TRUE)
#   CTRL$dat <- coverageGR(CTRL$dat)    # calculate coverage
#   VCP$dat <- coverageGR(VCP$dat)
#   
#   gene <- geneModelFromTxdb(txdb, org.Hs.eg.db, gr = gr, strand = as.character(strand))   # Build Gene model
#   # gene[(elementMetadata(gene)[, "feature"] == "utr3")]
#   # entrezIDforFMR1 <- BiocGenerics::get("METTL22", org.Hs.egSYMBOL2EG)
#   # three_prime_utr <- geneTrack(ids = entrezIDforFMR1, txdb = txdb)
#   # 
#   # geneTrack(ids = entrezIDforFMR1, txdb = txdb)
#   # , symbols = "METTL22", type = "transcript"
#   # columns(org.Hs.eg.db)
#   # keytypes(org.Hs.eg.db)
#   # txby <- transcriptsBy(txdb, by="three_prime_utr")
#   # txby
#   # threeUTRsByTranscript(txdb)
#   
#   # Set track arguments
#   optSty <- optimizeStyle(trackList(VCP, CTRL, gene))
#   trackList <- optSty$tracks
#   viewerStyle <- optSty$style  # viewerStyle <- trackViewerStyle()
#   setTrackViewerStyleParam(viewerStyle, "xaxis", FALSE)
#   setTrackViewerStyleParam(viewerStyle, "margin", c(.01, .10, -0.1, .01))
#   
#   # optimise coverage plots
#   for(i in 1:2) {
#     setTrackXscaleParam(trackList[[i]], "gp", list(cex=0.9))
#     setTrackStyleParam(trackList[[i]], "height", .36)
#     setTrackStyleParam(trackList[[i]], "color", c("firebrick2", "dodgerblue2"))
#     setTrackStyleParam(trackList[[i]], "ylabgp", list(cex=0.8))
#   }
#   if (length(trackList) > (plot_transcript_number + 2) ){
#     print(paste("There are ", length(trackList), " ensembl transcripts. Plotting only the first ", plot_transcript_number)) 
#     for(i in (plot_transcript_number + 3):length(trackList)) {
#       print(trackList[[i]]$name) 
#     }
#     idx = c((plot_transcript_number + 3):length(trackList))
#     trackList <- trackList[-idx] # only plot first selected transcript tracks - remove remaining transcript tracks
#   }
#   # modify the transcript tracks
#   for(i in 3:length(trackList)) {
#     setTrackStyleParam(trackList[[i]], "height", .28/(plot_transcript_number+2))
#     setTrackStyleParam(trackList[[i]], "ylabgp", list(cex=0))
#     setTrackStyleParam(trackList[[i]], "color", "black")
#   }
#   
#   if (limit.max.ctrl != "NA"){ 
#     setTrackStyleParam(trackList[[2]], "ylim", c(limit.min.ctrl, limit.max.ctrl))  ## Adjust the limit of CTRL y-axis
#   }
#   if (limit.max.vcp != "NA"){
#     setTrackStyleParam(trackList[[1]], "ylim", c(limit.min.vcp, limit.max.vcp))  ## Adjust the limit of VCP y-axis
#   }
#   
#   if (print2screen == TRUE ){
#     print("Printing to screen.")
#     vp <- viewTracks(trackList, gr=gr, viewerStyle=viewerStyle)
#     grid.text("Nuclear reads", x=.85, y=.76, just="bottom", gp=gpar(fontsize=6, col = "firebrick2"))
#     grid.text("Cytoplasm reads", x=.85, y=.41, just="bottom", gp=gpar(fontsize=6, col = "dodgerblue2"))
#     grid.text("Nuclear reads", x=.85, y=.36, just="bottom", gp=gpar(fontsize=6, col = "firebrick2"))
#     grid.text("Cytoplasm reads", x=0.85, y=.01, just="bottom", gp=gpar(fontsize=6, col = "dodgerblue2"))
#     grid.text(gene_name, x=.06, y=.88, just="bottom", gp=gpar(fontsize=10))#, rot=90)
#     dev.off()
#   }
#   if (save == TRUE ){
#     print(paste0("Saving figure as /camp/lab/luscomben/home/users/ziffo/projects/motor-neuron-mislocalisation/apa/coverage/", gene_name, ".nuc_cyt_coverage_plot.png"))
#     png(file = paste0("/camp/lab/luscomben/home/users/ziffo/projects/motor-neuron-mislocalisation/apa/coverage/", gene_name, ".nuc_cyt_coverage_plot.png"), height = 3.5, width = 6.5, units ="in", res = 300, type="cairo")
#     vp <- viewTracks(trackList, gr=gr, viewerStyle=viewerStyle)
#     grid.text("Nuclear reads", x=.85, y=.76, just="bottom", gp=gpar(fontsize=6, col = "firebrick2"))
#     grid.text("Cytoplasm reads", x=.85, y=.41, just="bottom", gp=gpar(fontsize=6, col = "dodgerblue2"))
#     grid.text("Nuclear reads", x=.85, y=.36, just="bottom", gp=gpar(fontsize=6, col = "firebrick2"))
#     grid.text("Cytoplasm reads", x=0.85, y=.01, just="bottom", gp=gpar(fontsize=6, col = "dodgerblue2"))
#     grid.text(gene_name, x=.06, y=.88, just="bottom", gp=gpar(fontsize=10))#, rot=90)
#     dev.off()
#   }}




track_viewer_nuc_cyt_ir <- function(Intron.GeneName.GeneID.Coords = "SAP25/ENSG00000205307/known-exon/7:100573379-100573596:-", IRFinder = mn.nuc_vs_cyt.vcp_vs_ctrl$irfinder, mut.nuc.bam = mn_nuc_vcp.bam, ctrl.nuc.bam = mn_nuc_ctrl.bam, mut.cyt.bam = mn_cyt_vcp.bam, ctrl.cyt.bam = mn_cyt_ctrl.bam, print2screen = TRUE, save = TRUE){
  if (!Intron.GeneName.GeneID.Coords %in% IRFinder$Intron.GeneName.GeneID.Coords){
    stop("Intron.GeneName.GeneID.Coords not found in IRFinder output")
  }
  
  intron_type = str_split_fixed(Intron.GeneName.GeneID.Coords, "/", 4)[,3] 
  gene_name <-  str_split_fixed(Intron.GeneName.GeneID.Coords, "/", 4)[,1]
  gene_id <-  str_split_fixed(Intron.GeneName.GeneID.Coords, "/", 4)[,2]
  coords <- str_split_fixed(Intron.GeneName.GeneID.Coords, "/",4)[,4]
  thechr <-  str_split_fixed(coords, ":", 3)[,1]  # without chr before number
  st <- as.numeric(str_split_fixed(coords, ":|-", 4)[,2])
  en <- as.numeric(str_split_fixed(coords, ":|-", 4)[,3])
  strand <- str_split_fixed(coords, ":|-", 4)[,4]
  intron_length = GenomicRanges::width(as.GenomicRange(Intron.GeneName.GeneID.Coords))
  xpad <- intron_length[1]*3.8
  gr <- GRanges(thechr, IRanges(st - xpad, en + xpad), strand = strand) # Granges for entire gene (not just the retained intron)
  # gr <- GRanges(thechr, IRanges(st - 100, en + 100), strand = strand) # Granges for retained intron with padding either side
  
  CTRL <- importBam(ctrl.nuc.bam, ctrl.cyt.bam, ranges=gr, pairs = TRUE) # import BAMs
  VCP <- importBam(mut.nuc.bam, mut.cyt.bam, ranges=gr, pairs = TRUE)
  CTRL$dat <- coverageGR(CTRL$dat)    # calculate coverage
  VCP$dat <- coverageGR(VCP$dat)
  
  gene <- geneModelFromTxdb(txdb, org.Hs.eg.db, gr = gr, strand = as.character(strand))   # Build Gene model
  
  # Set track arguments
  optSty <- optimizeStyle(trackList(VCP, CTRL, gene))
  trackList <- optSty$tracks
  viewerStyle <- optSty$style  # viewerStyle <- trackViewerStyle()
  setTrackViewerStyleParam(viewerStyle, "xaxis", FALSE)
  setTrackViewerStyleParam(viewerStyle, "margin", c(.01, .10, -0.1, .01))
  
  # optimise coverage plots
  for(i in 1:2) {
    setTrackXscaleParam(trackList[[i]], "gp", list(cex=0.9))
    setTrackStyleParam(trackList[[i]], "height", .36)
    setTrackStyleParam(trackList[[i]], "color", c("firebrick2", "dodgerblue2"))
    setTrackStyleParam(trackList[[i]], "ylabgp", list(cex=0.8))
  }
  
  if (length(trackList) > 8 ){
    print(paste("There are ", length(trackList) - 2, " ensembl transcripts. Plotting only the first 6.")) 
    for(i in 9:length(trackList)) {
      print(trackList[[i]]$name) 
    }
    idx = c(9:length(trackList))
    trackList <- trackList[-idx] # only plot first 6 transcript tracks - remove remaining transcript tracks
  }
  # modify the transcript tracks
  for(i in 3:length(trackList)) {
    setTrackStyleParam(trackList[[i]], "height", .03)
    setTrackStyleParam(trackList[[i]], "ylabgp", list(cex=0))
    setTrackStyleParam(trackList[[i]], "color", "black")
  }
  if (print2screen == TRUE ){
    print("Printing to screen.")
    vp <- viewTracks(trackList, gr=gr, viewerStyle=viewerStyle)
    addGuideLine(c(st, en), vp=vp, col = "grey", lwd = 3, lty = "dashed")
    addGuideLine(c(st, en), vp=vp, col = "grey", lwd = 3, lty = "dashed")
    grid.text("Nuclear reads", x=.85, y=.76, just="bottom", gp=gpar(fontsize=6, col = "firebrick2"))
    grid.text("Cytoplasm reads", x=.85, y=.41, just="bottom", gp=gpar(fontsize=6, col = "dodgerblue2"))
    grid.text("Nuclear reads", x=.85, y=.36, just="bottom", gp=gpar(fontsize=6, col = "firebrick2"))
    grid.text("Cytoplasm reads", x=0.85, y=.01, just="bottom", gp=gpar(fontsize=6, col = "dodgerblue2"))
    grid.text(gene_name, x=.06, y=.88, just="bottom", gp=gpar(fontsize=10))#, rot=90)
  }
  if (save == TRUE ){
    print(paste0("Saving figure as /camp/lab/luscomben/home/users/ziffo/projects/motor-neuron-mislocalisation/splicing/IRFinder/coverage/", gene_name, ".", thechr, ".", st, ".", en, ".nuc_cyt_IR_coverage_plot.png"))
    png(file = paste0("/camp/lab/luscomben/home/users/ziffo/projects/motor-neuron-mislocalisation/splicing/IRFinder/coverage/", gene_name, ".", thechr, ".", st, ".", en, ".nuc_cyt_IR_coverage_plot.png"), height = 3.5, width = 6.5, units ="in", res = 300)
    vp <- viewTracks(trackList, gr=gr, viewerStyle=viewerStyle)
    addGuideLine(c(st, en), vp=vp, col = "grey", lwd = 3, lty = "dashed")
    addGuideLine(c(st, en), vp=vp, col = "grey", lwd = 3, lty = "dashed")
    grid.text("Nuclear reads", x=.85, y=.76, just="bottom", gp=gpar(fontsize=6, col = "firebrick2"))
    grid.text("Cytoplasm reads", x=.85, y=.41, just="bottom", gp=gpar(fontsize=6, col = "dodgerblue2"))
    grid.text("Nuclear reads", x=.85, y=.36, just="bottom", gp=gpar(fontsize=6, col = "firebrick2"))
    grid.text("Cytoplasm reads", x=0.85, y=.01, just="bottom", gp=gpar(fontsize=6, col = "dodgerblue2"))
    grid.text(gene_name, x=.06, y=.88, just="bottom", gp=gpar(fontsize=10))#, rot=90)
    dev.off()
  }
}


track_viewer_nuc_cyt_clip <- function(gene.name = "FMR1", mut.nuc.bam = mn_nuc_vcp.bam, ctrl.nuc.bam = mn_nuc_ctrl.bam, mut.cyt.bam = mn_cyt_vcp.bam, ctrl.cyt.bam = mn_cyt_ctrl.bam,
                                      plot_transcript_number = 5, print2screen = TRUE, save = TRUE, limit.max.ctrl = "NA", limit.min.ctrl = "NA", limit.max.vcp = "NA", limit.min.vcp = "NA"){
  gtf.gene <- gtf %>% filter(gene_name == gene.name) %>% group_by(ensemblID, gene_name, strand, seqnames) %>% summarise(start = min(start), end = max(end), width = max(width))
  print(gtf.gene)
  if (nrow(gtf.gene) == 0){
    stop("ensembl ID not found in GTF")
  }
  gene_name <-  gtf.gene$gene_name
  thechr <-  gtf.gene$seqnames  # without chr before number
  st <- gtf.gene$start
  en <- gtf.gene$end
  strand <- gtf.gene$strand
  gr <- GRanges(thechr, IRanges(st - 100, en + 100), strand = strand) # Granges for retained intron with padding either side
  
  CTRL <- importBam(ctrl.nuc.bam, ctrl.cyt.bam, ranges=gr, pairs = TRUE) # import BAMs
  VCP <- importBam(mut.nuc.bam, mut.cyt.bam, ranges=gr, pairs = TRUE)
  
  CTRL$dat <- coverageGR(CTRL$dat)    # calculate coverage
  VCP$dat <- coverageGR(VCP$dat)
  
  gene <- geneModelFromTxdb(txdb, org.Hs.eg.db, gr = gr, strand = as.character(strand))   # Build Gene model
  
  # Import iCLIP peaks
  if(file.exists(paste0(here(camp_path,"/genomes/rbp-databases/clip-bedgraphs/HepG2/"), gene.name, ".merged.bed.gz"))){
    print(paste0(gene.name, " clip bed found in HepG2"))
    CLIP <- import(paste0(here(camp_path,"/genomes/rbp-databases/clip-bedgraphs/HepG2/"), gene.name, ".merged.bed.gz"), "BED")
  }else if(file.exists(paste0(here(camp_path,"/genomes/rbp-databases/clip-bedgraphs/K562/"), gene.name, ".merged.bed.gz"))){
    print(paste0(gene.name, " clip bed found in K562"))
    CLIP <- import(paste0(here(camp_path,"/genomes/rbp-databases/clip-bedgraphs/K562/"), gene.name, ".merged.bed.gz"), "BED")
  } else if(file.exists(paste0(here(camp_path,"/genomes/rbp-databases/clip-bedgraphs/Hela_hg19/"), gene.name, ".merged.bed.gz"))){
    print(paste0(gene.name, " clip bed found in Hela_hg19"))
    CLIP <- import(paste0(here(camp_path,"/genomes/rbp-databases/clip-bedgraphs/Hela_hg19/"), gene.name, ".merged.bed.gz"), "BED")
  }
  seqlevelsStyle(CLIP) <- "NCBI"
  CLIP <- new("track", dat=CLIP, type="data", format="BED")
  
  # Set track arguments
  optSty <- optimizeStyle(trackList(CLIP, VCP, CTRL, gene))
  trackList <- optSty$tracks
  viewerStyle <- optSty$style  # viewerStyle <- trackViewerStyle()
  setTrackViewerStyleParam(viewerStyle, "xaxis", FALSE)
  setTrackViewerStyleParam(viewerStyle, "margin", c(.01, .10, -0.1, .01))
  
  # optimise CLIP track
  setTrackXscaleParam(trackList[[1]], "gp", list(cex=0.9))
  setTrackStyleParam(trackList[[1]], "height", .1)
  setTrackStyleParam(trackList[[1]], "color", "forestgreen")
  setTrackStyleParam(trackList[[1]], "ylabgp", list(cex=0.8))
  setTrackStyleParam(trackList[[1]], "ylim", c(0, 1))  ## Adjust the limit of y-axis
  
  # optimise coverage tracks
  for(i in 2:3) {
    setTrackXscaleParam(trackList[[i]], "gp", list(cex=0.9))
    setTrackStyleParam(trackList[[i]], "height", .35)
    setTrackStyleParam(trackList[[i]], "color", c("firebrick2", "dodgerblue2"))
    setTrackStyleParam(trackList[[i]], "ylabgp", list(cex=0.8))
  }
  
  # modify the transcript tracks
  if (length(trackList) > (plot_transcript_number + 3) ){
    print(paste("There are ", length(trackList), " ensembl transcripts. Plotting only the first ", plot_transcript_number)) 
    for(i in (plot_transcript_number + 4):length(trackList)) {
      print(trackList[[i]]$name) 
    }
    idx = c((plot_transcript_number + 4):length(trackList))
    trackList <- trackList[-idx] # only plot first selected transcript tracks - remove remaining transcript tracks
  }
  for(i in 4:length(trackList)) {
    setTrackStyleParam(trackList[[i]], "height", .2/(plot_transcript_number+2))
    setTrackStyleParam(trackList[[i]], "ylabgp", list(cex=0))
    setTrackStyleParam(trackList[[i]], "color", "black")
  }
  
  if (limit.max.ctrl != "NA"){ 
    setTrackStyleParam(trackList[[2]], "ylim", c(limit.min.ctrl, limit.max.ctrl))  ## Adjust the limit of CTRL y-axis
  }
  if (limit.max.vcp != "NA"){
    setTrackStyleParam(trackList[[1]], "ylim", c(limit.min.vcp, limit.max.vcp))  ## Adjust the limit of VCP y-axis
  }
  
  if (print2screen == TRUE ){
    print("Printing to screen.")
    vp <- viewTracks(trackList, gr=gr, viewerStyle=viewerStyle)
    grid.text("Nuclear reads", x=.85, y=.80, just="bottom", gp=gpar(fontsize=6, col = "firebrick2"))
    grid.text("Cytoplasm reads", x=.85, y=.51, just="bottom", gp=gpar(fontsize=6, col = "dodgerblue2"))
    grid.text("Nuclear reads", x=.85, y=.45, just="bottom", gp=gpar(fontsize=6, col = "firebrick2"))
    grid.text("Cytoplasm reads", x=0.85, y=.16, just="bottom", gp=gpar(fontsize=6, col = "dodgerblue2"))
    grid.text(gene_name, x=.06, y=.90, just="bottom", gp=gpar(fontsize=10))#, rot=90)
    dev.off()
  }
  if (save == TRUE ){
    print(paste0("Saving figure as /camp/lab/luscomben/home/users/ziffo/projects/motor-neuron-mislocalisation/expression/coverage/", gene_name, ".nuc_cyt_CLIP_coverage_plot.png"))
    png(file = paste0("/camp/lab/luscomben/home/users/ziffo/projects/motor-neuron-mislocalisation/expression/coverage/", gene_name, ".nuc_cyt_CLIP_coverage_plot.png"), height = 3.5, width = 6.5, units ="in", res = 300)
    vp <- viewTracks(trackList, gr=gr, viewerStyle=viewerStyle)
    grid.text("Nuclear reads", x=.85, y=.80, just="bottom", gp=gpar(fontsize=6, col = "firebrick2"))
    grid.text("Cytoplasm reads", x=.85, y=.51, just="bottom", gp=gpar(fontsize=6, col = "dodgerblue2"))
    grid.text("Nuclear reads", x=.85, y=.45, just="bottom", gp=gpar(fontsize=6, col = "firebrick2"))
    grid.text("Cytoplasm reads", x=0.85, y=.16, just="bottom", gp=gpar(fontsize=6, col = "dodgerblue2"))
    grid.text(gene_name, x=.06, y=.90, just="bottom", gp=gpar(fontsize=10))#, rot=90)
    dev.off()
  }}


track_viewer <- function(gene_test = "SFPQ", txt = ac_who_vcp_vs_ctrl, mut.bam = ac_who_vcp, ctrl.bam = ac_who_ctrl, xpad = 1000){
  top_event <- txt %>% dplyr::filter(p0.05 == "significant") %>% dplyr::filter(gene_name == gene_test) %>% top_n(., 1, abs(IRratio.diff))
  IRratio_ctrl <- top_event %>% dplyr::pull(B.IRratio) # B = ctrl
  IRratio_ctrl <- format(round(IRratio_ctrl,2),nsmall = 2)
  IRratio_vcp <- top_event %>% dplyr::pull(A.IRratio) # A = vcp
  IRratio_vcp <- format(round(IRratio_vcp,2),nsmall = 2)
  # intron coords
  thechr <-  top_event %>% dplyr::pull(Chr) # without chr before number
  st <-  top_event %>% dplyr::pull(Start) # 35184593
  en <-  top_event %>% dplyr::pull(End) # 35187000
  strand <- top_event %>% dplyr::pull(Direction) 
  # whole gene coords
  # gene_test_gtf <- gtf %>% filter(gene_name %in% gene_test)
  # gene_chr <- gene_test_gtf %>% dplyr::pull(seqnames) %>% unique
  # gene_st <- gene_test_gtf %>% dplyr::pull(start) %>% min
  # gene_en <- gene_test_gtf %>% dplyr::pull(end) %>% max
  # gr <- GRanges(gene_chr, IRanges(gene_st, gene_en), strand= strand) # set width of tracks as entire gene
  gr <- GRanges(thechr, IRanges(st - xpad, en + xpad), strand= strand)
  CTRL <- importBam(ctrl.bam, ranges=gr, pairs = TRUE) # import BAM
  VCP <- importBam(mut.bam, ranges=gr, pairs = TRUE) # import BAM
  CTRL$dat <- coverageGR(CTRL$dat) # calculate coverage
  VCP$dat <- coverageGR(VCP$dat) # calculate coverage
  # Build Gene model
  # trs <- geneModelFromTxdb(txdb, org.Hs.eg.db, gr = gr) # build gene model showing transcripts
  # ids <- getGeneIDsFromTxDb(gr, txdb) # get ENSG of gene_test
  ids <- top_event %>% dplyr::pull(gene_id) # get ENSG of gene_test
  gene_track <- geneTrack(ids,txdb)[[1]] # gene all exons on single line
  
  # View the tracks
  optSty <- optimizeStyle(trackList(VCP, CTRL, gene_track))
  trackList <- optSty$tracks
  names(trackList)[3] <- gene_test # rename gene_track with gene_test name - shows as the gene label
  viewerStyle <- optSty$style
  setTrackViewerStyleParam(viewerStyle, "xaxis", FALSE)
  setTrackViewerStyleParam(viewerStyle, "margin", c(.01, .1, -0.1, .01))
  setTrackXscaleParam(trackList[[1]], "draw", TRUE)
  # setTrackXscaleParam(trackList[[1]], "gp", list(cex=0.9))
  setTrackStyleParam(trackList[[1]], "height", .4)
  setTrackStyleParam(trackList[[2]], "height", .4)
  setTrackStyleParam(trackList[[3]], "height", .1)
  setTrackStyleParam(trackList[[3]], "ylabgp", list(cex=1.5))
  setTrackStyleParam(trackList[[1]], "color", c("darkgrey", "black"))
  setTrackStyleParam(trackList[[2]], "color", c("darkgrey", "black"))
  setTrackStyleParam(trackList[[3]], "color", c("dodgerblue2", "black"))
  
  # save figure
  png(file = paste(here(camp_path,"/projects/astrocyte-ir-als/splicing/IRFinder/figures/coverage/"), gene_test, "_coverage_plot.png", sep = ""), height = 3.5, width = 10, units ="in", res = 300)
  vp <- viewTracks(trackList, gr=gr, viewerStyle=viewerStyle)
  addGuideLine(c(st, en), vp=vp, col = "firebrick2", lwd = 3, lty = "dashed")
  grid.text(IRratio_ctrl, x=.545, y=.6, just="bottom")
  grid.text(IRratio_vcp, x=.545, y=.2, just="bottom")
  dev.off()
}

# track_viewer(gene_test = "SFPQ", txt = ac_who_vcp_vs_ctrl, mut.bam = ac_who_vcp, ctrl.bam = ac_who_ctrl, xpad = 7000)

track_viewer_tx <- function(gene_test = "RECK", txt = ac_who_vcp_vs_ctrl, mut.bam = ac_who_vcp, ctrl.bam = ac_who_ctrl, path = NULL){
  txt %>% dplyr::filter(p0.05 == "significant" & gene_name == gene_test) %>% print
  top_event <- txt %>% dplyr::filter(p0.05 == "significant" & gene_name == gene_test) %>% top_n(., 1, abs(IRratio.diff))
  txt_nrow <- txt  %>% dplyr::filter(p0.05_reliable == "significant" & gene_name == gene_test) %>% nrow
  if(txt_nrow > 1){
    print(paste(gene_test, " has ", txt_nrow, " significant reliable IR events. Only the top deltaIR event will be plotted.", sep = ""))
  }
  IRratio_ctrl <- top_event %>% dplyr::pull(B.IRratio) # B = ctrl
  IRratio_ctrl <- format(round(IRratio_ctrl,2),nsmall = 2)
  IRratio_vcp <- top_event %>% dplyr::pull(A.IRratio) # A = vcp
  IRratio_vcp <- format(round(IRratio_vcp,2),nsmall = 2)
  thechr <-  top_event %>% dplyr::pull(Chr) # without chr before number
  st <-  top_event %>% dplyr::pull(Start) # 35184593
  en <-  top_event %>% dplyr::pull(End) # 35187000
  strand <- top_event %>% dplyr::pull(Direction) 
  intron_length <- top_event %>% dplyr::pull(intron_length)
  xpad <- intron_length[1]*3.8
  gr <- GRanges(thechr, IRanges(st - xpad, en + xpad), strand = strand) # Granges for entire gene (not just the retained intron)
  CTRL <- importBam(ctrl.bam, ranges=gr, pairs = TRUE) # import BAM
  VCP <- importBam(mut.bam, ranges=gr, pairs = TRUE) # import BAM
  CTRL$dat <- coverageGR(CTRL$dat) # calculate coverage
  VCP$dat <- coverageGR(VCP$dat) # calculate coverage
  # Build Gene model
  trs <- geneModelFromTxdb(txdb, org.Hs.eg.db, gr = gr) # build gene model showing transcripts
  # ids <- getGeneIDsFromTxDb(gr, txdb) # get ENSG of gene_test
  ids <- top_event %>% dplyr::pull(gene_id) # get ENSG of gene_test
  # gene_track <- geneTrack(ids,txdb)[[1]] # gene all exons on single line
  
  print("Creating the coverage and gene model tracks.")
  optSty <- optimizeStyle(trackList(VCP, CTRL, trs))
  trackList <- optSty$tracks
  viewerStyle <- optSty$style
  setTrackViewerStyleParam(viewerStyle, "xaxis", FALSE)
  setTrackViewerStyleParam(viewerStyle, "margin", c(.01, .15, -0.1, .01))
  setTrackXscaleParam(trackList[[1]], "draw", TRUE)
  # setTrackXscaleParam(trackList[[1]], "gp", list(cex=0.9))
  setTrackStyleParam(trackList[[1]], "height", .36)
  setTrackStyleParam(trackList[[2]], "height", .36)
  setTrackXscaleParam(trackList[[1]], "draw", TRUE)
  setTrackXscaleParam(trackList[[1]], "gp", list(cex=0.9))
  setTrackStyleParam(trackList[[1]], "color", c("darkgrey", "black"))
  setTrackStyleParam(trackList[[2]], "color", c("darkgrey", "black"))
  
  if (length(trackList) > 8 ){
    print(paste("There are ", length(trackList), " ensembl transcripts. Plotting only the first 6.")) 
    for(i in 9:length(trackList)) {
      print(trackList[[i]]$name) 
    }
    idx = c(9:length(trackList))
    trackList <- trackList[-idx] # only plot first 6 transcript tracks - remove remaining transcript tracks
  }
  # modify the transcript tracks
  for(i in 3:length(trackList)) {
    setTrackStyleParam(trackList[[i]], "height", .03)
    setTrackStyleParam(trackList[[i]], "ylabgp", list(cex=0))
    setTrackStyleParam(trackList[[i]], "color", c("dodgerblue2", "black"))
  }
  
  vp <- viewTracks(trackList, gr=gr, viewerStyle=viewerStyle)
  addGuideLine(c(st, en), vp=vp, col = "firebrick2", lwd = 3, lty = "dashed")
  addGuideLine(c(st, en), vp=vp, col = "firebrick2", lwd = 3, lty = "dashed")
  grid.text(IRratio_ctrl, x=.57, y=.6, just="bottom")
  grid.text(IRratio_vcp, x=.57, y=.2, just="bottom")
  grid.text(gene_test, x=.09, y=.87, just="bottom", gp=gpar(cex=1.5))
  
  print("Saving figure.")
  png(file = paste(here(camp_path,"/projects/astrocyte-ir-als/splicing/IRFinder/figures/coverage/", path, gene_test, "_tx_coverage_plot.png", sep = "")), height = 2.2, width = 7, units ="in", res = 300)
  vp <- viewTracks(trackList, gr=gr, viewerStyle=viewerStyle)
  addGuideLine(c(st, en), vp=vp, col = "firebrick2", lwd = 2, lty = "dashed")
  grid.text(IRratio_ctrl, x=.57, y=.65, just="bottom")
  grid.text(IRratio_vcp, x=.57, y=.2, just="bottom")
  grid.text(gene_test, x=.09, y=.87, just="bottom", gp=gpar(cex=1.5))
  dev.off()
}

# plot top delta IR ratio event irrespective of p.diff from IRFinder
track_viewer_tx_top_IRratio_event <- function(gene_test = "RECK", txt = ac_who_vcp_vs_ctrl, mut.bam = ac_who_vcp, ctrl.bam = ac_who_ctrl, path = NULL){
  txt %>% dplyr::filter(gene_name == gene_test) %>% print
  top_event <- txt %>% dplyr::filter(gene_name == gene_test) %>% top_n(., 1, abs(IRratio.diff))
  txt_nrow <- txt  %>% dplyr::filter(p0.05_reliable == "significant" & gene_name == gene_test) %>% nrow
  cat(blue(gene_test, " has ", txt_nrow, " significant reliable IR events", sep = ""))
  txt_all_nrow <- txt  %>% dplyr::filter(gene_name == gene_test) %>% nrow
  if (txt_all_nrow == 0){
    stop("Gene has no events in IRFinder differential output")
  }
  IRratio_ctrl <- top_event %>% dplyr::pull(B.IRratio) # B = ctrl
  IRratio_ctrl <- format(round(IRratio_ctrl,2),nsmall = 2)
  IRratio_vcp <- top_event %>% dplyr::pull(A.IRratio) # A = vcp
  IRratio_vcp <- format(round(IRratio_vcp,2),nsmall = 2)
  thechr <-  top_event %>% dplyr::pull(Chr) # without chr before number
  st <-  top_event %>% dplyr::pull(Start) # 35184593
  en <-  top_event %>% dplyr::pull(End) # 35187000
  strand <- top_event %>% dplyr::pull(Direction) 
  intron_length <- top_event %>% dplyr::pull(intron_length)
  xpad <- intron_length[1]*3.8
  gr <- GRanges(thechr, IRanges(st - xpad, en + xpad), strand = strand) # Granges for entire gene (not just the retained intron)
  CTRL <- importBam(ctrl.bam, ranges=gr, pairs = TRUE) # import BAM
  VCP <- importBam(mut.bam, ranges=gr, pairs = TRUE) # import BAM
  CTRL$dat <- coverageGR(CTRL$dat) # calculate coverage
  VCP$dat <- coverageGR(VCP$dat) # calculate coverage
  # Build Gene model
  trs <- geneModelFromTxdb(txdb, org.Hs.eg.db, gr = gr) # build gene model showing transcripts
  # ids <- getGeneIDsFromTxDb(gr, txdb) # get ENSG of gene_test
  ids <- top_event %>% dplyr::pull(gene_id) # get ENSG of gene_test
  # gene_track <- geneTrack(ids,txdb)[[1]] # gene all exons on single line
  
  print("Creating the coverage and gene model tracks.")
  optSty <- optimizeStyle(trackList(VCP, CTRL, trs))
  trackList <- optSty$tracks
  viewerStyle <- optSty$style
  setTrackViewerStyleParam(viewerStyle, "xaxis", FALSE)
  setTrackViewerStyleParam(viewerStyle, "margin", c(.01, .15, -0.1, .01))
  setTrackXscaleParam(trackList[[1]], "draw", TRUE)
  # setTrackXscaleParam(trackList[[1]], "gp", list(cex=0.9))
  setTrackStyleParam(trackList[[1]], "height", .36)
  setTrackStyleParam(trackList[[2]], "height", .36)
  setTrackXscaleParam(trackList[[1]], "draw", TRUE)
  setTrackXscaleParam(trackList[[1]], "gp", list(cex=0.9))
  setTrackStyleParam(trackList[[1]], "color", c("darkgrey", "black"))
  setTrackStyleParam(trackList[[2]], "color", c("darkgrey", "black"))
  
  if (length(trackList) > 8 ){
    print(paste("There are ", length(trackList), " ensembl transcripts. Plotting only the first 6.")) 
    for(i in 9:length(trackList)) {
      print(trackList[[i]]$name) 
    }
    idx = c(9:length(trackList))
    trackList <- trackList[-idx] # only plot first 6 transcript tracks - remove remaining transcript tracks
  }
  # modify the transcript tracks
  for(i in 3:length(trackList)) {
    setTrackStyleParam(trackList[[i]], "height", .03)
    setTrackStyleParam(trackList[[i]], "ylabgp", list(cex=0))
    setTrackStyleParam(trackList[[i]], "color", c("dodgerblue2", "black"))
  }
  
  vp <- viewTracks(trackList, gr=gr, viewerStyle=viewerStyle)
  addGuideLine(c(st, en), vp=vp, col = "firebrick2", lwd = 3, lty = "dashed")
  addGuideLine(c(st, en), vp=vp, col = "firebrick2", lwd = 3, lty = "dashed")
  grid.text(IRratio_ctrl, x=.57, y=.6, just="bottom")
  grid.text(IRratio_vcp, x=.57, y=.2, just="bottom")
  grid.text(gene_test, x=.09, y=.87, just="bottom", gp=gpar(cex=1.5))
  
  print("Saving figure.")
  png(file = paste(here(camp_path,"/projects/astrocyte-ir-als/splicing/IRFinder/figures/coverage/"), path, gene_test, "_topDeltaIRratio_tx_coverage_plot.png", sep = ""), height = 2.2, width = 7, units ="in", res = 300)
  vp <- viewTracks(trackList, gr=gr, viewerStyle=viewerStyle)
  addGuideLine(c(st, en), vp=vp, col = "firebrick2", lwd = 2, lty = "dashed")
  grid.text(IRratio_ctrl, x=.57, y=.65, just="bottom")
  grid.text(IRratio_vcp, x=.57, y=.2, just="bottom")
  grid.text(gene_test, x=.09, y=.87, just="bottom", gp=gpar(cex=1.5))
  dev.off()
}


# Plot multiple IR events on a single gene coverage plot
track_viewer_tx_multiple <- function(gene_test = "SFPQ", txt = ac_who_vcp_vs_ctrl, mut.bam = ac_who_vcp, ctrl.bam = ac_who_ctrl, path = NULL){
  txt_nrow <- txt  %>% dplyr::filter(p0.05_reliable == "significant" & gene_name == gene_test) %>% nrow
  if(txt_nrow > 1){
    print(paste(gene_test, " has ", txt_nrow, " significant reliable IR events.", sep = ""))
  }
  events <- txt %>% dplyr::filter(p0.05_reliable == "significant" & gene_name == gene_test) 
  
  # set GRanges for entire gene
  # whole gene coords
  gene_test_gtf <- gtf %>% filter(gene_name %in% gene_test)
  gene_chr <- gene_test_gtf %>% dplyr::pull(seqnames) %>% unique %>% as.numeric
  gene_st <- gene_test_gtf %>% dplyr::pull(start) %>% min %>% as.numeric
  gene_en <- gene_test_gtf %>% dplyr::pull(end) %>% max %>% as.numeric
  strand <- gene_test_gtf %>% dplyr::pull(strand) %>% unique %>% as.character
  gr <- GRanges(gene_chr, IRanges(gene_st, gene_en), strand= strand) # set width of tracks as entire gene
  CTRL <- importBam(ctrl.bam, ranges=gr, pairs = TRUE) # import BAM
  VCP <- importBam(mut.bam, ranges=gr, pairs = TRUE) # import BAM
  CTRL$dat <- coverageGR(CTRL$dat) # calculate coverage
  VCP$dat <- coverageGR(VCP$dat) # calculate coverage
  trs <- geneModelFromTxdb(txdb, org.Hs.eg.db, gr = gr) # build gene model showing transcripts
  
  print("Creating the coverage and gene model tracks.")
  optSty <- optimizeStyle(trackList(VCP, CTRL, trs))
  trackList <- optSty$tracks
  viewerStyle <- optSty$style
  # setTrackViewerStyleParam(viewerStyle, "xaxis", FALSE)
  # setTrackViewerStyleParam(viewerStyle, "margin", c(.01, .1, 0.03, .01))
  # setTrackXscaleParam(trackList[[1]], "draw", TRUE)
  # setTrackXscaleParam(trackList[[1]], "gp", list(cex=0.9))
  setTrackViewerStyleParam(viewerStyle, "xaxis", FALSE)
  setTrackViewerStyleParam(viewerStyle, "margin", c(.01, .15, -0.1, .01))
  setTrackStyleParam(trackList[[1]], "height", .36)
  setTrackStyleParam(trackList[[2]], "height", .36)
  setTrackXscaleParam(trackList[[1]], "draw", TRUE)
  setTrackXscaleParam(trackList[[1]], "gp", list(cex=0.9))
  setTrackStyleParam(trackList[[1]], "color", c("darkgrey", "black"))
  setTrackStyleParam(trackList[[2]], "color", c("darkgrey", "black"))
  
  # only show first 6 transcript tracks - remove remaining transcript tracks
  if (length(trackList) > 8 ){
    print("removing transcript tracks > 6:")
    for(i in 9:length(trackList)) {
      print(trackList[[i]]$name)
    }
    idx = c(9:length(trackList))
    trackList <- trackList[-idx]
    # [[9:length(trackList)]] <- NULL
  }
  # modify the transcript tracks
  for(i in 3:length(trackList)) {
    setTrackStyleParam(trackList[[i]], "height", .03)
    setTrackStyleParam(trackList[[i]], "ylabgp", list(cex=0))
    setTrackStyleParam(trackList[[i]], "color", c("dodgerblue2", "black"))
  }
  
  print("Extracting each IR event coords.")
  event <- list()
  IRratio_ctrl <- list()
  IRratio_vcp <- list()
  thechr <- list()
  st <- list()
  en <- list()
  strand <- list()
  for (i in 1:nrow(events)) {
    event[[i]] <- events[i, ] # iterate through rows of events df  
    IRratio_ctrl[[i]] <- event[[i]] %>% dplyr::pull(B.IRratio) # B = ctrl
    IRratio_ctrl[[i]] <- format(round(IRratio_ctrl[[i]],2),nsmall = 2)
    IRratio_vcp[[i]] <- event[[i]] %>% dplyr::pull(A.IRratio) # A = vcp
    IRratio_vcp[[i]] <- format(round(IRratio_vcp[[i]],2),nsmall = 2)
    thechr[[i]] <-  event[[i]] %>% dplyr::pull(Chr) # without chr before number
    st[[i]] <-  event[[i]] %>% dplyr::pull(Start) # 35184593
    en[[i]] <-  event[[i]] %>% dplyr::pull(End) # 35187000
    strand[[i]] <- event[[i]] %>% dplyr::pull(Direction) 
  }
  
  vp <- viewTracks(trackList, gr=gr, viewerStyle=viewerStyle)
  for (i in 1:nrow(events)) {
    addGuideLine(c(st[[i]], en[[i]]), vp=vp, col = "firebrick2", lwd = 1, lty = "dashed")
    # grid.text(IRratio_ctrl[[i]], x=.545, y=.6, just="bottom")
    # grid.text(IRratio_vcp[[i]], x=.545, y=.2, just="bottom")
  }
  grid.text(gene_test, x=.09, y=.87, just="bottom", gp=gpar(cex=1.5))
  
  print("Saving figure.")
  png(file = paste(here(camp_path,"/projects/astrocyte-ir-als/splicing/IRFinder/figures/coverage/"), path, gene_test, "_tx_multiple_coverage_plot.png", sep = ""), height = 2.5, width = 8, units ="in", res = 300)
  vp <- viewTracks(trackList, gr=gr, viewerStyle=viewerStyle)
  for (i in 1:nrow(events)) {
    addGuideLine(c(st[[i]], en[[i]]), vp=vp, col = "firebrick2", lwd = 1, lty = "dashed")
    # grid.text(IRratio_ctrl[[i]], x=.545, y=.6, just="bottom")
    # grid.text(IRratio_vcp[[i]], x=.545, y=.2, just="bottom")
  }
  grid.text(gene_test, x=.09, y=.87, just="bottom", gp=gpar(cex=1.5))
  dev.off()
}


track_viewer_ir_fractions <- function(gene_test = "PRPF4", txt = irfinder.ac_vcp_vs_ctrl.filt, mut.nuc.bam = ac_nuc_vcp.bam, ctrl.nuc.bam = ac_nuc_ctrl.bam, mut.cyt.bam = ac_cyt_vcp.bam, ctrl.cyt.bam = ac_cyt_ctrl.bam, path = NULL){
  txt %>% drop_na %>% arrange(IRratio.diff.nuc) %>% dplyr::filter(gene_name == gene_test) %>% print
  top_event <- txt %>% drop_na %>% arrange(IRratio.diff.nuc) %>% dplyr::filter(gene_name == gene_test) %>% top_n(., 1, abs(IRratio.diff.nuc))
  txt_nrow <- txt %>% arrange(IRratio.diff.nuc) %>% dplyr::filter(gene_name == gene_test) %>% nrow
  if(txt_nrow > 1){
    print(paste(gene_test, " has ", txt_nrow, " significant reliable nuclear IR events. Only the top deltaIR event will be plotted.", sep = ""))
  }
  IRratio_nuc_ctrl <- top_event %>% dplyr::pull(IRratio.ctrl.nuc) # B = ctrl
  IRratio_nuc_ctrl <- format(round(IRratio_nuc_ctrl,2),nsmall = 2)
  IRratio_nuc_vcp <- top_event %>% dplyr::pull(IRratio.vcp.nuc) # A = vcp
  IRratio_nuc_vcp <- format(round(IRratio_nuc_vcp,2),nsmall = 2)
  IRratio_cyt_ctrl <- top_event %>% dplyr::pull(IRratio.ctrl.cyt) # B = ctrl
  IRratio_cyt_ctrl <- format(round(IRratio_cyt_ctrl,2),nsmall = 2)
  IRratio_cyt_vcp <- top_event %>% dplyr::pull(IRratio.vcp.cyt) # A = vcp
  IRratio_cyt_vcp <- format(round(IRratio_cyt_vcp,2),nsmall = 2)
  thechr <-  top_event %>% dplyr::pull(Chr) # without chr before number
  st <-  top_event %>% dplyr::pull(Start) # 35184593
  en <-  top_event %>% dplyr::pull(End) # 35187000
  strand <- top_event %>% dplyr::pull(Direction) 
  intron_length <- top_event %>% dplyr::pull(intron_length)
  xpad <- intron_length[1]*3.8
  gr <- GRanges(thechr, IRanges(st - xpad, en + xpad), strand= strand) # Granges for retained intron with padding either side
  nuc.CTRL <- importBam(ctrl.nuc.bam, ranges=gr, pairs = TRUE) # import BAMs
  nuc.VCP <- importBam(mut.nuc.bam, ranges=gr, pairs = TRUE)
  cyt.CTRL <- importBam(ctrl.cyt.bam, ranges=gr, pairs = TRUE)
  cyt.VCP <- importBam(mut.cyt.bam, ranges=gr, pairs = TRUE)
  nuc.CTRL$dat <- coverageGR(nuc.CTRL$dat)    # calculate coverage
  nuc.VCP$dat <- coverageGR(nuc.VCP$dat)
  cyt.CTRL$dat <- coverageGR(cyt.CTRL$dat)
  cyt.VCP$dat <- coverageGR(cyt.VCP$dat)
  # Build Gene model
  trs <- geneModelFromTxdb(txdb, org.Hs.eg.db, gr = gr) # build gene model showing transcripts
  # ids <- getGeneIDsFromTxDb(gr, txdb) # get ENSG of gene_test
  ids <- top_event %>% dplyr::pull(gene_id) # get ENSG of gene_test
  # gene_track <- geneTrack(ids,txdb)[[1]] # gene all exons on single line
  
  print("Creating the coverage and gene model tracks.")
  optSty <- optimizeStyle(trackList(cyt.VCP, cyt.CTRL, nuc.VCP, nuc.CTRL, trs))
  trackList <- optSty$tracks
  viewerStyle <- optSty$style
  setTrackViewerStyleParam(viewerStyle, "xaxis", FALSE)
  setTrackViewerStyleParam(viewerStyle, "margin", c(.01, .10, -0.1, .01))
  # optimise coverage plots
  for(i in 1:4) {
    # setTrackStyleParam(trackList[[i]], "draw", TRUE)
    setTrackXscaleParam(trackList[[i]], "gp", list(cex=0.9))
    setTrackStyleParam(trackList[[i]], "height", .18)
    setTrackStyleParam(trackList[[i]], "color", c("darkgrey", "black"))
    setTrackStyleParam(trackList[[i]], "ylabgp", list(cex=0.8))
  }
  # optimise the transcript tracks
  if (length(trackList) > 10 ){
    print(paste("There are ", length(trackList), " ensembl transcripts. Plotting only the first 6.")) 
    for(i in 11:length(trackList)) {
      print(trackList[[i]]$name) 
    }
    idx = c(11:length(trackList))
    trackList <- trackList[-idx] # only plot first 6 transcript tracks - remove remaining transcript tracks
  }
  for(i in 5:length(trackList)) {
    setTrackStyleParam(trackList[[i]], "height", .03)
    setTrackStyleParam(trackList[[i]], "ylabgp", list(cex=0))
    setTrackStyleParam(trackList[[i]], "color", c("dodgerblue2", "black"))
  }
  
  vp <- viewTracks(trackList, gr=gr, viewerStyle=viewerStyle)
  addGuideLine(c(st, en), vp=vp, col = "firebrick2", lwd = 3, lty = "dashed")
  grid.text(IRratio_cyt_vcp, x=.545, y=.15, just="bottom")
  grid.text(IRratio_cyt_ctrl, x=.545, y=.35, just="bottom")
  grid.text(IRratio_nuc_vcp, x=.545, y=.55, just="bottom")
  grid.text(IRratio_nuc_ctrl, x=.545, y=.75, just="bottom")
  grid.text(gene_test, x=.05, y=.87, just="bottom", gp=gpar(cex=1.0))
  
  print("Saving figure.")
  png(file = paste("/camp/lab/luscomben/home/users/ziffo/projects/motor-neuron-mislocalisation/splicing/IRFinder/figures/coverage/", path, gene_test, "_tx_coverage_plot.png", sep = ""), height = 4, width = 7, units ="in", res = 300)
  vp <- viewTracks(trackList, gr=gr, viewerStyle=viewerStyle)
  addGuideLine(c(st, en), vp=vp, col = "firebrick2", lwd = 3, lty = "dashed")
  grid.text(IRratio_cyt_vcp, x=.545, y=.15, just="bottom")
  grid.text(IRratio_cyt_ctrl, x=.545, y=.35, just="bottom")
  grid.text(IRratio_nuc_vcp, x=.545, y=.55, just="bottom")
  grid.text(IRratio_nuc_ctrl, x=.545, y=.75, just="bottom")
  grid.text(gene_test, x=.06, y=.87, just="bottom", gp=gpar(cex=1.0))
  dev.off()
}



# DEP ---------------------------------------------------------------------

DEP.analysis <- function(protein_groups = proteinGroups_vcp_vs_ctrl, experimental_design = experimental_design_vcp_vs_ctrl, show_proteins_per_threshold = TRUE, threshold = 4, numerator = NULL, denomenator = "ctrl", 
                         design_formula = ~0 + condition, model_matrix = NULL, contrast = NULL, # formula, or model matrix if combining factors
                         pca_groups = c("condition"), #, "treatment"),
                         ir.ge.res = NULL, plot_dep = TRUE, plot_heatmaps = FALSE, file_path = "/camp/lab/luscomben/home/users/ziffo/projects/microglia-bulk-rnaseq/mass-spectrometry", file_prefix = NULL){
  DEP.se = make_se(protein_groups, grep("lfq_", colnames(protein_groups)), experimental_design)
  if(show_proteins_per_threshold == TRUE){ # this is  slow so turn FALSE once optimal threshold identified
    proteins_detected <- list()
    for (i in seq_along(unique(experimental_design$replicate))) {
      proteins_detected[[i]] = nrow(filter_missval(DEP.se, thr = unique(experimental_design$replicate)[[i]]))
    }
    filtering_proteins = tibble(filtering_threshold = c(1:max(experimental_design$replicate)), proteins_detected = unlist(proteins_detected))
    print(filtering_proteins, n=Inf)
    ggplot(filtering_proteins, aes(x = filtering_threshold, y = proteins_detected)) + geom_point() + geom_line() + theme_bw() + theme(panel.grid = element_blank(), axis.text=element_text(size=8), axis.title=element_text(size=10))
  }
  cat(blue(paste("filtering with threshold:", threshold, "\n")))
  DEP.filt = filter_missval(DEP.se, thr = threshold)
  DEP.norm = DEP.filt %>% normalize_vsn() 
  # conflict_prefer("expand", "S4Vectors")
  # unload("conflicted")
  if(is.null(model_matrix) == TRUE){ DEP.res = DEP.norm %>% impute(., fun = "MinProb", q = 0.01) %>% test_diff_oz(., type = "control", control = denomenator, design_formula = design_formula) 
  } else if (!is.null(model_matrix)) { DEP.res = DEP.norm %>% impute(., fun = "MinProb", q = 0.01) %>% test_diff_oz(., type = "manual", test = contrast, model_matrix = model_matrix) }
  write_csv(DEP.res, paste0(file_path, "/results/", file_prefix, ".dep_res_filter.thr.", threshold, ".csv"))
  
  if(!is.null(ir.ge.res)) {DEP.res.ge.ir <- DEP.res %>% left_join(select(ir.ge.res, -direction, -significant)) %>% rename(MS.log2FoldChange = log2FoldChange, MS.pvalue = pvalue, MS.padj = padj, MS.direction = direction) %>%
    mutate(significant = case_when(GE.pvalue < 0.05 & IR.pvalue < 0.05 & MS.pvalue < 0.05 ~ "all", GE.pvalue < 0.05 & MS.pvalue < 0.05 ~ "GE & MS", IR.pvalue < 0.05 & MS.pvalue < 0.05 ~ "IR & MS", GE.pvalue < 0.05 & MS.pvalue < 0.05 ~ "GE & MS",  GE.pvalue < 0.05 ~ "GE only", IR.pvalue < 0.05 ~ "IR only", MS.pvalue < 0.05 ~ "MS only", TRUE ~ "None"), 
           direction = case_when(IR.direction == "IR up" & GE.direction == "GE down" & MS.direction == "down" ~ "IR up, GE & MS down", IR.direction == "IR down" & GE.direction == "GE up" & MS.direction == "up" ~ "IR down, GE & MS up", IR.direction == "IR up" & GE.direction == "GE up" & MS.direction == "down" ~ "IR & GE up, MS down", IR.direction == "IR down" & GE.direction == "GE down" & MS.direction == "up" ~ "IR & GE down, MS up",
                                 IR.direction == "IR up" & GE.direction == "GE down" & MS.direction == "up" ~ "IR & MS up, GE down", IR.direction == "IR down" & GE.direction == "GE up" & MS.direction == "down" ~ "IR & MS down, GE up", IR.direction == "IR up" & GE.direction == "GE up" & MS.direction == "up" ~ "IR & GE & MS up", IR.direction == "IR down" & GE.direction == "GE down" & MS.direction == "down" ~ "IR & GE & MS down")) %>%
    dplyr::select(gene_name, MS.log2FoldChange, MS.padj, MS.pvalue, GE.log2FoldChange, GE.padj, GE.pvalue, GE.baseMean, IR.log2FoldChange, IR.pvalue, IR.baseMean, direction, significant, MS.direction, IR.direction, GE.direction)
  }
  
  if(plot_dep == TRUE){
    cat(blue("Analysis Complete. Now plotting...\n"))
    plot.pca = plot_pca_oz(DEP.norm, point_size = 1, indicate = pca_groups) + scale_colour_manual(values = c("firebrick2", "dodgerblue2", "forestgreen", "gold2")) + 
      theme_oz() + theme(panel.grid = element_blank(), axis.text=element_text(size=8), axis.title=element_text(size=10), legend.position = "top", legend.spacing = unit(0.01, 'cm')) + 
      guides(color=guide_legend(title="", keywidth=0.1, keyheight=0.15, default.unit="inch", reverse = TRUE, nrow = 1)) + ggtitle("") +
      ggforce::geom_mark_ellipse(aes(color = condition), show.legend = FALSE) + geom_text_repel(aes(label=colnames(DEP.norm)), size = 2.4, max.overlaps = Inf)
    if(!is.null(numerator)) { plot.volcano = volcano_plot(DEP.res, numerator = numerator, denomenator = denomenator, arrow_labels = TRUE) } else { plot.volcano = volcano_plot(DEP.res, arrow_labels = FALSE, title = contrast) }
    plot.res = ggarrange(plot.pca, plot.volcano, nrow = 2, heights = c(1,1))
    ggsave(plot.res, filename = paste0(file_path, "/figures/", file_prefix, ".dep_results_filter.thr.", threshold, ".png"), height = 10, width = 8) # as.ggplot(plot.cor), as.ggplot(plot.heatmap), 
    plot.protein.numbers = plot_numbers(DEP.filt) + scale_colour_manual(values = c("firebrick2", "dodgerblue2", "forestgreen", "gold2")) +  scale_fill_manual(values = c("firebrick2", "dodgerblue2", "forestgreen", "gold2")) + 
      theme_bw() + theme(panel.grid = element_blank(), axis.text=element_text(size=8), axis.title=element_text(size=10), axis.text.x= element_text(angle = 90)) #  barplot of the number of identified proteins per samples after filtering
    # plot.missingvalues = plot_missval(DEP.filt) # explore proteins with missing values with heatmap - bias of missing proteins in nuclear samples.
    plot.protein.intensities = plot_detect(DEP.filt) # missing proteins have lower intense proteins and slightly higher densities
    plot.normalisation = plot_normalization(DEP.filt , DEP.norm) + scale_colour_manual(values = c("firebrick2", "dodgerblue2", "forestgreen", "gold2")) +  scale_fill_manual(values = c("firebrick2", "dodgerblue2", "forestgreen", "gold2")) + 
      theme_bw() + theme(panel.grid = element_blank(), axis.text=element_text(size=8), axis.title=element_text(size=10)) # # Visualize normalization by boxplots for all samples before and after normalization
    ggsave(ggarrange(plot.protein.numbers, plot.protein.intensities, plot.normalisation, ncol = 3), filename = paste0(file_path, "/figures/", file_prefix, ".dep_qc_filter.thr.", threshold, ".png"), height = 10, width = 14)
    if(plot_heatmaps == TRUE){
      conflict_prefer("plot_cor", "DEP")
      plot.cor = grid.grabExpr(draw(plot_cor(DEP.norm, significant = TRUE, lower = 0, upper = 1, pal = "Reds"))) # plot.cor = plot_cor(DEP.norm, significant = TRUE, lower = 0, upper = 1, pal = "Reds")
      plot.heatmap = grid.grabExpr(draw(plot_heatmap(DEP.norm, type = "centered", kmeans = TRUE, show_row_names = FALSE, indicate = c("condition", "replicate"))))
      ggsave(ggarrange(plot.cor, plot.heatmap, ncol = 2), filename = paste0(file_path, "/figures/", file_prefix, ".dep_qc_heatmaps.thr.", threshold, ".png"), height = 10, width = 16) # plot.missingvalues, 
    } 
  }
  
  cat(blue("Returning results...\n"))
  if(!is.null(ir.ge.res)){ return(list(dep=DEP.norm,res=DEP.res, res.ge.ir = DEP.res.ge.ir)) #assay=DEP.assay,
  } else if(is.null(ir.ge.res) == TRUE){ return(list(dep=DEP.norm,res=DEP.res)) }#assay=DEP.assay,
  
}


test_diff_oz <- function(se, type = c("control", "all", "manual"), control = NULL, test = NULL, design_formula = ~ 0 + condition, model_matrix = NULL) {
  # assertthat::assert_that(inherits(se, "SummarizedExperiment"), is.character(type), class(design_formula) %in% c("formula","matrix","array")  # Show error if inputs are not the required classes
  type <- match.arg(type)  # Show error if inputs do not contain required columns
  col_data <- colData(se)
  raw <- assay(se)
  if(any(!c("name", "ID") %in% colnames(rowData(se, use.names = FALSE)))) {stop("'name' and/or 'ID' columns are not present in '",deparse(substitute(se)),"'\nRun make_unique() and make_se() to obtain the required columns",call. = FALSE) }
  if(any(!c("label", "condition", "replicate") %in% colnames(col_data))) { stop("'label', 'condition' and/or 'replicate' columns are not present in '", deparse(substitute(se)),  "'\nRun make_se() or make_se_parse() to obtain the required columns",  call. = FALSE) }
  if(any(is.na(raw))) { warning("Missing values in '", deparse(substitute(se)), "'") }
  if(!is.null(control)) { assertthat::assert_that(is.character(control), length(control) == 1)   # Show error if control input is not valid
    if(!control %in% unique(col_data$condition)) { stop("run test_diff() with a valid control.\nValid controls are: '", paste0(unique(col_data$condition), collapse = "', '"), "'", call. = FALSE) }
  }
  if(is.null(model_matrix)) {
    variables <- terms.formula(formula(design_formula)) %>% attr(., "variables") %>% as.character() %>%  .[-1]  # extract variables in formula
    if(any(!variables %in% colnames(col_data))) { stop("run make_diff() with an appropriate 'design_formula'") }  # Throw error if variables are not col_data columns
    if(variables[1] != "condition") { stop("first factor of 'design_formula' should be 'condition'") }
    for(var in variables) {  # Obtain variable factors
      temp <- factor(col_data[[var]])
      assign(var, temp)
    }
    model_matrix <- model.matrix(formula(design_formula), data = environment()) # Make model_matrix matrix from design formula
    conditions <- as.character(unique(condition))
  }
  colnames(model_matrix) <- gsub("condition", "", colnames(model_matrix)) 
  
  # Generate contrasts to be tested: Either make all possible combinations ("all"), only the contrasts versus the control sample ("control") or use manual contrasts
  if(type == "all") { # All possible combinations
    cntrst <- apply(utils::combn(conditions, 2), 2, paste, collapse = " - ")
    if(!is.null(control)) { # Make sure that contrast containing the control sample have the control as denominator
      flip <- grep(paste("^", control, sep = ""), cntrst)
      if(length(flip) >= 1) { cntrst[flip] <- cntrst[flip] %>% gsub(paste(control, "- ", sep = " "), "", .) %>% paste(" - ", control, sep = "") }
    }
  }
  if(type == "control") { if(is.null(control)) stop("run test_diff(type = 'control') with a 'control' argument") # Throw error if no control argument is present
    cntrst <- paste(conditions[!conditions %in% control], control, sep = " - ")  # Make contrasts
  }
  if(type == "manual") { if(is.null(test)) { stop("run test_diff(type = 'manual') with a 'test' argument") } # Throw error if no test argument is present
    assertthat::assert_that(is.character(test))
    # if(any(!unlist(strsplit(gsub("\\(|\\)","", test), "_vs_| - ")) %in% conditions)) { stop("run test_diff() with valid contrasts", ".\nValid contrasts should contain combinations of: '", paste0(conditions, collapse = "', '"), "', for example '", paste0(conditions[1], "_vs_", conditions[2]), "'.", call. = FALSE) }
    cntrst <- gsub("_vs_", " - ", test)
  }
  message("Tested contrast: ", cntrst)  # Print tested contrasts
  print(model_matrix)
  # Test for differential expression by empirical Bayes moderation of a linear model on the predefined contrasts
  fit <- lmFit(raw, design = model_matrix)
  made_contrasts <- makeContrasts(contrasts = cntrst, levels = model_matrix)
  contrast_fit <- contrasts.fit(fit, made_contrasts)
  
  if(any(is.na(raw))) {
    for(i in cntrst) {
      covariates <- strsplit(i, " - ") %>% unlist
      single_contrast <- makeContrasts(contrasts = i, levels = model_matrix[, covariates])
      single_contrast_fit <- contrasts.fit(fit[, covariates], single_contrast)
      contrast_fit$coefficients[, i] <- single_contrast_fit$coefficients[, 1]
      contrast_fit$stdev.unscaled[, i] <- single_contrast_fit$stdev.unscaled[, 1]
    }
  }
  
  eB_fit <- eBayes(contrast_fit)
  # Retrieve the differential expression test results
  res <- topTable(eB_fit, adjust="BH", sort.by = "t", coef = cntrst,  number = Inf, confint = TRUE) %>% drop_na(t) %>% as_tibble(rownames = "gene_name") %>%
    rename(log2FoldChange = contains("logFC"), padj = contains("adj.P.Val"), pvalue = contains("P.Value"), stat = t, baseMean = contains("AveExpr"))
  res$qval <- fdrtool::fdrtool(res$stat, plot = FALSE, verbose = FALSE)$qval
  # table <- res %>% mutate(comparison = gsub(" - ", "_vs_", comparison)) %>% gather(variable, value, -c(rowname,comparison)) %>%  mutate(variable = recode(variable, logFC = "diff", P.Value = "p.val", qval = "p.adj")) %>%
  # unite(temp, comparison, variable) %>%  spread(temp, value)
  return(res)
}

# test_diff_multifactor <- function(se, combined_design_factor = factor(paste(sample_replicates.untreated.tardbp_vs_ctrl.nuc_vs_cyt$fraction, sample_replicates.untreated.tardbp_vs_ctrl.nuc_vs_cyt$mutation, sep="."), levels = c("cytoplasm.ctrl","nuclear.ctrl","cytoplasm.tardbp","nuclear.tardbp")), 
#                                   contrast_formula = "(nuclear.tardbp - cytoplasm.tardbp) - (nuclear.ctrl - cytoplasm.ctrl)") {
#   assertthat::assert_that(inherits(se, "SummarizedExperiment"))   # Show error if inputs are not the required classes
#   if(any(!c("name", "ID") %in% colnames(rowData(se, use.names = FALSE)))) { stop("'name' and/or 'ID' columns are not present in '", deparse(substitute(se)), "'\nRun make_unique() and make_se() to obtain the required columns", call. = FALSE) }
#   if(any(!c("label", "replicate") %in% colnames(colData(se)))) { stop("'label' and/or 'replicate' columns are not present in '", deparse(substitute(se)), "'\nRun make_se() or make_se_parse() to obtain the required columns", call. = FALSE) }
#   if(any(is.na(assay(se)))) { warning("Missing values in '", deparse(substitute(se)), "'") }
#   
#   message("Levels: ", paste(levels(combined_design_factor)))
#   composite_design <- model.matrix(~0+combined_design_factor, data = environment()) # Make an appropriate design matrix
#   colnames(composite_design) <- levels(combined_design_factor)
#   # Test for differential expression by empirical Bayes moderation of a linear model on the predefined contrasts
#   message("Tested contrasts: ", paste(contrast_formula))  # Print tested contrasts
#   fit <- lmFit(assay(se), design = composite_design)
#   made_contrasts <- makeContrasts(contrasts = {{ contrast_formula }}, levels = composite_design)
#   contrast_fit <- contrasts.fit(fit, made_contrasts)
#   eB_fit <- eBayes(contrast_fit)
#   res <- topTable(eB_fit, adjust="BH", sort.by = "t", coef = contrast_formula, number = Inf, confint = TRUE) %>% drop_na(t) %>% as_tibble(rownames = "gene_name") %>%
#     rename(log2FoldChange = contains("logFC"), padj = contains("adj.P.Val"), pvalue = contains("P.Value"), stat = t, baseMean = contains("AveExpr"))
#   return(res)
# }


plot_pca_oz <- function(dep, x = 1, y = 2, indicate = c("condition", "replicate"),
                        label = FALSE, n = 500, point_size = 4, label_size = 3, plot = TRUE) {
  if(is.integer(x)) x <- as.numeric(x)
  if(is.integer(y)) y <- as.numeric(y)
  if(is.integer(n)) n <- as.numeric(n)
  if(is.integer(point_size)) point_size <- as.numeric(point_size)
  if(is.integer(label_size)) label_size <- as.numeric(label_size)
  # Show error if inputs are not the required classes
  assertthat::assert_that(inherits(dep, "SummarizedExperiment"),
                          is.numeric(x),
                          length(x) == 1,
                          is.numeric(y),
                          length(y) == 1,
                          is.numeric(n),
                          length(n) == 1,
                          is.character(indicate),
                          is.logical(label),
                          is.numeric(point_size),
                          length(point_size) == 1,
                          is.numeric(label_size),
                          length(label_size) == 1,
                          is.logical(plot),
                          length(plot) == 1)
  
  # Check for valid x and y values
  if(x > ncol(dep) | y > ncol(dep)) {
    stop(paste0("'x' and/or 'y' arguments are not valid\n",
                "Run plot_pca() with 'x' and 'y' <= ",
                ncol(dep), "."),
         call. = FALSE)
  }
  
  # Check for valid 'n' value
  if(n > nrow(dep)) {
    stop(paste0("'n' argument is not valid.\n",
                "Run plot_pca() with 'n' <= ",
                nrow(dep),
                "."),
         call. = FALSE)
  }
  
  # Check for valid 'indicate'
  columns <- colnames(colData(dep))
  if(!is.null(indicate)) {
    if(length(indicate) > 3) {
      stop("Too many features in 'indicate'
        Run plot_pca() with a maximum of 3 indicate features")
    }
    if(any(!indicate %in% columns)) {
      stop(paste0("'",
                  paste0(indicate, collapse = "' and/or '"),
                  "' column(s) is/are not present in ",
                  deparse(substitute(dep)),
                  ".\nValid columns are: '",
                  paste(columns, collapse = "', '"),
                  "'."),
           call. = FALSE)
    }
  }
  
  # Get the variance per protein and take the top n variable proteins
  var <- apply(assay(dep), 1, sd)
  df <- assay(dep)[order(var, decreasing = TRUE)[seq_len(n)],]
  df[which(!is.finite(df))] <- 0
  
  # Calculate PCA
  pca <- prcomp(t(df), scale = FALSE)
  pca_df <- pca$x %>% data.frame() %>% rownames_to_column() %>% left_join(., data.frame(colData(dep)), by = c("rowname" = "ID"))
  # Calculate the percentage of variance explained
  percent <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
  # Make factors of indicate features
  for(feat in indicate) { pca_df[[feat]] <- as.factor(pca_df[[feat]]) }
  
  # Plot the PCA plot
  p <- ggplot(pca_df, aes(get(paste0("PC", x)), get(paste0("PC", y)))) +
    labs(title = paste0("PCA plot - top ", n, " variable proteins"),
         x = paste0("PC", x, ": ", percent[x], "%"),
         y = paste0("PC", y, ": ", percent[y], "%")) +
    # coord_fixed() +
    theme_DEP1()
  
  if(length(indicate) == 0) {    p <- p + geom_point(size = point_size)  }
  if(length(indicate) == 1) { p <- p + geom_point(aes(col = pca_df[[indicate[1]]]), size = point_size) + labs(col = indicate[1])  }
  if(length(indicate) == 2) {    p <- p + geom_point(aes(col = pca_df[[indicate[1]]], shape = pca_df[[indicate[2]]]), size = point_size) + labs(col = indicate[1], shape = indicate[2])  }
  if(length(indicate) == 3) {    p <- p + geom_point(aes(col = pca_df[[indicate[1]]], shape = pca_df[[indicate[2]]]), size = point_size) +  facet_wrap(~pca_df[[indicate[3]]])
  labs(col = indicate[1], shape = indicate[2])
  }
  if(label) { p <- p + geom_text(aes(label = rowname), size = label_size) }
  if(plot) {
    return(p)
  } else {
    df <- pca_df %>%bselect(rowname, paste0("PC", c(x, y)), match(indicate, colnames(pca_df)))
    colnames(df)[1] <- "sample"
    return(df)
  }
}

plot_single_oz <- function(dep, proteins, type = c("contrast", "centered"), plot = TRUE, facet_nrow = NULL, facet_ncol = NULL) {
  # Show error if inputs are not the required classes
  assertthat::assert_that(inherits(dep, "SummarizedExperiment"), is.character(proteins), is.character(type), is.logical(plot), length(plot) == 1)
  
  # Show error if inputs do not contain required columns
  type <- match.arg(type)
  row_data <- rowData(dep, use.names = FALSE)
  if(any(!c("label", "condition", "replicate") %in% colnames(colData(dep)))) {
    stop("'label', 'condition' and/or 'replicate' columns are not present in '",
         deparse(substitute(dep)),
         "'\nRun make_se() or make_se_parse() to obtain the required columns",
         call. = FALSE)
  }
  if(length(grep("_p.adj|_diff", colnames(row_data))) < 1) {
    stop("'[contrast]_diff' and '[contrast]_p.adj' columns are not present in '",
         deparse(substitute(dep)),
         "'\nRun test_diff() to obtain the required columns",
         call. = FALSE)
  }
  if(!"name" %in% colnames(row_data)) {
    stop("'name' column not present in '",
         deparse(substitute(dep)),
         "'\nRun make_se() or make_se_parse() to obtain the required columns",
         call. = FALSE)
  }
  
  # Show error if an unvalid protein name is given
  if(all(!proteins %in% row_data$name)) {
    if(length(proteins) == 1) {
      rows <- grep(substr(proteins, 1, nchar(proteins) - 1),row_data$name)
      possibilities <- row_data$name[rows]
    } else {
      rows <- lapply(proteins, function(x)
        grep(substr(x, 1, nchar(x) - 1),row_data$name))
      possibilities <- row_data$name[unlist(rows)]
    }
    
    if(length(possibilities) > 0) {
      possibilities_msg <- paste0(
        "Do you mean: '",
        paste0(possibilities, collapse = "', '"),
        "'")
    } else {
      possibilities_msg <- NULL
    }
    stop("please run `plot_single()` with a valid protein names in the 'proteins' argument\n",
         possibilities_msg,
         call. = FALSE)
  }
  if(any(!proteins %in% row_data$name)) {
    proteins <- proteins[proteins %in% row_data$name]
    warning("Only used the following protein(s): '",
            paste0(proteins, collapse = "', '"),
            "'")
  }
  
  # Single protein
  subset <- dep[proteins]
  
  # Plot either the centered log-intensity values
  # per condition ('centered') or the average fold change of conditions
  # versus the control condition ('contrast') for a single protein
  if(type == "centered") {
    # Obtain protein-centered fold change values
    means <- rowMeans(assay(subset), na.rm = TRUE)
    df_reps <- data.frame(assay(subset) - means) %>%
      rownames_to_column() %>%
      gather(ID, val, -rowname) %>%
      left_join(., data.frame(colData(subset)), by = "ID")
    df_reps$replicate <- as.factor(df_reps$replicate)
    df <- df_reps %>%
      group_by(condition, rowname) %>%
      summarize(mean = mean(val, na.rm = TRUE),
                sd = sd(val, na.rm = TRUE),
                n = n()) %>%
      mutate(error = qnorm(0.975) * sd / sqrt(n),
             CI.L = mean - error,
             CI.R = mean + error) %>%
      as.data.frame()
    df$rowname <- parse_factor(df$rowname, levels = proteins)
    
    # Plot the centered intensity values for the replicates and the mean
    p <- ggplot(df, aes(condition, mean)) +
      geom_hline(yintercept = 0) +
      ggrastr::rasterise(geom_sina(data = df_reps, aes(condition, val, col = condition), adjust = 1, size = 1), dpi = 300) + 
      ggrastr::rasterise(geom_boxplot(aes(fill = condition), colour = "black", width=0.2, alpha = 0.5, outlier.shape = NA), dpi = 300) +  
      geom_errorbar(aes(ymin = CI.L, ymax = CI.R), width = 0.3) +
      labs(x = "", y = "log2 centered intensity", col = "condition") +
      scale_colour_manual(values = c("dodgerblue2","firebrick2","forestgreen","gold2","purple")) + scale_fill_manual(values =c("dodgerblue2","firebrick2","forestgreen","gold2","purple")) +
      facet_wrap(~rowname, scales = "free", nrow = {{ facet_nrow }}, ncol = {{ facet_ncol }}) + theme_oz(xaxis_arrow = FALSE) + theme(legend.position = "none")
  }
  if(type == "contrast") {
    # Select values for a single protein
    df <- rowData(subset, use.names = FALSE) %>%
      data.frame() %>%
      select(name,
             ends_with("_diff"),
             ends_with("_CI.L"),
             ends_with("_CI.R")) %>%
      gather(var, val, -name) %>%
      mutate(contrast = gsub("_diff|_CI.L|_CI.R", "", var),
             var = gsub(".*_", "", var)) %>%
      spread(var, val)
    df$name <- parse_factor(df$name, levels = proteins)
    suffix <- get_suffix(df$contrast)
    if(length(suffix)) {df$contrast <- delete_suffix(df$contrast)}
    # Plot the average fold change of conditions versus the control condition
    p <- ggplot(df, aes(contrast, diff)) +
      geom_hline(yintercept = 0) +
      geom_col(colour = "black", fill = "grey") +
      geom_errorbar(aes(ymin = CI.L, ymax = CI.R), width = 0.3) +
      labs(x = suffix, y = expression(log[2]~"Fold change"~"(\u00B195% CI)")) +
      facet_wrap(~name) +
      theme_oz(xaxis_arrow = FALSE)
  }
  if(plot) {
    return(p)
  } else {
    if(type == "centered") {
      df <- df %>%
        select(rowname, condition, mean, CI.L, CI.R)
      colnames(df) <- c("protein", "condition",
                        "log2_intensity", "CI.L", "CI.R")
    }
    if(type == "contrast") {
      df <- df %>%
        select(name, contrast, diff, CI.L, CI.R) %>%
        mutate(contrast = paste0(contrast, suffix))
      colnames(df) <- c("protein", "contrast",
                        "log2_fold_change", "CI.L", "CI.R")
    }
    return(df)
  }
}


# library(biomaRt)
# ensembl99 <- useEnsembl(biomart = 'genes', dataset = 'hsapiens_gene_ensembl', version = 99)
# filters = listFilters(ensembl99)
# attributes = listAttributes(ensembl99)
# attributes %>% filter(grepl("uniprot",name))
# uniprot2gene <- getBM(attributes=c("uniprot_gn_symbol", "uniprot_gn_id", "uniprotswissprot", "ensembl_gene_id",  "external_gene_name"), mart= ensembl99, useCache = FALSE) %>% as_tibble() #%>%
# select(uniprot_accession = uniprot_gn_id, gene_id = ensembl_gene_id, gene_name = external_gene_name) %>% distinct()
# uniprot2gene %>% filter(external_gene_name == "TARDBP")
# saveRDS(uniprot2gene, here(camp_path, "genomes/annotation/uniprot2gene.rds"))
# uniprot2gene = readRDS(here(camp_path, "genomes/annotation/uniprot2gene.rds"))


# Conflicts --------------
library(conflicted) # for package conflicts
conflict_prefer("select", "dplyr") # let dplyr::select win package conflict. can also do: select <- dplyr::select
conflict_prefer("filter", "dplyr")
conflict_prefer("mutate", "dplyr")
conflict_prefer("count", "dplyr")
conflict_prefer("n", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("list", "base")
conflict_prefer("first", "dplyr")
conflict_prefer("last", "dplyr")
conflict_prefer("match", "base")
conflict_prefer("which.max", "base")
conflict_prefer("genes", "GenomicFeatures")
conflict_prefer("as.data.frame", "base")
conflict_prefer("desc", "dplyr")
conflict_prefer("pheatmap", "ComplexHeatmap")
conflict_prefer("rowVars", "matrixStats")
conflict_prefer("name", "GenomicScores")
conflict_prefer("unique", "base")
conflict_prefer("make_clean_names", "janitor")
conflict_prefer("problems", "readr")
conflict_prefer("melt", "reshape2")
conflict_prefer("set", "dendextend")
conflict_prefer("fisher.test", "stats")
conflict_prefer("intersect", "base")
conflict_prefer("expand", "S4Vectors") # for DEP function
conflict_prefer("results", "DESeq2")

