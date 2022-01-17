# load("/camp/home/ziffo/home/projects/ipsc-mn-als-meta/expression/deseq2/ipsc_mn_als_meta.RData")
source("/camp/lab/luscomben/home/users/ziffo/scripts/functions/OnDemand_R_functions.R")
# source("/camp/lab/luscomben/home/users/ziffo/scripts/functions/all_r_functions.R")


# # Metadata ----------------------------------------------------------------

### hiPSC datasets
catanese.metadata <- read_csv("/camp/lab/luscomben/home/shared/projects/patani-collab/public-data/ipsc-mn-c9orf72-fus-catanese-2021/sample-details/samplesheet.csv") %>% distinct(sample, .keep_all = TRUE) %>%
  mutate(file_salmon = file.path("/camp/lab/luscomben/home/shared/projects/patani-collab/public-data/ipsc-mn-c9orf72-fus-catanese-2021/nfcore/star_salmon", sample, "quant.sf"), 
         file_irfinder = file.path("/camp/lab/luscomben/home/shared/projects/patani-collab/public-data/ipsc-mn-c9orf72-fus-catanese-2021/splicing/irfinder", sample), dataset = "catanese", age = as.numeric(gsub("Day in Vitro ","",age)))
lee.mutants.metadata <- read_csv("/camp/lab/luscomben/home/shared/projects/patani-collab/public-data/ipsc-occular-spinal-mn-lee-2021/sample-details/samplesheet.csv") %>% distinct(sample, .keep_all = TRUE) %>%
  mutate(file_salmon = file.path("/camp/lab/luscomben/home/shared/projects/patani-collab/public-data/ipsc-occular-spinal-mn-lee-2021/nfcore/star_salmon", sample, "quant.sf"), 
         file_irfinder = file.path("/camp/lab/luscomben/home/shared/projects/patani-collab/public-data/ipsc-occular-spinal-mn-lee-2021/splicing/irfinder", sample), dataset = "lee", condition = "als", replicate = as.numeric(gsub("REP","",replicate))) %>% filter(location == "spinal")
lee.controls.metadata <- read_csv("/camp/lab/luscomben/home/shared/projects/patani-collab/public-data/ipsc-occular-spinal-mn-lee-2021/controls/sample-details/samplesheet.csv") %>% distinct(sample, .keep_all = TRUE) %>%
  mutate(file_salmon = file.path("/camp/lab/luscomben/home/shared/projects/patani-collab/public-data/ipsc-occular-spinal-mn-lee-2021/controls/nfcore/star_salmon", sample, "quant.sf"), 
         file_irfinder = file.path("/camp/lab/luscomben/home/shared/projects/patani-collab/public-data/ipsc-occular-spinal-mn-lee-2021/splicing/irfinder", sample), dataset = "lee", condition = "ctrl", instrument = "Illumina NovaSeq 6000", total_size = read_count) %>% 
  filter(location == "spinal")
lee.metadata = bind_rows(lee.mutants.metadata, lee.controls.metadata)
dafianca.c9orf72.metadata <- read_csv("/camp/lab/luscomben/home/shared/projects/patani-collab/public-data/ipsc-mn-c9orf72-dafinca-2020/sample-details/samplesheet.csv") %>% distinct(sample, .keep_all = TRUE) %>%
  mutate(file_salmon = file.path("/camp/lab/luscomben/home/shared/projects/patani-collab/public-data/ipsc-mn-c9orf72-dafinca-2020/nfcore/star_salmon", sample, "quant.sf"), 
         file_irfinder = file.path("/camp/lab/luscomben/home/shared/projects/patani-collab/public-data/ipsc-mn-c9orf72-dafinca-2020/splicing/irfinder", sample), dataset = "dafianca.c9orf72")
dafianca.tardbp.metadata <- read_csv("/camp/lab/luscomben/home/shared/projects/patani-collab/public-data/ipsc-mn-tardbp-dafinca-2020/sample-details/samplesheet.csv") %>% distinct(sample, .keep_all = TRUE) %>%
  mutate(file_salmon = file.path("/camp/lab/luscomben/home/shared/projects/patani-collab/public-data/ipsc-mn-tardbp-dafinca-2020/nfcore/star_salmon", sample, "quant.sf"), 
         file_irfinder = file.path("/camp/lab/luscomben/home/shared/projects/patani-collab/public-data/ipsc-mn-tardbp-dafinca-2020/splicing/irfinder", sample), dataset = "dafianca.tardbp")
luisier.metadata <- read_csv("/camp/lab/luscomben/home/shared/projects/patani-collab/motor-neuron-vcp-luisier-2018/sample-details/samplesheet.csv") %>% distinct(sample, .keep_all = TRUE) %>%
  mutate(file_salmon = file.path("/camp/lab/luscomben/home/shared/projects/patani-collab/motor-neuron-vcp-luisier-2018/nfcore/star_salmon", sample, "quant.sf"), 
         file_irfinder = file.path("/camp/lab/luscomben/home/shared/projects/patani-collab/motor-neuron-vcp-luisier-2018/splicing/irfinder", sample), dataset = "luisier", total_size = bases)
mitchell.metadata <- read_csv("/camp/lab/luscomben/home/shared/projects/patani-collab/inter-neuron-bulk-rnaseq/sample-details/samplesheet.csv") %>% distinct(sample, .keep_all = TRUE) %>%
  mutate(file_salmon = file.path("/camp/lab/luscomben/home/shared/projects/patani-collab/inter-neuron-bulk-rnaseq/nfcore/star_salmon", sample, "quant.sf"), 
         file_irfinder = file.path("/camp/lab/luscomben/home/shared/projects/patani-collab/inter-neuron-bulk-rnaseq/splicing/irfinder", sample), dataset = "mitchell",
         mutation = case_when(condition == "als" ~ "vcp", TRUE ~ "ctrl"), days_of_differentiation = 25, age = 25, instrument = "Illumina HiSeq 4000",library_layout = "PAIRED", total_size = 5000000000, avg_spot_len = 100) %>% filter(day == 25, celltype == "mn")
wang.metadata <- read_csv("/camp/lab/luscomben/home/shared/projects/patani-collab/public-data/ipsc-mn-sod1-wang-2017/sample-details/samplesheet.csv") %>% distinct(sample, .keep_all = TRUE) %>%
  mutate(file_salmon = file.path("/camp/lab/luscomben/home/shared/projects/patani-collab/public-data/ipsc-mn-sod1-wang-2017/nfcore/star_salmon", sample, "quant.sf"), 
         file_irfinder = file.path("/camp/lab/luscomben/home/shared/projects/patani-collab/public-data/ipsc-mn-sod1-wang-2017/splicing/irfinder", sample), dataset = "wang")
moccia.metadata <- read_csv("/camp/lab/luscomben/home/shared/projects/patani-collab/public-data/ipsc-mn-sod1-moccia-2014/sample-details/samplesheet.csv") %>% distinct(sample, .keep_all = TRUE) %>%
  mutate(file_salmon = file.path("/camp/lab/luscomben/home/shared/projects/patani-collab/public-data/ipsc-mn-sod1-moccia-2014/nfcore/star_salmon", sample, "quant.sf"), 
         file_irfinder = file.path("/camp/lab/luscomben/home/shared/projects/patani-collab/public-data/ipsc-mn-sod1-moccia-2014/splicing/irfinder", sample), dataset = "moccia")
sareen.metadata <- read_csv("/camp/lab/luscomben/home/shared/projects/patani-collab/public-data/ipsc-mn-c9orf72-sareen-2013/sample-details/samplesheet.csv") %>% distinct(sample, .keep_all = TRUE) %>%
  mutate(file_salmon = file.path("/camp/lab/luscomben/home/shared/projects/patani-collab/public-data/ipsc-mn-c9orf72-sareen-2013/nfcore/star_salmon", sample, "quant.sf"), 
         file_irfinder = file.path("/camp/lab/luscomben/home/shared/projects/patani-collab/public-data/ipsc-mn-c9orf72-sareen-2013/splicing/irfinder", sample), dataset = "sareen")
sterneckert.metadata <- read_csv("/camp/lab/luscomben/home/shared/projects/patani-collab/public-data/ipsc-mn-c9orf72-sterneckert-2020/sample-details/samplesheet.csv") %>% distinct(sample, .keep_all = TRUE) %>%
  mutate(file_salmon = file.path("/camp/lab/luscomben/home/shared/projects/patani-collab/public-data/ipsc-mn-c9orf72-sterneckert-2020/nfcore/star_salmon", sample, "quant.sf"), 
         file_irfinder = file.path("/camp/lab/luscomben/home/shared/projects/patani-collab/public-data/ipsc-mn-c9orf72-sterneckert-2020/splicing/irfinder", sample), dataset = "sterneckert")
kapeli.metadata <- read_csv("/camp/lab/luscomben/home/shared/projects/patani-collab/public-data/ipsc-mn-fus-kapeli-2016/sample-details/samplesheet.csv") %>% distinct(sample, .keep_all = TRUE) %>%
  mutate(file_salmon = file.path("/camp/lab/luscomben/home/shared/projects/patani-collab/public-data/ipsc-mn-fus-kapeli-2016/nfcore/star_salmon", sample, "quant.sf"), 
         file_irfinder = file.path("/camp/lab/luscomben/home/shared/projects/patani-collab/public-data/ipsc-mn-fus-kapeli-2016/splicing/irfinder", sample), dataset = "kapeli")
melamed.metadata <- read_csv("/camp/lab/luscomben/home/shared/projects/patani-collab/public-data/ipsc-mn-tardbp-melamed-2019/sample-details/samplesheet.csv") %>% distinct(sample, .keep_all = TRUE) %>%
  mutate(file_salmon = file.path("/camp/lab/luscomben/home/shared/projects/patani-collab/public-data/ipsc-mn-tardbp-melamed-2019/nfcore/star_salmon", sample, "quant.sf"), 
         file_irfinder = file.path("/camp/lab/luscomben/home/shared/projects/patani-collab/public-data/ipsc-mn-tardbp-melamed-2019/splicing/irfinder", sample), dataset = "melamed")
desantis.metadata <- read_csv("/camp/lab/luscomben/home/shared/projects/patani-collab/public-data/ipsc-mn-fus-desantis-2017/sample-details/samplesheet.csv") %>% distinct(sample, .keep_all = TRUE) %>%
  mutate(file_salmon = file.path("/camp/lab/luscomben/home/shared/projects/patani-collab/public-data/ipsc-mn-fus-desantis-2017/nfcore/star_salmon", sample, "quant.sf"), 
         file_irfinder = file.path("/camp/lab/luscomben/home/shared/projects/patani-collab/public-data/ipsc-mn-fus-desantis-2017/splicing/irfinder", sample), days_of_differentiation = 19, dataset = "desantis")
smith.metadata <- read_csv("/camp/lab/luscomben/home/shared/projects/patani-collab/public-data/ipsc-mn-tardbp-smith-2021/sample-details/samplesheet.csv") %>% distinct(sample, .keep_all = TRUE) %>%
  mutate(file_salmon = file.path("/camp/lab/luscomben/home/shared/projects/patani-collab/public-data/ipsc-mn-tardbp-smith-2021/nfcore/star_salmon", sample, "quant.sf"), 
         file_irfinder = file.path("/camp/lab/luscomben/home/shared/projects/patani-collab/public-data/ipsc-mn-tardbp-smith-2021/splicing/irfinder", sample), dataset = "smith")
singapore.metadata <- read_csv("/camp/lab/luscomben/home/shared/projects/patani-collab/public-data/ipsc-mn-sod1-singapore-2017/sample-details/samplesheet.csv") %>% distinct(sample, .keep_all = TRUE) %>%
  mutate(file_salmon = file.path("/camp/lab/luscomben/home/shared/projects/patani-collab/public-data/ipsc-mn-sod1-singapore-2017/nfcore/star_salmon", sample, "quant.sf"), 
         file_irfinder = file.path("/camp/lab/luscomben/home/shared/projects/patani-collab/public-data/ipsc-mn-sod1-singapore-2017/splicing/irfinder", sample), dataset = "singapore") %>% select(-age)

# ALS hiPSC meta-analysis all datasets
ipsc_mn_als_datasets.metadata <- bind_rows(catanese.metadata, lee.metadata, dafianca.c9orf72.metadata, dafianca.tardbp.metadata, luisier.metadata, mitchell.metadata, wang.metadata, moccia.metadata, sareen.metadata, sterneckert.metadata, kapeli.metadata, melamed.metadata, desantis.metadata, smith.metadata, singapore.metadata) %>% 
  mutate(sample_name = sample, sample = paste0(dataset,"_",sample), condition = factor(condition, levels = c("ctrl", "als")), dataset = as.factor(dataset), mutation = as.factor(mutation),
         DIV = as.factor(case_when(dataset %in% c("catanese","lee","luisier","sareen","kapeli","smith","singapore") ~ 35, dataset %in% c("dafianca.c9orf72", "dafianca.tardbp") ~ 30, dataset == "wang" ~ 12, dataset == "melamed" ~ 28, dataset == "desantis" ~ 19, dataset == "moccia" ~ 30, 
                                   dataset == "sterneckert" ~ 40, dataset == "mitchell" ~ 25)), instrument = gsub("Illumina ","",instrument)) %>%
  select(sample, dataset, condition, mutation, sample_name, replicate, file_salmon, file_irfinder, everything())
# FUS
ipsc_mn_fus_datasets.metadata <- bind_rows(filter(catanese.metadata, mutation %in% c("fus","iso","ctrl")), kapeli.metadata, desantis.metadata) %>% 
  mutate(sample_name = sample, sample = paste0(dataset,"_",sample), condition = factor(condition, levels = c("ctrl", "als")), dataset = as.factor(dataset), mutation = factor(mutation, levels = c("ctrl", "iso", "fus"))) %>%
  select(sample, dataset, condition, mutation, sample_name, replicate, file_salmon, file_irfinder, everything())
# SOD1
ipsc_mn_sod1_datasets.metadata <- bind_rows(filter(lee.metadata, mutation %in% c("sod1","iso","ctrl")), wang.metadata, moccia.metadata, singapore.metadata) %>% 
  mutate(sample_name = sample, sample = paste0(dataset,"_",sample), condition = factor(condition, levels = c("ctrl", "als")), dataset = as.factor(dataset), mutation = factor(mutation, levels = c("ctrl", "iso", "sod1"))) %>%
  select(sample, dataset, condition, mutation, sample_name, replicate, file_salmon, file_irfinder, everything())
# C9orf72
ipsc_mn_c9orf72_datasets.metadata <- bind_rows(filter(catanese.metadata, mutation %in% c("c9orf72","iso","ctrl")), filter(lee.metadata, mutation %in% c("c9orf72","iso","ctrl")), dafianca.c9orf72.metadata, sareen.metadata, sterneckert.metadata) %>% 
  mutate(sample_name = sample, sample = paste0(dataset,"_",sample), condition = factor(condition, levels = c("ctrl", "als")), dataset = as.factor(dataset), mutation = factor(mutation, levels = c("ctrl", "iso", "c9orf72"))) %>%
  select(sample, dataset, condition, mutation, sample_name, replicate, file_salmon, file_irfinder, everything())
# TARDBP
ipsc_mn_tardbp_datasets.metadata <- bind_rows(dafianca.tardbp.metadata, melamed.metadata, smith.metadata) %>% 
  mutate(sample_name = sample, sample = paste0(dataset,"_",sample), condition = factor(condition, levels = c("ctrl", "als")), dataset = as.factor(dataset), mutation = factor(mutation, levels = c("ctrl", "iso", "tardbp"))) %>%
  select(sample, dataset, condition, mutation, sample_name, replicate, file_salmon, file_irfinder, everything())
# VCP
ipsc_mn_vcp_datasets.metadata <- bind_rows(luisier.metadata, mitchell.metadata) %>% mutate(sample_name = sample, sample = paste0(dataset,"_",sample), condition = factor(condition, levels = c("ctrl", "als")), dataset = as.factor(dataset), mutation = factor(mutation, levels = c("ctrl", "iso", "vcp"))) %>%
  select(sample, dataset, condition, mutation, sample_name, replicate, file_salmon, file_irfinder, everything())


# # DESeq2 ------------------------------------------------------------------
# 
# ### hiPSC
# Individual mutations
ipsc_mn_fus_datasets = DESeq.analysis(metadata = ipsc_mn_fus_datasets.metadata, design = ~ dataset + condition, contrast = "condition_als_vs_ctrl", species = "human", transcript.level = FALSE, IRFinder = TRUE, IRFinder_design = ~ dataset + condition + condition:IRFinder)
ipsc_mn_sod1_datasets = DESeq.analysis(metadata = ipsc_mn_sod1_datasets.metadata, design = ~ dataset + condition, contrast = "condition_als_vs_ctrl", species = "human", transcript.level = FALSE, IRFinder = TRUE, IRFinder_design = ~ dataset + condition + condition:IRFinder)
ipsc_mn_c9orf72_datasets = DESeq.analysis(metadata = ipsc_mn_c9orf72_datasets.metadata, design = ~ dataset + condition, contrast = "condition_als_vs_ctrl", species = "human", transcript.level = FALSE, IRFinder = TRUE, IRFinder_design = ~ dataset + condition + condition:IRFinder)
ipsc_mn_tardbp_datasets = DESeq.analysis(metadata = ipsc_mn_tardbp_datasets.metadata, design = ~ dataset + condition, contrast = "condition_als_vs_ctrl", species = "human", transcript.level = FALSE, IRFinder = TRUE, IRFinder_design = ~ dataset + condition + condition:IRFinder)
ipsc_mn_vcp_datasets = DESeq.analysis(metadata = ipsc_mn_vcp_datasets.metadata, design = ~ dataset + condition, contrast = "condition_als_vs_ctrl", species = "human", transcript.level = FALSE, IRFinder = TRUE, IRFinder_design = ~ dataset + condition + condition:IRFinder)

# all mutations
ipsc_mn_als_datasets = DESeq.analysis(metadata = ipsc_mn_als_datasets.metadata, design = ~ dataset + condition, contrast = "condition_als_vs_ctrl", species = "human", transcript.level = FALSE, IRFinder = TRUE, IRFinder_design = ~ dataset + condition + condition:IRFinder)

print("saving")
save.image("/camp/home/ziffo/home/projects/ipsc-mn-als-meta/expression/deseq2/ipsc_mn_als_meta.RData")
print("complete")
