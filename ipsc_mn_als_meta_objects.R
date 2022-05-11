# run this script in Rstudio server OnDemand Terminal Tab (not with conda env r4.0.3) with:
# cd /camp/home/ziffo/home/projects/ipsc-mn-als-meta/scripts
# sbatch -N 1 -c 12 --mem=72G -t 72:00:00 --wrap="Rscript /camp/lab/luscomben/home/users/ziffo/projects/ipsc-mn-als-meta/scripts/ipsc_mn_als_meta/ipsc_mn_als_meta_objects.R" --mail-type=ALL,ARRAY_TASKS --mail-user=ziffo@crick.ac.uk --job-name=DESeq2objects
# sbatch -N 1 -c 6 --mem=50G -t 8:00:00 --wrap="Rscript /camp/lab/luscomben/home/users/ziffo/projects/ipsc-mn-als-meta/scripts/ipsc_mn_als_meta/ipsc_mn_als_meta_objects.R" --mail-type=ALL,ARRAY_TASKS --mail-user=ziffo@crick.ac.uk --job-name=DESeq2objects

# load("/camp/home/ziffo/home/projects/ipsc-mn-als-meta/expression/deseq2/ipsc_mn_als_meta.RData")
source("/camp/lab/luscomben/home/users/ziffo/scripts/functions/OnDemand_R_functions.R")

# Samplesheets ----------------------------------------------------------------
answerals.metadata = read_csv("/camp/project/proj-luscombn-patani/working/answerals/sample-details/samplesheet.csv") %>% distinct(sample, .keep_all = TRUE) %>%
  mutate(database_dir = "/camp/project/proj-luscombn-patani/working/answerals", dataset = "answerals", dataset_sample = paste0(dataset,"_",sample), study_accession = dataset) %>%
  filter(condition %in% c("als", "ctrl")) %>% # remove 1 C9orf72 carrier and 4 other mnd
  filter(!(condition == "ctrl" & mutation %in% c("optn", "als2"))) %>% # remove the 2 controls that have mutations in OPTN and ALS2
  mutate(condition = factor(condition, levels = c("ctrl", "als")))
answerals.als_only.metadata = answerals.metadata %>% filter(condition == "als")
answerals.progression.metadata = answerals.als_only.metadata %>% filter(progression %in% c("fast","slow")) %>% mutate(progression = factor(progression, levels = c("slow","fast")))
answerals.age_onset.metadata = answerals.als_only.metadata %>% drop_na(age_onset) %>% mutate(early_late_onset = factor(early_late_onset, levels = c("late","early")))
answerals.mortality.metadata = answerals.als_only.metadata %>% drop_na(survival) %>% mutate(early_late_death = factor(early_late_death, levels = c("late","early")))
answerals.onset_site.metadata = answerals.metadata %>% mutate(onset_site_detailed = replace_na(onset_site_detailed, "ctrl"), onset_site = replace_na(onset_site, "ctrl"), onset_site = factor(onset_site, levels = c("ctrl","other","limb","bulbar"))) %>% #onset_site = factor(onset_site, levels = c("ctrl", "limb","axial","mixed","bulbar")), 
  filter(!(condition =="als" & onset_site =="ctrl"))
answerals.als_only.onset_site.metadata = answerals.onset_site.metadata %>% filter(condition == "als") %>% mutate(onset_site = fct_drop(onset_site)) #onset_site = fct_drop(onset_site), 
neurolincs.metadata = read_csv("/camp/project/proj-luscombn-patani/working/neurolincs/sample-details/samplesheet.csv") %>% distinct(sample, .keep_all = TRUE) %>% filter(DIV != 0) %>% # diMN = 0007..., iMN A-042-... labels
  mutate(database_dir = "/camp/project/proj-luscombn-patani/working/neurolincs", condition = factor(condition, levels = c("ctrl", "als")), dataset = case_when(dataset == "neurolincs0" ~ "neurolincs.diMN",dataset == "neurolincsA" ~ "neurolincs.iMN"), study_accession = dataset)
catanese.metadata <- read_csv("/camp/project/proj-luscombn-patani/working/public-data/ipsc-mn-c9orf72-fus-catanese-2021/sample-details/samplesheet.csv") %>% distinct(sample, .keep_all = TRUE) %>%
  mutate(database_dir = "/camp/project/proj-luscombn-patani/working/public-data/ipsc-mn-c9orf72-fus-catanese-2021", dataset = "catanese", DIV = 35)%>% select(-age)
dafianca.c9orf72.metadata <- read_csv("/camp/project/proj-luscombn-patani/working/public-data/ipsc-mn-c9orf72-dafinca-2020/sample-details/samplesheet.csv") %>% distinct(sample, .keep_all = TRUE) %>%
  mutate(database_dir = "/camp/project/proj-luscombn-patani/working/public-data/ipsc-mn-c9orf72-dafinca-2020", dataset = "dafianca.c9orf72", DIV = 19)
dafianca.tardbp.metadata <- read_csv("/camp/project/proj-luscombn-patani/working/public-data/ipsc-mn-tardbp-dafinca-2020/sample-details/samplesheet.csv") %>% distinct(sample, .keep_all = TRUE) %>%
  mutate(database_dir = "/camp/project/proj-luscombn-patani/working/public-data/ipsc-mn-tardbp-dafinca-2020", dataset = "dafianca.tardbp", DIV = 19)
luisier.metadata <- read_csv("/camp/lab/luscomben/home/shared/projects/patani-collab/motor-neuron-vcp-luisier-2018/sample-details/samplesheet.csv") %>% distinct(sample, .keep_all = TRUE) %>%
  filter(day == 35) %>%
  mutate(database_dir = "/camp/lab/luscomben/home/shared/projects/patani-collab/motor-neuron-vcp-luisier-2018", dataset = "luisier", total_size = bases, name = gsub("d35_","",sample), DIV = 35)
mitchell.metadata <- read_csv("/camp/lab/luscomben/home/shared/projects/patani-collab/inter-neuron-bulk-rnaseq/sample-details/samplesheet.csv") %>% distinct(sample, .keep_all = TRUE) %>%
  filter(day == 25, celltype == "mn") %>%
  mutate(database_dir = "/camp/lab/luscomben/home/shared/projects/patani-collab/inter-neuron-bulk-rnaseq", dataset = "mitchell", DIV = 25, instrument = "Illumina HiSeq 4000",library_layout = "PAIRED", mutation = case_when(condition == "als" ~ "vcp", TRUE ~ "ctrl"), total_size = 5000000000, avg_spot_len = 100)
wang.metadata <- read_csv("/camp/project/proj-luscombn-patani/working/public-data/ipsc-mn-sod1-wang-2017/sample-details/samplesheet.csv") %>% distinct(sample, .keep_all = TRUE) %>%
  mutate(database_dir = "/camp/project/proj-luscombn-patani/working/public-data/ipsc-mn-sod1-wang-2017", dataset = "wang", DIV = 12)
kiskinis.metadata <- read_csv("/camp/project/proj-luscombn-patani/working/public-data/ipsc-mn-sod1-moccia-2014/sample-details/samplesheet.csv") %>% distinct(sample, .keep_all = TRUE) %>%
  mutate(database_dir = "/camp/project/proj-luscombn-patani/working/public-data/ipsc-mn-sod1-moccia-2014", dataset = "kiskinis", DIV = 45)
sareen.metadata <- read_csv("/camp/project/proj-luscombn-patani/working/public-data/ipsc-mn-c9orf72-sareen-2013/sample-details/samplesheet.csv") %>% distinct(sample, .keep_all = TRUE) %>%
  mutate(database_dir = "/camp/project/proj-luscombn-patani/working/public-data/ipsc-mn-c9orf72-sareen-2013", dataset = "sareen", DIV = 35)
sterneckert.metadata <- read_csv("/camp/project/proj-luscombn-patani/working/public-data/ipsc-mn-c9orf72-sterneckert-2020/sample-details/samplesheet.csv") %>% distinct(sample, .keep_all = TRUE) %>%
  mutate(database_dir = "/camp/project/proj-luscombn-patani/working/public-data/ipsc-mn-c9orf72-sterneckert-2020", dataset = "sterneckert", DIV = 40)
kapeli.metadata <- read_csv("/camp/project/proj-luscombn-patani/working/public-data/ipsc-mn-fus-kapeli-2016/sample-details/samplesheet.csv") %>% distinct(sample, .keep_all = TRUE) %>%
  filter(sh_rna != "scramble") %>%
  mutate(database_dir = "/camp/project/proj-luscombn-patani/working/public-data/ipsc-mn-fus-kapeli-2016", dataset = "kapeli", DIV = 35)
desantis.metadata <- read_csv("/camp/project/proj-luscombn-patani/working/public-data/ipsc-mn-fus-desantis-2017/sample-details/samplesheet.csv") %>% distinct(sample, .keep_all = TRUE) %>%
  mutate(database_dir = "/camp/project/proj-luscombn-patani/working/public-data/ipsc-mn-fus-desantis-2017", DIV = 19, dataset = "desantis") %>% select(-days_of_differentiation)
smith.metadata <- read_csv("/camp/project/proj-luscombn-patani/working/public-data/ipsc-mn-tardbp-smith-2021/sample-details/samplesheet.csv") %>% distinct(sample, .keep_all = TRUE) %>%
  mutate(database_dir = "/camp/project/proj-luscombn-patani/working/public-data/ipsc-mn-tardbp-smith-2021", dataset = "smith", DIV = 45)
bhinge.metadata <- read_csv("/camp/project/proj-luscombn-patani/working/public-data/ipsc-mn-sod1-bhinge-2017/sample-details/samplesheet.csv") %>% distinct(sample, .keep_all = TRUE) %>%
  mutate(database_dir = "/camp/project/proj-luscombn-patani/working/public-data/ipsc-mn-sod1-bhinge-2017", dataset = "bhinge", DIV = 32) %>% select(-age)
# lee.metadata <- read_csv("/camp/project/proj-luscombn-patani/working/public-data/ipsc-occular-spinal-mn-lee-2021/sample-details/samplesheet.csv") %>% distinct(sample, .keep_all = TRUE) %>%
#   mutate(database_dir = "/camp/project/proj-luscombn-patani/working/public-data/ipsc-occular-spinal-mn-lee-2021", dataset = "lee", condition = case_when(mutation == "ctrl" ~ "ctrl", TRUE ~ "als")) %>% filter(location == "spinal")

# Metadata merge ----------------------------------------------------------------
ipsc_mn_als_datasets.metadata <- bind_rows(answerals.metadata, neurolincs.metadata, catanese.metadata, dafianca.c9orf72.metadata, dafianca.tardbp.metadata, luisier.metadata, mitchell.metadata,
                                           wang.metadata, kiskinis.metadata, sareen.metadata, sterneckert.metadata, kapeli.metadata, desantis.metadata, smith.metadata, bhinge.metadata) %>%
  mutate(dataset_sample = paste0(dataset,"_",sample), condition = factor(condition, levels = c("ctrl", "als")), dataset = as.factor(dataset), mutation = as.factor(mutation),
         DIV = as.factor(DIV), instrument = gsub("Illumina ","",instrument),
         library_type = case_when(dataset %in% c("catanese","dafinca","dafianca","kapeli","kiskinis","luisier","mitchell","sareen","smith","sterneckert","wang") ~ "poly(A)", TRUE ~ "Total")) %>%
  select(sample, dataset, condition, mutation, dataset_sample, replicate, database_dir, everything())
# add intron_read_pct
multiqc_stats.tsv = ipsc_mn_als_datasets.metadata %>% distinct(database_dir) %>% mutate(multiqc_stat_paths = paste0(database_dir,"/nfcore/multiqc/star_salmon/multiqc_data/multiqc_rseqc_read_distribution.txt")) %>% pull(multiqc_stat_paths) %>% map(read_tsv, show_col_types = FALSE)
names(multiqc_stats.tsv) =  ipsc_mn_als_datasets.metadata %>% distinct(database_dir,.keep_all=TRUE) %>% pull(database_dir)
multiqc_rseqc_read_distribution <- map_dfr(multiqc_stats.tsv, bind_rows, .id = "database_dir") %>% clean_names()
ipsc_mn_als_datasets.metadata = ipsc_mn_als_datasets.metadata %>% left_join(select(multiqc_rseqc_read_distribution, database_dir, sample, introns_tag_pct)) # need unique sample names - only works if sample & dataset are unique. If multiple dataset names within a single multiqc_general_stats.txt then will only join the first dataset name
ipsc_mn_als_datasets.noiso.metadata <- ipsc_mn_als_datasets.metadata %>% filter(mutation != "iso") # remove 15 isogenic corrected samples
ipsc_mn_familial_datasets.metadata = ipsc_mn_als_datasets.metadata %>% filter(mutation != "sporadic")

ipsc_mn_sporadic_datasets.metadata <- bind_rows(filter(answerals.metadata, mutation %in% c("sporadic","ctrl")), filter(neurolincs.metadata, mutation %in% c("sporadic","ctrl"))) %>%
  mutate(dataset_sample = paste0(dataset,"_",sample), condition = factor(condition, levels = c("ctrl", "als")), dataset = as.factor(dataset), mutation = factor(mutation, levels = c("ctrl", "sporadic"))) %>%
  select(sample, dataset, condition, mutation, dataset_sample, database_dir, everything())

ipsc_mn_c9orf72_datasets.metadata <- bind_rows(filter(answerals.metadata, mutation %in% c("c9orf72","ctrl")), filter(neurolincs.metadata, mutation %in% c("c9orf72","ctrl")), filter(catanese.metadata, mutation %in% c("c9orf72","iso","ctrl")), dafianca.c9orf72.metadata, sareen.metadata, sterneckert.metadata) %>%
  mutate(dataset_sample = paste0(dataset,"_",sample), condition = factor(condition, levels = c("ctrl", "als")), dataset = as.factor(dataset), mutation = factor(mutation, levels = c("ctrl", "iso", "c9orf72")),
         library_type = case_when(dataset %in% c("catanese","dafinca","dafianca","kapeli","kiskinis","luisier","mitchell","sareen","smith","sterneckert","wang") ~ "poly(A)", TRUE ~ "Total")) %>%
  select(sample, dataset, condition, mutation, dataset_sample, replicate, database_dir, everything())

ipsc_mn_fus_datasets.metadata <- bind_rows(filter(answerals.metadata, mutation %in% c("fus","ctrl")), filter(catanese.metadata, mutation %in% c("fus","iso","ctrl")), kapeli.metadata, desantis.metadata) %>%
  mutate(dataset_sample = paste0(dataset,"_",sample), condition = factor(condition, levels = c("ctrl", "als")), dataset = as.factor(dataset), mutation = factor(mutation, levels = c("ctrl", "iso", "fus")),
         library_type = case_when(dataset %in% c("catanese","dafinca","dafianca","kapeli","kiskinis","luisier","mitchell","sareen","smith","sterneckert","wang") ~ "poly(A)", TRUE ~ "Total")) %>%
  select(sample, dataset, condition, mutation, dataset_sample, replicate, database_dir, everything())

ipsc_mn_sod1_datasets.metadata <- bind_rows(filter(answerals.metadata, mutation %in% c("sod1","ctrl")), filter(neurolincs.metadata, mutation %in% c("sod1","ctrl")), wang.metadata, kiskinis.metadata, bhinge.metadata) %>%
  mutate(dataset_sample = paste0(dataset,"_",sample), condition = factor(condition, levels = c("ctrl", "als")), dataset = as.factor(dataset), mutation = factor(mutation, levels = c("ctrl", "iso", "sod1")),
         library_type = case_when(dataset %in% c("catanese","dafinca","dafianca","kapeli","kiskinis","luisier","mitchell","sareen","smith","sterneckert","wang") ~ "poly(A)", TRUE ~ "Total")) %>%
  select(sample, dataset, condition, mutation, dataset_sample, replicate, database_dir, everything())

ipsc_mn_tardbp_datasets.metadata <- bind_rows(filter(answerals.metadata, mutation %in% c("tardbp","ctrl")), dafianca.tardbp.metadata, smith.metadata) %>%
  mutate(dataset_sample = paste0(dataset,"_",sample), condition = factor(condition, levels = c("ctrl", "als")), dataset = as.factor(dataset), mutation = factor(mutation, levels = c("ctrl", "iso", "tardbp")),
         library_type = case_when(dataset %in% c("catanese","dafinca","dafianca","kapeli","kiskinis","luisier","mitchell","sareen","smith","sterneckert","wang") ~ "poly(A)", TRUE ~ "Total")) %>%
  select(sample, dataset, condition, mutation, dataset_sample, replicate, database_dir, everything())

ipsc_mn_vcp_datasets.metadata <- bind_rows(luisier.metadata, mitchell.metadata) %>% mutate(dataset_sample = paste0(dataset,"_",sample), condition = factor(condition, levels = c("ctrl", "als")), dataset = as.factor(dataset), mutation = factor(mutation, levels = c("ctrl", "iso", "vcp"))) %>%
  select(sample, dataset, condition, mutation, dataset_sample, replicate, database_dir, everything())

ipsc_mn_answerals_neurolincs.metadata <- bind_rows(answerals.metadata, neurolincs.metadata) %>%
  mutate(dataset_sample = paste0(dataset,"_",sample), condition = factor(condition, levels = c("ctrl", "als")), dataset = as.factor(dataset), mutation = as.factor(mutation), DIV = as.factor(DIV), instrument = gsub("Illumina ","",instrument)) %>%
  select(sample, dataset, condition, mutation, dataset_sample, database_dir, everything())

ipsc_mn_c9orf72_sporadic.metadata <- bind_rows(ipsc_mn_c9orf72_datasets.metadata, ipsc_mn_sporadic_datasets.metadata) %>%  filter(condition != "ctrl") %>%
  mutate(dataset_sample = paste0(dataset,"_",sample), mutation = factor(mutation, levels = c("sporadic", "c9orf72")), dataset = as.factor(dataset), DIV = as.factor(DIV), instrument = gsub("Illumina ","",instrument)) %>%
  select(sample, dataset, condition, mutation, dataset_sample, database_dir, everything())
ipsc_mn_tdp43_positive_datasets.metadata = ipsc_mn_als_datasets.metadata %>% filter(!mutation %in% c("sod1","fus"), !dataset %in% c("kapeli","desantis","wang","kiskinis","bhinge"))
ipsc_mn_tdp43_negative_datasets.metadata = bind_rows(ipsc_mn_fus_datasets.metadata, ipsc_mn_sod1_datasets.metadata) %>% distinct(dataset_sample, .keep_all = TRUE)
ipsc_mn_tdp_pathology.metadata = ipsc_mn_als_datasets.metadata %>% filter(condition == "als") %>% mutate(tdp_pathology = case_when(mutation %in% c("sod1","fus") ~ "no", TRUE ~ "yes"))
ipsc_mn_mutant_sporadic.metadata = ipsc_mn_als_datasets.metadata %>% filter(condition == "als") %>% mutate(mutant_sporadic = factor(case_when(mutation == "sporadic" ~ "sporadic", TRUE ~ "mutant"), levels = c("sporadic", "mutant")))

save.image("/camp/home/ziffo/home/projects/ipsc-mn-als-meta/expression/deseq2/ipsc_mn_als_meta.metadata.RData")

# # DESeq2 ------------------------------------------------------------------

### ALS ------------------------------------------------------------------
# ipsc_mn_als_datasets = DESeq.analysis(metadata = ipsc_mn_als_datasets.metadata, unique_names = "dataset_sample", design = ~ dataset + condition, contrast = "condition_als_vs_ctrl", run_irfinder = TRUE, irfinder_design = ~ dataset + condition + condition:IRFinder)
# saveRDS(ipsc_mn_als_datasets, "/camp/home/ziffo/home/projects/ipsc-mn-als-meta/expression/deseq2/ipsc_mn_als_datasets.rds")
# 
# ipsc_mn_sporadic_datasets = DESeq.analysis(metadata = ipsc_mn_sporadic_datasets.metadata, unique_names = "dataset_sample", design = ~ dataset + condition, contrast = "condition_als_vs_ctrl", species = "human", run_irfinder = TRUE, irfinder_design = ~ dataset + condition + condition:IRFinder)
# saveRDS(ipsc_mn_sporadic_datasets, "/camp/home/ziffo/home/projects/ipsc-mn-als-meta/expression/deseq2/ipsc_mn_sporadic_datasets.rds")
# ipsc_mn_c9orf72_datasets = DESeq.analysis(metadata = ipsc_mn_c9orf72_datasets.metadata, unique_names = "dataset_sample", design = ~ dataset + condition, contrast = "condition_als_vs_ctrl", species = "human", run_irfinder = TRUE, irfinder_design = ~ dataset + condition + condition:IRFinder)
# saveRDS(ipsc_mn_c9orf72_datasets, "/camp/home/ziffo/home/projects/ipsc-mn-als-meta/expression/deseq2/ipsc_mn_c9orf72_datasets.rds")
# ipsc_mn_fus_datasets = DESeq.analysis(metadata = ipsc_mn_fus_datasets.metadata, unique_names = "dataset_sample", design = ~ dataset + condition, contrast = "condition_als_vs_ctrl", species = "human", run_irfinder = TRUE, irfinder_design = ~ dataset + condition + condition:IRFinder)
# saveRDS(ipsc_mn_fus_datasets, "/camp/home/ziffo/home/projects/ipsc-mn-als-meta/expression/deseq2/ipsc_mn_fus_datasets.rds")
# ipsc_mn_sod1_datasets = DESeq.analysis(metadata = ipsc_mn_sod1_datasets.metadata, unique_names = "dataset_sample", design = ~ dataset + condition, contrast = "condition_als_vs_ctrl", species = "human", run_irfinder = TRUE, irfinder_design = ~ dataset + condition + condition:IRFinder)
# saveRDS(ipsc_mn_sod1_datasets, "/camp/home/ziffo/home/projects/ipsc-mn-als-meta/expression/deseq2/ipsc_mn_sod1_datasets.rds")
# ipsc_mn_tardbp_datasets = DESeq.analysis(metadata = ipsc_mn_tardbp_datasets.metadata, unique_names = "dataset_sample", design = ~ dataset + condition, contrast = "condition_als_vs_ctrl", species = "human", run_irfinder = TRUE, irfinder_design = ~ dataset + condition + condition:IRFinder)
# saveRDS(ipsc_mn_tardbp_datasets, "/camp/home/ziffo/home/projects/ipsc-mn-als-meta/expression/deseq2/ipsc_mn_tardbp_datasets.rds")
# ipsc_mn_vcp_datasets = DESeq.analysis(metadata = ipsc_mn_vcp_datasets.metadata, unique_names = "dataset_sample", design = ~ dataset + condition, contrast = "condition_als_vs_ctrl", species = "human", run_irfinder = TRUE, irfinder_design = ~ dataset + condition + condition:IRFinder)
# saveRDS(ipsc_mn_vcp_datasets, "/camp/home/ziffo/home/projects/ipsc-mn-als-meta/expression/deseq2/ipsc_mn_vcp_datasets.rds")

### ALS mutation groupings ------------------------------------------------------------------
# ipsc_mn_familial_datasets = DESeq.analysis(metadata = ipsc_mn_familial_datasets.metadata, unique_names = "dataset_sample", design = ~ dataset + condition, contrast = "condition_als_vs_ctrl", species = "human", run_irfinder = TRUE, irfinder_design = ~ dataset + condition + condition:IRFinder)
# saveRDS(ipsc_mn_familial_datasets, "/camp/home/ziffo/home/projects/ipsc-mn-als-meta/expression/deseq2/ipsc_mn_familial_datasets.rds")
# ipsc_mn_c9orf72_sporadic = DESeq.analysis(metadata = ipsc_mn_c9orf72_sporadic.metadata, unique_names = "dataset_sample", design = ~ dataset + mutation, contrast = "mutation_c9orf72_vs_sporadic", species = "human")
# saveRDS(ipsc_mn_c9orf72_sporadic, "/camp/home/ziffo/home/projects/ipsc-mn-als-meta/expression/deseq2/ipsc_mn_c9orf72_sporadic.rds")
# ipsc_mn_tdp43_positive_datasets = DESeq.analysis(metadata = ipsc_mn_tdp43_positive_datasets.metadata, unique_names = "dataset_sample", design = ~ dataset + condition, contrast = "condition_als_vs_ctrl", species = "human")
# saveRDS(ipsc_mn_tdp43_positive_datasets, "/camp/home/ziffo/home/projects/ipsc-mn-als-meta/expression/deseq2/ipsc_mn_tdp43_positive_datasets.rds")
# ipsc_mn_tdp43_negative_datasets = DESeq.analysis(metadata = ipsc_mn_tdp43_negative_datasets.metadata, unique_names = "dataset_sample", design = ~ dataset + condition, contrast = "condition_als_vs_ctrl", species = "human")
# saveRDS(ipsc_mn_tdp43_negative_datasets, "/camp/home/ziffo/home/projects/ipsc-mn-als-meta/expression/deseq2/ipsc_mn_tdp43_negative_datasets.rds")
# ipsc_mn_tdp_pathology = DESeq.analysis(metadata = ipsc_mn_tdp_pathology.metadata, unique_names = "dataset_sample", design = ~ dataset + tdp_pathology, contrast = "tdp_pathology_yes_vs_no", species = "human")
# saveRDS(ipsc_mn_tdp_pathology, "/camp/home/ziffo/home/projects/ipsc-mn-als-meta/expression/deseq2/ipsc_mn_tdp_pathology.rds")
# ipsc_mn_mutant_sporadic = DESeq.analysis(metadata = ipsc_mn_mutant_sporadic.metadata, unique_names = "dataset_sample", design = ~ dataset + mutant_sporadic, contrast = "mutant_sporadic_mutant_vs_sporadic", species = "human")
# saveRDS(ipsc_mn_mutant_sporadic, "/camp/home/ziffo/home/projects/ipsc-mn-als-meta/expression/deseq2/ipsc_mn_mutant_sporadic.rds")
# ipsc_mn_als_datasets.noiso = DESeq.analysis(metadata = ipsc_mn_als_datasets.noiso.metadata, unique_names = "dataset_sample", design = ~ dataset + condition, contrast = "condition_als_vs_ctrl", run_irfinder = TRUE, irfinder_design = ~ dataset + condition + condition:IRFinder)
# saveRDS(ipsc_mn_als_datasets.noiso, "/camp/home/ziffo/home/projects/ipsc-mn-als-meta/expression/deseq2/ipsc_mn_als_datasets.noiso.rds")
# ipsc_mn_answerals.mutation = DESeq.analysis(metadata = ipsc_mn_als_datasets.metadata, unique_names = "dataset_sample", LRT_full = ~dataset + mutation, LRT_reduced = ~dataset, species = "human") # Mutation LRT
# saveRDS(ipsc_mn_answerals.mutation, "/camp/home/ziffo/home/projects/ipsc-mn-als-meta/expression/deseq2/ipsc_mn_answerals.mutation.rds")

### Library type polyA vs Total ------------------------------------------------------------------
# ipsc_mn_als_datasets.polyA = DESeq.analysis(metadata = filter(ipsc_mn_als_datasets.metadata, library_type == "poly(A)"), unique_names = "dataset_sample", design = ~ dataset + condition, contrast = "condition_als_vs_ctrl", run_irfinder = TRUE, irfinder_design = ~ dataset + condition + condition:IRFinder)
# saveRDS(ipsc_mn_als_datasets.polyA, "/camp/home/ziffo/home/projects/ipsc-mn-als-meta/expression/deseq2/ipsc_mn_als_datasets.polyA.rds")
# ipsc_mn_als_datasets.total = DESeq.analysis(metadata = filter(ipsc_mn_als_datasets.metadata, library_type == "Total"), unique_names = "dataset_sample", design = ~ dataset + condition, contrast = "condition_als_vs_ctrl", run_irfinder = TRUE, irfinder_design = ~ dataset + condition + condition:IRFinder)
# saveRDS(ipsc_mn_als_datasets.total, "/camp/home/ziffo/home/projects/ipsc-mn-als-meta/expression/deseq2/ipsc_mn_als_datasets.total.rds")
# ipsc_mn_c9orf72_datasets.total = DESeq.analysis(metadata = filter(ipsc_mn_c9orf72_datasets.metadata, library_type == "Total"), unique_names = "dataset_sample", design = ~ dataset + condition, contrast = "condition_als_vs_ctrl", run_irfinder = TRUE, irfinder_design = ~ dataset + condition + condition:IRFinder)
# saveRDS(ipsc_mn_c9orf72_datasets.total, "/camp/home/ziffo/home/projects/ipsc-mn-als-meta/expression/deseq2/ipsc_mn_c9orf72_datasets.total.rds")
ipsc_mn_fus_datasets.polyA = DESeq.analysis(metadata = filter(ipsc_mn_fus_datasets.metadata, library_type == "poly(A)"), unique_names = "dataset_sample", design = ~ dataset + condition, contrast = "condition_als_vs_ctrl", run_irfinder = TRUE, irfinder_design = ~ dataset + condition + condition:IRFinder)
saveRDS(ipsc_mn_fus_datasets.polyA, "/camp/home/ziffo/home/projects/ipsc-mn-als-meta/expression/deseq2/ipsc_mn_fus_datasets.polyA.rds")
# ipsc_mn_sod1_datasets.total = DESeq.analysis(metadata = filter(ipsc_mn_sod1_datasets.metadata, library_type == "Total"), unique_names = "dataset_sample", design = ~ dataset + condition, contrast = "condition_als_vs_ctrl", run_irfinder = TRUE, irfinder_design = ~ dataset + condition + condition:IRFinder)
# saveRDS(ipsc_mn_sod1_datasets.total, "/camp/home/ziffo/home/projects/ipsc-mn-als-meta/expression/deseq2/ipsc_mn_sod1_datasets.total.rds")
ipsc_mn_tardbp_datasets.total = DESeq.analysis(metadata = filter(ipsc_mn_tardbp_datasets.metadata, library_type == "Total"), unique_names = "dataset_sample", design = ~ dataset + condition, contrast = "condition_als_vs_ctrl", run_irfinder = TRUE, irfinder_design = ~ dataset + condition + condition:IRFinder)
saveRDS(ipsc_mn_tardbp_datasets.total, "/camp/home/ziffo/home/projects/ipsc-mn-als-meta/expression/deseq2/ipsc_mn_tardbp_datasets.total.rds")

### AnswerALS Clinical Correlation ------------------------------------------------------------------
# ipsc_mn_answerals = DESeq.analysis(metadata = answerals.metadata, unique_names = "sample", design = ~ condition, contrast = "condition_als_vs_ctrl", species = "human", run_irfinder = TRUE, irfinder_design = ~ condition + condition:IRFinder)
# saveRDS(ipsc_mn_answerals, "/camp/home/ziffo/home/projects/ipsc-mn-als-meta/expression/deseq2/ipsc_mn_answerals.rds")
# ipsc_mn_answerals.alsfrs_slope = DESeq.analysis(metadata = answerals.progression.metadata, unique_names = "sample", LRT_full = ~ns(alsfrs_slope, 3), LRT_reduced = ~1, species = "human")
# saveRDS(ipsc_mn_answerals.alsfrs_slope, "/camp/home/ziffo/home/projects/ipsc-mn-als-meta/expression/deseq2/ipsc_mn_answerals.alsfrs_slope.rds")
# ipsc_mn_answerals.alsfrs_r_progression_slope = DESeq.analysis(metadata = drop_na(answerals.als_only.metadata,alsfrs_r_progression_slope), unique_names = "sample",  LRT_full = ~ns(alsfrs_r_progression_slope, 3), LRT_reduced = ~1, species = "human")
# saveRDS(ipsc_mn_answerals.alsfrs_r_progression_slope, "/camp/home/ziffo/home/projects/ipsc-mn-als-meta/expression/deseq2/ipsc_mn_answerals.alsfrs_r_progression_slope.rds")
# ipsc_mn_answerals.age_onset = DESeq.analysis(metadata = answerals.age_onset.metadata, unique_names = "sample", LRT_full = ~ns(age_onset, 3), LRT_reduced = ~1, species = "human")
# saveRDS(ipsc_mn_answerals.age_onset, "/camp/home/ziffo/home/projects/ipsc-mn-als-meta/expression/deseq2/ipsc_mn_answerals.age_onset.rds")
# ipsc_mn_answerals.survival = DESeq.analysis(metadata = answerals.mortality.metadata, unique_names = "sample", LRT_full = ~ns(survival, 3), LRT_reduced = ~1, species = "human")
# saveRDS(ipsc_mn_answerals.survival, "/camp/home/ziffo/home/projects/ipsc-mn-als-meta/expression/deseq2/ipsc_mn_answerals.survival.rds")
# ipsc_mn_answerals.ck = DESeq.analysis(metadata = drop_na(answerals.als_only.metadata, ck), unique_names = "sample", LRT_full = ~ns(ck, 3), LRT_reduced = ~1, species = "human")
# saveRDS(ipsc_mn_answerals.ck, "/camp/home/ziffo/home/projects/ipsc-mn-als-meta/expression/deseq2/ipsc_mn_answerals.ck.rds")
# ipsc_mn_answerals.el_escorial = DESeq.analysis(metadata = drop_na(answerals.als_only.metadata, elescrlr), unique_names = "sample", LRT_full = ~ns(elescrlr, 3), LRT_reduced = ~1, species = "human")
# saveRDS(ipsc_mn_answerals.el_escorial, "/camp/home/ziffo/home/projects/ipsc-mn-als-meta/expression/deseq2/ipsc_mn_answerals.el_escorial.rds")

# ipsc_mn_answerals.onset_site = DESeq.analysis(metadata = answerals.onset_site.metadata, unique_names = "sample", LRT_full = ~onset_site, LRT_reduced = ~1, species = "human")
# saveRDS(ipsc_mn_answerals.onset_site, "/camp/home/ziffo/home/projects/ipsc-mn-als-meta/expression/deseq2/ipsc_mn_answerals.onset_site.rds")
# ipsc_mn_answerals.als_only.onset_site = DESeq.analysis(metadata = answerals.als_only.onset_site.metadata, unique_names = "sample", LRT_full = ~onset_site, LRT_reduced = ~1, species = "human")
# saveRDS(ipsc_mn_answerals.als_only.onset_site, "/camp/home/ziffo/home/projects/ipsc-mn-als-meta/expression/deseq2/ipsc_mn_answerals.als_only.onset_site.rds")
# ipsc_mn_answerals.bulbar_ctrl = DESeq.analysis(metadata = filter(answerals.onset_site.metadata, onset_site %in% c("ctrl","bulbar")), unique_names = "sample", design = ~ onset_site, contrast = "onset_site_bulbar_vs_ctrl", species = "human")
# saveRDS(ipsc_mn_answerals.bulbar_ctrl, "/camp/home/ziffo/home/projects/ipsc-mn-als-meta/expression/deseq2/ipsc_mn_answerals.bulbar_ctrl.rds")
# ipsc_mn_answerals.limb_ctrl = DESeq.analysis(metadata = filter(answerals.onset_site.metadata, onset_site %in% c("ctrl","limb")), unique_names = "sample", design = ~ onset_site, contrast = "onset_site_limb_vs_ctrl", species = "human")
# saveRDS(ipsc_mn_answerals.limb_ctrl, "/camp/home/ziffo/home/projects/ipsc-mn-als-meta/expression/deseq2/ipsc_mn_answerals.limb_ctrl.rds")
# ipsc_mn_answerals.other_ctrl = DESeq.analysis(metadata = filter(answerals.onset_site.metadata, onset_site %in% c("ctrl","other")), unique_names = "sample", design = ~ onset_site, contrast = "onset_site_other_vs_ctrl", species = "human")
# saveRDS(ipsc_mn_answerals.other_ctrl, "/camp/home/ziffo/home/projects/ipsc-mn-als-meta/expression/deseq2/ipsc_mn_answerals.other_ctrl.rds")
# ipsc_mn_answerals.mixed_onset = DESeq.analysis(metadata = filter(answerals.onset_site.metadata, onset_site %in% c("ctrl","mixed")), unique_names = "sample", design = ~ onset_site, contrast = "onset_site_mixed_vs_ctrl", species = "human")
# saveRDS(ipsc_mn_answerals.mixed_onset, "/camp/home/ziffo/home/projects/ipsc-mn-als-meta/expression/deseq2/ipsc_mn_answerals.mixed_onset.rds")
# ipsc_mn_answerals.limb_bulbar = DESeq.analysis(metadata = filter(answerals.onset_site.metadata, onset_site %in% c("bulbar","limb")), unique_names = "sample", design = ~ onset_site, contrast = "onset_site_bulbar_vs_limb", species = "human")
# saveRDS(ipsc_mn_answerals.limb_bulbar, "/camp/home/ziffo/home/projects/ipsc-mn-als-meta/expression/deseq2/ipsc_mn_answerals.limb_bulbar.rds")

# ipsc_mn_answerals.mutant = DESeq.analysis(metadata = answerals.als_only.metadata, unique_names = "sample", design = ~ mutant, contrast = "mutant_yes_vs_no", species = "human")
# saveRDS(ipsc_mn_answerals.mutant, "/camp/home/ziffo/home/projects/ipsc-mn-als-meta/expression/deseq2/ipsc_mn_answerals.mutant.rds")
# ipsc_mn_answerals.family_history = DESeq.analysis(metadata = drop_na(answerals.als_only.metadata, family_history), unique_names = "sample", design = ~ family_history, contrast = "family_history_yes_vs_no", species = "human")
# saveRDS(ipsc_mn_answerals.family_history, "/camp/home/ziffo/home/projects/ipsc-mn-als-meta/expression/deseq2/ipsc_mn_answerals.family_history.rds")

# ### NYGC Post Mortem ------------------------------------------------------------------
# # counts
# GSE153960_Gene_counts_matrix_RSEM_Prudencio_et_al_2020.txt.col_names = read_tsv("/camp/project/proj-luscombn-patani/working/nygc-als-consortium/counts/GSE153960/GSE153960_Gene_counts_matrix_RSEM_Prudencio_et_al_2020.txt") %>% clean_names()
# GSE153960_Gene_counts_matrix_RSEM_Prudencio_et_al_2020.txt = read_tsv("/camp/project/proj-luscombn-patani/working/nygc-als-consortium/counts/GSE153960/GSE153960_Gene_counts_matrix_RSEM_Prudencio_et_al_2020.txt", skip = 1, col_names = FALSE) %>% select(-X1)
# colnames(GSE153960_Gene_counts_matrix_RSEM_Prudencio_et_al_2020.txt) = colnames(GSE153960_Gene_counts_matrix_RSEM_Prudencio_et_al_2020.txt.col_names)
# GSE153960_Gene_counts_matrix_RSEM_Prudencio_et_al_2020 = GSE153960_Gene_counts_matrix_RSEM_Prudencio_et_al_2020.txt %>% mutate(gene_id = str_split_fixed(ensembl_id, "\\.", 2)[,1]) %>% select(gene_id, everything(), -ensembl_id) %>% distinct(gene_id, .keep_all = TRUE)
# # metadata
# samples = colnames(GSE153960_Gene_counts_matrix_RSEM_Prudencio_et_al_2020.txt.col_names)
# GSE153960_metadata.txt = read_csv("/camp/project/proj-luscombn-patani/working/nygc-als-consortium/counts/GSE153960/GSE153960_metadata.txt") %>% clean_names() %>% mutate(sample = gsub("-","_",sample_id_alt), sample = tolower(sample))
# GSE153960_metadata.txt %>% filter(sample %in% samples) # 2,038 / 2,256
# GSE153960_metadata.txt %>% get_dupes(sample) # 792
# GSE153960_metadata = GSE153960_metadata.txt %>% filter(sample %in% samples) %>% select(sample, run, geo_accession_exp, instrument, library_preparation_method, project, sample_id_alt, subject_id, tissue, group) %>%
#   mutate(condition = case_when(group == "Non-Neurological Control" ~ "ctrl", grepl("ALS Spectrum", group) ~ "als", group == "Other MND" ~ "other_mnd", grepl("Other Neurological D", group) ~ "other_neuro"), tissue.detailed = to_any_case(tissue),
#          tissue.simplified = case_when(grepl("spinal",tissue.detailed) ~ "spinal", grepl("motor",tissue.detailed) ~ "motor_cortex", TRUE ~ "other"), library_preparation_method = to_any_case(library_preparation_method), ) %>%
#   distinct(sample, .keep_all = TRUE)
# GSE153960_metadata %>% count(condition)
# GSE153960_metadata %>% count(library_preparation_method, instrument)
# GSE153960_metadata %>% get_dupes(sample) # 0
# GSE153960_metadata.als_ctrl = GSE153960_metadata %>% filter(condition %in% c("als","ctrl")) %>% mutate(condition = factor(condition, levels = c("ctrl","als")), sample_clean = sample) %>% column_to_rownames(var="sample")
# GSE153960_metadata.als_ctrl %>% count(condition)
# GSE153960_Gene_counts_matrix.als_ctrl.mat =  GSE153960_Gene_counts_matrix_RSEM_Prudencio_et_al_2020 %>% select(gene_id, contains(rownames(GSE153960_metadata.als_ctrl))) %>% column_to_rownames(var = "gene_id") %>% as.matrix() %>% round() # keep counts only for samples in metadata
# GSE153960_metadata.als_ctrl.spinal_cord = GSE153960_metadata.als_ctrl %>% filter(tissue.detailed %in% c("spinal_cord_cervical","spinal_cord_thoracic","spinal_cord_lumbar")) # 436 ALS, 94 CTRL biopsy samples
# 
# # all spinalcord
# GSE153960_metadata.als_ctrl.spinal_cord %>% distinct(subject_id, .keep_all = TRUE) %>% count(condition) # ALS 206 CTRL 56 unique patients
# GSE153960_Gene_counts_matrix.als_ctrl.spinal_cord.mat =  GSE153960_Gene_counts_matrix_RSEM_Prudencio_et_al_2020 %>% select(gene_id, contains(rownames(GSE153960_metadata.als_ctrl.spinal_cord))) %>% column_to_rownames(var = "gene_id") %>% as.matrix() %>% round() # keep counts only for samples in metadata
# nygc_postmortem_spinal_cord.als_vs_ctrl.tissue = list()
# nygc_postmortem_spinal_cord.als_vs_ctrl.tissue$dds <- DESeqDataSetFromMatrix(countData = GSE153960_Gene_counts_matrix.als_ctrl.spinal_cord.mat, colData = GSE153960_metadata.als_ctrl.spinal_cord, design = ~ library_preparation_method + tissue + condition) %>% DESeq() # adjust for cevical / thoracic / lumbar spinal cord site
# nygc_postmortem_spinal_cord.als_vs_ctrl.tissue$vsd <- vst(nygc_postmortem_spinal_cord.als_vs_ctrl.tissue$dds, blind=TRUE)
# nygc_postmortem_spinal_cord.als_vs_ctrl.tissue$vsd.counts <- as_tibble(assay(nygc_postmortem_spinal_cord.als_vs_ctrl.tissue$vsd), rownames = "gene_id")  %>% left_join(gene2ens)
# nygc_postmortem_spinal_cord.als_vs_ctrl.tissue$res <- DESeq2::results(nygc_postmortem_spinal_cord.als_vs_ctrl.tissue$dds, name = "condition_als_vs_ctrl") %>% as_tibble(rownames = "gene_id") %>% left_join(gene2ens) %>% arrange(-abs(stat))
# saveRDS(nygc_postmortem_spinal_cord.als_vs_ctrl.tissue, "/camp/lab/luscomben/home/users/ziffo/projects/ipsc-mn-als-meta/expression/deseq2/nygc_postmortem_spinal_cord.als_vs_ctrl.tissue.rds")
# 
# # NYGC mutation metadata
# nygc_metadata_mutation= read_csv("/camp/project/proj-luscombn-patani/working/nygc-als-consortium/sample-details/samplesheet.csv") %>% mutate(sample_id_alt = sample, subject_id = external_subject_id) %>% select(subject_id, mutation) %>% distinct(across(everything())) %>% 
#   filter(mutation %in% c("c9orf72","fus","tardbp","vcp","sod1","ctrl","sporadic","none")) %>%
#   filter(!(mutation %in% c("none","sporadic") & subject_id %in% c("JHU79","NEUDG727NTW", "NEUKN179PH8", "NEUTB398FC9","NEUYY221HZA", "PF-UCL-111", "PF-UCL-113", "PF-UCL-115","PF-UCL-46","PF-UCL-52","PF-UCL-78")))
# nygc_metadata_mutation %>% get_dupes(subject_id)
# GSE153960_metadata.als_ctrl.spinal_cord.mutation = GSE153960_metadata.als_ctrl.spinal_cord %>% left_join(nygc_metadata_mutation) %>% drop_na(mutation) %>% mutate(mutation = case_when(condition == "ctrl" ~ "ctrl", TRUE ~ mutation)) %>% column_to_rownames(var="sample_clean")
# GSE153960_metadata.als_ctrl.spinal_cord.mutation %>% count(mutation) %>% arrange(-n)
# # mutation   n
# # 1   tardbp 277
# # 2     ctrl  94
# # 3  c9orf72  61
# # 4 sporadic  58
# # 5     sod1   8
# # 6      fus   4
# # 7      vcp   2
# 
# # C9orf72
# GSE153960_metadata.als_ctrl.spinal_cord.c9orf72 = GSE153960_metadata.als_ctrl.spinal_cord.mutation %>% filter(mutation %in% c("c9orf72","ctrl")) # 436 ALS, 94 CTRL biopsy samples
# GSE153960_Gene_counts_matrix.als_ctrl.spinal_cord.c9orf72.mat =  GSE153960_Gene_counts_matrix_RSEM_Prudencio_et_al_2020 %>% select(gene_id, contains(rownames(GSE153960_metadata.als_ctrl.spinal_cord.c9orf72))) %>% column_to_rownames(var = "gene_id") %>% as.matrix() %>% round() # keep counts only for samples in metadata
# nygc_postmortem_spinal_cord.c9orf72_vs_ctrl = list()
# nygc_postmortem_spinal_cord.c9orf72_vs_ctrl$dds <- DESeqDataSetFromMatrix(countData = GSE153960_Gene_counts_matrix.als_ctrl.spinal_cord.c9orf72.mat, colData = GSE153960_metadata.als_ctrl.spinal_cord.c9orf72, design = ~ library_preparation_method + tissue + condition) %>% DESeq()
# # nygc_postmortem_spinal_cord.c9orf72_vs_ctrl$vsd <- vst(nygc_postmortem_spinal_cord.c9orf72_vs_ctrl$dds, blind=TRUE)
# # nygc_postmortem_spinal_cord.c9orf72_vs_ctrl$vsd.counts <- as_tibble(assay(nygc_postmortem_spinal_cord.c9orf72_vs_ctrl$vsd), rownames = "gene_id")  %>% left_join(gene2ens)
# nygc_postmortem_spinal_cord.c9orf72_vs_ctrl$res <- DESeq2::results(nygc_postmortem_spinal_cord.c9orf72_vs_ctrl$dds, name = "condition_als_vs_ctrl") %>% as_tibble(rownames = "gene_id") %>% left_join(gene2ens) %>% arrange(-abs(stat))
# saveRDS(nygc_postmortem_spinal_cord.c9orf72_vs_ctrl, "/camp/lab/luscomben/home/users/ziffo/projects/ipsc-mn-als-meta/expression/deseq2/nygc_postmortem_spinal_cord.c9orf72_vs_ctrl.rds")
# 
# # tardbp
# GSE153960_metadata.als_ctrl.spinal_cord.tardbp = GSE153960_metadata.als_ctrl.spinal_cord.mutation %>% filter(mutation %in% c("tardbp","ctrl")) # 436 ALS, 94 CTRL biopsy samples
# GSE153960_Gene_counts_matrix.als_ctrl.spinal_cord.tardbp.mat =  GSE153960_Gene_counts_matrix_RSEM_Prudencio_et_al_2020 %>% select(gene_id, contains(rownames(GSE153960_metadata.als_ctrl.spinal_cord.tardbp))) %>% column_to_rownames(var = "gene_id") %>% as.matrix() %>% round() # keep counts only for samples in metadata
# nygc_postmortem_spinal_cord.tardbp_vs_ctrl = list()
# nygc_postmortem_spinal_cord.tardbp_vs_ctrl$dds <- DESeqDataSetFromMatrix(countData = GSE153960_Gene_counts_matrix.als_ctrl.spinal_cord.tardbp.mat, colData = GSE153960_metadata.als_ctrl.spinal_cord.tardbp, design = ~ library_preparation_method + tissue + condition) %>% DESeq()
# # nygc_postmortem_spinal_cord.tardbp_vs_ctrl$vsd <- vst(nygc_postmortem_spinal_cord.tardbp_vs_ctrl$dds, blind=TRUE)
# # nygc_postmortem_spinal_cord.tardbp_vs_ctrl$vsd.counts <- as_tibble(assay(nygc_postmortem_spinal_cord.tardbp_vs_ctrl$vsd), rownames = "gene_id")  %>% left_join(gene2ens)
# nygc_postmortem_spinal_cord.tardbp_vs_ctrl$res <- DESeq2::results(nygc_postmortem_spinal_cord.tardbp_vs_ctrl$dds, name = "condition_als_vs_ctrl") %>% as_tibble(rownames = "gene_id") %>% left_join(gene2ens) %>% arrange(-abs(stat))
# saveRDS(nygc_postmortem_spinal_cord.tardbp_vs_ctrl, "/camp/lab/luscomben/home/users/ziffo/projects/ipsc-mn-als-meta/expression/deseq2/nygc_postmortem_spinal_cord.tardbp_vs_ctrl.rds")
# 
# # sporadic
# GSE153960_metadata.als_ctrl.spinal_cord.sporadic = GSE153960_metadata.als_ctrl.spinal_cord.mutation %>% filter(mutation %in% c("sporadic","ctrl")) # 436 ALS, 94 CTRL biopsy samples
# GSE153960_Gene_counts_matrix.als_ctrl.spinal_cord.sporadic.mat =  GSE153960_Gene_counts_matrix_RSEM_Prudencio_et_al_2020 %>% select(gene_id, contains(rownames(GSE153960_metadata.als_ctrl.spinal_cord.sporadic))) %>% column_to_rownames(var = "gene_id") %>% as.matrix() %>% round() # keep counts only for samples in metadata
# nygc_postmortem_spinal_cord.sporadic_vs_ctrl = list()
# nygc_postmortem_spinal_cord.sporadic_vs_ctrl$dds <- DESeqDataSetFromMatrix(countData = GSE153960_Gene_counts_matrix.als_ctrl.spinal_cord.sporadic.mat, colData = GSE153960_metadata.als_ctrl.spinal_cord.sporadic, design = ~ library_preparation_method + tissue + condition) %>% DESeq()
# # nygc_postmortem_spinal_cord.sporadic_vs_ctrl$vsd <- vst(nygc_postmortem_spinal_cord.sporadic_vs_ctrl$dds, blind=TRUE)
# # nygc_postmortem_spinal_cord.sporadic_vs_ctrl$vsd.counts <- as_tibble(assay(nygc_postmortem_spinal_cord.sporadic_vs_ctrl$vsd), rownames = "gene_id")  %>% left_join(gene2ens)
# nygc_postmortem_spinal_cord.sporadic_vs_ctrl$res <- DESeq2::results(nygc_postmortem_spinal_cord.sporadic_vs_ctrl$dds, name = "condition_als_vs_ctrl") %>% as_tibble(rownames = "gene_id") %>% left_join(gene2ens) %>% arrange(-abs(stat))
# saveRDS(nygc_postmortem_spinal_cord.sporadic_vs_ctrl, "/camp/lab/luscomben/home/users/ziffo/projects/ipsc-mn-als-meta/expression/deseq2/nygc_postmortem_spinal_cord.sporadic_vs_ctrl.rds")
# 
# # fus
# GSE153960_metadata.als_ctrl.spinal_cord.fus = GSE153960_metadata.als_ctrl.spinal_cord.mutation %>% filter(mutation %in% c("fus","ctrl")) # 436 ALS, 94 CTRL biopsy samples
# GSE153960_Gene_counts_matrix.als_ctrl.spinal_cord.fus.mat =  GSE153960_Gene_counts_matrix_RSEM_Prudencio_et_al_2020 %>% select(gene_id, contains(rownames(GSE153960_metadata.als_ctrl.spinal_cord.fus))) %>% column_to_rownames(var = "gene_id") %>% as.matrix() %>% round() # keep counts only for samples in metadata
# nygc_postmortem_spinal_cord.fus_vs_ctrl = list()
# nygc_postmortem_spinal_cord.fus_vs_ctrl$dds <- DESeqDataSetFromMatrix(countData = GSE153960_Gene_counts_matrix.als_ctrl.spinal_cord.fus.mat, colData = GSE153960_metadata.als_ctrl.spinal_cord.fus, design = ~ library_preparation_method + tissue + condition) %>% DESeq()
# # nygc_postmortem_spinal_cord.fus_vs_ctrl$vsd <- vst(nygc_postmortem_spinal_cord.fus_vs_ctrl$dds, blind=TRUE)
# # nygc_postmortem_spinal_cord.fus_vs_ctrl$vsd.counts <- as_tibble(assay(nygc_postmortem_spinal_cord.fus_vs_ctrl$vsd), rownames = "gene_id")  %>% left_join(gene2ens)
# nygc_postmortem_spinal_cord.fus_vs_ctrl$res <- DESeq2::results(nygc_postmortem_spinal_cord.fus_vs_ctrl$dds, name = "condition_als_vs_ctrl") %>% as_tibble(rownames = "gene_id") %>% left_join(gene2ens) %>% arrange(-abs(stat))
# saveRDS(nygc_postmortem_spinal_cord.fus_vs_ctrl, "/camp/lab/luscomben/home/users/ziffo/projects/ipsc-mn-als-meta/expression/deseq2/nygc_postmortem_spinal_cord.fus_vs_ctrl.rds")
# 
# # sod1
# GSE153960_metadata.als_ctrl.spinal_cord.sod1 = GSE153960_metadata.als_ctrl.spinal_cord.mutation %>% filter(mutation %in% c("sod1","ctrl")) # 436 ALS, 94 CTRL biopsy samples
# GSE153960_Gene_counts_matrix.als_ctrl.spinal_cord.sod1.mat =  GSE153960_Gene_counts_matrix_RSEM_Prudencio_et_al_2020 %>% select(gene_id, contains(rownames(GSE153960_metadata.als_ctrl.spinal_cord.sod1))) %>% column_to_rownames(var = "gene_id") %>% as.matrix() %>% round() # keep counts only for samples in metadata
# nygc_postmortem_spinal_cord.sod1_vs_ctrl = list()
# nygc_postmortem_spinal_cord.sod1_vs_ctrl$dds <- DESeqDataSetFromMatrix(countData = GSE153960_Gene_counts_matrix.als_ctrl.spinal_cord.sod1.mat, colData = GSE153960_metadata.als_ctrl.spinal_cord.sod1, design = ~ library_preparation_method + tissue + condition) %>% DESeq()
# # nygc_postmortem_spinal_cord.sod1_vs_ctrl$vsd <- vst(nygc_postmortem_spinal_cord.sod1_vs_ctrl$dds, blind=TRUE)
# # nygc_postmortem_spinal_cord.sod1_vs_ctrl$vsd.counts <- as_tibble(assay(nygc_postmortem_spinal_cord.sod1_vs_ctrl$vsd), rownames = "gene_id")  %>% left_join(gene2ens)
# nygc_postmortem_spinal_cord.sod1_vs_ctrl$res <- DESeq2::results(nygc_postmortem_spinal_cord.sod1_vs_ctrl$dds, name = "condition_als_vs_ctrl") %>% as_tibble(rownames = "gene_id") %>% left_join(gene2ens) %>% arrange(-abs(stat))
# saveRDS(nygc_postmortem_spinal_cord.sod1_vs_ctrl, "/camp/lab/luscomben/home/users/ziffo/projects/ipsc-mn-als-meta/expression/deseq2/nygc_postmortem_spinal_cord.sod1_vs_ctrl.rds")

# # als_vs_ctrl across all tissues: design = ~ library_preparation_method + tissue.detailed + conditions
# nygc_postmortem_als_vs_ctrl = list()
# nygc_postmortem_als_vs_ctrl$dds <- DESeqDataSetFromMatrix(countData = GSE153960_Gene_counts_matrix.als_ctrl.mat, colData = GSE153960_metadata.als_ctrl, design = ~ library_preparation_method + tissue.detailed + condition) %>% DESeq()
# nygc_postmortem_als_vs_ctrl$vsd <- vst(nygc_postmortem_als_vs_ctrl$dds, blind=TRUE)
# nygc_postmortem_als_vs_ctrl$vsd.counts <- as_tibble(assay(nygc_postmortem_als_vs_ctrl$vsd), rownames = "gene_id")  %>% left_join(gene2ens)
# nygc_postmortem_als_vs_ctrl$res <- DESeq2::results(nygc_postmortem_als_vs_ctrl$dds, name = "condition_als_vs_ctrl") %>% as_tibble(rownames = "gene_id") %>% left_join(gene2ens) %>% arrange(-abs(stat))
# saveRDS(nygc_postmortem_als_vs_ctrl, "/camp/lab/luscomben/home/users/ziffo/projects/ipsc-mn-als-meta/expression/deseq2/nygc_postmortem_als_vs_ctrl.rds")
# 
# # tissue.detailed: using LRT: design = ~ library_preparation_method + tissue.detailed + condition + tissue.detailed:condition
# nygc_postmortem_tissue.detailed.als_vs_ctrl = list()
# nygc_postmortem_tissue.detailed.als_vs_ctrl$dds <- DESeqDataSetFromMatrix(countData = GSE153960_Gene_counts_matrix.als_ctrl.mat, colData = GSE153960_metadata.als_ctrl, design = ~ library_preparation_method + tissue.detailed + condition + tissue.detailed:condition) %>%
#   DESeq(test = "LRT", reduced = ~ library_preparation_method + tissue.detailed + condition)
# resultsNames(nygc_postmortem_tissue.detailed.als_vs_ctrl$dds)
# nygc_postmortem_tissue.detailed.als_vs_ctrl$res <- DESeq2::results(nygc_postmortem_tissue.detailed.als_vs_ctrl$dds) %>% as_tibble(rownames = "gene_id") %>% left_join(gene2ens) %>% arrange(-abs(stat))
# nygc_postmortem_tissue.detailed.als_vs_ctrl$vsd <- vst(nygc_postmortem_tissue.detailed.als_vs_ctrl$dds, blind=TRUE)
# nygc_postmortem_tissue.detailed.als_vs_ctrl$vsd.counts <- as_tibble(assay(nygc_postmortem_tissue.detailed.als_vs_ctrl$vsd), rownames = "gene_id") %>% left_join(gene2ens)
# print(resultsNames(nygc_postmortem_tissue.detailed.als_vs_ctrl$dds))
# # compare individual tissues using names argument
# saveRDS(nygc_postmortem_tissue.detailed.als_vs_ctrl, "/camp/lab/luscomben/home/users/ziffo/projects/ipsc-mn-als-meta/expression/deseq2/nygc_postmortem_tissue.detailed.als_vs_ctrl.rds")
# # 
# # tissue.simplified: using LRT: design = ~ library_preparation_method + tissue.simplified + condition + tissue.simplified:condition
# nygc_postmortem_tissue.simplified.als_vs_ctrl = list()
# nygc_postmortem_tissue.simplified.als_vs_ctrl$dds <- DESeqDataSetFromMatrix(countData = GSE153960_Gene_counts_matrix.als_ctrl.mat, colData = GSE153960_metadata.als_ctrl, design = ~ library_preparation_method + tissue.simplified + condition + tissue.simplified:condition) %>%
#   DESeq(test = "LRT", reduced = ~ library_preparation_method + tissue.simplified + condition)
# resultsNames(nygc_postmortem_tissue.simplified.als_vs_ctrl$dds)
# nygc_postmortem_tissue.simplified.als_vs_ctrl$res <- DESeq2::results(nygc_postmortem_tissue.simplified.als_vs_ctrl$dds) %>% as_tibble(rownames = "gene_id") %>% left_join(gene2ens) %>% arrange(-abs(stat))
# nygc_postmortem_tissue.simplified.als_vs_ctrl$vsd <- vst(nygc_postmortem_tissue.simplified.als_vs_ctrl$dds, blind=TRUE)
# nygc_postmortem_tissue.simplified.als_vs_ctrl$vsd.counts <- as_tibble(assay(nygc_postmortem_tissue.simplified.als_vs_ctrl$vsd), rownames = "gene_id") %>% left_join(gene2ens)
# print(resultsNames(nygc_postmortem_tissue.simplified.als_vs_ctrl$dds))
# # compare individual tissues using names argument: spinal_vs_motor_cortex
# saveRDS(nygc_postmortem_tissue.simplified.als_vs_ctrl, "/camp/lab/luscomben/home/users/ziffo/projects/ipsc-mn-als-meta/expression/deseq2/nygc_postmortem_tissue.simplified.als_vs_ctrl.rds")
# 
# # spinal_cord_cervical
# GSE153960_metadata.als_ctrl.cervical = GSE153960_metadata.als_ctrl %>% filter(tissue.detailed %in% c("spinal_cord_cervical"))
# GSE153960_Gene_counts_matrix.als_ctrl.cervical.mat =  GSE153960_Gene_counts_matrix_RSEM_Prudencio_et_al_2020 %>% select(gene_id, contains(rownames(GSE153960_metadata.als_ctrl.cervical))) %>% column_to_rownames(var = "gene_id") %>% as.matrix() %>% round() # keep counts only for samples in metadata
# nygc_postmortem_cervical.als_vs_ctrl = list()
# nygc_postmortem_cervical.als_vs_ctrl$dds <- DESeqDataSetFromMatrix(countData = GSE153960_Gene_counts_matrix.als_ctrl.cervical.mat, colData = GSE153960_metadata.als_ctrl.cervical, design = ~ library_preparation_method + condition) %>% DESeq()
# nygc_postmortem_cervical.als_vs_ctrl$vsd <- vst(nygc_postmortem_cervical.als_vs_ctrl$dds, blind=TRUE)
# nygc_postmortem_cervical.als_vs_ctrl$vsd.counts <- as_tibble(assay(nygc_postmortem_cervical.als_vs_ctrl$vsd), rownames = "gene_id")  %>% left_join(gene2ens)
# nygc_postmortem_cervical.als_vs_ctrl$res <- DESeq2::results(nygc_postmortem_cervical.als_vs_ctrl$dds, name = "condition_als_vs_ctrl") %>% as_tibble(rownames = "gene_id") %>% left_join(gene2ens) %>% arrange(-abs(stat))
# saveRDS(nygc_postmortem_cervical.als_vs_ctrl, "/camp/lab/luscomben/home/users/ziffo/projects/ipsc-mn-als-meta/expression/deseq2/nygc_postmortem_cervical.als_vs_ctrl.rds")
# 
# # spinal_cord_thoracic
# GSE153960_metadata.als_ctrl.thoracic = GSE153960_metadata.als_ctrl %>% filter(tissue.detailed %in% c("spinal_cord_thoracic"))
# GSE153960_Gene_counts_matrix.als_ctrl.thoracic.mat =  GSE153960_Gene_counts_matrix_RSEM_Prudencio_et_al_2020 %>% select(gene_id, contains(rownames(GSE153960_metadata.als_ctrl.thoracic))) %>% column_to_rownames(var = "gene_id") %>% as.matrix() %>% round() # keep counts only for samples in metadata
# nygc_postmortem_thoracic.als_vs_ctrl = list()
# nygc_postmortem_thoracic.als_vs_ctrl$dds <- DESeqDataSetFromMatrix(countData = GSE153960_Gene_counts_matrix.als_ctrl.thoracic.mat, colData = GSE153960_metadata.als_ctrl.thoracic, design = ~ library_preparation_method + condition) %>% DESeq()
# nygc_postmortem_thoracic.als_vs_ctrl$vsd <- vst(nygc_postmortem_thoracic.als_vs_ctrl$dds, blind=TRUE)
# nygc_postmortem_thoracic.als_vs_ctrl$vsd.counts <- as_tibble(assay(nygc_postmortem_thoracic.als_vs_ctrl$vsd), rownames = "gene_id")  %>% left_join(gene2ens)
# nygc_postmortem_thoracic.als_vs_ctrl$res <- DESeq2::results(nygc_postmortem_thoracic.als_vs_ctrl$dds, name = "condition_als_vs_ctrl") %>% as_tibble(rownames = "gene_id") %>% left_join(gene2ens) %>% arrange(-abs(stat)) 
# saveRDS(nygc_postmortem_thoracic.als_vs_ctrl, "/camp/lab/luscomben/home/users/ziffo/projects/ipsc-mn-als-meta/expression/deseq2/nygc_postmortem_thoracic.als_vs_ctrl.rds")
# 
# # spinal_cord_lumbar
# GSE153960_metadata.als_ctrl.lumbar = GSE153960_metadata.als_ctrl %>% filter(tissue.detailed %in% c("spinal_cord_lumbar"))
# GSE153960_Gene_counts_matrix.als_ctrl.lumbar.mat =  GSE153960_Gene_counts_matrix_RSEM_Prudencio_et_al_2020 %>% select(gene_id, contains(rownames(GSE153960_metadata.als_ctrl.lumbar))) %>% column_to_rownames(var = "gene_id") %>% as.matrix() %>% round() # keep counts only for samples in metadata
# nygc_postmortem_lumbar.als_vs_ctrl = list()
# nygc_postmortem_lumbar.als_vs_ctrl$dds <- DESeqDataSetFromMatrix(countData = GSE153960_Gene_counts_matrix.als_ctrl.lumbar.mat, colData = GSE153960_metadata.als_ctrl.lumbar, design = ~ library_preparation_method + condition) %>% DESeq()
# nygc_postmortem_lumbar.als_vs_ctrl$vsd <- vst(nygc_postmortem_lumbar.als_vs_ctrl$dds, blind=TRUE)
# nygc_postmortem_lumbar.als_vs_ctrl$vsd.counts <- as_tibble(assay(nygc_postmortem_lumbar.als_vs_ctrl$vsd), rownames = "gene_id")  %>% left_join(gene2ens)
# nygc_postmortem_lumbar.als_vs_ctrl$res <- DESeq2::results(nygc_postmortem_lumbar.als_vs_ctrl$dds, name = "condition_als_vs_ctrl") %>% as_tibble(rownames = "gene_id") %>% left_join(gene2ens) %>% arrange(-abs(stat)) 
# saveRDS(nygc_postmortem_lumbar.als_vs_ctrl, "/camp/lab/luscomben/home/users/ziffo/projects/ipsc-mn-als-meta/expression/deseq2/nygc_postmortem_lumbar.als_vs_ctrl.rds")
# 
# # motor_cortex: cortex_motor_lateral, cortex_motor_medial, cortex_motor_unspecified
# GSE153960_metadata.als_ctrl.motor_cortex = GSE153960_metadata.als_ctrl %>% filter(tissue.detailed %in% c("cortex_motor_lateral","cortex_motor_medial","cortex_motor_unspecified"))
# GSE153960_Gene_counts_matrix.als_ctrl.motor_cortex.mat =  GSE153960_Gene_counts_matrix_RSEM_Prudencio_et_al_2020 %>% select(gene_id, contains(rownames(GSE153960_metadata.als_ctrl.motor_cortex))) %>% column_to_rownames(var = "gene_id") %>% as.matrix() %>% round() # keep counts only for samples in metadata
# nygc_postmortem_motor_cortex.als_vs_ctrl = list()
# nygc_postmortem_motor_cortex.als_vs_ctrl$dds <- DESeqDataSetFromMatrix(countData = GSE153960_Gene_counts_matrix.als_ctrl.motor_cortex.mat, colData = GSE153960_metadata.als_ctrl.motor_cortex, design = ~ library_preparation_method + tissue + condition) %>% DESeq()
# nygc_postmortem_motor_cortex.als_vs_ctrl$vsd <- vst(nygc_postmortem_motor_cortex.als_vs_ctrl$dds, blind=TRUE)
# nygc_postmortem_motor_cortex.als_vs_ctrl$vsd.counts <- as_tibble(assay(nygc_postmortem_motor_cortex.als_vs_ctrl$vsd), rownames = "gene_id")  %>% left_join(gene2ens)
# nygc_postmortem_motor_cortex.als_vs_ctrl$res <- DESeq2::results(nygc_postmortem_motor_cortex.als_vs_ctrl$dds, name = "condition_als_vs_ctrl") %>% as_tibble(rownames = "gene_id") %>% left_join(gene2ens) %>% arrange(-abs(stat)) 
# saveRDS(nygc_postmortem_motor_cortex.als_vs_ctrl, "/camp/lab/luscomben/home/users/ziffo/projects/ipsc-mn-als-meta/expression/deseq2/nygc_postmortem_motor_cortex.als_vs_ctrl.rds")

## Knockdown ------------------------------------------------------------------
# tdp43 knockdown datasets
# tdp43_kd_brown.metadata <- read_csv("/camp/project/proj-luscombn-patani/working/public-data/ipsc-mn-tardbp-kd-brown-2022/sample-details/samplesheet.csv") %>% distinct(sample, .keep_all = TRUE) %>%
#   mutate(database_dir = "/camp/project/proj-luscombn-patani/working/public-data/ipsc-mn-tardbp-kd-brown-2022", condition = factor(condition, levels = c("ctrl", "tdp43kd")), dataset = "brown")
# tdp43_kd_klim.metadata <- read_csv("/camp/project/proj-luscombn-patani/working/public-data/ipsc-mn-tardbp-kd-klim-2019/sample-details/samplesheet.csv") %>% distinct(sample, .keep_all = TRUE) %>%
#   mutate(database_dir = "/camp/project/proj-luscombn-patani/working/public-data/ipsc-mn-tardbp-kd-klim-2019", condition = factor(condition, levels = c("ctrl", "tdp43kd")), dataset = "klim")
# tdp43_kd_appocher.metadata <- read_csv("/camp/project/proj-luscombn-patani/working/public-data/shsy5y-tardbp-kd-appocher-2017/sample-details/samplesheet.csv") %>% distinct(sample, .keep_all = TRUE) %>%
#   mutate(database_dir = "/camp/project/proj-luscombn-patani/working/public-data/shsy5y-tardbp-kd-appocher-2017", condition = factor(condition, levels = c("ctrl", "tdp43kd")), dataset = "appocher")
# neuronal_nuclei_tdp43_liu.metadata <- read_csv("/camp/project/proj-luscombn-patani/working/public-data/neuronal-nuclei-tdp43-liu-2019/sample-details/samplesheet.csv") %>% distinct(sample, .keep_all = TRUE) %>%
#   mutate(database_dir = "/camp/project/proj-luscombn-patani/working/public-data/neuronal-nuclei-tdp43-liu-2019", group = condition, condition = case_when(condition == "tdp43pos" ~ "ctrl", condition == "tdp43neg" ~ "tdp43kd"), condition = factor(condition, levels = c("ctrl", "tdp43kd")), dataset = "liu")
# kd_tyzack.metadata <- read_csv("/camp/lab/luscomben/home/shared/projects/patani-collab/motor-neuron-rbp-knockdown-rnaseq/sample-details/samplesheet.csv") %>% filter(fraction == "whole", group != "untreated") %>%
#   mutate(database_dir = "/camp/lab/luscomben/home/shared/projects/patani-collab/motor-neuron-rbp-knockdown-rnaseq", dataset = "tyzack", condition = group) %>%
#   distinct(sample, .keep_all = TRUE) %>% select(sample, dataset, condition, group, database_dir, everything())
# tdp43_kd_tyzack.metadata <- kd_tyzack.metadata %>% filter(group %in% c("scramble","tdp43kd")) %>% mutate(condition = case_when(group == "scramble" ~ "ctrl", TRUE ~ "tdp43kd"), condition = factor(condition, levels = c("ctrl", "tdp43kd")))
# fus_kd_tyzack.metadata <- kd_tyzack.metadata %>% filter(group %in% c("scramble","fuskd")) %>% mutate(condition = case_when(group == "scramble" ~ "ctrl", TRUE ~ "fuskd"), condition = factor(condition, levels = c("ctrl", "fuskd")))
# sfpq_kd_tyzack.metadata <- kd_tyzack.metadata %>% filter(group %in% c("scramble","sfpqkd")) %>% mutate(condition = case_when(group == "scramble" ~ "ctrl", TRUE ~ "sfpqkd"), condition = factor(condition, levels = c("ctrl", "sfpqkd")))
# tdp43_kd_datasets.metadata <- bind_rows(tdp43_kd_brown.metadata, tdp43_kd_klim.metadata, neuronal_nuclei_tdp43_liu.metadata, tdp43_kd_tyzack.metadata, tdp43_kd_appocher.metadata) %>%
#   mutate(dataset_sample = paste0(dataset,"_",sample), condition = factor(condition, levels = c("ctrl", "tdp43kd")), dataset = as.factor(dataset)) %>%
#   select(sample, dataset, condition, dataset_sample, replicate, database_dir, everything())
# tdp43_kd_datasets_notyzack.metadata <- bind_rows(tdp43_kd_brown.metadata, tdp43_kd_klim.metadata, neuronal_nuclei_tdp43_liu.metadata, tdp43_kd_appocher.metadata) %>%
#   mutate(dataset_sample = paste0(dataset,"_",sample), condition = factor(condition, levels = c("ctrl", "tdp43kd")), dataset = as.factor(dataset)) %>%
#   select(sample, dataset, condition, dataset_sample, database_dir, everything())

# kd_tyzack = DESeq.analysis(metadata = kd_tyzack.metadata, unique_names = "sample", design = ~ 1)
# saveRDS(kd_tyzack, "/camp/home/ziffo/home/projects/ipsc-mn-als-meta/expression/deseq2/kd_tyzack.rds")
# tdp43_kd_tyzack = DESeq.analysis(metadata = tdp43_kd_tyzack.metadata, unique_names = "sample", design = ~ cellline + condition, contrast = "condition_tdp43kd_vs_ctrl", run_irfinder = TRUE, irfinder_design = ~ cellline + condition + condition:IRFinder)
# saveRDS(tdp43_kd_tyzack, "/camp/home/ziffo/home/projects/ipsc-mn-als-meta/expression/deseq2/tdp43_kd_tyzack.rds")
# fus_kd_tyzack = DESeq.analysis(metadata = fus_kd_tyzack.metadata, unique_names = "sample", design = ~ cellline + condition, contrast = "condition_fuskd_vs_ctrl", run_irfinder = TRUE, irfinder_design = ~ cellline + condition + condition:IRFinder)
# saveRDS(fus_kd_tyzack, "/camp/home/ziffo/home/projects/ipsc-mn-als-meta/expression/deseq2/fus_kd_tyzack.rds")
# sfpq_kd_tyzack = DESeq.analysis(metadata = sfpq_kd_tyzack.metadata, unique_names = "sample", design = ~ cellline + condition, contrast = "condition_sfpqkd_vs_ctrl", run_irfinder = TRUE, irfinder_design = ~ cellline + condition + condition:IRFinder)
# saveRDS(sfpq_kd_tyzack, "/camp/home/ziffo/home/projects/ipsc-mn-als-meta/expression/deseq2/sfpq_kd_tyzack.rds")
# tdp43_kd_brown = DESeq.analysis(metadata = tdp43_kd_brown.metadata, unique_names = "sample", design = ~ condition, contrast = "condition_tdp43kd_vs_ctrl", run_irfinder = TRUE, irfinder_design = ~ condition + condition:IRFinder)
# saveRDS(tdp43_kd_brown, "/camp/home/ziffo/home/projects/ipsc-mn-als-meta/expression/deseq2/tdp43_kd_brown.rds")
# tdp43_kd_klim = DESeq.analysis(metadata = tdp43_kd_klim.metadata, unique_names = "sample", design = ~ condition, contrast = "condition_tdp43kd_vs_ctrl", run_irfinder = TRUE, irfinder_design = ~ condition + condition:IRFinder)
# saveRDS(tdp43_kd_klim, "/camp/home/ziffo/home/projects/ipsc-mn-als-meta/expression/deseq2/tdp43_kd_klim.rds")
# neuronal_nuclei_tdp43_liu = DESeq.analysis(metadata = neuronal_nuclei_tdp43_liu.metadata, unique_names = "sample", design = ~ condition, contrast = "condition_tdp43kd_vs_ctrl", run_irfinder = TRUE, irfinder_design = ~ condition + condition:IRFinder)
# saveRDS(neuronal_nuclei_tdp43_liu, "/camp/home/ziffo/home/projects/ipsc-mn-als-meta/expression/deseq2/neuronal_nuclei_tdp43_liu.rds")
# tdp43_kd_appocher = DESeq.analysis(metadata = tdp43_kd_appocher.metadata, unique_names = "sample", design = ~ condition, contrast = "condition_tdp43kd_vs_ctrl", run_irfinder = TRUE, irfinder_design = ~ condition + condition:IRFinder)
# saveRDS(tdp43_kd_appocher, "/camp/home/ziffo/home/projects/ipsc-mn-als-meta/expression/deseq2/tdp43_kd_appocher.rds")
# tdp43_kd_datasets_notyzack = DESeq.analysis(metadata = tdp43_kd_datasets_notyzack.metadata, unique_names = "dataset_sample", design = ~ dataset + condition, contrast = "condition_tdp43kd_vs_ctrl", run_irfinder = FALSE, irfinder_design = ~ dataset + condition + condition:IRFinder)
# saveRDS(tdp43_kd_datasets_notyzack, "/camp/home/ziffo/home/projects/ipsc-mn-als-meta/expression/deseq2/tdp43_kd_datasets_notyzack.rds")
# tdp43_kd_datasets = DESeq.analysis(metadata = tdp43_kd_datasets.metadata, unique_names = "dataset_sample", design = ~ dataset + condition, contrast = "condition_tdp43kd_vs_ctrl", run_irfinder = TRUE, irfinder_design = ~ dataset + condition + condition:IRFinder)
# saveRDS(tdp43_kd_datasets, "/camp/home/ziffo/home/projects/ipsc-mn-als-meta/expression/deseq2/tdp43_kd_datasets.rds")

# print("saving")
# save.image("/camp/home/ziffo/home/projects/ipsc-mn-als-meta/expression/deseq2/ipsc_mn_als_meta.RData")
print("complete")
