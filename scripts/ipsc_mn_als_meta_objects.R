# run this script in Rstudio server OnDemand Terminal Tab (not with conda env r4.0.3) with:
# cd /camp/home/ziffo/home/projects/ipsc-mn-als-meta/scripts
# sbatch -N 1 -c 12 --mem=0 -t 72:00:00 --wrap="Rscript /camp/lab/luscomben/home/users/ziffo/projects/ipsc-mn-als-meta/scripts/ipsc_mn_als_meta/ipsc_mn_als_meta_objects.R" --mail-type=ALL,ARRAY_TASKS --mail-user=oliver.ziff@crick.ac.uk --job-name=DESeq2objects
# sbatch -N 1 -c 12 --mem=72G -t 72:00:00 --wrap="Rscript /camp/lab/luscomben/home/users/ziffo/projects/ipsc-mn-als-meta/scripts/ipsc_mn_als_meta/ipsc_mn_als_meta_objects.R" --mail-type=ALL,ARRAY_TASKS --mail-user=oliver.ziff@crick.ac.uk --job-name=DESeq2objects --dependency=afterany:38886199
# sbatch -N 1 -c 6 --mem=50G -t 8:00:00 --wrap="Rscript /camp/lab/luscomben/home/users/ziffo/projects/ipsc-mn-als-meta/scripts/ipsc_mn_als_meta/ipsc_mn_als_meta_objects.R" --mail-type=ALL,ARRAY_TASKS --mail-user=oliver.ziff@crick.ac.uk --job-name=DESeq2objects 
# sbatch --cpus-per-task=10 --mem=1500G -N 1 --partition=hmem -t 72:00:00 --wrap="Rscript /camp/lab/luscomben/home/users/ziffo/projects/ipsc-mn-als-meta/scripts/ipsc_mn_als_meta/ipsc_mn_als_meta_objects.R" --mail-type=ALL,ARRAY_TASKS --mail-user=oliver.ziff@crick.ac.uk --job-name=DESeq2objects

# source("/Volumes/lab-luscomben/home/users/ziffo/scripts/functions/mac_R_functions.R") 
source("/camp/lab/luscomben/home/users/ziffo/scripts/functions/OnDemand_R_functions.R")
# camp_path = here("/Volumes/lab-luscomben/home/users/ziffo")
camp_path = here("/camp/lab/luscomben/home/users/ziffo")
# shared_path = here("/Volumes/lab-luscomben/home/users/ziffo/proj-luscombn-patani/working")
shared_path = here("/camp/project/proj-luscombn-patani/working")
# collab_path = here("/Volumes/lab-luscomben/home/users/ziffo/patani-collab")
collab_path = here("/camp/lab/luscomben/home/users/ziffo/patani-collab")
# proj_path = here("/Volumes/lab-luscomben/home/users/ziffo/projects/ipsc-mn-als-meta")
proj_path = here("/camp/lab/luscomben/home/users/ziffo/projects/ipsc-mn-als-meta")
# load(here(proj_path,"expression/deseq2/ipsc_mn_als_meta.metadata.RData")) # overwrites with old functions


          
# Samplesheets ----------------------------------------------------------------
answerals.gene_tpm.male_counts = read_tsv(here(shared_path,"answerals/nfcore/star_salmon/salmon.merged.gene_tpm.tsv")) %>% filter(gene_name %in% c('RPS4Y1', 'EIFAY', 'DDX3Y', 'KDM5D')) %>% select(-gene_id) %>% adorn_totals("row", name = "male_counts") %>% as_tibble %>% 
  filter(gene_name == "male_counts") %>% transpose_tibble(old_rownames = "gene_name", new_rownames = "sample") %>% mutate(gender = case_when(male_counts > 10 ~ "male", TRUE ~ "female"))
answerals.metadata = read_csv(here(shared_path,"answerals/sample-details/samplesheet.csv")) %>% distinct(sample, .keep_all = TRUE) %>%
  mutate(database_dir = here(shared_path,"answerals"), dataset = "answerals", dataset_sample = paste0(dataset,"_",sample), study_accession = dataset, library_layout = "PAIRED") %>%
  filter(condition %in% c("als", "ctrl")) %>% # remove 1 C9orf72 carrier and 4 other mnd
  filter(!(condition == "ctrl" & mutation %in% c("optn", "als2"))) %>% # remove the 2 controls that have mutations in OPTN and ALS2
  mutate(condition = factor(condition, levels = c("ctrl", "als"))) %>% left_join(select(answerals.gene_tpm.male_counts, sample, gender))
neurolincs.gene_tpm.male_counts = read_tsv(here(shared_path,"neurolincs/nfcore/star_salmon/salmon.merged.gene_tpm.tsv")) %>% filter(gene_name %in% c('RPS4Y1', 'EIFAY', 'DDX3Y', 'KDM5D')) %>% select(-gene_id) %>% adorn_totals("row", name = "male_counts") %>% as_tibble %>% 
  filter(gene_name == "male_counts") %>% transpose_tibble(old_rownames = "gene_name", new_rownames = "sample") %>% mutate(gender = case_when(male_counts > 10 ~ "male", TRUE ~ "female"), sample = gsub("^X","",sample), sample = gsub("\\.","-",sample))
neurolincs.metadata = read_csv(here(shared_path,"neurolincs/sample-details/samplesheet.csv")) %>% distinct(sample, .keep_all = TRUE) %>% filter(DIV != 0) %>% # diMN = 0007..., iMN A-042-... labels
  mutate(database_dir = here(shared_path,"neurolincs"), condition = factor(condition, levels = c("ctrl", "als")), dataset = case_when(dataset == "neurolincs0" ~ "neurolincs.diMN",dataset == "neurolincsA" ~ "neurolincs.iMN"), study_accession = dataset) %>% left_join(select(neurolincs.gene_tpm.male_counts, sample, gender))
neurolincs.ipsc.metadata = read_csv(here(shared_path,"neurolincs/sample-details/samplesheet.csv")) %>% distinct(sample, .keep_all = TRUE) %>% filter(DIV == 0) %>% # only iPSCs
  mutate(database_dir = here(shared_path,"neurolincs"), condition = factor(condition, levels = c("ctrl", "als")), dataset = case_when(dataset == "neurolincsI" ~ "neurolincs.iPSC"), study_accession = dataset) %>% left_join(select(neurolincs.gene_tpm.male_counts, sample, gender))
catanese.gene_tpm.male_counts = read_tsv(here(shared_path,"public-data/ipsc-mn-c9orf72-fus-catanese-2021/nfcore/star_salmon/salmon.merged.gene_tpm.tsv")) %>% filter(gene_name %in% c('RPS4Y1', 'EIFAY', 'DDX3Y', 'KDM5D')) %>% select(-gene_id) %>% adorn_totals("row", name = "male_counts") %>% as_tibble %>% 
  filter(gene_name == "male_counts") %>% transpose_tibble(old_rownames = "gene_name", new_rownames = "sample") %>% mutate(gender = case_when(male_counts > 10 ~ "male", TRUE ~ "female"))
catanese.metadata <- read_csv(here(shared_path,"public-data/ipsc-mn-c9orf72-fus-catanese-2021/sample-details/samplesheet.csv")) %>% distinct(sample, .keep_all = TRUE) %>%
  mutate(database_dir = here(shared_path,"public-data/ipsc-mn-c9orf72-fus-catanese-2021"), dataset = "catanese", DIV = 35)%>% select(-age) %>% left_join(select(catanese.gene_tpm.male_counts, sample, gender))
dafinca.c9orf72.gene_tpm.male_counts = read_tsv(here(shared_path,"public-data/ipsc-mn-c9orf72-dafinca-2020/nfcore/star_salmon/salmon.merged.gene_tpm.tsv")) %>% filter(gene_name %in% c('RPS4Y1', 'EIFAY', 'DDX3Y', 'KDM5D')) %>% select(-gene_id) %>% adorn_totals("row", name = "male_counts") %>% as_tibble %>% 
  filter(gene_name == "male_counts") %>% transpose_tibble(old_rownames = "gene_name", new_rownames = "sample") %>% mutate(gender = case_when(male_counts > 10 ~ "male", TRUE ~ "female"))
dafinca.c9orf72.metadata <- read_csv(here(shared_path,"public-data/ipsc-mn-c9orf72-dafinca-2020/sample-details/samplesheet.csv")) %>% distinct(sample, .keep_all = TRUE) %>%
  mutate(database_dir = here(shared_path,"public-data/ipsc-mn-c9orf72-dafinca-2020"), dataset = "dafinca.c9orf72", DIV = 19) %>% left_join(select(dafinca.c9orf72.gene_tpm.male_counts, sample, gender))
dafinca.tardbp.gene_tpm.male_counts = read_tsv(here(shared_path,"public-data/ipsc-mn-tardbp-dafinca-2020/nfcore/star_salmon/salmon.merged.gene_tpm.tsv")) %>% filter(gene_name %in% c('RPS4Y1', 'EIFAY', 'DDX3Y', 'KDM5D')) %>% select(-gene_id) %>% adorn_totals("row", name = "male_counts") %>% as_tibble %>% 
  filter(gene_name == "male_counts") %>% transpose_tibble(old_rownames = "gene_name", new_rownames = "sample") %>% mutate(gender = case_when(male_counts > 10 ~ "male", TRUE ~ "female"))
dafinca.tardbp.metadata <- read_csv(here(shared_path,"public-data/ipsc-mn-tardbp-dafinca-2020/sample-details/samplesheet.csv")) %>% distinct(sample, .keep_all = TRUE) %>%
  mutate(database_dir = here(shared_path,"public-data/ipsc-mn-tardbp-dafinca-2020"), dataset = "dafinca.tardbp", DIV = 19) %>% left_join(select(dafinca.tardbp.gene_tpm.male_counts, sample, gender))
luisier.gene_tpm.male_counts = read_tsv(here(collab_path,"motor-neuron-vcp-luisier-2018/nfcore/star_salmon/salmon.merged.gene_tpm.tsv")) %>% filter(gene_name %in% c('RPS4Y1', 'EIFAY', 'DDX3Y', 'KDM5D')) %>% select(-gene_id) %>% adorn_totals("row", name = "male_counts") %>% as_tibble %>% 
  filter(gene_name == "male_counts") %>% transpose_tibble(old_rownames = "gene_name", new_rownames = "sample") %>% mutate(gender = case_when(male_counts > 10 ~ "male", TRUE ~ "female"))
luisier.metadata <- read_csv(here(collab_path,"motor-neuron-vcp-luisier-2018/sample-details/samplesheet.csv")) %>% distinct(sample, .keep_all = TRUE) %>%  filter(day == 35) %>%
  mutate(database_dir = here(collab_path,"motor-neuron-vcp-luisier-2018"), dataset = "luisier", total_size = bases, name = gsub("d35_","",sample), DIV = 35) %>% left_join(select(luisier.gene_tpm.male_counts, sample, gender))
wang.gene_tpm.male_counts = read_tsv(here(shared_path,"public-data/ipsc-mn-sod1-wang-2017/nfcore/star_salmon/salmon.merged.gene_tpm.tsv")) %>% filter(gene_name %in% c('RPS4Y1', 'EIFAY', 'DDX3Y', 'KDM5D')) %>% select(-gene_id) %>% adorn_totals("row", name = "male_counts") %>% as_tibble %>% 
  filter(gene_name == "male_counts") %>% transpose_tibble(old_rownames = "gene_name", new_rownames = "sample") %>% mutate(gender = case_when(male_counts > 10 ~ "male", TRUE ~ "female"))
wang.metadata <- read_csv(here(shared_path,"public-data/ipsc-mn-sod1-wang-2017/sample-details/samplesheet.csv")) %>% distinct(sample, .keep_all = TRUE) %>%
  mutate(database_dir = here(shared_path,"public-data/ipsc-mn-sod1-wang-2017"), dataset = "wang", DIV = 12) %>% left_join(select(wang.gene_tpm.male_counts, sample, gender))
kiskinis.gene_tpm.male_counts = read_tsv(here(shared_path,"public-data/ipsc-mn-sod1-moccia-2014/nfcore/star_salmon/salmon.merged.gene_tpm.tsv")) %>% filter(gene_name %in% c('RPS4Y1', 'EIFAY', 'DDX3Y', 'KDM5D')) %>% select(-gene_id) %>% adorn_totals("row", name = "male_counts") %>% as_tibble %>% 
  filter(gene_name == "male_counts") %>% transpose_tibble(old_rownames = "gene_name", new_rownames = "sample") %>% mutate(gender = case_when(male_counts > 10 ~ "male", TRUE ~ "female"))
kiskinis.metadata <- read_csv(here(shared_path,"public-data/ipsc-mn-sod1-moccia-2014/sample-details/samplesheet.csv")) %>% distinct(sample, .keep_all = TRUE) %>%
  mutate(database_dir = here(shared_path,"public-data/ipsc-mn-sod1-moccia-2014"), dataset = "kiskinis", DIV = 45) %>% left_join(select(kiskinis.gene_tpm.male_counts, sample, gender))
sareen.gene_tpm.male_counts = read_tsv(here(shared_path,"public-data/ipsc-mn-c9orf72-sareen-2013/nfcore/star_salmon/salmon.merged.gene_tpm.tsv")) %>% filter(gene_name %in% c('RPS4Y1', 'EIFAY', 'DDX3Y', 'KDM5D')) %>% select(-gene_id) %>% adorn_totals("row", name = "male_counts") %>% as_tibble %>% 
  filter(gene_name == "male_counts") %>% transpose_tibble(old_rownames = "gene_name", new_rownames = "sample") %>% mutate(gender = case_when(male_counts > 10 ~ "male", TRUE ~ "female"))
sareen.metadata <- read_csv(here(shared_path,"public-data/ipsc-mn-c9orf72-sareen-2013/sample-details/samplesheet.csv")) %>% distinct(sample, .keep_all = TRUE) %>%
  mutate(database_dir = here(shared_path,"public-data/ipsc-mn-c9orf72-sareen-2013"), dataset = "sareen", DIV = 35) %>% left_join(select(sareen.gene_tpm.male_counts, sample, gender))
sommer.gene_tpm.male_counts = read_tsv(here(shared_path,"public-data/ipsc-mn-c9orf72-sommer-2022/nfcore/star_salmon/salmon.merged.gene_tpm.tsv")) %>% filter(gene_name %in% c('RPS4Y1', 'EIFAY', 'DDX3Y', 'KDM5D')) %>% select(-gene_id) %>% adorn_totals("row", name = "male_counts") %>% as_tibble %>% 
  filter(gene_name == "male_counts") %>% transpose_tibble(old_rownames = "gene_name", new_rownames = "sample") %>% mutate(gender = case_when(male_counts > 10 ~ "male", TRUE ~ "female"))
sommer.metadata <- read_csv(here(shared_path,"public-data/ipsc-mn-c9orf72-sommer-2022/sample-details/samplesheet.csv")) %>% distinct(sample, .keep_all = TRUE) %>%
  mutate(database_dir = here(shared_path,"public-data/ipsc-mn-c9orf72-sommer-2022"), dataset = "sommer") %>% left_join(select(sommer.gene_tpm.male_counts, sample, gender))
sterneckert.gene_tpm.male_counts = read_tsv(here(shared_path,"public-data/ipsc-mn-c9orf72-sterneckert-2020/nfcore/star_salmon/salmon.merged.gene_tpm.tsv")) %>% filter(gene_name %in% c('RPS4Y1', 'EIFAY', 'DDX3Y', 'KDM5D')) %>% select(-gene_id) %>% adorn_totals("row", name = "male_counts") %>% as_tibble %>% 
  filter(gene_name == "male_counts") %>% transpose_tibble(old_rownames = "gene_name", new_rownames = "sample") %>% mutate(gender = case_when(male_counts > 10 ~ "male", TRUE ~ "female"))
sterneckert.metadata <- read_csv(here(shared_path,"public-data/ipsc-mn-c9orf72-sterneckert-2020/sample-details/samplesheet.csv")) %>% distinct(sample, .keep_all = TRUE) %>%
  mutate(database_dir = here(shared_path,"public-data/ipsc-mn-c9orf72-sterneckert-2020"), dataset = "sterneckert", DIV = 40) %>% left_join(select(sterneckert.gene_tpm.male_counts, sample, gender))
kapeli.gene_tpm.male_counts = read_tsv(here(shared_path,"public-data/ipsc-mn-fus-kapeli-2016/nfcore/star_salmon/salmon.merged.gene_tpm.tsv")) %>% filter(gene_name %in% c('RPS4Y1', 'EIFAY', 'DDX3Y', 'KDM5D')) %>% select(-gene_id) %>% adorn_totals("row", name = "male_counts") %>% as_tibble %>% 
  filter(gene_name == "male_counts") %>% transpose_tibble(old_rownames = "gene_name", new_rownames = "sample") %>% mutate(gender = case_when(male_counts > 10 ~ "male", TRUE ~ "female"))
kapeli.metadata <- read_csv(here(shared_path,"public-data/ipsc-mn-fus-kapeli-2016/sample-details/samplesheet.csv")) %>% distinct(sample, .keep_all = TRUE) %>%  filter(sh_rna != "scramble") %>%
  mutate(database_dir = here(shared_path,"public-data/ipsc-mn-fus-kapeli-2016"), dataset = "kapeli", DIV = 35) %>% left_join(select(kapeli.gene_tpm.male_counts, sample, gender))
desantis.gene_tpm.male_counts = read_tsv(here(shared_path,"public-data/ipsc-mn-fus-desantis-2017/nfcore/star_salmon/salmon.merged.gene_tpm.tsv")) %>% filter(gene_name %in% c('RPS4Y1', 'EIFAY', 'DDX3Y', 'KDM5D')) %>% select(-gene_id) %>% adorn_totals("row", name = "male_counts") %>% as_tibble %>% 
  filter(gene_name == "male_counts") %>% transpose_tibble(old_rownames = "gene_name", new_rownames = "sample") %>% mutate(gender = case_when(male_counts > 10 ~ "male", TRUE ~ "female"))
desantis.metadata <- read_csv(here(shared_path,"public-data/ipsc-mn-fus-desantis-2017/sample-details/samplesheet.csv")) %>% distinct(sample, .keep_all = TRUE) %>%
  mutate(database_dir = here(shared_path,"public-data/ipsc-mn-fus-desantis-2017"), DIV = 19, dataset = "desantis") %>% select(-days_of_differentiation) %>% left_join(select(desantis.gene_tpm.male_counts, sample, gender))
smith.gene_tpm.male_counts = read_tsv(here(shared_path,"public-data/ipsc-mn-tardbp-smith-2021/nfcore/star_salmon/salmon.merged.gene_tpm.tsv")) %>% filter(gene_name %in% c('RPS4Y1', 'EIFAY', 'DDX3Y', 'KDM5D')) %>% select(-gene_id) %>% adorn_totals("row", name = "male_counts") %>% as_tibble %>% 
  filter(gene_name == "male_counts") %>% transpose_tibble(old_rownames = "gene_name", new_rownames = "sample") %>% mutate(gender = case_when(male_counts > 10 ~ "male", TRUE ~ "female"))
smith.metadata <- read_csv(here(shared_path,"public-data/ipsc-mn-tardbp-smith-2021/sample-details/samplesheet.csv")) %>% distinct(sample, .keep_all = TRUE) %>%
  mutate(database_dir = here(shared_path,"public-data/ipsc-mn-tardbp-smith-2021"), dataset = "smith", DIV = 45) %>% left_join(select(smith.gene_tpm.male_counts, sample, gender))
bhinge.gene_tpm.male_counts = read_tsv(here(shared_path,"public-data/ipsc-mn-sod1-bhinge-2017/nfcore/star_salmon/salmon.merged.gene_tpm.tsv")) %>% filter(gene_name %in% c('RPS4Y1', 'EIFAY', 'DDX3Y', 'KDM5D')) %>% select(-gene_id) %>% adorn_totals("row", name = "male_counts") %>% as_tibble %>% 
  filter(gene_name == "male_counts") %>% transpose_tibble(old_rownames = "gene_name", new_rownames = "sample") %>% mutate(gender = case_when(male_counts > 10 ~ "male", TRUE ~ "female"))
bhinge.metadata <- read_csv(here(shared_path,"public-data/ipsc-mn-sod1-bhinge-2017/sample-details/samplesheet.csv")) %>% distinct(sample, .keep_all = TRUE) %>%
  mutate(database_dir = here(shared_path,"public-data/ipsc-mn-sod1-bhinge-2017"), dataset = "bhinge", DIV = 32) %>% select(-age) %>% left_join(select(bhinge.gene_tpm.male_counts, sample, gender))
hawkins.gene_tpm.male_counts = read_tsv(here(shared_path,"public-data/ipsc-mn-fus-hawkins-2022/nfcore/star_salmon/salmon.merged.gene_tpm.tsv")) %>% filter(gene_name %in% c('RPS4Y1', 'EIFAY', 'DDX3Y', 'KDM5D')) %>% select(-gene_id) %>% adorn_totals("row", name = "male_counts") %>% as_tibble %>% 
  filter(gene_name == "male_counts") %>% transpose_tibble(old_rownames = "gene_name", new_rownames = "sample") %>% mutate(gender = case_when(male_counts > 10 ~ "male", TRUE ~ "female"))
hawkins.metadata <- read_csv(here(shared_path,"public-data/ipsc-mn-fus-hawkins-2022/sample-details/samplesheet.csv")) %>% distinct(sample, .keep_all = TRUE) %>%
  mutate(database_dir = here(shared_path,"public-data/ipsc-mn-fus-hawkins-2022")) %>% left_join(select(hawkins.gene_tpm.male_counts, sample, gender))
lee.gene_tpm.male_counts = read_tsv(here(shared_path,"public-data/ipsc-occular-spinal-mn-lee-2021/nfcore/star_salmon/salmon.merged.gene_tpm.tsv")) %>% filter(gene_name %in% c('RPS4Y1', 'EIFAY', 'DDX3Y', 'KDM5D')) %>% select(-gene_id) %>% adorn_totals("row", name = "male_counts") %>% as_tibble %>% 
  filter(gene_name == "male_counts") %>% transpose_tibble(old_rownames = "gene_name", new_rownames = "sample") %>% mutate(gender = case_when(male_counts > 10 ~ "male", TRUE ~ "female"))
lee.metadata <- read_csv(here(shared_path,"public-data/ipsc-occular-spinal-mn-lee-2021/sample-details/samplesheet.csv")) %>% distinct(sample, .keep_all = TRUE) %>%
  mutate(database_dir = here(shared_path,"public-data/ipsc-occular-spinal-mn-lee-2021"), dataset = case_when(mutation == "ctrl" ~ "lee.CTRL", TRUE ~ "lee.ALS"), condition = case_when(mutation == "ctrl" ~ "ctrl", TRUE ~ "als"),
         strandedness = case_when(condition == "als" ~ "reverse", TRUE ~ strandedness)) %>% filter(location == "spinal") %>% left_join(select(lee.gene_tpm.male_counts, sample, gender))


# Metadata merge ----------------------------------------------------------------
ipsc_mn_als_datasets_with_lee.metadata <- bind_rows(answerals.metadata, neurolincs.metadata, catanese.metadata, dafinca.c9orf72.metadata, dafinca.tardbp.metadata, luisier.metadata, #mitchell.metadata,
                                           wang.metadata, kiskinis.metadata, sareen.metadata, sommer.metadata, sterneckert.metadata, kapeli.metadata, desantis.metadata, smith.metadata, bhinge.metadata, hawkins.metadata, lee.metadata) %>%
  mutate(dataset_sample = paste0(dataset,"_",sample), condition = factor(condition, levels = c("ctrl", "als")), dataset_simplified = case_when(dataset %in% c("dafinca.c9orf72","dafinca.tardbp") ~ "dafinca", TRUE ~ dataset), 
         dataset = as.factor(dataset), mutation = as.factor(mutation), DIV = as.factor(DIV), instrument = gsub("Illumina ","",instrument), 
         library_type = case_when(dataset %in% c("catanese","dafinca.c9orf72","dafinca.tardbp","kapeli","kiskinis","luisier","mitchell","sareen","sommer","smith","sterneckert","wang","lee.CTRL") ~ "poly(A)", TRUE ~ "Total")) %>%
  select(sample, dataset, condition, mutation, dataset_sample, replicate, database_dir, everything())
ipsc_mn_als_datasets.metadata <- ipsc_mn_als_datasets_with_lee.metadata %>% filter(!dataset %in% c("lee.CTRL","lee.ALS")) %>% mutate(dataset = fct_drop(dataset), dataset_simplified = fct_drop(dataset_simplified))
ipsc_mn_als_datasets.polyA.metadata <- filter(ipsc_mn_als_datasets.metadata, library_type == "poly(A)")
ipsc_mn_als_datasets.total.metadata <- filter(ipsc_mn_als_datasets.metadata, library_type == "Total")
ipsc_mn_als_datasets.single_end.metadata <- filter(ipsc_mn_als_datasets.metadata, library_layout == "SINGLE")
ipsc_mn_als_datasets.paired_end.metadata <- filter(ipsc_mn_als_datasets.metadata, library_layout == "PAIRED")
ipsc_mn_als_datasets.noiso.metadata <- ipsc_mn_als_datasets.metadata %>% filter(!dataset %in% c("bhinge","kiskinis","wang"), mutation != "iso") # remove 3 datasets with only isogenic corrected controls as well as other isogenic controls in mixed datasets - "dafinca.c9orf72", "catanese"
ipsc_mn_familial_datasets.metadata = ipsc_mn_als_datasets.metadata %>% filter(mutation != "sporadic")
ipsc_mn_sporadic_datasets.metadata <- bind_rows(filter(answerals.metadata, mutation %in% c("sporadic","ctrl")), filter(neurolincs.metadata, mutation %in% c("sporadic","ctrl"))) %>%
  mutate(dataset_sample = paste0(dataset,"_",sample), condition = factor(condition, levels = c("ctrl", "als")), dataset = as.factor(dataset), mutation = factor(mutation, levels = c("ctrl", "sporadic"))) %>%
  select(sample, dataset, condition, mutation, dataset_sample, database_dir, everything())

ipsc_mn_c9orf72_datasets.metadata <- bind_rows(filter(answerals.metadata, mutation %in% c("c9orf72","ctrl")), filter(neurolincs.metadata, mutation %in% c("c9orf72","ctrl")), filter(catanese.metadata, mutation %in% c("c9orf72","iso","ctrl")), 
                                               dafinca.c9orf72.metadata, sareen.metadata, sommer.metadata, sterneckert.metadata) %>%
  mutate(dataset_sample = paste0(dataset,"_",sample), condition = factor(condition, levels = c("ctrl", "als")), dataset = as.factor(dataset), mutation = factor(mutation, levels = c("ctrl", "iso", "c9orf72")),
         library_type = case_when(dataset %in% c("catanese","dafinca.c9orf72","dafinca.tardbp","kapeli","kiskinis","luisier","mitchell","sareen","sommer","smith","sterneckert","wang") ~ "poly(A)", TRUE ~ "Total")) %>%
  select(sample, dataset, condition, mutation, dataset_sample, replicate, database_dir, everything())

ipsc_mn_fus_datasets.metadata <- bind_rows(filter(answerals.metadata, mutation %in% c("fus","ctrl")), filter(catanese.metadata, mutation %in% c("fus","iso","ctrl")), kapeli.metadata, desantis.metadata, hawkins.metadata) %>%
  mutate(dataset_sample = paste0(dataset,"_",sample), condition = factor(condition, levels = c("ctrl", "als")), dataset = as.factor(dataset), mutation = factor(mutation, levels = c("ctrl", "iso", "fus")),
         library_type = case_when(dataset %in% c("catanese","dafinca.c9orf72","dafinca.tardbp","kapeli","kiskinis","luisier","mitchell","sareen","sommer","smith","sterneckert","wang") ~ "poly(A)", TRUE ~ "Total")) %>%
  select(sample, dataset, condition, mutation, dataset_sample, replicate, database_dir, everything())

ipsc_mn_sod1_datasets.metadata <- bind_rows(filter(answerals.metadata, mutation %in% c("sod1","ctrl")), filter(neurolincs.metadata, mutation %in% c("sod1","ctrl")), wang.metadata, kiskinis.metadata, bhinge.metadata) %>%
  mutate(dataset_sample = paste0(dataset,"_",sample), condition = factor(condition, levels = c("ctrl", "als")), dataset = as.factor(dataset), mutation = factor(mutation, levels = c("ctrl", "iso", "sod1")),
         library_type = case_when(dataset %in% c("catanese","dafinca.c9orf72","dafinca.tardbp","kapeli","kiskinis","luisier","mitchell","sareen","sommer","smith","sterneckert","wang") ~ "poly(A)", TRUE ~ "Total")) %>%
  select(sample, dataset, condition, mutation, dataset_sample, replicate, database_dir, everything())

ipsc_mn_tardbp_datasets.metadata <- bind_rows(filter(answerals.metadata, mutation %in% c("tardbp","ctrl")), dafinca.tardbp.metadata, smith.metadata) %>%
  mutate(dataset_sample = paste0(dataset,"_",sample), condition = factor(condition, levels = c("ctrl", "als")), dataset = as.factor(dataset), mutation = factor(mutation, levels = c("ctrl", "iso", "tardbp")),
         library_type = case_when(dataset %in% c("catanese","dafinca.c9orf72","dafinca.tardbp","kapeli","kiskinis","luisier","mitchell","sareen","sommer","smith","sterneckert","wang") ~ "poly(A)", TRUE ~ "Total")) %>%
  select(sample, dataset, condition, mutation, dataset_sample, replicate, database_dir, everything())

ipsc_mn_mutant_sporadic.metadata = ipsc_mn_als_datasets.metadata %>% filter(condition == "als") %>% mutate(mutant_sporadic = factor(case_when(mutation == "sporadic" ~ "sporadic", TRUE ~ "mutant"), levels = c("sporadic", "mutant")))

write_csv(select(ipsc_mn_als_datasets_with_lee.metadata, sample = dataset_sample, fastq_1, fastq_2, strandedness, dataset, condition, mutation, gender), here(proj_path,"sample-details/ipsc_mn_als_datasets_with_lee.samplesheet.csv"), na = "")
write_csv(select(ipsc_mn_als_datasets.metadata, sample = dataset_sample, fastq_1, fastq_2, strandedness, dataset, condition, mutation, gender), here(proj_path,"sample-details/ipsc_mn_als_datasets.samplesheet.csv"), na = "")
ipsc_mn_als_datasets.metadata.rnavar = ipsc_mn_als_datasets.metadata %>% mutate(dataset_sample = gsub("_",".",dataset_sample))
write_csv(select(ipsc_mn_als_datasets.metadata.rnavar, sample = dataset_sample, fastq_1, fastq_2, strandedness, dataset, condition, mutation, gender), here(proj_path,"sample-details/ipsc_mn_als_datasets.rnavar.samplesheet.csv"), na = "")
write_csv(select(ipsc_mn_als_datasets.paired_end.metadata, sample = dataset_sample, fastq_1, fastq_2, strandedness, dataset, condition, mutation, gender), here(proj_path,"sample-details/ipsc_mn_als_datasets.paired_end.samplesheet.csv"), na = "")

rm(answerals.gene_tpm.male_counts, neurolincs.gene_tpm.male_counts, catanese.gene_tpm.male_counts, dafinca.c9orf72.gene_tpm.male_counts, dafinca.tardbp.gene_tpm.male_counts,  luisier.gene_tpm.male_counts, lee.gene_tpm.male_counts,
   wang.gene_tpm.male_counts, kiskinis.gene_tpm.male_counts, sareen.gene_tpm.male_counts, sommer.gene_tpm.male_counts, sterneckert.gene_tpm.male_counts, kapeli.gene_tpm.male_counts, desantis.gene_tpm.male_counts, smith.gene_tpm.male_counts, bhinge.gene_tpm.male_counts, 
   hawkins.gene_tpm.male_counts) 

save.image(here(proj_path,"expression/deseq2/ipsc_mn_als_meta.metadata.RData"))

# DESeq2 ------------------------------------------------------------------

# ### ALS ------------------------------------------------------------------
# ipsc_mn_als_datasets_with_lee = DESeq.analysis(metadata = ipsc_mn_als_datasets_with_lee.metadata, unique_names = "dataset_sample")
# saveRDS(ipsc_mn_als_datasets_with_lee, here(proj_path,"expression/deseq2/ipsc_mn_als_datasets_with_lee.rds"))

# # estimate technical variation with RUVg method
# ipsc_mn_als_datasets_with_lee.metadata.salmon_files = ipsc_mn_als_datasets_with_lee.metadata %>% mutate(file_salmon = file.path(database_dir, "nfcore/star_salmon", ipsc_mn_als_datasets_with_lee.metadata[["sample"]], "quant.sf")) %>% pull(file_salmon)
# names(ipsc_mn_als_datasets_with_lee.metadata.salmon_files) = ipsc_mn_als_datasets_with_lee.metadata %>% pull(dataset_sample) # sample filename in nfcore outdir
# rownames(ipsc_mn_als_datasets_with_lee.metadata) <- ipsc_mn_als_datasets_with_lee.metadata %>% pull(dataset_sample)  # unique sample name
# ipsc_mn_als_datasets_with_lee.dds.ruv <- tximport(ipsc_mn_als_datasets_with_lee.metadata.salmon_files, type="salmon", tx2gene=tx2gene) %>% DESeqDataSetFromTximport(colData = ipsc_mn_als_datasets_with_lee.metadata, design = ~condition)
# ipsc_mn_als_datasets_with_lee.vsd.ruv.lrt <- vst(ipsc_mn_als_datasets_with_lee.dds.ruv, blind=FALSE)
# ipsc_mn_als_datasets_with_lee.dds.ruv.lrt = ipsc_mn_als_datasets_with_lee.dds.ruv %>% DESeq(test = "LRT", reduced = ~1) #, fitType="glmGamPoi")
# ipsc_mn_als_datasets_with_lee.res.ruv.lrt <- DESeq2::results(ipsc_mn_als_datasets_with_lee.dds.ruv.lrt)
# table(ipsc_mn_als_datasets_with_lee.res.ruv.lrt$padj < .1) # FALSE  TRUE # 13385 20117
# set <- newSeqExpressionSet(counts(ipsc_mn_als_datasets_with_lee.dds.ruv.lrt))
# set <- betweenLaneNormalization(set, which="upper")
# not_sig <- rownames(ipsc_mn_als_datasets_with_lee.res.ruv.lrt)[which(ipsc_mn_als_datasets_with_lee.res.ruv.lrt$pvalue > .1)]
# empirical <- rownames(set)[ rownames(set) %in% not_sig ]
# set <- RUVg(set, empirical, k=5)
# pdat <- pData(set) # %>% mutate(condition = case_when) # W_1 and W_2 etc are unwanted variation factors
# pdat$condition <- ipsc_mn_als_datasets_with_lee.metadata$condition # visualize how the factors of unwanted variation describe the samples in the PC1 and PC2 space:
# ipsc_mn_als_datasets_with_lee.vsd.ruv.lrt$W1 <- pdat$W_1
# ipsc_mn_als_datasets_with_lee.vsd.ruv.lrt$W2 <- pdat$W_2
# ipsc_mn_als_datasets_with_lee.vsd.ruv.lrt$W3 <- pdat$W_3
# ipsc_mn_als_datasets_with_lee.vsd.ruv.lrt$W4 <- pdat$W_4
# ipsc_mn_als_datasets_with_lee.vsd.ruv.lrt$W5 <- pdat$W_5
# saveRDS(ipsc_mn_als_datasets_with_lee.vsd.ruv.lrt, here(proj_path,"expression/deseq2/ipsc_mn_als_datasets_with_lee.vsd.ruv.lrt.rds"))

# saveRDS(ipsc_mn_als_datasets, here(proj_path,"expression/deseq2/ipsc_mn_als_datasets.rds"))
# ipsc_mn_als_datasets.noiso = DESeq.analysis(metadata = ipsc_mn_als_datasets.noiso.metadata, unique_names = "dataset_sample", design = ~ gender + dataset + condition, contrast = "condition_als_vs_ctrl")
# saveRDS(ipsc_mn_als_datasets.noiso, here(proj_path,"expression/deseq2/ipsc_mn_als_datasets.noiso.rds"))

# ipsc_mn_sporadic_datasets = DESeq.analysis(metadata = ipsc_mn_sporadic_datasets.metadata, unique_names = "dataset_sample", design = ~ gender + dataset + condition, contrast = "condition_als_vs_ctrl", species = "human")
# saveRDS(ipsc_mn_sporadic_datasets, here(proj_path,"expression/deseq2/ipsc_mn_sporadic_datasets.rds"))
# ipsc_mn_c9orf72_datasets = DESeq.analysis(metadata = ipsc_mn_c9orf72_datasets.metadata, unique_names = "dataset_sample", design = ~ gender + dataset + condition, contrast = "condition_als_vs_ctrl", species = "human")
# saveRDS(ipsc_mn_c9orf72_datasets, here(proj_path,"expression/deseq2/ipsc_mn_c9orf72_datasets.rds"))
# ipsc_mn_fus_datasets = DESeq.analysis(metadata = ipsc_mn_fus_datasets.metadata, unique_names = "dataset_sample", design = ~ gender + dataset + condition, contrast = "condition_als_vs_ctrl", species = "human")
# saveRDS(ipsc_mn_fus_datasets, here(proj_path,"expression/deseq2/ipsc_mn_fus_datasets.rds"))
# ipsc_mn_sod1_datasets = DESeq.analysis(metadata = ipsc_mn_sod1_datasets.metadata, unique_names = "dataset_sample", design = ~ gender + dataset + condition, contrast = "condition_als_vs_ctrl", species = "human")
# saveRDS(ipsc_mn_sod1_datasets, here(proj_path,"expression/deseq2/ipsc_mn_sod1_datasets.rds"))
# ipsc_mn_tardbp_datasets = DESeq.analysis(metadata = ipsc_mn_tardbp_datasets.metadata, unique_names = "dataset_sample", design = ~ gender + dataset + condition, contrast = "condition_als_vs_ctrl", species = "human")
# saveRDS(ipsc_mn_tardbp_datasets, here(proj_path,"expression/deseq2/ipsc_mn_tardbp_datasets.rds"))

# ### ALS mutation groupings ------------------------------------------------------------------
# ipsc_mn_familial_datasets = DESeq.analysis(metadata = ipsc_mn_familial_datasets.metadata, unique_names = "dataset_sample", design = ~ gender + dataset + condition, contrast = "condition_als_vs_ctrl", species = "human")
# saveRDS(ipsc_mn_familial_datasets, here(proj_path,"expression/deseq2/ipsc_mn_familial_datasets.rds"))
# ipsc_mn_mutant_sporadic = DESeq.analysis(metadata = ipsc_mn_mutant_sporadic.metadata, unique_names = "dataset_sample", design = ~ gender +  dataset + mutant_sporadic, contrast = "mutant_sporadic_mutant_vs_sporadic", species = "human")
# saveRDS(ipsc_mn_mutant_sporadic, here(proj_path,"expression/deseq2/ipsc_mn_mutant_sporadic.rds"))

# ### Library type polyA vs Total ------------------------------------------------------------------
# ipsc_mn_als_datasets.polyA = DESeq.analysis(metadata = ipsc_mn_als_datasets.polyA.metadata, unique_names = "dataset_sample", design = ~ gender + dataset + condition, contrast = "condition_als_vs_ctrl", run_irfinder = TRUE, irfinder_design = ~ gender + dataset + condition + condition:IRFinder)
# saveRDS(ipsc_mn_als_datasets.polyA, here(proj_path,"expression/deseq2/ipsc_mn_als_datasets.polyA.rds"))
# ipsc_mn_als_datasets.total = DESeq.analysis(metadata = ipsc_mn_als_datasets.total.metadata, unique_names = "dataset_sample", design = ~ gender + dataset + condition, contrast = "condition_als_vs_ctrl", run_irfinder = TRUE, irfinder_design = ~ gender + dataset + condition + condition:IRFinder)
# saveRDS(ipsc_mn_als_datasets.total, here(proj_path,"expression/deseq2/ipsc_mn_als_datasets.total.rds"))
# ipsc_mn_c9orf72_datasets.polyA = DESeq.analysis(metadata = filter(ipsc_mn_c9orf72_datasets.metadata, library_type == "poly(A)"), unique_names = "dataset_sample", design = ~ gender + dataset + condition, contrast = "condition_als_vs_ctrl", run_irfinder = TRUE, irfinder_design = ~ gender + dataset + condition + condition:IRFinder)
# saveRDS(ipsc_mn_c9orf72_datasets.polyA, here(proj_path,"expression/deseq2/ipsc_mn_c9orf72_datasets.polyA.rds"))
# ipsc_mn_c9orf72_datasets.total = DESeq.analysis(metadata = filter(ipsc_mn_c9orf72_datasets.metadata, library_type == "Total"), unique_names = "dataset_sample", design = ~ gender + dataset + condition, contrast = "condition_als_vs_ctrl", run_irfinder = TRUE, irfinder_design = ~ gender + dataset + condition + condition:IRFinder)
# saveRDS(ipsc_mn_c9orf72_datasets.total, here(proj_path,"expression/deseq2/ipsc_mn_c9orf72_datasets.total.rds"))
# ipsc_mn_fus_datasets.polyA = DESeq.analysis(metadata = filter(ipsc_mn_fus_datasets.metadata, library_type == "poly(A)"), unique_names = "dataset_sample", design = ~ gender + dataset + condition, contrast = "condition_als_vs_ctrl", run_irfinder = TRUE, irfinder_design = ~ gender + dataset + condition + condition:IRFinder)
# saveRDS(ipsc_mn_fus_datasets.polyA, here(proj_path,"expression/deseq2/ipsc_mn_fus_datasets.polyA.rds"))
# ipsc_mn_fus_datasets.total = DESeq.analysis(metadata = filter(ipsc_mn_fus_datasets.metadata, library_type == "Total"), unique_names = "dataset_sample", design = ~ gender + dataset + condition, contrast = "condition_als_vs_ctrl", run_irfinder = TRUE, irfinder_design = ~ gender + dataset + condition + condition:IRFinder)
# saveRDS(ipsc_mn_fus_datasets.total, here(proj_path,"expression/deseq2/ipsc_mn_fus_datasets.total.rds"))
# ipsc_mn_sod1_datasets.polyA = DESeq.analysis(metadata = filter(ipsc_mn_sod1_datasets.metadata, library_type == "poly(A)"), unique_names = "dataset_sample", design = ~ dataset + condition, contrast = "condition_als_vs_ctrl", run_irfinder = TRUE, irfinder_design = ~ dataset + condition + condition:IRFinder)
# saveRDS(ipsc_mn_sod1_datasets.polyA, here(proj_path,"expression/deseq2/ipsc_mn_sod1_datasets.polyA.rds"))
# ipsc_mn_sod1_datasets.total = DESeq.analysis(metadata = filter(ipsc_mn_sod1_datasets.metadata, library_type == "Total"), unique_names = "dataset_sample", design = ~ gender + dataset + condition, contrast = "condition_als_vs_ctrl", run_irfinder = TRUE, irfinder_design = ~ gender + dataset + condition + condition:IRFinder)
# saveRDS(ipsc_mn_sod1_datasets.total, here(proj_path,"expression/deseq2/ipsc_mn_sod1_datasets.total.rds"))
# ipsc_mn_tardbp_datasets.polyA = DESeq.analysis(metadata = filter(ipsc_mn_tardbp_datasets.metadata, library_type == "poly(A)"), unique_names = "dataset_sample", design = ~ gender + dataset + condition, contrast = "condition_als_vs_ctrl", run_irfinder = TRUE, irfinder_design = ~ gender + dataset + condition + condition:IRFinder)
# saveRDS(ipsc_mn_tardbp_datasets.polyA, here(proj_path,"expression/deseq2/ipsc_mn_tardbp_datasets.polyA.rds"))
# ipsc_mn_tardbp_datasets.total = DESeq.analysis(metadata = filter(ipsc_mn_tardbp_datasets.metadata, library_type == "Total"), unique_names = "dataset_sample", design = ~ gender + condition, contrast = "condition_als_vs_ctrl", run_irfinder = TRUE, irfinder_design = ~ gender + condition + condition:IRFinder)
# saveRDS(ipsc_mn_tardbp_datasets.total, here(proj_path,"expression/deseq2/ipsc_mn_tardbp_datasets.total.rds"))

# # ### NYGC Post Mortem ------------------------------------------------------------------
# # spinal cord counts from Jack Humphrey https://zenodo.org/record/6385747#.Yod3Y5PMJqs
# nygc_cervical_cord_rsem_counts = read_tsv(here(shared_path,"nygc-als-consortium/counts/spinal_cord_jack_humphrey/Cervical_Spinal_Cord/Cervical_Spinal_Cord_gene_counts.tsv")) # 14
# nygc_thoracic_cord_rsem_counts = read_tsv(here(shared_path,"nygc-als-consortium/counts/spinal_cord_jack_humphrey/Thoracic_Spinal_Cord/Thoracic_Spinal_Cord_gene_counts.tsv")) # 52
# nygc_lumbar_cord_rsem_counts = read_tsv(here(shared_path,"nygc-als-consortium/counts/spinal_cord_jack_humphrey/Lumbar_Spinal_Cord/Lumbar_Spinal_Cord_gene_counts.tsv")) # 154
# nygc_spinal_cord_rsem_counts.join = full_join(nygc_cervical_cord_rsem_counts, nygc_thoracic_cord_rsem_counts, by = c("ensembl_id","gene_name")) %>% full_join(nygc_lumbar_cord_rsem_counts, by = c("ensembl_id","gene_name")) # 228
# 
# nygc_cervical_cord_metadata = read_tsv(here(shared_path,"nygc-als-consortium/counts/spinal_cord_jack_humphrey/Cervical_Spinal_Cord/Cervical_Spinal_Cord_metadata.tsv")) # 174
# nygc_thoracic_cord_metadata = read_tsv(here(shared_path,"nygc-als-consortium/counts/spinal_cord_jack_humphrey/Thoracic_Spinal_Cord/Thoracic_Spinal_Cord_metadata.tsv")) # 52
# nygc_lumbar_cord_metadata = read_tsv(here(shared_path,"nygc-als-consortium/counts/spinal_cord_jack_humphrey/Lumbar_Spinal_Cord/Lumbar_Spinal_Cord_metadata.tsv")) # 154
# nygc_spinal_cord_metadata.bind = bind_rows(nygc_cervical_cord_metadata, nygc_thoracic_cord_metadata, nygc_lumbar_cord_metadata) %>% # 380
#   mutate(condition = case_when(disease == "Control" ~ "ctrl", grepl("ALS",disease) ~ "als"), condition = factor(condition, levels = c("ctrl","als")), library_prep = as.factor(library_prep),
#          gender = tolower(sex), mutation = tolower(mutations), mutation = case_when(condition == "als" & mutation == "none" ~ "sporadic", TRUE ~ mutation))
# nygc_spinal_cord_metadata.bind %>% count(mutation)
# multiple_donor_dna_id = nygc_spinal_cord_metadata.bind %>% get_dupes(dna_id) %>% pull(dna_id) # 325 individuals have multiple spinal cord samples. For these we take only their Cervical sample
# nygc_spinal_cord_metadata = nygc_spinal_cord_metadata.bind %>% mutate(tissue = factor(tissue, levels = c("Cervical_Spinal_Cord", "Thoracic_Spinal_Cord", "Lumbar_Spinal_Cord")), sample = rna_id) %>% arrange(tissue) %>% distinct(dna_id, .keep_all = TRUE) %>% column_to_rownames(var = "rna_id") # 203 / 380
# nygc_spinal_cord_rsem_counts = nygc_spinal_cord_rsem_counts.join %>% select(gene_id = ensembl_id, gene_name, all_of(nygc_spinal_cord_metadata$sample)) # 205 / 228
# nygc_spinal_cord_rsem_counts.mat =  nygc_spinal_cord_rsem_counts %>% select(gene_id, all_of(rownames(nygc_spinal_cord_metadata))) %>% column_to_rownames(var = "gene_id") %>% as.matrix() %>% round() # keep counts only for samples in metadata
# 
# nygc_postmortem_spinal_cord.als_vs_ctrl = list()
# nygc_postmortem_spinal_cord.als_vs_ctrl$dds <- DESeqDataSetFromMatrix(countData = nygc_spinal_cord_rsem_counts.mat, colData = nygc_spinal_cord_metadata, design = ~ gender + library_prep + tissue + condition) %>% DESeq() # adjust for spinal cord site
# nygc_postmortem_spinal_cord.als_vs_ctrl$vsd <- vst(nygc_postmortem_spinal_cord.als_vs_ctrl$dds, blind=TRUE)
# nygc_postmortem_spinal_cord.als_vs_ctrl$vsd.counts <- as_tibble(assay(nygc_postmortem_spinal_cord.als_vs_ctrl$vsd), rownames = "gene_id")  %>% left_join(gene2ens)
# nygc_postmortem_spinal_cord.als_vs_ctrl$res <- DESeq2::results(nygc_postmortem_spinal_cord.als_vs_ctrl$dds, name = "condition_als_vs_ctrl") %>% as_tibble(rownames = "gene_id") %>% left_join(gene2ens) %>% arrange(-abs(stat))
# saveRDS(nygc_postmortem_spinal_cord.als_vs_ctrl, here(proj_path,"expression/deseq2/nygc_postmortem_spinal_cord.als_vs_ctrl.rds"))
# 
# # NYGC mutation metadata
# nygc_spinal_cord_metadata %>% count(mutation) %>% arrange(-n)
# # mutation   n
# # 1 sporadic 115
# # 2     none  50
# # 3  c9orf72  29
# # 4     sod1   5
# # 5      fus   2
# # 6      ang   1
# # 7     optn   1
# 
# # Sporadic
# nygc_spinal_cord_metadata.sporadic = nygc_spinal_cord_metadata %>% filter(mutation %in% c("sporadic","none"))
# nygc_spinal_cord_rsem_counts.sporadic.mat =  nygc_spinal_cord_rsem_counts %>% select(gene_id, all_of(rownames(nygc_spinal_cord_metadata.sporadic))) %>% column_to_rownames(var = "gene_id") %>% as.matrix() %>% round()
# nygc_postmortem_spinal_cord.sporadic_vs_ctrl = list()
# nygc_postmortem_spinal_cord.sporadic_vs_ctrl$dds <- DESeqDataSetFromMatrix(countData = nygc_spinal_cord_rsem_counts.sporadic.mat, colData = nygc_spinal_cord_metadata.sporadic, design = ~ gender + library_prep + tissue + condition) %>% DESeq()
# nygc_postmortem_spinal_cord.sporadic_vs_ctrl$res <- DESeq2::results(nygc_postmortem_spinal_cord.sporadic_vs_ctrl$dds, name = "condition_als_vs_ctrl") %>% as_tibble(rownames = "gene_id") %>% left_join(gene2ens) %>% arrange(-abs(stat))
# saveRDS(nygc_postmortem_spinal_cord.sporadic_vs_ctrl, here(proj_path,"expression/deseq2/nygc_postmortem_spinal_cord.sporadic_vs_ctrl.rds"))
# 
# # C9orf72
# nygc_spinal_cord_metadata.c9orf72 = nygc_spinal_cord_metadata %>% filter(mutation %in% c("c9orf72","none"))
# nygc_spinal_cord_rsem_counts.c9orf72.mat =  nygc_spinal_cord_rsem_counts %>% select(gene_id, all_of(rownames(nygc_spinal_cord_metadata.c9orf72))) %>% column_to_rownames(var = "gene_id") %>% as.matrix() %>% round()
# nygc_postmortem_spinal_cord.c9orf72_vs_ctrl = list()
# nygc_postmortem_spinal_cord.c9orf72_vs_ctrl$dds <- DESeqDataSetFromMatrix(countData = nygc_spinal_cord_rsem_counts.c9orf72.mat, colData = nygc_spinal_cord_metadata.c9orf72, design = ~ gender + library_prep + tissue + condition) %>% DESeq()
# nygc_postmortem_spinal_cord.c9orf72_vs_ctrl$res <- DESeq2::results(nygc_postmortem_spinal_cord.c9orf72_vs_ctrl$dds, name = "condition_als_vs_ctrl") %>% as_tibble(rownames = "gene_id") %>% left_join(gene2ens) %>% arrange(-abs(stat))
# saveRDS(nygc_postmortem_spinal_cord.c9orf72_vs_ctrl, here(proj_path,"expression/deseq2/nygc_postmortem_spinal_cord.c9orf72_vs_ctrl.rds"))
# 
# # SOD1
# nygc_spinal_cord_metadata.sod1 = nygc_spinal_cord_metadata %>% filter(mutation %in% c("sod1","none"))
# nygc_spinal_cord_metadata.sod1 %>% count(condition,mutation) # 29 ALS, 50 CTRL biopsy samples
# nygc_spinal_cord_rsem_counts.sod1.mat =  nygc_spinal_cord_rsem_counts %>% select(gene_id, all_of(rownames(nygc_spinal_cord_metadata.sod1))) %>% column_to_rownames(var = "gene_id") %>% as.matrix() %>% round()
# nygc_postmortem_spinal_cord.sod1_vs_ctrl = list()
# nygc_postmortem_spinal_cord.sod1_vs_ctrl$dds <- DESeqDataSetFromMatrix(countData = nygc_spinal_cord_rsem_counts.sod1.mat, colData = nygc_spinal_cord_metadata.sod1, design = ~ gender + library_prep + tissue + condition) %>% DESeq()
# nygc_postmortem_spinal_cord.sod1_vs_ctrl$res <- DESeq2::results(nygc_postmortem_spinal_cord.sod1_vs_ctrl$dds, name = "condition_als_vs_ctrl") %>% as_tibble(rownames = "gene_id") %>% left_join(gene2ens) %>% arrange(-abs(stat))
# saveRDS(nygc_postmortem_spinal_cord.sod1_vs_ctrl, here(proj_path,"expression/deseq2/nygc_postmortem_spinal_cord.sod1_vs_ctrl.rds"))
# 
# # FUS
# nygc_spinal_cord_metadata.fus = nygc_spinal_cord_metadata %>% filter(mutation %in% c("fus","none"))
# nygc_spinal_cord_rsem_counts.fus.mat =  nygc_spinal_cord_rsem_counts %>% select(gene_id, all_of(rownames(nygc_spinal_cord_metadata.fus))) %>% column_to_rownames(var = "gene_id") %>% as.matrix() %>% round()
# nygc_postmortem_spinal_cord.fus_vs_ctrl = list()
# nygc_postmortem_spinal_cord.fus_vs_ctrl$dds <- DESeqDataSetFromMatrix(countData = nygc_spinal_cord_rsem_counts.fus.mat, colData = nygc_spinal_cord_metadata.fus, design = ~ gender + library_prep + tissue + condition) %>% DESeq()
# nygc_postmortem_spinal_cord.fus_vs_ctrl$res <- DESeq2::results(nygc_postmortem_spinal_cord.fus_vs_ctrl$dds, name = "condition_als_vs_ctrl") %>% as_tibble(rownames = "gene_id") %>% left_join(gene2ens) %>% arrange(-abs(stat))
# saveRDS(nygc_postmortem_spinal_cord.fus_vs_ctrl, here(proj_path,"expression/deseq2/nygc_postmortem_spinal_cord.fus_vs_ctrl.rds"))


# ## Knockdown datasets ------------------------------------------------------------------
# tdp43_kd_brown.metadata <- read_csv(here(shared_path,"public-data/ipsc-mn-tardbp-kd-brown-2022/sample-details/samplesheet.csv")) %>% distinct(sample, .keep_all = TRUE) %>%
#   mutate(database_dir = here(shared_path,"public-data/ipsc-mn-tardbp-kd-brown-2022"), condition = factor(condition, levels = c("ctrl", "tdp43kd")), dataset = "brown")
# tdp43_kd_klim.metadata <- read_csv(here(shared_path,"public-data/ipsc-mn-tardbp-kd-klim-2019/sample-details/samplesheet.csv")) %>% distinct(sample, .keep_all = TRUE) %>%
#   mutate(database_dir = here(shared_path,"public-data/ipsc-mn-tardbp-kd-klim-2019"), condition = factor(condition, levels = c("ctrl", "tdp43kd")), dataset = "klim")
# tdp43_kd_appocher.metadata <- read_csv(here(shared_path,"public-data/shsy5y-tardbp-kd-appocher-2017/sample-details/samplesheet.csv")) %>% distinct(sample, .keep_all = TRUE) %>%
#   mutate(database_dir = here(shared_path,"public-data/shsy5y-tardbp-kd-appocher-2017"), condition = factor(condition, levels = c("ctrl", "tdp43kd")), dataset = "appocher")
# neuronal_nuclei_tdp43_liu.metadata <- read_csv(here(shared_path,"public-data/neuronal-nuclei-tdp43-liu-2019/sample-details/samplesheet.csv")) %>% distinct(sample, .keep_all = TRUE) %>%
#   mutate(database_dir = here(shared_path,"public-data/neuronal-nuclei-tdp43-liu-2019"), group = condition, condition = case_when(condition == "tdp43pos" ~ "ctrl", condition == "tdp43neg" ~ "tdp43kd"), condition = factor(condition, levels = c("ctrl", "tdp43kd")), dataset = "liu")
# kd_tyzack.metadata <- read_csv(here(shared_path,"patani-data/motor-neuron-rbp-knockdown-rnaseq/sample-details/samplesheet.csv")) %>% filter(fraction == "whole", group != "untreated") %>%
#   mutate(database_dir = here(shared_path,"patani-data/motor-neuron-rbp-knockdown-rnaseq"), dataset = "tyzack", condition = group) %>%
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
# 
# kd_tyzack = DESeq.analysis(metadata = kd_tyzack.metadata, unique_names = "sample", design = ~ 1)
# saveRDS(kd_tyzack, here(proj_path,"expression/deseq2/kd_tyzack.rds"))
# tdp43_kd_tyzack = DESeq.analysis(metadata = tdp43_kd_tyzack.metadata, unique_names = "sample", design = ~ cellline + condition, contrast = "condition_tdp43kd_vs_ctrl", run_irfinder = TRUE, irfinder_design = ~ cellline + condition + condition:IRFinder)
# saveRDS(tdp43_kd_tyzack, here(proj_path,"expression/deseq2/tdp43_kd_tyzack.rds"))
# fus_kd_tyzack = DESeq.analysis(metadata = fus_kd_tyzack.metadata, unique_names = "sample", design = ~ cellline + condition, contrast = "condition_fuskd_vs_ctrl", run_irfinder = TRUE, irfinder_design = ~ cellline + condition + condition:IRFinder)
# saveRDS(fus_kd_tyzack, here(proj_path,"expression/deseq2/fus_kd_tyzack.rds"))
# sfpq_kd_tyzack = DESeq.analysis(metadata = sfpq_kd_tyzack.metadata, unique_names = "sample", design = ~ cellline + condition, contrast = "condition_sfpqkd_vs_ctrl", run_irfinder = TRUE, irfinder_design = ~ cellline + condition + condition:IRFinder)
# saveRDS(sfpq_kd_tyzack, here(proj_path,"expression/deseq2/sfpq_kd_tyzack.rds"))
# tdp43_kd_brown = DESeq.analysis(metadata = tdp43_kd_brown.metadata, unique_names = "sample", design = ~ condition, contrast = "condition_tdp43kd_vs_ctrl", run_irfinder = TRUE, irfinder_design = ~ condition + condition:IRFinder)
# saveRDS(tdp43_kd_brown, here(proj_path,"expression/deseq2/tdp43_kd_brown.rds"))
# tdp43_kd_klim = DESeq.analysis(metadata = tdp43_kd_klim.metadata, unique_names = "sample", design = ~ condition, contrast = "condition_tdp43kd_vs_ctrl", run_irfinder = TRUE, irfinder_design = ~ condition + condition:IRFinder)
# saveRDS(tdp43_kd_klim, here(proj_path,"expression/deseq2/tdp43_kd_klim.rds"))
# neuronal_nuclei_tdp43_liu = DESeq.analysis(metadata = neuronal_nuclei_tdp43_liu.metadata, unique_names = "sample", design = ~ condition, contrast = "condition_tdp43kd_vs_ctrl", run_irfinder = TRUE, irfinder_design = ~ condition + condition:IRFinder)
# saveRDS(neuronal_nuclei_tdp43_liu, here(proj_path,"expression/deseq2/neuronal_nuclei_tdp43_liu.rds"))
# tdp43_kd_appocher = DESeq.analysis(metadata = tdp43_kd_appocher.metadata, unique_names = "sample", design = ~ condition, contrast = "condition_tdp43kd_vs_ctrl", run_irfinder = TRUE, irfinder_design = ~ condition + condition:IRFinder)
# saveRDS(tdp43_kd_appocher, here(proj_path,"expression/deseq2/tdp43_kd_appocher.rds"))
# tdp43_kd_datasets_notyzack = DESeq.analysis(metadata = tdp43_kd_datasets_notyzack.metadata, unique_names = "dataset_sample", design = ~ dataset + condition, contrast = "condition_tdp43kd_vs_ctrl", irfinder_design = ~ dataset + condition + condition:IRFinder)
# saveRDS(tdp43_kd_datasets_notyzack, here(proj_path,"expression/deseq2/tdp43_kd_datasets_notyzack.rds"))
# tdp43_kd_datasets = DESeq.analysis(metadata = tdp43_kd_datasets.metadata, unique_names = "dataset_sample", design = ~ dataset + condition, contrast = "condition_tdp43kd_vs_ctrl", run_irfinder = TRUE, irfinder_design = ~ dataset + condition + condition:IRFinder)
# saveRDS(tdp43_kd_datasets, here(proj_path,"expression/deseq2/tdp43_kd_datasets.rds"))


# print("saving")
# save.image(here(proj_path,"expression/deseq2/ipsc_mn_als_meta.RData"))
print("complete")

