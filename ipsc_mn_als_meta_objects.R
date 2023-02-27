# run this script in Rstudio server OnDemand Terminal Tab (not with conda env r4.0.3) with:
# cd /camp/home/ziffo/home/projects/ipsc-mn-als-meta/scripts
# sbatch -N 1 -c 12 --mem=0 -t 72:00:00 --wrap="Rscript /camp/lab/luscomben/home/users/ziffo/projects/ipsc-mn-als-meta/scripts/ipsc_mn_als_meta/ipsc_mn_als_meta_objects.R" --mail-type=ALL,ARRAY_TASKS --mail-user=oliver.ziff@crick.ac.uk --job-name=DESeq2objects
# sbatch -N 1 -c 12 --mem=72G -t 72:00:00 --wrap="Rscript /camp/lab/luscomben/home/users/ziffo/projects/ipsc-mn-als-meta/scripts/ipsc_mn_als_meta/ipsc_mn_als_meta_objects.R" --mail-type=ALL,ARRAY_TASKS --mail-user=oliver.ziff@crick.ac.uk --job-name=DESeq2objects --dependency=afterany:44191817,44191802
# sbatch -N 1 -c 6 --mem=50G -t 8:00:00 --wrap="Rscript /camp/lab/luscomben/home/users/ziffo/projects/ipsc-mn-als-meta/scripts/ipsc_mn_als_meta/ipsc_mn_als_meta_objects.R" --mail-type=ALL,ARRAY_TASKS --mail-user=oliver.ziff@crick.ac.uk --job-name=DESeq2objects 
# sbatch --cpus-per-task=10 --mem=1500G -N 1 --partition=hmem -t 72:00:00 --wrap="Rscript /camp/lab/luscomben/home/users/ziffo/projects/ipsc-mn-als-meta/scripts/ipsc_mn_als_meta/ipsc_mn_als_meta_objects.R" --mail-type=ALL,ARRAY_TASKS --mail-user=oliver.ziff@crick.ac.uk --job-name=DESeq2objects

library(here)
# camp_path = here("/Volumes/lab-luscomben/home/users/ziffo")
camp_path = here("/camp/lab/luscomben/home/users/ziffo")
# shared_path = here("/Volumes/lab-luscomben/home/users/ziffo/proj-luscombn-patani/working")
shared_path = here("/camp/project/proj-luscombn-patani/working")
# collab_path = here("/Volumes/lab-luscomben/home/users/ziffo/patani-collab")
collab_path = here("/camp/lab/luscomben/home/users/ziffo/patani-collab")
# proj_path = here("/Volumes/lab-luscomben/home/users/ziffo/projects/ipsc-mn-als-meta")
proj_path = here("/camp/lab/luscomben/home/users/ziffo/projects/ipsc-mn-als-meta")
# load(here(proj_path,"expression/deseq2/ipsc_mn_als_meta.metadata.RData")) # overwrites with old functions
source(here(camp_path,"scripts/functions/R_workspace.R"))
# load(here(proj_path,"scripts/untreated_deseq_objects.RData"))
# load(here(proj_path,"expression/deseq2/ipsc_mn_als_meta.metadata.RData")) # load metadata before functions to avoid overwriting with old functions

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
  mutate(database_dir = here(shared_path,"neurolincs"), condition = factor(condition, levels = c("ctrl", "als")), dataset = case_when(dataset == "neurolincsI" ~ "neurolincs.iPSC"), study_accession = dataset, sample2 = gsub("\\.","-",sample)) %>% left_join(select(neurolincs.gene_tpm.male_counts, sample2=sample, gender))
catanese.gene_tpm.male_counts = read_tsv(here(shared_path,"public-data/ipsc-mn-c9orf72-fus-catanese-2021/nfcore/star_salmon/salmon.merged.gene_tpm.tsv")) %>% filter(gene_name %in% c('RPS4Y1', 'EIFAY', 'DDX3Y', 'KDM5D')) %>% select(-gene_id) %>% adorn_totals("row", name = "male_counts") %>% as_tibble %>% 
  filter(gene_name == "male_counts") %>% transpose_tibble(old_rownames = "gene_name", new_rownames = "sample") %>% mutate(gender = case_when(male_counts > 10 ~ "male", TRUE ~ "female"))
catanese.metadata <- read_csv(here(shared_path,"public-data/ipsc-mn-c9orf72-fus-catanese-2021/sample-details/samplesheet.csv")) %>% distinct(sample, .keep_all = TRUE) %>%
  mutate(database_dir = here(shared_path,"public-data/ipsc-mn-c9orf72-fus-catanese-2021"), dataset = "catanese", DIV = 35, condition = factor(condition, levels = c("ctrl", "als"))) %>% select(-age) %>% left_join(select(catanese.gene_tpm.male_counts, sample, gender))
dafinca.c9orf72.gene_tpm.male_counts = read_tsv(here(shared_path,"public-data/ipsc-mn-c9orf72-dafinca-2020/nfcore/star_salmon/salmon.merged.gene_tpm.tsv")) %>% filter(gene_name %in% c('RPS4Y1', 'EIFAY', 'DDX3Y', 'KDM5D')) %>% select(-gene_id) %>% adorn_totals("row", name = "male_counts") %>% as_tibble %>% 
  filter(gene_name == "male_counts") %>% transpose_tibble(old_rownames = "gene_name", new_rownames = "sample") %>% mutate(gender = case_when(male_counts > 10 ~ "male", TRUE ~ "female"))
dafinca.c9orf72.metadata <- read_csv(here(shared_path,"public-data/ipsc-mn-c9orf72-dafinca-2020/sample-details/samplesheet.csv")) %>% distinct(sample, .keep_all = TRUE) %>%
  mutate(database_dir = here(shared_path,"public-data/ipsc-mn-c9orf72-dafinca-2020"), dataset = "dafinca.c9orf72", DIV = 19, condition = factor(condition, levels = c("ctrl", "als"))) %>% left_join(select(dafinca.c9orf72.gene_tpm.male_counts, sample, gender))
dafinca.tardbp.gene_tpm.male_counts = read_tsv(here(shared_path,"public-data/ipsc-mn-tardbp-dafinca-2020/nfcore/star_salmon/salmon.merged.gene_tpm.tsv")) %>% filter(gene_name %in% c('RPS4Y1', 'EIFAY', 'DDX3Y', 'KDM5D')) %>% select(-gene_id) %>% adorn_totals("row", name = "male_counts") %>% as_tibble %>% 
  filter(gene_name == "male_counts") %>% transpose_tibble(old_rownames = "gene_name", new_rownames = "sample") %>% mutate(gender = case_when(male_counts > 10 ~ "male", TRUE ~ "female"))
dafinca.tardbp.metadata <- read_csv(here(shared_path,"public-data/ipsc-mn-tardbp-dafinca-2020/sample-details/samplesheet.csv")) %>% distinct(sample, .keep_all = TRUE) %>%
  mutate(database_dir = here(shared_path,"public-data/ipsc-mn-tardbp-dafinca-2020"), dataset = "dafinca.tardbp", DIV = 19, condition = factor(condition, levels = c("ctrl", "als"))) %>% left_join(select(dafinca.tardbp.gene_tpm.male_counts, sample, gender))
luisier.gene_tpm.male_counts = read_tsv(here(collab_path,"motor-neuron-vcp-luisier-2018/nfcore/star_salmon/salmon.merged.gene_tpm.tsv")) %>% filter(gene_name %in% c('RPS4Y1', 'EIFAY', 'DDX3Y', 'KDM5D')) %>% select(-gene_id) %>% adorn_totals("row", name = "male_counts") %>% as_tibble %>% 
  filter(gene_name == "male_counts") %>% transpose_tibble(old_rownames = "gene_name", new_rownames = "sample") %>% mutate(gender = case_when(male_counts > 10 ~ "male", TRUE ~ "female"))
luisier.metadata <- read_csv(here(collab_path,"motor-neuron-vcp-luisier-2018/sample-details/samplesheet.csv")) %>% distinct(sample, .keep_all = TRUE) %>%  filter(day == 35) %>%
  mutate(database_dir = here(collab_path,"motor-neuron-vcp-luisier-2018"), dataset = "luisier", total_size = bases, name = gsub("d35_","",sample), DIV = 35, condition = factor(condition, levels = c("ctrl", "als"))) %>% left_join(select(luisier.gene_tpm.male_counts, sample, gender))
wang.gene_tpm.male_counts = read_tsv(here(shared_path,"public-data/ipsc-mn-sod1-wang-2017/nfcore/star_salmon/salmon.merged.gene_tpm.tsv")) %>% filter(gene_name %in% c('RPS4Y1', 'EIFAY', 'DDX3Y', 'KDM5D')) %>% select(-gene_id) %>% adorn_totals("row", name = "male_counts") %>% as_tibble %>% 
  filter(gene_name == "male_counts") %>% transpose_tibble(old_rownames = "gene_name", new_rownames = "sample") %>% mutate(gender = case_when(male_counts > 10 ~ "male", TRUE ~ "female"))
wang.metadata <- read_csv(here(shared_path,"public-data/ipsc-mn-sod1-wang-2017/sample-details/samplesheet.csv")) %>% distinct(sample, .keep_all = TRUE) %>%
  mutate(database_dir = here(shared_path,"public-data/ipsc-mn-sod1-wang-2017"), dataset = "wang", DIV = 12, condition = factor(condition, levels = c("ctrl", "als"))) %>% left_join(select(wang.gene_tpm.male_counts, sample, gender))
kiskinis.gene_tpm.male_counts = read_tsv(here(shared_path,"public-data/ipsc-mn-sod1-moccia-2014/nfcore/star_salmon/salmon.merged.gene_tpm.tsv")) %>% filter(gene_name %in% c('RPS4Y1', 'EIFAY', 'DDX3Y', 'KDM5D')) %>% select(-gene_id) %>% adorn_totals("row", name = "male_counts") %>% as_tibble %>% 
  filter(gene_name == "male_counts") %>% transpose_tibble(old_rownames = "gene_name", new_rownames = "sample") %>% mutate(gender = case_when(male_counts > 10 ~ "male", TRUE ~ "female"))
kiskinis.metadata <- read_csv(here(shared_path,"public-data/ipsc-mn-sod1-moccia-2014/sample-details/samplesheet.csv")) %>% distinct(sample, .keep_all = TRUE) %>%
  mutate(database_dir = here(shared_path,"public-data/ipsc-mn-sod1-moccia-2014"), dataset = "kiskinis", DIV = 45, condition = factor(condition, levels = c("ctrl", "als"))) %>% left_join(select(kiskinis.gene_tpm.male_counts, sample, gender))
sareen.gene_tpm.male_counts = read_tsv(here(shared_path,"public-data/ipsc-mn-c9orf72-sareen-2013/nfcore/star_salmon/salmon.merged.gene_tpm.tsv")) %>% filter(gene_name %in% c('RPS4Y1', 'EIFAY', 'DDX3Y', 'KDM5D')) %>% select(-gene_id) %>% adorn_totals("row", name = "male_counts") %>% as_tibble %>% 
  filter(gene_name == "male_counts") %>% transpose_tibble(old_rownames = "gene_name", new_rownames = "sample") %>% mutate(gender = case_when(male_counts > 10 ~ "male", TRUE ~ "female"))
sareen.metadata <- read_csv(here(shared_path,"public-data/ipsc-mn-c9orf72-sareen-2013/sample-details/samplesheet.csv")) %>% distinct(sample, .keep_all = TRUE) %>%
  mutate(database_dir = here(shared_path,"public-data/ipsc-mn-c9orf72-sareen-2013"), dataset = "sareen", DIV = 35, condition = factor(condition, levels = c("ctrl", "als"))) %>% left_join(select(sareen.gene_tpm.male_counts, sample, gender))
sommer.gene_tpm.male_counts = read_tsv(here(shared_path,"public-data/ipsc-mn-c9orf72-sommer-2022/nfcore/star_salmon/salmon.merged.gene_tpm.tsv")) %>% filter(gene_name %in% c('RPS4Y1', 'EIFAY', 'DDX3Y', 'KDM5D')) %>% select(-gene_id) %>% adorn_totals("row", name = "male_counts") %>% as_tibble %>% 
  filter(gene_name == "male_counts") %>% transpose_tibble(old_rownames = "gene_name", new_rownames = "sample") %>% mutate(gender = case_when(male_counts > 10 ~ "male", TRUE ~ "female"))
sommer.metadata <- read_csv(here(shared_path,"public-data/ipsc-mn-c9orf72-sommer-2022/sample-details/samplesheet.csv")) %>% distinct(sample, .keep_all = TRUE) %>%
  mutate(database_dir = here(shared_path,"public-data/ipsc-mn-c9orf72-sommer-2022"), dataset = "sommer", condition = factor(condition, levels = c("ctrl", "als"))) %>% left_join(select(sommer.gene_tpm.male_counts, sample, gender))
sterneckert.gene_tpm.male_counts = read_tsv(here(shared_path,"public-data/ipsc-mn-c9orf72-sterneckert-2020/nfcore/star_salmon/salmon.merged.gene_tpm.tsv")) %>% filter(gene_name %in% c('RPS4Y1', 'EIFAY', 'DDX3Y', 'KDM5D')) %>% select(-gene_id) %>% adorn_totals("row", name = "male_counts") %>% as_tibble %>% 
  filter(gene_name == "male_counts") %>% transpose_tibble(old_rownames = "gene_name", new_rownames = "sample") %>% mutate(gender = case_when(male_counts > 10 ~ "male", TRUE ~ "female"))
sterneckert.metadata <- read_csv(here(shared_path,"public-data/ipsc-mn-c9orf72-sterneckert-2020/sample-details/samplesheet.csv")) %>% distinct(sample, .keep_all = TRUE) %>%
  mutate(database_dir = here(shared_path,"public-data/ipsc-mn-c9orf72-sterneckert-2020"), dataset = "sterneckert", DIV = 40, condition = factor(condition, levels = c("ctrl", "als"))) %>% left_join(select(sterneckert.gene_tpm.male_counts, sample, gender))
kapeli.gene_tpm.male_counts = read_tsv(here(shared_path,"public-data/ipsc-mn-fus-kapeli-2016/nfcore/star_salmon/salmon.merged.gene_tpm.tsv")) %>% filter(gene_name %in% c('RPS4Y1', 'EIFAY', 'DDX3Y', 'KDM5D')) %>% select(-gene_id) %>% adorn_totals("row", name = "male_counts") %>% as_tibble %>% 
  filter(gene_name == "male_counts") %>% transpose_tibble(old_rownames = "gene_name", new_rownames = "sample") %>% mutate(gender = case_when(male_counts > 10 ~ "male", TRUE ~ "female"))
kapeli.metadata <- read_csv(here(shared_path,"public-data/ipsc-mn-fus-kapeli-2016/sample-details/samplesheet.csv")) %>% distinct(sample, .keep_all = TRUE) %>%  filter(sh_rna != "scramble") %>%
  mutate(database_dir = here(shared_path,"public-data/ipsc-mn-fus-kapeli-2016"), dataset = "kapeli", DIV = 35, condition = factor(condition, levels = c("ctrl", "als"))) %>% left_join(select(kapeli.gene_tpm.male_counts, sample, gender))
desantis.gene_tpm.male_counts = read_tsv(here(shared_path,"public-data/ipsc-mn-fus-desantis-2017/nfcore/star_salmon/salmon.merged.gene_tpm.tsv")) %>% filter(gene_name %in% c('RPS4Y1', 'EIFAY', 'DDX3Y', 'KDM5D')) %>% select(-gene_id) %>% adorn_totals("row", name = "male_counts") %>% as_tibble %>% 
  filter(gene_name == "male_counts") %>% transpose_tibble(old_rownames = "gene_name", new_rownames = "sample") %>% mutate(gender = case_when(male_counts > 10 ~ "male", TRUE ~ "female"))
desantis.metadata <- read_csv(here(shared_path,"public-data/ipsc-mn-fus-desantis-2017/sample-details/samplesheet.csv")) %>% distinct(sample, .keep_all = TRUE) %>%
  mutate(database_dir = here(shared_path,"public-data/ipsc-mn-fus-desantis-2017"), DIV = 19, dataset = "desantis", condition = factor(condition, levels = c("ctrl", "als"))) %>% select(-days_of_differentiation) %>% left_join(select(desantis.gene_tpm.male_counts, sample, gender))
smith.gene_tpm.male_counts = read_tsv(here(shared_path,"public-data/ipsc-mn-tardbp-smith-2021/nfcore/star_salmon/salmon.merged.gene_tpm.tsv")) %>% filter(gene_name %in% c('RPS4Y1', 'EIFAY', 'DDX3Y', 'KDM5D')) %>% select(-gene_id) %>% adorn_totals("row", name = "male_counts") %>% as_tibble %>% 
  filter(gene_name == "male_counts") %>% transpose_tibble(old_rownames = "gene_name", new_rownames = "sample") %>% mutate(gender = case_when(male_counts > 10 ~ "male", TRUE ~ "female"))
smith.metadata <- read_csv(here(shared_path,"public-data/ipsc-mn-tardbp-smith-2021/sample-details/samplesheet.csv")) %>% distinct(sample, .keep_all = TRUE) %>%
  mutate(database_dir = here(shared_path,"public-data/ipsc-mn-tardbp-smith-2021"), dataset = "smith", DIV = 45, condition = factor(condition, levels = c("ctrl", "als"))) %>% left_join(select(smith.gene_tpm.male_counts, sample, gender))
bhinge.gene_tpm.male_counts = read_tsv(here(shared_path,"public-data/ipsc-mn-sod1-bhinge-2017/nfcore/star_salmon/salmon.merged.gene_tpm.tsv")) %>% filter(gene_name %in% c('RPS4Y1', 'EIFAY', 'DDX3Y', 'KDM5D')) %>% select(-gene_id) %>% adorn_totals("row", name = "male_counts") %>% as_tibble %>% 
  filter(gene_name == "male_counts") %>% transpose_tibble(old_rownames = "gene_name", new_rownames = "sample") %>% mutate(gender = case_when(male_counts > 10 ~ "male", TRUE ~ "female"))
bhinge.metadata <- read_csv(here(shared_path,"public-data/ipsc-mn-sod1-bhinge-2017/sample-details/samplesheet.csv")) %>% distinct(sample, .keep_all = TRUE) %>%
  mutate(database_dir = here(shared_path,"public-data/ipsc-mn-sod1-bhinge-2017"), dataset = "bhinge", DIV = 32, condition = factor(condition, levels = c("ctrl", "als"))) %>% select(-age) %>% left_join(select(bhinge.gene_tpm.male_counts, sample, gender))
hawkins.gene_tpm.male_counts = read_tsv(here(shared_path,"public-data/ipsc-mn-fus-hawkins-2022/nfcore/star_salmon/salmon.merged.gene_tpm.tsv")) %>% filter(gene_name %in% c('RPS4Y1', 'EIFAY', 'DDX3Y', 'KDM5D')) %>% select(-gene_id) %>% adorn_totals("row", name = "male_counts") %>% as_tibble %>% 
  filter(gene_name == "male_counts") %>% transpose_tibble(old_rownames = "gene_name", new_rownames = "sample") %>% mutate(gender = case_when(male_counts > 10 ~ "male", TRUE ~ "female"))
hawkins.metadata <- read_csv(here(shared_path,"public-data/ipsc-mn-fus-hawkins-2022/sample-details/samplesheet.csv")) %>% distinct(sample, .keep_all = TRUE) %>%
  mutate(database_dir = here(shared_path,"public-data/ipsc-mn-fus-hawkins-2022"), condition = factor(condition, levels = c("ctrl", "als"))) %>% left_join(select(hawkins.gene_tpm.male_counts, sample, gender))
lee.gene_tpm.male_counts = read_tsv(here(shared_path,"public-data/ipsc-occular-spinal-mn-lee-2021/nfcore/star_salmon/salmon.merged.gene_tpm.tsv")) %>% filter(gene_name %in% c('RPS4Y1', 'EIFAY', 'DDX3Y', 'KDM5D')) %>% select(-gene_id) %>% adorn_totals("row", name = "male_counts") %>% as_tibble %>% 
  filter(gene_name == "male_counts") %>% transpose_tibble(old_rownames = "gene_name", new_rownames = "sample") %>% mutate(gender = case_when(male_counts > 10 ~ "male", TRUE ~ "female"))
lee.metadata <- read_csv(here(shared_path,"public-data/ipsc-occular-spinal-mn-lee-2021/sample-details/samplesheet.csv")) %>% distinct(sample, .keep_all = TRUE) %>%
  mutate(database_dir = here(shared_path,"public-data/ipsc-occular-spinal-mn-lee-2021"), dataset = case_when(mutation == "ctrl" ~ "lee.CTRL", TRUE ~ "lee.ALS"), condition = case_when(mutation == "ctrl" ~ "ctrl", TRUE ~ "als"),
         strandedness = case_when(condition == "als" ~ "reverse", TRUE ~ strandedness), condition = factor(condition, levels = c("ctrl", "als"))) %>% filter(location == "spinal") %>% left_join(select(lee.gene_tpm.male_counts, sample, gender))
# markmiller.gene_tpm.male_counts = read_tsv(here(shared_path,"public-data/ipsc-mn-hnrnpa2b1-markmiller-2021/nfcore/star_salmon/salmon.merged.gene_tpm.tsv")) %>% filter(gene_name %in% c('RPS4Y1', 'EIFAY', 'DDX3Y', 'KDM5D')) %>% select(-gene_id) %>% adorn_totals("row", name = "male_counts") %>% as_tibble %>% 
#   filter(gene_name == "male_counts") %>% transpose_tibble(old_rownames = "gene_name", new_rownames = "sample") %>% mutate(gender = case_when(male_counts > 10 ~ "male", TRUE ~ "female"))
# markmiller.metadata <- read_csv(here(shared_path,"public-data/ipsc-mn-hnrnpa2b1-markmiller-2021/sample-details/samplesheet.csv")) %>% distinct(sample, .keep_all = TRUE) %>%
#   mutate(database_dir = here(shared_path,"public-data/ipsc-mn-hnrnpa2b1-markmiller-2021"), condition = factor(condition, levels = c("ctrl", "als"))) %>% left_join(select(markmiller.gene_tpm.male_counts, sample, gender))
# dash.gene_tpm.male_counts = read_tsv(here(shared_path,"public-data/ipsc-mn-sod1-tardbp-dash-2022/nfcore/star_salmon/salmon.merged.gene_tpm.tsv")) %>% filter(gene_name %in% c('RPS4Y1', 'EIFAY', 'DDX3Y', 'KDM5D')) %>% select(-gene_id) %>% adorn_totals("row", name = "male_counts") %>% as_tibble %>%
#   filter(gene_name == "male_counts") %>% transpose_tibble(old_rownames = "gene_name", new_rownames = "sample") %>% mutate(gender = case_when(male_counts > 10 ~ "male", TRUE ~ "female"))
# dash.metadata <- read_csv(here(shared_path,"public-data/ipsc-mn-sod1-tardbp-dash-2022/sample-details/samplesheet.csv")) %>% distinct(sample, .keep_all = TRUE) %>%
#   mutate(database_dir = here(shared_path,"public-data/ipsc-mn-sod1-tardbp-dash-2022"), condition = factor(condition, levels = c("ctrl", "als"))) %>% left_join(select(dash.gene_tpm.male_counts, sample, gender))

# counts_all = read_tsv("/camp/lab/luscomben/home/users/ziffo/projects/ipsc-mn-als-meta/nfcore/star_salmon/salmon.merged.gene_counts.tsv")

nygc.gene_tpm.male_counts = read_tsv(here(shared_path,"nygc-als-consortium/nfcore/star_salmon/salmon.merged.gene_tpm.tsv")) %>% filter(gene_name %in% c('RPS4Y1', 'EIFAY', 'DDX3Y', 'KDM5D')) %>% select(-gene_id) %>% adorn_totals("row", name = "male_counts") %>% as_tibble %>% 
  filter(gene_name == "male_counts") %>% transpose_tibble(old_rownames = "gene_name", new_rownames = "sample") %>% mutate(gender = case_when(male_counts > 10 ~ "male", TRUE ~ "female"))
nygc_spinal_cord_metadata = read_csv(here(shared_path, "nygc-als-consortium/sample-details/nygc_spinal_cord_metadata.csv")) %>% distinct(sample, .keep_all = TRUE) %>%
  mutate(database_dir = here(shared_path,"nygc-als-consortium")) %>% mutate(condition = factor(condition, levels = c("ctrl", "als"))) %>% left_join(select(nygc.gene_tpm.male_counts, sample, gender))


# Metadata merge ----------------------------------------------------------------
ipsc_mn_als_datasets_with_lee_ipsc.metadata <- bind_rows(neurolincs.ipsc.metadata, answerals.metadata, neurolincs.metadata, catanese.metadata, dafinca.c9orf72.metadata, dafinca.tardbp.metadata, luisier.metadata, 
                                                         wang.metadata, kiskinis.metadata, sareen.metadata, sommer.metadata, sterneckert.metadata, kapeli.metadata, desantis.metadata, smith.metadata, bhinge.metadata, hawkins.metadata, lee.metadata) %>%
  mutate(dataset_sample = paste0(dataset,"_",sample), condition = factor(condition, levels = c("ctrl", "als")), dataset_simplified = case_when(dataset %in% c("dafinca.c9orf72","dafinca.tardbp") ~ "dafinca", TRUE ~ dataset), 
         dataset = as.factor(dataset), mutation = as.factor(mutation), DIV = as.factor(DIV), instrument = gsub("Illumina ","",instrument), 
         library_type = case_when(dataset %in% c("catanese","dafinca.c9orf72","dafinca.tardbp","kapeli","kiskinis","luisier","mitchell","sareen","sommer","smith","sterneckert","wang","lee.CTRL") ~ "poly(A)", TRUE ~ "Total")) %>%
  select(sample, dataset, condition, mutation, dataset_sample, replicate, database_dir, everything())
ipsc_mn_als_datasets_with_lee.metadata <- ipsc_mn_als_datasets_with_lee_ipsc.metadata %>% filter(!dataset %in% c("neurolincs.iPSC")) %>% mutate(dataset = fct_drop(dataset), dataset_simplified = fct_drop(dataset_simplified))
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
ipsc_mn_tdp_pathology.metadata = filter(ipsc_mn_als_datasets.metadata, mutation != "ctrl") %>% mutate(tdp_pathology = case_when(mutation %in% c("sod1","fus") ~ "no", TRUE ~ "yes"))

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
   hawkins.gene_tpm.male_counts, nygc.gene_tpm.male_counts) 

save.image(here(proj_path,"expression/deseq2/ipsc_mn_als_meta.metadata.RData"))

# DESeq2 ------------------------------------------------------------------

# ### ALS ------------------------------------------------------------------
# ipsc_mn_als_datasets_with_lee_ipsc = DESeq.analysis(metadata = ipsc_mn_als_datasets_with_lee_ipsc.metadata, unique_names = "dataset_sample")
# saveRDS(ipsc_mn_als_datasets_with_lee_ipsc, here(proj_path,"expression/deseq2/ipsc_mn_als_datasets_with_lee_ipsc.rds"))
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

# ipsc_mn_als_datasets = DESeq.analysis(metadata = ipsc_mn_als_datasets.metadata, unique_names = "dataset_sample", design = ~ gender + dataset + condition, contrast = "condition_als_vs_ctrl")
# saveRDS(ipsc_mn_als_datasets, here(proj_path,"expression/deseq2/ipsc_mn_als_datasets.rds"))
# ipsc_mn_als_datasets.noiso = DESeq.analysis(metadata = ipsc_mn_als_datasets.noiso.metadata, unique_names = "dataset_sample", design = ~ gender + dataset + condition, contrast = "condition_als_vs_ctrl")
# saveRDS(ipsc_mn_als_datasets.noiso, here(proj_path,"expression/deseq2/ipsc_mn_als_datasets.noiso.rds"))

# POU5F1_gene_counts = read_tsv(here(proj_path,"nfcore/star_salmon/salmon.merged.gene_tpm.tsv")) %>% filter(gene_name %in% c("POU5F1"))
# POU5F1_gene_counts.metadata = POU5F1_gene_counts %>% select(-gene_id) %>% pivot_longer(!gene_name, names_to = "dataset_sample", values_to = "value") %>% left_join(ipsc_mn_als_datasets_with_lee.metadata)
# POU5F1_expressing_ipsns = filter(POU5F1_gene_counts.metadata, value > 5) 
# ipsc_mn_als_datasets.noPOU5F1 = DESeq.analysis(metadata = filter(ipsc_mn_als_datasets.metadata, !dataset_sample %in% POU5F1_expressing_ipsns$dataset_sample), unique_names = "dataset_sample", design = ~ gender + dataset + condition, contrast = "condition_als_vs_ctrl")
# saveRDS(ipsc_mn_als_datasets.noPOU5F1, here(proj_path,"expression/deseq2/ipsc_mn_als_datasets.noPOU5F1.rds"))

# ### ALS mutation groupings ------------------------------------------------------------------
# ipsc_mn_sporadic_datasets = DESeq.analysis(metadata = ipsc_mn_sporadic_datasets.metadata, unique_names = "dataset_sample", design = ~ gender + dataset + condition, contrast = "condition_als_vs_ctrl", species = "human")
# ipsc_mn_c9orf72_datasets = DESeq.analysis(metadata = ipsc_mn_c9orf72_datasets.metadata, unique_names = "dataset_sample", design = ~ gender + dataset + condition, contrast = "condition_als_vs_ctrl", species = "human")
# ipsc_mn_fus_datasets = DESeq.analysis(metadata = ipsc_mn_fus_datasets.metadata, unique_names = "dataset_sample", design = ~ gender + dataset + condition, contrast = "condition_als_vs_ctrl", species = "human")
# ipsc_mn_sod1_datasets = DESeq.analysis(metadata = ipsc_mn_sod1_datasets.metadata, unique_names = "dataset_sample", design = ~ gender + dataset + condition, contrast = "condition_als_vs_ctrl", species = "human")
# ipsc_mn_tardbp_datasets = DESeq.analysis(metadata = ipsc_mn_tardbp_datasets.metadata, unique_names = "dataset_sample", design = ~ gender + dataset + condition, contrast = "condition_als_vs_ctrl", species = "human")
# save(ipsc_mn_sporadic_datasets, ipsc_mn_c9orf72_datasets, ipsc_mn_fus_datasets, ipsc_mn_sod1_datasets,ipsc_mn_tardbp_datasets,
#      file = here(proj_path, "expression/deseq2/ipsc_mn_genetic_groups.RData"))
# 
# ipsc_mn_tdp43_positive_datasets = DESeq.analysis(metadata = filter(ipsc_mn_als_datasets.metadata, !mutation %in% c("sod1","fus")), unique_names = "dataset_sample", design = ~ gender + dataset + condition, contrast = "condition_als_vs_ctrl", species = "human")
# ipsc_mn_tdp43_negative_datasets = DESeq.analysis(metadata = filter(ipsc_mn_als_datasets.metadata, mutation %in% c("sod1","fus","ctrl"), dataset %in% c("answerals","bhinge","catanese", "desantis", "hawkins", "kapeli", "kiskinis", "neurolincs.diMN", "wang")), unique_names = "dataset_sample", design = ~ gender + dataset + condition, contrast = "condition_als_vs_ctrl", species = "human")
# ipsc_mn_tdp_pathology = DESeq.analysis(metadata = ipsc_mn_tdp_pathology.metadata, unique_names = "dataset_sample", design = ~ gender + dataset + tdp_pathology, contrast = "tdp_pathology_yes_vs_no", species = "human")
# save(ipsc_mn_tdp43_positive_datasets, ipsc_mn_tdp43_negative_datasets, ipsc_mn_tdp_pathology,
#      file = here(proj_path, "expression/deseq2/ipsc_mn_tdp43_status.RData"))
# 
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

# ### Datasets independently ------------------------------------------------------------------
# answerals = DESeq.analysis(metadata = answerals.metadata, design = ~ gender + condition, contrast = "condition_als_vs_ctrl")
# neurolincs.diMN = DESeq.analysis(metadata = filter(neurolincs.metadata, grepl("^0007",sample)), design = ~ gender + condition, contrast = "condition_als_vs_ctrl")
# neurolincs.iMN = DESeq.analysis(metadata = filter(neurolincs.metadata, grepl("^A-042",sample)), design = ~ gender + condition, contrast = "condition_als_vs_ctrl")
# catanese = DESeq.analysis(metadata = catanese.metadata, design = ~ condition, contrast = "condition_als_vs_ctrl")
# dafinca.c9orf72 = DESeq.analysis(metadata = dafinca.c9orf72.metadata, design = ~ condition, contrast = "condition_als_vs_ctrl")
# dafinca.tardbp = DESeq.analysis(metadata = dafinca.tardbp.metadata, design = ~ condition, contrast = "condition_als_vs_ctrl")
# luisier = DESeq.analysis(metadata = luisier.metadata, design = ~ condition, contrast = "condition_als_vs_ctrl")
# wang = DESeq.analysis(metadata = wang.metadata, design = ~ condition, contrast = "condition_als_vs_ctrl")
# kiskinis = DESeq.analysis(metadata = kiskinis.metadata, design = ~ condition, contrast = "condition_als_vs_ctrl")
# sareen = DESeq.analysis(metadata = sareen.metadata, design = ~ gender + condition, contrast = "condition_als_vs_ctrl")
# sommer = DESeq.analysis(metadata = sommer.metadata, design = ~ gender + condition, contrast = "condition_als_vs_ctrl")
# sterneckert = DESeq.analysis(metadata = sterneckert.metadata, design = ~ condition, contrast = "condition_als_vs_ctrl")
# kapeli = DESeq.analysis(metadata = kapeli.metadata, design = ~ condition, contrast = "condition_als_vs_ctrl")
# desantis = DESeq.analysis(metadata = desantis.metadata, design = ~ condition, contrast = "condition_als_vs_ctrl")
# smith = DESeq.analysis(metadata = smith.metadata, design = ~ condition, contrast = "condition_als_vs_ctrl")
# bhinge = DESeq.analysis(metadata = bhinge.metadata, design = ~ condition, contrast = "condition_als_vs_ctrl")
# hawkins = DESeq.analysis(metadata = hawkins.metadata, design = ~ condition, contrast = "condition_als_vs_ctrl")
# #markmiller = DESeq.analysis(metadata = markmiller.metadata, design = ~ condition, contrast = "condition_als_vs_ctrl")
# # dash.tardbp = DESeq.analysis(metadata = filter(dash.metadata,mutation %in% c("ctrl","tardbp")), design = ~ gender + condition, contrast = "condition_als_vs_ctrl")
# # dash.sod1 = DESeq.analysis(metadata = filter(dash.metadata,mutation %in% c("ctrl","sod1")), design = ~ gender + condition, contrast = "condition_als_vs_ctrl")

# save(answerals, neurolincs.diMN, neurolincs.iMN, catanese, dafinca.c9orf72, dafinca.tardbp, luisier, wang, kiskinis, sareen, sommer, sterneckert, kapeli, desantis, smith, bhinge, hawkins, #markmiller, dash.tardbp, dash.sod1,
#   file = here(proj_path, "expression/deseq2/datasets_separated.RData"))

# neurolincs.iPSC = DESeq.analysis(metadata = neurolincs.ipsc.metadata)
# saveRDS(neurolincs.iPSC, file = here(proj_path, "expression/deseq2/neurolincs.ipsc.rds"))

# ### NYGC Post Mortem ------------------------------------------------------------------
# nygc_postmortem_spinal_cord = DESeq.analysis(metadata = nygc_spinal_cord_metadata, design = ~ gender + library_prep + sample_source + condition, contrast = "condition_als_vs_ctrl")
# nygc_postmortem_spinal_cord.tx = DESeq.analysis(metadata = nygc_spinal_cord_metadata, design = ~ gender + library_prep + sample_source + condition, transcript.level = TRUE, gene.level = FALSE, contrast = "condition_als_vs_ctrl")

# 
# nygc_spinal_cord_metadata %>% count(mutation) %>% arrange(-n)
# mutation     n
# <chr>    <int>
# 1 sporadic   161
# 2 ctrl        57
# 3 c9orf72     36
# 4 sod1         5
# 5 fus          2
# 6 nefh         2
# 7 setx         2
# 8 ang          1
# 9 dctn1        1
# 10 fig4         1
# 11 optn         1
# 12 ubqln2       1
# 13 vcp          1
# nygc_postmortem_spinal_cord.sporadic = DESeq.analysis(metadata = filter(nygc_spinal_cord_metadata, mutation %in% c("sporadic","ctrl")), design = ~ gender + library_prep + sample_source + condition, contrast = "condition_als_vs_ctrl")
# nygc_postmortem_spinal_cord.c9orf72 = DESeq.analysis(metadata = filter(nygc_spinal_cord_metadata, mutation %in% c("c9orf72","ctrl")), design = ~ gender + library_prep + sample_source + condition, contrast = "condition_als_vs_ctrl")
# nygc_postmortem_spinal_cord.sod1 = DESeq.analysis(metadata = filter(nygc_spinal_cord_metadata, mutation %in% c("sod1","ctrl")), design = ~ gender + library_prep + sample_source + condition, contrast = "condition_als_vs_ctrl")
# nygc_postmortem_spinal_cord.fus = DESeq.analysis(metadata = filter(nygc_spinal_cord_metadata, mutation %in% c("fus","ctrl")), design = ~ gender + library_prep + sample_source + condition, contrast = "condition_als_vs_ctrl")
# save(nygc_postmortem_spinal_cord, nygc_postmortem_spinal_cord.sporadic, nygc_postmortem_spinal_cord.c9orf72, nygc_postmortem_spinal_cord.sod1, nygc_postmortem_spinal_cord.fus,
#      file = here(proj_path, "expression/deseq2/nygc_postmortem_spinal_cord.RData"))
# 
# nygc_postmortem_spinal_cord_tdp43_positive_datasets = DESeq.analysis(metadata = filter(nygc_spinal_cord_metadata, !mutation %in% c("sod1","fus")), design = ~ gender + library_prep + sample_source + condition, contrast = "condition_als_vs_ctrl")
# nygc_postmortem_spinal_cord_tdp43_negative_datasets = DESeq.analysis(metadata = filter(nygc_spinal_cord_metadata, mutation %in% c("sod1","fus","ctrl")), design = ~ gender + library_prep + sample_source + condition, contrast = "condition_als_vs_ctrl")
# nygc_postmortem_spinal_cord_tdp_pathology.metadata = filter(nygc_spinal_cord_metadata, mutation != "ctrl") %>% mutate(tdp_pathology = case_when(mutation %in% c("sod1","fus") ~ "no", TRUE ~ "yes"))
# nygc_postmortem_spinal_cord_tdp_pathology = DESeq.analysis(metadata = nygc_postmortem_spinal_cord_tdp_pathology.metadata, design = ~ gender + library_prep + sample_source + tdp_pathology, contrast = "tdp_pathology_yes_vs_no")
# save(nygc_postmortem_spinal_cord_tdp43_positive_datasets, nygc_postmortem_spinal_cord_tdp43_negative_datasets, nygc_postmortem_spinal_cord_tdp_pathology,
#      file = here(proj_path, "expression/deseq2/nygc_postmortem_spinal_cord_tdp43_status.RData"))
# 
# ## TDP-43 Knockdown datasets ------------------------------------------------------------------
neuronal_nuclei_tdp43_liu.metadata <- read_csv(here(shared_path,"public-data/neuronal-nuclei-tdp43-liu-2019/sample-details/samplesheet.csv")) %>% distinct(sample, .keep_all = TRUE) %>%
  mutate(database_dir = here(shared_path,"public-data/neuronal-nuclei-tdp43-liu-2019"), group = condition, condition = case_when(condition == "tdp43pos" ~ "ctrl", condition == "tdp43neg" ~ "tdp43kd"), 
         condition = factor(condition, levels = c("ctrl", "tdp43kd")), dataset = "liu")
tdp43_kd_melamed.metadata <- read_csv(here(shared_path,"public-data/shsy5y-tdp43-kd-melamed-2019/sample-details/samplesheet.csv")) %>% distinct(sample, .keep_all = TRUE) %>%
  mutate(database_dir = here(shared_path,"public-data/shsy5y-tdp43-kd-melamed-2019"), condition = factor(condition, levels = c("ctrl", "tdp43kd")), dataset = "melamed")
# tdp43_kd_kapeli.metadata <- read_csv(here(shared_path,"public-data/ipsc-mn-tdp43-kd-kapeli-2016/sample-details/samplesheet.csv")) %>% distinct(sample, .keep_all = TRUE) %>%
#   mutate(database_dir = here(shared_path,"public-data/ipsc-mn-tdp43-kd-kapeli-2016"), condition = factor(condition, levels = c("ctrl", "tdp43kd")), dataset = "kapeli")
tdp43_kd_ferguson.metadata <- read_csv(here(shared_path,"public-data/hela-tdp43-kd-ferguson-2019/sample-details/samplesheet.csv")) %>% distinct(sample, .keep_all = TRUE) %>%
  mutate(database_dir = here(shared_path,"public-data/hela-tdp43-kd-ferguson-2019"), condition = factor(condition, levels = c("ctrl", "tdp43kd")), dataset = "ferguson")
tdp43_kd_tam.metadata <- read_csv(here(shared_path,"public-data/shsy5y-tdp43-kd-tam-2019/sample-details/samplesheet.csv")) %>% distinct(sample, .keep_all = TRUE) %>%
  mutate(database_dir = here(shared_path,"public-data/shsy5y-tdp43-kd-tam-2019"), condition = factor(condition, levels = c("ctrl", "tdp43kd")), dataset = "tam")
tdp43_kd_klim.metadata <- read_csv(here(shared_path,"public-data/ipsc-mn-tardbp-kd-klim-2019/sample-details/samplesheet.csv")) %>% distinct(sample, .keep_all = TRUE) %>%
  mutate(database_dir = here(shared_path,"public-data/ipsc-mn-tardbp-kd-klim-2019"), condition = factor(condition, levels = c("ctrl", "tdp43kd")), dataset = "klim")
tdp43_kd_appocher.metadata <- read_csv(here(shared_path,"public-data/shsy5y-tardbp-kd-appocher-2017/sample-details/samplesheet.csv")) %>% distinct(sample, .keep_all = TRUE) %>%
  mutate(database_dir = here(shared_path,"public-data/shsy5y-tardbp-kd-appocher-2017"), condition = factor(condition, levels = c("ctrl", "tdp43kd")), dataset = "appocher")
tdp43_kd_dunker.metadata <- read_csv(here(shared_path,"public-data/shsy5y-tdp43-kd-dunker-2021/sample-details/samplesheet.csv")) %>% distinct(sample, .keep_all = TRUE) %>%
  mutate(database_dir = here(shared_path,"public-data/shsy5y-tdp43-kd-dunker-2021"), condition = factor(condition, levels = c("ctrl", "tdp43kd")), dataset = "dunker")
tdp43_kd_brown_ipsmn.metadata <- read_csv(here(shared_path,"public-data/ipsc-mn-tardbp-kd-brown-2022/sample-details/samplesheet.csv")) %>% distinct(sample, .keep_all = TRUE) %>%
  mutate(database_dir = here(shared_path,"public-data/ipsc-mn-tardbp-kd-brown-2022"), condition = factor(condition, levels = c("ctrl", "tdp43kd")), dataset = "brown_ipsmn")
# tdp43_kd_brown_shsy5y.metadata <- read_csv(here(shared_path,"public-data/shsy5y-tdp43-kd-brown-2022/sample-details/samplesheet.csv")) %>% distinct(sample, .keep_all = TRUE) %>% filter(celltype == "shsy5y") %>%
#   mutate(database_dir = here(shared_path,"public-data/shsy5y-tdp43-kd-brown-2022"), condition = factor(condition, levels = c("ctrl", "tdp43kd")), dataset = "brown_shsy5y")
tdp43_kd_datasets.metadata <- bind_rows(tdp43_kd_klim.metadata, tdp43_kd_appocher.metadata, tdp43_kd_melamed.metadata, tdp43_kd_ferguson.metadata, tdp43_kd_tam.metadata, tdp43_kd_dunker.metadata, tdp43_kd_brown_ipsmn.metadata) %>% #tdp43_kd_kapeli.metadata, tdp43_kd_brown_shsy5y.metadata, tdp43_kd_brown_skndz.metadata, 
  mutate(dataset_sample = paste0(dataset,"_",sample), condition = factor(condition, levels = c("ctrl", "tdp43kd")), dataset = as.factor(dataset)) %>%
  select(sample, dataset, condition, dataset_sample, database_dir, everything())

neuronal_nuclei_tdp43_liu = DESeq.analysis(metadata = neuronal_nuclei_tdp43_liu.metadata, design = ~ condition, contrast = "condition_tdp43kd_vs_ctrl")
saveRDS(neuronal_nuclei_tdp43_liu, here(proj_path,"expression/deseq2/neuronal_nuclei_tdp43_liu.rds"))
tdp43_kd_melamed = DESeq.analysis(metadata = tdp43_kd_melamed.metadata, design = ~ condition, contrast = "condition_tdp43kd_vs_ctrl")
saveRDS(tdp43_kd_melamed, here(proj_path,"expression/deseq2/tdp43_kd_melamed.rds"))
tdp43_kd_kapeli = DESeq.analysis(metadata = tdp43_kd_kapeli.metadata, design = ~ condition, contrast = "condition_tdp43kd_vs_ctrl")
saveRDS(tdp43_kd_kapeli, here(proj_path,"expression/deseq2/tdp43_kd_kapeli.rds"))
tdp43_kd_ferguson = DESeq.analysis(metadata = tdp43_kd_ferguson.metadata, design = ~ condition, contrast = "condition_tdp43kd_vs_ctrl")
saveRDS(tdp43_kd_ferguson, here(proj_path,"expression/deseq2/tdp43_kd_ferguson.rds"))
tdp43_kd_tam = DESeq.analysis(metadata = tdp43_kd_tam.metadata, design = ~ condition, contrast = "condition_tdp43kd_vs_ctrl")
saveRDS(tdp43_kd_tam, here(proj_path,"expression/deseq2/tdp43_kd_tam.rds"))
tdp43_kd_klim = DESeq.analysis(metadata = tdp43_kd_klim.metadata, design = ~ condition, contrast = "condition_tdp43kd_vs_ctrl")
saveRDS(tdp43_kd_klim, here(proj_path,"expression/deseq2/tdp43_kd_klim.rds"))
tdp43_kd_appocher = DESeq.analysis(metadata = tdp43_kd_appocher.metadata, design = ~ condition, contrast = "condition_tdp43kd_vs_ctrl")
saveRDS(tdp43_kd_appocher, here(proj_path,"expression/deseq2/tdp43_kd_appocher.rds"))
tdp43_kd_dunker = DESeq.analysis(metadata = tdp43_kd_dunker.metadata, design = ~ condition, contrast = "condition_tdp43kd_vs_ctrl")
saveRDS(tdp43_kd_dunker, here(proj_path,"expression/deseq2/tdp43_kd_dunker.rds"))
tdp43_kd_brown_ipsmn = DESeq.analysis(metadata = tdp43_kd_brown_ipsmn.metadata, design = ~ condition, contrast = "condition_tdp43kd_vs_ctrl")
saveRDS(tdp43_kd_brown_ipsmn, here(proj_path,"expression/deseq2/tdp43_kd_brown_ipsmn.rds"))
tdp43_kd_brown_shsy5y = DESeq.analysis(metadata = tdp43_kd_brown_shsy5y.metadata, design = ~ condition, contrast = "condition_tdp43kd_vs_ctrl")
saveRDS(tdp43_kd_brown_shsy5y, here(proj_path,"expression/deseq2/tdp43_kd_brown_shsy5y.rds"))
tdp43_kd_brown_skndz = DESeq.analysis(metadata = tdp43_kd_brown_skndz.metadata, design = ~ condition, contrast = "condition_tdp43kd_vs_ctrl")
saveRDS(tdp43_kd_brown_skndz, here(proj_path,"expression/deseq2/tdp43_kd_brown_skndz.rds"))
tdp43_kd_datasets = DESeq.analysis(metadata = tdp43_kd_datasets.metadata, unique_names = "dataset_sample", design = ~ dataset + condition, contrast = "condition_tdp43kd_vs_ctrl")
saveRDS(tdp43_kd_datasets, here(proj_path,"expression/deseq2/tdp43_kd_datasets.rds"))

# print("saving")
# save.image(here(proj_path,"expression/deseq2/ipsc_mn_als_meta.RData"))
print("complete")

