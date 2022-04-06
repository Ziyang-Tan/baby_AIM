library(dplyr)
library(readr)
library(ggplot2)
library(ggpubr)
library(stringr)

proj_id <- c('P25158_1002')
sample_list = list(baby=c('BabyAIM_1', 'BabyAIM_2', 'BabyAIM_3', 
                          'BabyAIM_4', 'BabyAIM_5', 'BabyAIM_6'))
# load data
source('src/load_BD_scTCR.R')
glob_path <- '/Users/tan/OneDrive - KI.SE/TCR_processed_data/single cell'
raw_tcr <- lapply(proj_id, BD_load_VDJ, dir_path = glob_path) %>% do.call(what = rbind)
sample_tag <- lapply(proj_id, BD_load_sample_tag, dir_path = glob_path) %>% do.call(what = rbind)

cell_type <- read_csv('data/cell_types.csv')

raw_tcr_merge <- left_join(raw_tcr, sample_tag, by = 'unique_index') %>%
  left_join(cell_type, by = 'unique_index') %>%
  mutate(Sample_Name = case_when(
    is.na(Sample_Name) ~ proj_id,
    TRUE ~ Sample_Name
  ))

# summarise 
df <- raw_tcr_merge %>% filter(!Sample_Name %in% c('Multiplet', 'Undetermined'))
table_summary <- df %>% group_by(proj_id) %>% summarise(demultiplexed = n()) %>% left_join(
  df %>% filter(at_least_one_chain) %>% group_by(proj_id) %>% summarise(at_least_one_CDR3 = n())
) %>% left_join(
  df %>% filter(TCR_Paired_Chains) %>% group_by(proj_id) %>% summarise(paired = n())
)

data <- raw_tcr_merge %>%
  filter(!is.na(TCR_Beta_Delta_CDR3_Nucleotide_Dominant)) %>%
  filter(!is.na(TCR_Alpha_Gamma_CDR3_Nucleotide_Dominant)) %>%
  filter(!Sample_Name %in% c('Multiplet', 'Undetermined')) %>%
  mutate(CDR3_concat = paste0(TCR_Alpha_Gamma_CDR3_Nucleotide_Dominant, '_', 
                              TCR_Beta_Delta_CDR3_Nucleotide_Dominant),
         CDR3aa_concat = paste0(TCR_Alpha_Gamma_CDR3_Translation_Dominant, '_', 
                                TCR_Beta_Delta_CDR3_Translation_Dominant)) %>%
  mutate(clone_id = as.character(as.numeric(as.factor(CDR3_concat))))
write_csv(data, file = 'data/scTCR_data_merge.csv')
clone_id_map <- data %>% select(CDR3_concat, clone_id) %>% unique()
clone_exp <- data %>%
  group_by(CDR3_concat, Sample_Name) %>%
  summarise(clone_count = n()) %>%
  ungroup() %>%
  inner_join(clone_id_map, by='CDR3_concat')
  #inner_join(data %>% select(clone_id, CDR3aa_concat) %>% unique(), by='clone_id')
write_csv(clone_exp, file = 'data/clone_expansion.csv')

# donut chart 
source("src/clone_expansion_plots.R")
patient<- 'baby'
for (sub_name in c('CD4T', 'CD8T', 'gdT')){
  fig_dir <- file.path('figures', 'cartridge_wise', patient)
  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
  clone_exp_sub <- data %>%
    filter(cell_type == sub_name) %>%
    # mutate(Sample_Name = paste0(proj_id, '_', Sample_Name)) %>% # when proj_wise plots are needed
    group_by(CDR3_concat, Sample_Name) %>%
    summarise(clone_count = n()) %>%
    ungroup() %>%
    inner_join(clone_id_map, by='CDR3_concat')
  g_list1 <- lapply(sample_list[[patient]],function(x){clone_expansion_donut(x, clone_exp_sub)})
  ggarrange(plotlist = g_list1, ncol = 3, nrow = 5) %>%
    ggexport(filename = file.path(fig_dir, paste0(patient, '_clone_expansion_', sub_name, '.pdf')), 
             width = 10, height = 20)
  g <- clone_expansion_alluvium(patient, clone_exp_sub) + labs(title = paste0(patient, '_gdT'))
  ggsave(plot = g, filename = file.path(fig_dir, paste0(patient, '_top_clone_changes_', sub_name, '.pdf')))
}

# alluvium
#g_list2 <- lapply(c('ISAC99', 'ISAC35'), function(x){clone_expansion_alluvium(x, clone_exp)})
#ggarrange(plotlist = g_list2, ncol = 1, nrow = 2) %>%
#  ggexport(filename = 'shared_clones_ratio.pdf', width = 10, height = 20)
# expansion at different visit
data <- data %>%
  left_join(clone_exp, by = c('CDR3_concat', 'Sample_Name')) %>%
  group_by(Sample_Name) %>%
  mutate(clone_ratio = clone_count/n())

df_exp <- data %>%
  select(clone_id, clone_count, clone_ratio, Sample_Name) %>%
  filter(clone_count > 1) %>%
  mutate(clone_ratio_bin = case_when(
    clone_ratio < 0.003 ~ '<0.3%',
    clone_ratio >= 0.003 & clone_ratio < 0.005 ~ '0.3-0.5%',
    clone_ratio >= 0.005 & clone_ratio < 0.025 ~ '0.5%-2.5%',
    clone_ratio >= 0.025 ~ '>=2.5%'
  )) %>%
  mutate(clone_ratio_bin = factor(clone_ratio_bin, levels = c('<0.3%', '0.3-0.5%', '0.5%-2.5%', '>=2.5%'))) %>%
  unique()

ggplot(df_exp, aes(fill=clone_ratio_bin, x=Sample_Name)) + 
  geom_bar(position="fill", stat="count") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# typing
tmp = data %>%
  filter(TCR_Paired_Chains)%>%
  select(TCR_Alpha_Gamma_V_gene_Dominant, TCR_Beta_Delta_V_gene_Dominant, Sample_Name) %>% 
  filter(!is.na(TCR_Alpha_Gamma_V_gene_Dominant)) %>% 
  filter(!is.na(TCR_Beta_Delta_V_gene_Dominant)) %>% 
  mutate(type_a=str_extract(TCR_Alpha_Gamma_V_gene_Dominant, "^.{4}")) %>% 
  mutate(type_b=str_extract(TCR_Beta_Delta_V_gene_Dominant, "^.{4}")) %>% 
  mutate(type = paste0(type_a, '_', type_b))%>%
  group_by(type) %>% 
  count()