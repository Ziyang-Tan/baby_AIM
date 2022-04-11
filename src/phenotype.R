library(dplyr)
library(readr)
source("src/clone_expansion_plots.R")
source("src/customized plots.R")
library(ggplot2)
library(Seurat)
library(SeuratDisk)
library(wesanderson)

sc_data <- read_csv(file = 'data/scTCR_data_merge.csv') %>%
  mutate(seurat_clusters = as.character(seurat_clusters))
clone_exp <- read_csv(file = 'data/clone_expansion.csv')
obj <- LoadH5Seurat('data/seurat_results.h5Seurat')
cell_type <- read_csv('data/cell_types.csv')
sample_info <- read_csv('data/babyAIM_sample_info.csv')
rownames(cell_type) <- cell_type$unique_index
obj$cell_type <- cell_type$cell_type
obj$Sample_Name <- as.factor(obj$Sample_Name)


## temp
# add time point info

##

## all cells

# marker name mapping: 
# syntax: "common_name: name in the panel"
# CD69: CD69, CD25: IL2RA, 
# CD134 (OX40): TNFRSF4, CD137 (TNFRSF9): TNFRSF9, 
# CD274: CD274, CD279 (PD1): PDCD1, 
# CD70: CD70, LAG-3/CD223: LAG3, CD154: CD40LG
AIM_gene_list <- c('CD69', 'IL2RA', 'TNFRSF4', 'TNFRSF9', 'CD274', 'PDCD1', 'CD70', 'LAG3', 'CD40LG')
FeaturePlot(obj, features = AIM_gene_list)
ggsave('figures/AIM_markers.pdf')

FeatureScatter(obj, feature1 = 'TNFRSF9', feature2 = 'TNFRSF4', slot = 'data', group.by = 'cell_type') + 
  theme(plot.title = element_blank()) +
  xlab('TNFRSF9(CD137)') + ylab('TNFRSF4(CD134)')
ggsave('figures/CD134_CD137.pdf')

all_features <- FindAllMarkers(obj)
top_features <- all_features %>% 
  group_by(cluster) %>% 
  slice_max(n=10, order_by = avg_log2FC)
DoHeatmap(obj, features = top_features$gene)
ggsave('figures/feature_expression_all_cells.pdf', height = 15, width = 15)

cluster_stats <- obj[[]] %>% 
  group_by(seurat_clusters, Sample_Name) %>% 
  summarise(count=n()) %>%
  ungroup() %>%
  tidyr::complete(seurat_clusters, Sample_Name, fill = list(count=0)) %>%
  group_by(Sample_Name) %>%
  mutate(freq = round(count/sum(count), 3))

sub_stats <- cluster_stats %>% filter(seurat_clusters == 10)  # 10 is Tregs
ggplot(sub_stats, aes(x=Sample_Name, y=freq, label=count)) +
  geom_bar(stat = 'identity') + geom_text(position = position_stack(vjust = 0.5)) + 
  labs(title = 'all seurat cluster 10 (Tregs)') + ylab('frequency among all cells') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave('figures/all_Tregs_frequency_by_sample.pdf')
##


## AIMpos
# clone expansion of cells with AIM pos gene expression
obj_sub <- subset(obj, 
                  subset = CD69>1 | CD40LG>1 | TNFRSF9>1 | TNFRSF4>1 |
                    IL2RA>1 | CD274>1 | PDCD1>1 | CD70>1 | LAG3>1) # cells express any of the AIM markers
data <- sc_data %>% 
  filter(unique_index %in% obj_sub$unique_index) %>%
  group_by(clone_id) %>%
  tally() %>%
  filter(n>1)
target <- sc_data %>% filter(clone_id %in% data$clone_id,
                             unique_index %in% obj_sub$unique_index) %>%
  arrange(clone_id)
write.csv(target, file='data/tcr_clone_expansion_AIMpos.csv')

# gene expression of AIMpos cells
tmp <- sc_data$clone_id
names(tmp) <- sc_data$unique_index
obj_sub$clone_id = tmp
Idents(obj_sub) <- 'seurat_clusters'
all_features <- FindAllMarkers(obj_sub)
top_features <- all_features %>% 
  group_by(cluster) %>% 
  slice_max(n=10, order_by = avg_log2FC)
gene_list <- c((top_features %>% filter(!(gene %in% AIM_gene_list)))$gene, AIM_gene_list)
DoHeatmap(obj_sub, features = gene_list)
ggsave('figures/feature_expression_all_AIMpos_cells.pdf', height = 15, width = 15)

# gene expression of AIMpos and clonal expanded cells 
obj_sub_sub <- subset(obj_sub, subset = unique_index %in% target$unique_index)
tmp <- target$clone_id
names(tmp) <- target$unique_index
obj_sub_sub$clone_id <- tmp

DoHeatmap(obj_sub_sub, feature=gene_list,
          group.by = 'clone_id') +
  theme(legend.position = "none")
ggsave('figures/feature_expression_AIMpos_expanded_cells.pdf', height = 15, width = 15)

# freq of AIMpos and clonal expanded by sample
freq_barplot_AIMpos_and_expanded_by_sample<- function(d, type){
  # type can be either AIMpos or AIMpos_expanded
  d_stats <- d %>% group_by(Sample_Name, !!sym(type)) %>%
    summarise(count = n()) %>%
    ungroup() %>%
    tidyr::complete(Sample_Name, !!sym(type), fill = list(count=0)) %>%
    group_by(Sample_Name) %>%
    mutate(freq = round(count/sum(count), 3))
  
  g <- ggplot(d_stats %>% filter(!!sym(type)), aes(x=Sample_Name, y=freq, label=count)) +
    geom_bar(stat = 'identity') + geom_text(position = position_stack(vjust = 0.5)) + 
    labs(title = type) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  return(g)
}

get_AIMpos_expanded_stats <- function(all_meta, AIMpos_meta, AIMpos_expanded_meta, expanded_meta){
  d <- all_meta %>% 
    left_join(AIMpos_meta %>% 
                select(unique_index) %>% 
                tibble::add_column(AIMpos = TRUE), by='unique_index') %>%
    left_join(AIMpos_expanded_meta %>%
                select(unique_index) %>%
                tibble::add_column(AIMpos_expanded = TRUE), by='unique_index') %>%
    left_join(expanded_meta %>%
                select(unique_index) %>%
                tibble::add_column(expanded =TRUE), by='unique_index') %>%
    select(unique_index, AIMpos, expanded, AIMpos_expanded, Sample_Name, seurat_clusters) %>%
    tidyr::replace_na(list(AIMpos = FALSE, AIMpos_expanded=FALSE))
  return(d)
}
expanded_meta <- sc_data %>% filter(clone_id %in% (clone_exp %>% filter(clone_count>1))$clone_id)
d1 <- get_AIMpos_expanded_stats(obj[[]], obj_sub[[]], obj_sub_sub[[]], expanded_meta)
d2 <- get_AIMpos_expanded_stats(sc_data, obj_sub[[]], obj_sub_sub[[]], expanded_meta)
g1 <- freq_barplot_AIMpos_and_expanded_by_sample(d1, 'AIMpos') + ylab('frequency among all cel')
g2 <- freq_barplot_AIMpos_and_expanded_by_sample(d2, 'AIMpos_expanded') + ylab('frequency among cells with paired TCR info')
g1+g2
ggsave('figures/frequency changes of AIMpos.pdf')
# for Tregs (seurat cluster 10)
freq_barplot_AIMpos_and_expanded_by_sample(d1 %>% filter(seurat_clusters=='10'), 'AIMpos') + 
  ylab('frequency among all Tregs') + labs(title = 'AIMpos Tregs')
ggsave('figures/AIMpos_Tregs_frequency_by_sample.pdf')



## AIMpos per marker
df_AIMpos <- lapply(AIM_gene_list, function(AIM_marker){
  obj_tmp <- subset(obj, subset = !!sym(AIM_marker)>1)
  tmp <- tibble(unique_index = obj_tmp$unique_index)
  tmp[, AIM_marker] <- 1
  return(tmp)
}) %>% Reduce(f = full_join)
meta_AIMpos <- obj[[]] %>% left_join(df_AIMpos, by='unique_index') %>% 
  mutate(across(all_of(AIM_gene_list),~tidyr::replace_na(.x, 0)))

# for Tregs (cluster 10)
df <- lapply(AIM_gene_list, function(AIM_marker){
  meta_AIMpos %>% filter(seurat_clusters == '10') %>% 
    group_by(Sample_Name, !!sym(AIM_marker)) %>% 
    tally() %>%
    ungroup() %>%
    tidyr::complete(Sample_Name, !!sym(AIM_marker), fill = list(n=0)) %>%
    group_by(Sample_Name) %>%
    mutate(freq = n/sum(n)) %>% 
    filter(!!sym(AIM_marker)==1) %>%
    rename(AIM = !!sym(AIM_marker)) %>%
    mutate(AIM = AIM_marker)
}) %>% do.call(what = rbind)
tmp <- df %>% 
  filter(!(Sample_Name %in% c('BabyAIM_2', 'BabyAIM_4'))) %>%
  left_join(sample_info, by='Sample_Name') %>%
  mutate(label = if_else(Sample_Name == 'BabyAIM_1' | source=='Mother', AIM, NA_character_)) %>%
  filter(AIM != 'IL2RA') # CD25 is a Treg marker so it's not induced by activation
labeled_line_graph(data = tmp, 
                   aes_line = aes(x=timepoint, y=freq, group=AIM, color=AIM, label = label),
                   ref = list(source='Mother'),
                   padding = 0.15) +
  ylab('Frequency among Tregs') + 
  ggtitle('Frequency changes of each AIM')
ggsave('figures/eachAIM_Tregs_frequency_by_sample.pdf')

# for non-Treg CD4 T
df <- lapply(AIM_gene_list, function(AIM_marker){
  meta_AIMpos %>% 
    filter(seurat_clusters != '10',
           cell_type == 'CD4T') %>% 
    group_by(Sample_Name, !!sym(AIM_marker)) %>% 
    tally() %>%
    ungroup() %>%
    tidyr::complete(Sample_Name, !!sym(AIM_marker), fill = list(n=0)) %>%
    group_by(Sample_Name) %>%
    mutate(freq = n/sum(n)) %>% 
    filter(!!sym(AIM_marker)==1) %>%
    rename(AIM = !!sym(AIM_marker)) %>%
    mutate(AIM = AIM_marker)
}) %>% do.call(what = rbind)
tmp <- df %>% 
  filter(!(Sample_Name %in% c('BabyAIM_2', 'BabyAIM_4'))) %>%
  left_join(sample_info, by='Sample_Name') %>%
  mutate(label = if_else(timepoint == min(sample_info$timepoint) | source=='Mother', AIM, NA_character_))
labeled_line_graph(data = tmp, 
                   aes_line = aes(x=timepoint, y=freq, group=AIM, color=AIM, label = label),
                   ref = list(source='Mother'),
                   padding = 0.15) +
  ylab('Frequency among non-Treg CD4T cells') + 
  ggtitle('Frequency changes of each AIM')
ggsave('figures/eachAIM_nonTregsCD4_frequency_by_sample.pdf')











# 
# ##
# 
# # top clone phenotype
# top_clone_id <- get_top_expansion_id(clone_exp, 10)
# data <- sc_data %>% filter(clone_id %in% top_clone_id)
# #obj <- ScaleData(obj)
# 
# obj_exp <- subset(obj, subset = unique_index %in% data$unique_index)
# obj_exp <- ScaleData(obj_exp)
# 
# FeaturePlot(obj_exp, features = c('CD3E', 'CD4', 'CD8A', 'ITGAE'))
# 
# tmp <- obj_exp[[]] %>% left_join(data, by = 'unique_index')
# obj_exp$clone_id <- tmp$clone_id
# Idents(obj_exp) <- 'clone_id'
# all_features <- FindAllMarkers(obj_exp)
# top_features <- all_features %>% 
#   group_by(cluster) %>% 
#   slice_max(n=10, order_by = avg_log2FC)
# DimPlot(obj_exp, group.by = 'clone_id')
# 
# tmp <- GetAssayData(obj_exp, slot = 'scale.data')
# d <- hclust(dist(tmp[unique(top_features$gene),], method = 'euclidean'), method = 'average')
# DoHeatmap(obj_exp, group.by = 'clone_id', features = c(d$labels[d$order]))
# ggsave('figures/top_clone_phenotype.pdf', width = 15, height = 10, units = 'in')
# 
# # search for CD8 T with high CD39 and CD11b and their clone expansions
# 
# top_clone_id <- get_top_expansion_id(clone_exp %>% filter(grepl('ISAC99', Sample_Name)), 20)
# data <- sc_data %>% filter(clone_id %in% top_clone_id)
# obj_exp <- subset(obj, subset = unique_index %in% data$unique_index)
# FeaturePlot(obj_exp, features = c('CD3E', 'CD8A', 'ENTPD1', 'ITGAM'))
# FeaturePlot(obj, features = c('CD3E', 'CD8A', 'ENTPD1', 'ITGAM'))
# ggsave('figures/all_cells_CD39_CD11b.pdf')
# 
# tmp <- obj_exp[[]] %>% left_join(data, by = 'unique_index')
# obj_exp$clone_id <- tmp$clone_id
# Idents(obj_exp) <- 'clone_id'
# DoHeatmap(obj_exp, slot='data', group.by = 'clone_id', features = c('CD3E', 'CD8A', 'ENTPD1', 'ITGAM')) + 
#   scale_fill_gradientn(colors = c("black", "yellow"))
# ggsave('figures/ISAC99_CD39_CD11b_top_clones.pdf', width = 15, height = 5, units = 'in')







