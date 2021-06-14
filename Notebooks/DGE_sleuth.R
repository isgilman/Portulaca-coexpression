library("sleuth")
library("statmod")
library("stringr")
library("dplyr")
library("tidyr")
library("reshape2")
library("ggplot2")
library("splines")

#####################################
###~~~~~~~~~~ P. amilis ~~~~~~~~~~###
#####################################

kallisto_amilis_dir <- "~/Dropbox/Drought-experiments/Drought_Experiment_3/farnam-drought-diff-exp/Kallisto/Portulaca-amilis"

# Metadata
Pamilis_metadata = read.csv("~/Dropbox/Drought-experiments/Drought_Experiment_3/farnam-drought-diff-exp/P-amilis.metadata.csv")
Pamilis_metadata =  dplyr::mutate(Pamilis_metadata, path=file.path(kallisto_amilis_dir, sequencer_name, 'abundance.h5'))
exp_design = Pamilis_metadata[, c("sequencer_name", 'time_point', 'treatment', 'path')]
exp_design$time_point = exp_design$time_point %% 4
colnames(exp_design) = c('sample', 'time_point', 'treatment', 'path')

# Transcript-gene map
ttg = read.csv("~/Dropbox/Drought-experiments/Drought_Experiment_3/farnam-drought-diff-exp/transcript-gene-map.txt", sep = ' ')
colnames(ttg) = c('target_id', 'gene_id')

# Preparing the analysis
so = sleuth_prep(sample_to_covariates = exp_design, target_mapping = ttg, aggregation_column = 'gene_id',  extra_bootstrap_summary=TRUE, gene_mode=FALSE)
## Using gene_mode=TRUE to get gene-level TPM
so_gm = sleuth_prep(sample_to_covariates = exp_design, target_mapping = ttg, aggregation_column = 'gene_id', extra_bootstrap_summary=TRUE, gene_mode=TRUE)

so_norm_matrix = sleuth_to_matrix(so, which_df = "obs_norm", which_units = "est_counts")
so_gm_norm_matrix = sleuth_to_matrix(so_gm, which_df = "obs_norm", which_units = "scaled_reads_per_base")
# Save gene mode abundance data for other analyses
write.csv(so_gm_norm_matrix, "~/Dropbox/Drought-experiments/Drought_Experiment_3/Abundance-analysis/sleuth/P-amils-sleuth-gene-norm-filt.csv")

# The analysis
## Time only
so = sleuth_fit(obj = so, formula = ~time_point, fit_name = 'reduced')
so_gm = sleuth_fit(obj = so, formula = ~time_point, fit_name = 'reduced', gene_mode=TRUE)
## Full model
so = sleuth_fit(obj = so, formula = ~time_point + treatment, fit_name = 'full')
so_gm = sleuth_fit(obj = so, formula = ~time_point + treatment, fit_name = 'full', gene_mode=TRUE)
## Compare
so = sleuth_lrt(obj = so, null_model = 'reduced', alt_model = 'full')
so_gm = sleuth_lrt(obj = so, null_model = 'reduced', alt_model = 'full')

# Save objects
save(so, file="~/Dropbox/Drought-experiments/Drought_Experiment_3/Abundance-analysis/sleuth/P-amilis-sleuth-object.RData")
save(so_gm, file="~/Dropbox/Drought-experiments/Drought_Experiment_3/Abundance-analysis/sleuth/P-amilis-sleuth-object.RData")
load(file="~/Dropbox/Drought-experiments/Drought_Experiment_3/Abundance-analysis/sleuth/P-amilis-sleuth-object.RData")
load(file="~/Dropbox/Drought-experiments/Drought_Experiment_3/Abundance-analysis/sleuth/P-amilis-sleuth-object.RData")

# Obtaining gene-level DE results
sleuth_table_gene = sleuth_results(obj = so, test ='reduced:full', test_type = 'lrt', show_all = FALSE, pval_aggregate=TRUE)
sleuth_table_gene_gm = sleuth_results(obj = so_gm, test ='reduced:full', test_type = 'lrt', show_all = FALSE, pval_aggregate=TRUE)

# Save gene mode differential abundance statistics
write.csv(sleuth_table_gene, "~/Dropbox/Drought-experiments/Drought_Experiment_3/Abundance-analysis/sleuth/P-amilis-sleuth-DEgenes.csv")
# Look at significant genes
sleuth_table_gene = dplyr::filter(sleuth_table_gene, qval <= 0.01)

#####################################
##~~~~~~~~~~ P. oleracea ~~~~~~~~~~##
#####################################
kallisto_oleracea_dir <- "~/Dropbox/Drought-experiments/Drought_Experiment_3/farnam-drought-diff-exp/Kallisto/Portulaca-oleracea/no-filter-cluster98/"

# Metadata
Poleracea_metadata = read.csv("~/Dropbox/Drought-experiments/Drought_Experiment_3/farnam-drought-diff-exp/P-oleracea.metadata.csv")
Poleracea_metadata =  dplyr::mutate(Poleracea_metadata, path=file.path(kallisto_oleracea_dir, sequencer_name, 'abundance.h5'))
exp_design = Poleracea_metadata[, c("sequencer_name", 'time_point', 'treatment', 'path')]
exp_design$time_point = exp_design$time_point %% 4
colnames(exp_design) = c('sample', 'time_point', 'treatment', 'path')

# Transcript-gene map
ttg = read.csv("~/Dropbox/Drought-experiments/Drought_Experiment_3/farnam-drought-diff-exp/Kallisto/Portulaca-oleracea/no-filter-cluster98/Trinity.fasta.gene_trans_map", sep = '\t')
colnames(ttg) = c('gene_id', 'target_id')

# Preparing the analysis
so = sleuth_prep(sample_to_covariates = exp_design, target_mapping = ttg, aggregation_column = 'gene_id', extra_bootstrap_summary=TRUE, gene_mode=FALSE)
## Using gene_mode=TRUE to get gene-level TPM
so_gm = sleuth_prep(sample_to_covariates = exp_design, target_mapping = ttg, aggregation_column = 'gene_id', extra_bootstrap_summary=TRUE, gene_mode=TRUE)

so_norm_matrix = sleuth_to_matrix(so, which_df = "obs_norm", which_units = "est_counts")
so_gm_norm_matrix = sleuth_to_matrix(so_gm, which_df = "obs_norm", which_units = "scaled_reads_per_base")
# Save gene mode abundance data for other analyses
write.csv(so_gm_norm_matrix, "~/Dropbox/Drought-experiments/Drought_Experiment_3/Abundance-analysis/sleuth/P-oleracea-sleuth-gene-norm-filt.csv")

# The analysis
## Time only
so = sleuth_fit(obj = so, formula = ~time_point, fit_name = 'reduced')
so_gm = sleuth_fit(obj = so, formula = ~time_point, fit_name = 'reduced', gene_mode=TRUE)
## Full model
so = sleuth_fit(obj = so, formula = ~time_point + treatment, fit_name = 'full')
so_gm = sleuth_fit(obj = so, formula = ~time_point + treatment, fit_name = 'full', gene_mode=TRUE)
## Compare
so = sleuth_lrt(obj = so, null_model = 'reduced', alt_model = 'full')
so_gm = sleuth_lrt(obj = so, null_model = 'reduced', alt_model = 'full')

# Save objects
save(so, file="~/Dropbox/Drought-experiments/Drought_Experiment_3/Abundance-analysis/sleuth/P-amilis-oleracea-object.RData")
save(so_gm, file="~/Dropbox/Drought-experiments/Drought_Experiment_3/Abundance-analysis/sleuth/P-oleracea-sleuth-object.RData")
load(file="~/Dropbox/Drought-experiments/Drought_Experiment_3/Abundance-analysis/sleuth/P-oleracea-sleuth-object.RData")
load(file="~/Dropbox/Drought-experiments/Drought_Experiment_3/Abundance-analysis/sleuth/P-oleracea-sleuth-object.RData")

# Obtaining gene-level DE results
sleuth_table_gene = sleuth_results(obj = so, test ='reduced:full', test_type = 'lrt', show_all = FALSE, pval_aggregate=TRUE)
sleuth_table_gene_gm = sleuth_results(obj = so_gm, test ='reduced:full', test_type = 'lrt', show_all = FALSE, pval_aggregate=TRUE)

# Save gene mode differential abundance statistics
write.csv(sleuth_table_gene, "~/Dropbox/Drought-experiments/Drought_Experiment_3/Abundance-analysis/sleuth/P-oleracea-sleuth-DEgenes.csv")
# Look at significant genes
sleuth_table_gene = dplyr::filter(sleuth_table_gene, qval <= 0.01)

