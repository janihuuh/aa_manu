
# get data
df <- readRDS("data/OLINK/data.rds")
cytokine_df <- df %>% dplyr::select(OlinkID,Assay)
cytokine_df <- cytokine_df[!duplicated(cytokine_df), ]

# make df for PCA analysis
df_temp <- df %>%
  group_by(OlinkID) %>%
  mutate(assay_var = var(NPX, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(!(assay_var == 0 | is.na(assay_var))) %>%
  dplyr::select(-assay_var) %>%
  setDT()
df_wide <- tidyr::spread(df_temp[!is.na(NPX),c("SampleID", "OlinkID", "NPX")], key = OlinkID, value = NPX)

# count missing values percentages for all assays
percent_missingness <- colSums(is.na(df_wide[, -c(1:2)]))/nrow(df_wide)

#u se median imputing for assays with less than 5% missing values
PERCENT_CUTOFF <- 0.05
if(any(percent_missingness <= PERCENT_CUTOFF & percent_missingness > 0)){
  imputed_assays_index <- which(percent_missingness <= PERCENT_CUTOFF & percent_missingness > 0)
  percent_missingness <- percent_missingness[-imputed_assays_index]
  imputed_assays_index <- imputed_assays_index + 2
  imputed_assays <- colnames(df_wide)[imputed_assays_index]
  
  df_wide <- df_wide %>%
    mutate_at(tidyselect::all_of(imputed_assays),
              ~ifelse(is.na(.x), median(.x, na.rm = TRUE), .x))
  
}

# remove cytokines with more than 5% missing values
if(any(percent_missingness > 0.05)) {
  removed_assays_index <- which(percent_missingness < 0.05)
  percent_missingness <- percent_missingness[-removed_assays_index]
  removed_assays_index <- removed_assays_index + 2
  removed_assays <- colnames(df_wide)[removed_assays_index]
  df_wide <- df_wide[, -removed_assays_index]
}

df_wide_matrix <- df_wide %>%
  tibble::column_to_rownames('SampleID') %>%
  as.matrix()




#### Perform PCA and UMAP and kmeans
df_scaled  <- data.frame(scale(df_wide_matrix))
pca_df     <- prcomp(df_scaled, scale. = FALSE, center = TRUE)
pca_x_df   <- pca_df$x[,1:5] %>% as.data.frame()
umap_df    <- uwot::umap(pca_x_df) %>% as.data.frame()
colnames(umap_df) <- c("UMAP1", "UMAP2")

pca_imporantce_df <- data.frame(PC = paste0("PC", 1:76), importance = round(pca_df$sdev / sum(pca_df$sdev) * 100, 2))

pca_rot_df <- pca_df$rotation[,1:5] %>% as.data.frame() %>% add_rownames(var = "OlinkID") %>% left_join(cytokine_df, by = "OlinkID")
pca_rot_df <- pca_rot_df[!duplicated(pca_rot_df),]

## Pool all data
dimred_df <- bind_cols(pca_x_df, umap_df) %>% mutate(SampleID = rownames(pca_x_df))
dimred_df <- dimred_df %>% left_join(sample_info) %>% left_join(km_df)

## Cluster
Rphenograph_out              <- Rphenograph::Rphenograph(umap_df, k = 15)
dimred_df$phenograph_cluster <- factor(igraph::membership(Rphenograph_out[[2]]))
dimred_df$phenograph_cluster <- ifelse(dimred_df$phenograph_cluster %in% c(3,4), "3", dimred_df$phenograph_cluster)

fwrite(dimred_df, "results/OLINK/unsupervised/clusters_pca_umap.txt", sep = "\t", quote = F, row.names = F)

## Calcluate DE cytokines
cytokine_test_df <- df_scaled %>%
  add_rownames(var = "SampleID") %>% melt(id = "SampleID") %>% 
  left_join(cytokine_df, by = c("variable" = "OlinkID")) %>% 
  left_join(dplyr::select(dimred_df, SampleID, cluster), by = "SampleID") 

cytokine_test_df <- df_temp %>%
  left_join(dplyr::select(dimred_df, SampleID, cluster), by = "SampleID") 

krskl.df <- cytokine_test_df %>% dplyr::select(SampleID,Assay,NPX,cluster)
krskl.df <- krskl.df[!duplicated(krskl.df), ]

krskl.p.df <- lapply(unique(krskl.df$Assay), FUN = function(x){
  message(x)
  y <- krskl.df %>% filter(Assay == x)
  if(length(unique(y$cluster)) == 3){
    kruskal.test(NPX~cluster, data = y) %>% broom::tidy() %>% mutate(variable = x)
  }
}) %>% rbindlist() %>% mutate(p.adj = p.adjust(p.value, method = "BH")) %>% arrange(p.adj)# %>% left_join(var_df)
fwrite(krskl.p.df, "results/OLINK/unsupervised/cluster_krskl_p.txt", sep = "\t", quote = F, row.names = F)

wlcx.aa.like.p.df <- lapply(unique(krskl.df$Assay), FUN = function(x){
  message(x)
  y <- krskl.df %>% filter(Assay == x) %>% mutate(cluster = ifelse(cluster == "1 Healthy", "Healthy", "AA-like cytokine"))
  if(length(unique(y$cluster)) == 2){
    wilcox.test(NPX~cluster, data = y) %>% broom::tidy() %>% mutate(variable = x)
  }
}) %>% rbindlist() %>% mutate(p.adj = p.adjust(p.value, method = "BH")) %>% arrange(p.adj) # %>% left_join(var_df)
fwrite(wlcx.aa.like.p.df, "results/OLINK/unsupervised/wlcx.aa.like.p.df.txt", sep = "\t", quote = F, row.names = F)
wlcx.aa.like.sigf.p.df <- wlcx.aa.like.p.df %>% filter(p.adj < 0.05)


wlcx.aa.cl2v3.p.df <- lapply(unique(krskl.df$Assay), FUN = function(x){
  
  message(x)
  y <- krskl.df %>% filter(Assay == x) %>% filter(cluster != "1 Healthy")
  
  if(length(unique(y$cluster)) == 2){
    
    wilcox.test(NPX~cluster, data = y) %>% broom::tidy() %>% mutate(variable = x) %>% 
      mutate(median.x = median(subset(y, cluster == "2 AA")$NPX), median.y = median(subset(y, cluster != "2 AA")$NPX),
             mean.x   = mean(subset(y, cluster == "2 AA")$NPX), mean.y   = mean(subset(y, cluster != "2 AA")$NPX), 
             log2fc.median = log2(median.y/median.x), 
             log2fc.mean = log2(mean.y/mean.x))
  }
}) %>% rbindlist() %>% mutate(p.adj = p.adjust(p.value, method = "BH")) %>% arrange(p.adj) # %>% left_join(var_df)

colnames(wlcx.aa.cl2v3.p.df) <- gsub(pattern = "\\.x", ".2.AA", colnames(wlcx.aa.cl2v3.p.df))
colnames(wlcx.aa.cl2v3.p.df) <- gsub(pattern = "\\.y", ".3.AA", colnames(wlcx.aa.cl2v3.p.df))
fwrite(wlcx.aa.cl2v3.p.df, "results/OLINK/unsupervised/wlcx.aa.cl2v3.p.df.txt", sep = "\t", quote = F, row.names = F)
