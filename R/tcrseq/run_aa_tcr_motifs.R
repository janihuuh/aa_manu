## Read GLIPH2 files
aa_res       <- readGliphFile("src/bash/gliph2/oct21/local_folder/aa_df_cluster.csv") %>% filter(number_subject >= 3 & number_unique_cdr3 >= 3 & vb_score < 0.05)
itp_res      <- readGliphFile("src/bash/gliph2/oct21/local_folder3/itp_df_cluster.csv") %>% filter(number_subject >= 3 & number_unique_cdr3 >= 3 & vb_score < 0.05)
ra_res       <- readGliphFile("src/bash/gliph2/oct21/local_folder3/ra_df_cluster.csv") %>% filter(number_subject >= 3 & number_unique_cdr3 >= 3 & vb_score < 0.05)
mds_res      <- readGliphFile("src/bash/gliph2/oct21/local_folder3/mds_df_cluster.csv") %>% filter(number_subject >= 3 & number_unique_cdr3 >= 3 & vb_score < 0.05)
healthy_res  <- readGliphFile("src/bash/gliph2/oct21/ctrl_df_cluster.csv") %>% filter(number_subject >= 3 & number_unique_cdr3 >= 3 & vb_score < 0.05)
healthy4_res <- readGliphFile("src/bash/gliph2/oct21/ctrl4_df_cluster.csv") %>% filter(number_subject >= 3 & number_unique_cdr3 >= 3 & vb_score < 0.05)

## 
gliph_aa_res1  <- readGliphFile("src/bash/gliph2/oct21/aa_sample_df_seed1_cluster.csv") %>% mutate(seed = "1") %>% filter(number_subject >= 3 & number_unique_cdr3 >= 3 & vb_score < 0.05) %>% group_by(pattern) %>% summarise(n=n())
gliph_aa_res2  <- readGliphFile("src/bash/gliph2/oct21/aa_sample_df_seed2_cluster.csv") %>% mutate(seed = "2") %>% filter(number_subject >= 3 & number_unique_cdr3 >= 3 & vb_score < 0.05) %>% group_by(pattern) %>% summarise(n=n())
gliph_aa_res3  <- readGliphFile("src/bash/gliph2/oct21/aa_sample_df_seed3_cluster.csv") %>% mutate(seed = "3") %>% filter(number_subject >= 3 & number_unique_cdr3 >= 3 & vb_score < 0.05) %>% group_by(pattern) %>% summarise(n=n())
gliph_aa_res4  <- readGliphFile("src/bash/gliph2/oct21/aa_sample_df_seed4_cluster.csv") %>% mutate(seed = "4") %>% filter(number_subject >= 3 & number_unique_cdr3 >= 3 & vb_score < 0.05) %>% group_by(pattern) %>% summarise(n=n())
gliph_aa_res5  <- readGliphFile("src/bash/gliph2/oct21/aa_sample_df_seed5_cluster.csv") %>% mutate(seed = "5") %>% filter(number_subject >= 3 & number_unique_cdr3 >= 3 & vb_score < 0.05) %>% group_by(pattern) %>% summarise(n=n())
gliph_aa_res6  <- readGliphFile("src/bash/gliph2/oct21/aa_sample_df_seed6_cluster.csv") %>% mutate(seed = "6") %>% filter(number_subject >= 3 & number_unique_cdr3 >= 3 & vb_score < 0.05) %>% group_by(pattern) %>% summarise(n=n())
gliph_aa_res7  <- readGliphFile("src/bash/gliph2/oct21/aa_sample_df_seed7_cluster.csv") %>% mutate(seed = "7") %>% filter(number_subject >= 3 & number_unique_cdr3 >= 3 & vb_score < 0.05) %>% group_by(pattern) %>% summarise(n=n())
gliph_aa_res8  <- readGliphFile("src/bash/gliph2/oct21/aa_sample_df_seed8_cluster.csv") %>% mutate(seed = "8") %>% filter(number_subject >= 3 & number_unique_cdr3 >= 3 & vb_score < 0.05) %>% group_by(pattern) %>% summarise(n=n())
gliph_aa_res9  <- readGliphFile("src/bash/gliph2/oct21/aa_sample_df_seed9_cluster.csv") %>% mutate(seed = "9") %>% filter(number_subject >= 3 & number_unique_cdr3 >= 3 & vb_score < 0.05) %>% group_by(pattern) %>% summarise(n=n())
gliph_aa_res10 <- readGliphFile("src/bash/gliph2/oct21/aa_sample_df_seed10_cluster.csv") %>% mutate(seed = "10") %>% filter(number_subject >= 3 & number_unique_cdr3 >= 3 & vb_score < 0.05) %>% group_by(pattern) %>% summarise(n=n())

gliph_full <- gliph_aa_res1 %>% 
  full_join(gliph_aa_res2, by = "pattern") %>% 
  full_join(gliph_aa_res2, by = "pattern") %>% 
  full_join(gliph_aa_res3, by = "pattern") %>% 
  full_join(gliph_aa_res4, by = "pattern") %>% 
  # full_join(gliph_aa_res5, by = "pattern") %>% 
  full_join(gliph_aa_res6, by = "pattern") %>% 
  full_join(gliph_aa_res7, by = "pattern") %>% 
  full_join(gliph_aa_res8, by = "pattern") %>% 
  full_join(gliph_aa_res9, by = "pattern") %>% 
  full_join(gliph_aa_res10, by = "pattern")

gliph_full[is.na(gliph_full)] <- 0


#### Select AA motifs seen in at least 7 CV folds
df                     <- apply(gliph_full[,-1], 1, FUN = function(x) table(is.na(x)) %>% as.data.frame()) %>% rbindlist() %>% filter(Var1 == F)
aa_gliph_robust_motifs <- gliph_full[df$Freq >= 7, ] %>% pull(pattern)
aa_robust_res          <- gliph_aa_full %>% filter(pattern %in% aa_gliph_robust_motifs)


#### Select AA motifs not seen in healthy, RA, ITP, or MDS i) more than 5 people ii) more than 5 TCRs
aa_patterns        <- aa_robust_res %>% group_by(pattern, aa_number_subject = number_subject, aa_number_unique_cdr3 = number_unique_cdr3) %>% summarise(aa_max = max(Freq))
itp_patterns       <- itp_res %>% group_by(pattern, itp_number_subject = number_subject, itp_number_unique_cdr3 = number_unique_cdr3) %>% summarise(itp_max = max(Freq))
ra_patterns        <- ra_res %>% group_by(pattern, ra_number_subject = number_subject, ra_number_unique_cdr3 = number_unique_cdr3) %>% summarise(ra_max = max(Freq))
mds_patterns       <- mds_res %>% group_by(pattern, mds_number_subject = number_subject, mds_number_unique_cdr3 = number_unique_cdr3) %>% summarise(mds_max = max(Freq))
healthy_patterns   <- healthy_res %>% group_by(pattern, healthy_number_subject = number_subject, healthy_number_unique_cdr3 = number_unique_cdr3) %>% summarise(healthy_max = max(Freq))
healthyc4_patterns <- healthy4_res %>% group_by(pattern, healthyc4_number_subject = number_subject, healthyc4_number_unique_cdr3 = number_unique_cdr3) %>% summarise(healthyc4_max = max(Freq))

total_patterns <- aa_patterns %>% 
  left_join(itp_patterns, by = "pattern") %>% 
  left_join(ra_patterns, by = "pattern") %>% 
  left_join(mds_patterns, by = "pattern") %>% 
  left_join(healthy_patterns, by = "pattern") %>% 
  left_join(healthyc4_patterns, by = "pattern")

## Patterns must be at least three amino acids long
total_patterns[is.na(total_patterns)] <- 0
total_patterns <- total_patterns[nchar(total_patterns$pattern) > 3, ]

total_patterns <- total_patterns %>% 
  
  filter(itp_number_subject < 3) %>%
  filter(ra_number_subject < 3) %>%
  filter(mds_number_subject < 3) %>%
  filter(healthy_number_subject < 1) %>%
  filter(healthyc4_number_subject < 3) %>%
  
  filter(itp_number_unique_cdr3 < 5) %>% 
  filter(ra_number_unique_cdr3 < 5) %>% 
  filter(mds_number_unique_cdr3 < 5) %>% 
  filter(healthyc4_number_subject < 5) %>%
  filter(healthy_number_unique_cdr3 < 1) %>% 
  
  filter(itp_max < 0.001) %>% 
  filter(ra_max < 0.001) %>% 
  filter(mds_max < 0.001) %>% 
  filter(healthy_max < 0.001) %>% 
  filter(healthyc4_max < 0.001) %>% 
  
  filter(aa_max > 0.001) 
  

## Prune
aa_robust_res      <- aa_robust_res %>% filter(pattern %in% total_patterns$pattern)
aa_robust_res_uniq <- aa_robust_res[!duplicated(aa_robust_res$pattern), ]
fwrite(aa_robust_res_uniq, "/supplementary_tables/aa_robust_res_uniq.txt", sep = "\t", quote = F, row.names = F)
