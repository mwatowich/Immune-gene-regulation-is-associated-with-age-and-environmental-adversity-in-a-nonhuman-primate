### Longitudinal vs cross-sectional samples

DMR_prop <- DMR_counts_age/DMR_totalcounts_age
dnam_meta_repeats <- dnam_meta[dnam_meta$animal_id %in% dnam_meta[duplicated(dnam_meta$animal_id),]$animal_id,]
# 139 repeated samples, 67 unique individuals

# For each DMR, get the proportion methylated for each animal that has a repeated sample 
rep_meth <- do.call(rbind, lapply(rownames(DMR_prop), function(y) { 
  do.call(rbind, lapply(unique(dnam_meta_repeats$animal_id), function(x) {
    lids_tmp <- dnam_meta_repeats[dnam_meta_repeats$animal_id == x,]$updated_lid 
    prop_tmp <- as.data.frame(t(DMR_prop[rownames(DMR_prop) == y, 
                                         colnames(DMR_prop) %in% lids_tmp])) 
    colnames(prop_tmp) <- "prop"
    prop_tmp$dmr <- y
    prop_tmp$updated_lid <- rownames(prop_tmp)
    return(prop_tmp)
  }))
}))
# Join with metadata to get animal's age to find which LID is old/young sample
rep_meth_2 <- left_join(rep_meth, dnam_meta[,c("updated_lid","age","animal_id")], by="updated_lid")

# For each animal, get the oldest and youngest sample. Subtract old - young proportion 
rep_meth_3 <- left_join(rep_meth_2 %>% 
                    group_by(animal_id) %>% 
                    slice_max(age, n = 1) %>% 
                    mutate(prop_older = prop) %>% 
                    dplyr::select(prop_older, dmr, animal_id),
                  rep_meth_2 %>% 
                    group_by(animal_id) %>% 
                    slice_min(age, n = 1)  %>% 
                    mutate(prop_younger = prop) %>% 
                    dplyr::select(prop_younger, dmr, animal_id), 
                  by = c("dmr","animal_id")) %>% 
  mutate(prop_diff = prop_older-prop_younger)

# Remove DMRs with >75% of animals missing data 
na_check <- rep_meth_3[,c("dmr","prop_diff","animal_id")] %>% pivot_wider(names_from = animal_id, values_from = prop_diff)
na_check$n_NAs <- apply(na_check, 1, function(x) sum(is.na(x)))
rep_meth_3 <- rep_meth_3[!rep_meth_3$dmr %in% na_check[na_check$n_NAs > 50,]$dmr,]

# 
rep_meth_4 <- do.call(rbind, lapply(unique(rep_meth_3$dmr), function(x) { 
  # Proportion per site that are hypermethylated with age in the same individual 
  tmp<-data.frame(n_hypermeth = dim(rep_meth_3[rep_meth_3$dmr == x & complete.cases(rep_meth_3$prop_diff) & rep_meth_3$prop_diff>0,])[1]/
                    dim(rep_meth_3[rep_meth_3$dmr == x & complete.cases(rep_meth_3$prop_diff),])[1], 
                  dmr = x)
  # Median proportion difference for each DMR 
  tmp$median_prop_diff <- median(rep_meth_3[rep_meth_3$dmr == x & complete.cases(rep_meth_3$prop_diff),]$prop_diff)
  return(tmp)
}))
rep_meth_5 <- left_join(rep_meth_4, 
                  age_dmrs,
                  by = c("dmr" = "region"))

# Does the % hypermethylated (#animals per DMR) correlate with direction from cross-sectional model? 
cor.test(rep_meth_5$n_hypermeth, rep_meth_5$median_effect)
plot(rep_meth_5$n_hypermeth,rep_meth_5$median_effect)
# Does the median proportion methylated (% hypermeth) correlate with the std beta for age in cross-sectional model? 
cor.test(rep_meth_5$median_prop_diff,rep_meth_5$median_effect)
plot(rep_meth_5$median_prop_diff,rep_meth_5$median_effect)
