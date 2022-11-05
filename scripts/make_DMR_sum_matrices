### Make matrices with the sum of (i) methylated and (ii) total counts in DMR per sample

## Make file of all sites in each DMR 
# "dmrs2k_age_regionOnly.bed" is a bed file of the chr, start, stop, and region of each DMR
#bedtools intersect -a /path/all_sites.bed -b /path/dmrs2k_age_regionOnly.bed -wb > /path/sites_in_dmrs2k_age.bed


## Load in site-level counts matrices ("dnam_count" and dnam_totalcount)


## If needed, read in file with DMRs and row for each site in each DMR (made above)
#sites_age_dmrs <- read.delim("/path/sites_in_dmrs2k_age.bed", 
#                             header = F) %>% dplyr::select(!"V4")
#colnames(sites_age_dmrs) <- c("chr", "site_start","site_end","dmr_start","dmr_end","dmr_region")
#sites_age_dmrs <- sites_age_dmrs %>% 
#  mutate(site = paste0(chr,"_",site_start))

## Summarize methylated counts in age DMRs - list sites in each DMR, sum counts per individual across the sites (ie. 3 sites in a DMR would go to 1 row)
DMR_counts_age <- do.call(rbind, lapply(unique(sites_age_dmrs$dmr_region), function(x) { 
 dmr_sites <- sites_age_dmrs[sites_age_dmrs$dmr_region == x, ]$site
  tmp <- as.data.frame(t(apply(dnam_count[rownames(dnam_count) %in% dmr_sites,], MARGIN = 2, sum)))
  rownames(tmp) <- x
  return(tmp)
  }))

## Summarize total counts in age DMRs regions 
DMR_totalcounts_age <- do.call(rbind, lapply(unique(sites_age_dmrs$dmr_region), function(x) { 
  dmr_sites <- sites_age_dmrs[sites_age_dmrs$dmr_region == x, ]$site
  tmp <- as.data.frame(t(apply(dnam_totalcount[rownames(dnam_totalcount) %in% dmr_sites,], MARGIN = 2, sum)))
  rownames(tmp) <- x
  return(tmp)
}))
