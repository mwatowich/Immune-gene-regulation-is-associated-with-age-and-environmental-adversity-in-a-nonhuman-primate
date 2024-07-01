### Make matrices with the sum of (i) methylated and (ii) total counts in DMR per sample

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
