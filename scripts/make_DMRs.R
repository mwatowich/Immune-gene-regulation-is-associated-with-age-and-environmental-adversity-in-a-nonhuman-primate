### Code to make DMRs and background DMRs for RRBS data for Watowich 2022

dnam_out <- read.csv("path/to/model/output",header = T,sep = ",")

# Prep DMR df for variable
# get the beta, se, stdbeta, pvalue, and qvalue from the site-level model, as well as the chromosome and start and stop positions of the CpG
dnam_dss_fmt <- dnam_out[,c("beta_age", "se_age","stdbeta_age","pval_age","chr","site","site","qval_age")]
rownames(dnam_dss_fmt) <- paste0(dnam_dss_fmt$chr, "_", dnam_dss_fmt$site)
colnames(dnam_dss_fmt) <- c("beta", "se_beta","stat","raw_pval","chr","start","pos","pval")


### Call DMR function ---------
dmr_function <- function(input_df, fdr_cutoff_strict=0.05, fdr_cutoff_relaxed=0.1, fdr_col_name="pval", effect_size_name="stat", minCG=4, sep=1000, prop_sig=0.5) {
  input_df=input_df[order(input_df$chr,input_df$pos),]
  sig_df=input_df[input_df[,fdr_col_name]<fdr_cutoff_relaxed,]
  pos=sig_df$pos
  pos.diff <- abs(c(as.integer(0), diff(pos)))
  idx.jump <- which(pos.diff>sep|pos.diff<0)
  regions <- rbind(c(1, idx.jump), c(idx.jump-1, length(pos)))
  regions_bed=data.frame(chr = sig_df[regions[1,],"chr"],
                         start = sig_df[regions[1,],"pos"],
                         stop = sig_df[regions[2,],"pos"],
                         n_CG_relaxed = (regions[2,]-regions[1,])+1)
  regions_bed=regions_bed[regions_bed$n_CG_relaxed>=minCG,]
  for (i in 1:nrow(regions_bed)){
    chr=regions_bed[i,"chr"]
    start=regions_bed[i,"start"]
    stop=regions_bed[i,"stop"]
    tmp_df=input_df[input_df$chr==chr & input_df$pos>=start & input_df$pos<=stop,]
    n_CG=nrow(tmp_df)
    n_CG_strict=sum(as.numeric(tmp_df[,fdr_col_name])<fdr_cutoff_strict)
    median_effect=median(tmp_df[,effect_size_name])
    n_concordant=sum(sign(tmp_df[,effect_size_name])==sign(median_effect))
    n_sig_concordant=sum(sign(tmp_df[tmp_df[,fdr_col_name]<fdr_cutoff_relaxed,effect_size_name])==sign(median_effect))
    regions_bed[i,c("n_CG","n_CG_strict","median_effect","n_concordant","n_sig_concordant")]=c(n_CG,n_CG_strict,median_effect,n_concordant,n_sig_concordant)
  }
  regions_bed$p_sig=regions_bed$n_CG_relaxed/regions_bed$n_CG
  regions_bed$p_sig_concordant=regions_bed$n_sig_concordant/regions_bed$n_CG_relaxed
  regions_bed=regions_bed[regions_bed$p_sig>=prop_sig&regions_bed$n_CG_strict>0,]
  regions_bed$length=regions_bed$stop-regions_bed$start+1
  regions_bed
}


### Permutations to establish DMR CpG cutoff --------
# Use minCG=2 in  permutations: focal site and X number of other 10% sites.
perms_out <- do.call(rbind, lapply(X = 1:10, function(x) { 
  set.seed((1+1))
  dnam_dss_fmt_perm <- dnam_dss_fmt
  dnam_dss_fmt_perm$pval <- sample(x = dnam_dss_fmt_perm$pval, size = length(dnam_dss_fmt$pval), replace = F)
  dmrs2k_perm <- dmr_function(input_df=dnam_dss_fmt_perm, 
                              fdr_cutoff_strict=0.05, 
                              fdr_cutoff_relaxed=0.1, 
                              fdr_col_name="pval", 
                              effect_size_name="stat", 
                              minCG=2,
                              sep=1000, prop_sig=0.01)
  out <- data.frame(median_n_CG = median(dmrs2k_perm$n_CG), 
                    median_n_CG_relaxed = median(dmrs2k_perm$n_CG_relaxed), 
                    median_n_CG_strict = median(dmrs2k_perm$n_CG_strict), 
                    perm = x)
  return(out)
}))
perms_out$median_n_CG_relaxed


### Call DMRs --------
dmrs <- dmr_function(input_df=dnam_dss_fmt, 
                           fdr_cutoff_strict=0.05, 
                           fdr_cutoff_relaxed=0.1, 
                           fdr_col_name="pval", 
                           effect_size_name="stat", 
                           minCG=4, sep=1000, prop_sig=0.01)
median(dmrs2k$n_CG_relaxed)
dmrs2k$p_concordant <- (dmrs2k$n_concordant/dmrs2k$n_CG)

# Filter for 75% in same direction 
#hist(dmrs2k$p_sig_concordant,breaks = 50)
#hist(dmrs2k$p_concordant,breaks = 50)
dmrs2k <- dmrs2k[dmrs2k$p_concordant>=0.75 & dmrs2k$p_sig_concordant>=0.75,]

# Check length of DMRs 
hist(dmrs2k$length,breaks=100)

# write an identifier region column
dmrs2k$region <- paste0(dmrs2k$chr, "_", dmrs2k$start, "_", dmrs2k$stop)


### Call background DMRs ----------
dmrs2k_bkgd_tmp <- dmr_function(input_df=dnam_dss_fmt, 
                                    fdr_cutoff_strict=1, 
                                    fdr_cutoff_relaxed=1, 
                                    fdr_col_name="pval", 
                                    effect_size_name="stat", 
                                    minCG=4, sep=1000, prop_sig=0.01) 
# Keep minCG at the same level, and set pvalue to 1 (akin to permuting) 
dmrs2k_bkgd_tmp$p_concordant <- (dmrs2k_bkgd_tmp$n_concordant/dmrs2k_bkgd_tmp$n_CG)

## Use bedr intersect to find/remove the 'real' DMRs from the background set 
library(bedr)
# format chromosome 
dmrs2k$chr <- paste0("chr",dmrs2k$chr)
dmrs2k_bkgd_tmp$chr <- paste0("chr",dmrs2k_bkgd_tmp$chr)
# format region so have this after the merge 
dmrs2k$region <- paste0("chr",dmrs2k$region)
dmrs2k_bkgd_tmp$region_bkgd <- paste0(dmrs2k_bkgd_tmp$chr, "_", 
                                          dmrs2k_bkgd_tmp$start, "_", 
                                          dmrs2k_bkgd_tmp$stop)
# check if valid 
is.valid.region(dmrs2k)
is.valid.region(dmrs2k_bkgd_tmp)

# check if already sorted
is.sorted.region(dmrs2k)
is.sorted.region(dmrs2k_bkgd_tmp)

# sort lexographically
dmrs2k_bkgd_sort <- bedr.sort.region(dmrs2k_bkgd_tmp)

# intersect, then remove the intersected ones from background 
dmrs2k_intersect <- bedr.join.region(dmrs2k, dmrs2k_bkgd_sort)

# remove the 'real' DMRs from bkgd 
head(dmrs2k_intersect,2)
dmrs2k_bkgd <- dmrs2k_bkgd_sort[!dmrs2k_bkgd_sort$region_bkgd %in%
                                          dmrs2k_intersect$region_bkgd,]
dim(dmrs2k_bkgd)

## Format the chr column so it is usable by bedtools
head(dmrs2k_bkgd,2)
dmrs2k_bkgd$chr <- gsub(pattern = "chr",replacement = "",x = dmrs2k_bkgd$chr)
dmrs2k_bkgd$region_bkgd <- gsub(pattern = "chr",replacement = "",x = dmrs2k_bkgd$region_bkgd)
