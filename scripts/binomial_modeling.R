### CpG-site level binomial models

library(doParallel)
library(parallel)
library(Rcpp)
library(MuMIn)
library(Matrix)

dnam_out <- as.data.frame(do.call(rbind, lapply(X = row_start:row_end, FUN = function(x){
  dnam_count_new <- dnam_count[x, colSums(is.na(dnam_count[x,])) == 0]
  dnam_totalcount_new <- dnam_totalcount[x, colnames(dnam_count_new)]
  dnam_mod_mat <- as.data.frame(model.matrix(~ scale(sex) + scale(age) + scale(hurricane) + scale(batch), dnam_meta))
  model_DNA = pqlseq(RawCountDataSet=dnam_count_new, 
                         Covariates=dnam_mod_mat[,-2],
                         Phenotypes=dnam_mod_mat[,2],
                         RelatednessMatrix=dnam_kinship_new,
                         LibSize=dnam_totalcount_new,
                         fit.model="BMM")
    for (i in c(1:ncol(dnam_mod_mat[,-1]))) {
      colnames(model_DNA) <- gsub(pattern = paste0("covariates", i), 
                                  replacement = paste0(colnames(dnam_mod_mat)[-c(1:2)][i]),
                                  x = colnames(model_DNA))
    }
    colnames(model_DNA) <- gsub(pattern = "predictor", 
                                replacement = paste0(colnames(dnam_mod_mat)[2]), 
                                x = colnames(model_DNA))
    rownames(model_DNA) <- rownames(dnam_count_new)
    return(model_DNA)
})))
