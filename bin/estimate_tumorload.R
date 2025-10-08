library(data.table)
library(mclust)

#sample_cov_tab <- cov_tab[sample == "386"]
estimate_tumor_ratio_from_hist_one_sample <- function(sample_cov_tab,TL_vec,values_to_consider){
  densityMclust_res <- densityMclust(sample_cov_tab$cov,plot = F)
  
  hist <- hist(sample_cov_tab[cov < mean(cov) * 2 & cov > mean(cov) / 2]$cov,200,plot = F)
  normal_cov <- hist$mids[which.max(hist$counts)]
  TL_mat <- sapply(TL_vec,function(TL) TL * ((1:3)/2 * normal_cov) + normal_cov * (1 - TL))
  
  peak_centers_to_use <- densityMclust_res$parameters$mean[densityMclust_res$parameters$mean < 1.75 * normal_cov &
                                                             densityMclust_res$parameters$mean > 0.25 * normal_cov &
                                                             densityMclust_res$parameters$pro >  mean(densityMclust_res$parameters$pro) - sd(densityMclust_res$parameters$pro) &
                                                             densityMclust_res$parameters$pro <  mean(densityMclust_res$parameters$pro) + sd(densityMclust_res$parameters$pro)]
  
  
  est_scores_for_TL_vec <- apply(TL_mat,2,function(test_TL_vec){
    sum(sapply(peak_centers_to_use,function(peak_cov)  min(abs(test_TL_vec - peak_cov))))
  })
  
  res_tab <- data.table(TL = TL_vec,score = est_scores_for_TL_vec,tab_id = seq_along(TL_vec))
  res_tab[,norm_score := score / TL]
  res_tab <- cbind(res_tab,t(TL_mat))
  setorder(res_tab,score)
  
  return(res_tab)
}

estimate_tumor_ratio_from_hist <- function(cov_tab,TL_vec = seq(0.02,1,by = 0.01),values_to_consider = 1){
  
  est_ratio_score_tab <- cov_tab[,estimate_tumor_ratio_from_hist_one_sample(.SD,TL_vec,values_to_consider),by = sample]

  return(est_ratio_score_tab)
}


###
### my_test_job
###


# rpart_segmentation <- function(tab,minsplit=100, cp=.002){
#   tree <- rpart(cov ~ start, tab ,control=rpart.control(minsplit=minsplit, cp=cp))
#   enc_vec <- rle(predict(tree))
#   return(rep(seq_along(enc_vec$values),enc_vec$lengths))
# }

#sample_chr_cov_tab <- sample_cov_tab[chr == "1"]
estimate_tumor_ratio_from_hist_one_sample_one_chr <- function(sample_chr_cov_tab,TL_vec,values_to_consider,sd_multiplier){
  densityMclust_res <- densityMclust(sample_chr_cov_tab$cov,plot = F)
  
  hist <- hist(sample_chr_cov_tab[cov < mean(cov) * 2]$cov,200,plot = F)
  normal_cov <- hist$mids[which.max(hist$counts)]
  TL_mat <- sapply(TL_vec,function(TL) TL * ((1:3)/2 * normal_cov) + normal_cov * (1 - TL))
  
  peak_centers_to_use <- densityMclust_res$parameters$mean[densityMclust_res$parameters$mean < 1.75 * normal_cov &
                                                             densityMclust_res$parameters$mean > 0.25 * normal_cov &
                                                             densityMclust_res$parameters$pro >  mean(densityMclust_res$parameters$pro) - sd(densityMclust_res$parameters$pro) &
                                                             densityMclust_res$parameters$pro <  mean(densityMclust_res$parameters$pro) + sd(densityMclust_res$parameters$pro)]
  
  
  est_scores_for_TL_vec <- apply(TL_mat,2,function(test_TL_vec){
    sum(sapply(peak_centers_to_use,function(peak_cov)  min(abs(test_TL_vec - peak_cov))))
  })
  
  res_tab <- data.table(TL = TL_vec,score = est_scores_for_TL_vec,tab_id = seq_along(TL_vec))
  res_tab[,norm_score := score / TL]
  res_tab <- cbind(res_tab,t(TL_mat))
  setorder(res_tab,score)
  
}


# ## printing
# pdf("cov_hist_with_ctDNA_ratio_est.pdf",width = 9,height = 9)
# for(select_sample in unique(cov_tab$sample)){
#   est_ratio_score_tab <- est_ratio_score_tab[sample == select_sample]
#   normal_cov <- est_ratio_score_tab$V2[1]
#   TL <- est_tumor_ratio_tab[sample == select_sample]$TL
#   TL_vec <- TL * ((0:4)/2 * normal_cov) + normal_cov * (1 - TL)
#   names(TL_vec) <- 0:4
# 
#   max_hist_value <- median(cov_tab[sample == select_sample]$cov) + 3*sd(cov_tab[sample == select_sample]$cov)
#   min_hist_value <- max(median(cov_tab[sample == select_sample]$cov) - 3*sd(cov_tab[sample == select_sample]$cov),0)
#   hist <- hist(cov_tab[cov > min_hist_value & cov < max_hist_value & sample == select_sample]$cov,400,plot = F)
#   plot(hist,main = select_sample,xlab = "Coverage")
#   abline(v=TL_vec, col="red", lty=2, lwd=1)
#   # abline(v=peak_centers_to_use, col="green", lty=2, lwd=0.5)
#   text(TL_vec + 3*sd(cov_tab[sample == select_sample]$cov) / 20,max(hist$counts) , names(TL_vec),col="red")
# }



# run_all <- function(){
#   
#   
#   tab <- fread("/mnt/ssd/ssd_3/references/homo_sapiens/GRCh38-p10/other/GC_content_profile/GC_profile_50000.cnp")
#   # setnames(tab)
#   # tab[,V2 := V2 + 1]
#   fwrite(tab,file = "GC_profile_50000.cnp",sep = "\t",col.names = F)
#   
#   # cov_tab[,cn_id := rpart_segmentation(.SD),by = .(sample,chr)]
#   # cov_tab[,cnv_mean_cov := median(cov),by = .(sample,chr,cn_id)]
#   
#   # fwrite(cov_tab,file = "20_sample_norm_cov_table.tsv",sep = "\t")
#   cov_tab <- fread("20_sample_norm_cov_table.tsv")
#   
#   est_tumor_ratio_tab <- estimate_tumor_ratio_from_hist(cov_tab)
#   
#   write.xlsx(est_tumor_ratio_tab,"est_ctDNA_ratio_tab.xlsx")
#   
#   
#   control_freec_tab[,sum(V3 - V2), by = sample]
#   
#   
#   
#   cfreec_res_files <- list.files(paste0("input_files/CNV/variant_calls/",list.files("input_files/CNV/variant_calls"),"/control_freec/"),pattern = "bam_CNVs",full.names = T)
#   cfreec_tab <- lapply(cfreec_res_files,fread)
#   names(cfreec_tab) <- gsub(".bam_CNVs","",list.files(paste0("input_files/CNV/variant_calls/",list.files("input_files/CNV/variant_calls"),"/control_freec/"),pattern = "bam_CNVs"))
#   cfreec_tab <- rbindlist(cfreec_tab,use.names = T,idcol = "sample")
#   
#   
#   final_res <- sapply(res,function(tab) weighted.mean(tab$TL,1 - (tab$score / sum(tab$score))))
#   final_res <- sapply(res,function(tab) tab$TL[1])
#   
#   final_res2 <- sapply(res2,function(tab) weighted.mean(tab$TL,1 - (tab$score / sum(tab$score))))
#   
#   res_tab <- data.table(sample = names(final_res),TL1 = final_res,TL2 = final_res2[names(final_res)])
#   
#   
# 
# 
#   
#   test_cov_tab
#   
#   cn_count_to_normal <- cn_count / 2
#   
#   mu <- TL * ((0:6)/2 * normal_cov) + normal_cov * (1 - TL)
#   
#   
#   x <- c(191,	4.2,
#          211,	47.5,
#          334,	11.9,
#          338,	32.2,
#          372,	5.6,
#          347,	21.1,
#          386,	21.8,
#          404,	9.6,
#          432,	28.7)
#   x <- as.data.table(matrix(x,ncol = 2,byrow = T))
#   setnames(x,c("sample","Total_%"))
#   x[,sample := as.character(sample)]
#   x[,`Total_%` := `Total_%` / 100]
#   TL_est <- merge(TL_est,x,by = "sample",all.x = T)
#   TL_est[,x := pred_TL / `Total_%`]
# 
# }


#vykreslit hist + estimates
#nacist data z control freecu 
#vypocitat jabcontool varianty
#zkombinovat a mrknout na vysledek
#poslat killianovi



# 
# res2 <- lapply(unique(cov_tab$sample),function(select_sample){
#   hist <- hist(test_cov_tab[cov < 6000 & sample == select_sample]$cov,200,plot = F)
#   hist_counts <- hist$counts
#   hist_counts[hist_counts < 20] <- 20
#   peaks <- which(hist_counts > c(0,hist_counts[-length(hist_counts)]) & hist_counts > c(hist_counts[-1],0))
#   peaks <- peaks[order(hist$counts[peaks])]
#   
#   # res <- normalmixEM(test_cov_tab[sample == select_sample]$cov)
#   
#   for(peak in peaks){
#     merge_peak <- peak
#     for(dominant_peak in peaks[which(hist$counts[peak] < hist$counts[peaks])]){
#       if(all(hist$counts[peak:dominant_peak] > hist$counts[peak] * 0.9)){
#         if(hist$counts[dominant_peak] > hist$counts[merge_peak]){
#           merge_peak <- dominant_peak
#         }
#       }
#     }
#     if(peak != merge_peak){
#       peaks <- setdiff(peaks,peak)
#     }
#   }
#   
#   normal_cov <- hist$mids[which.max(hist$counts)]
#   
#   test_TL_vec <- seq(0.05,1,by = 0.01)
#   mat <- sapply(seq(0.05,1,by = 0.01),function(TL) TL * ((0:6)/2 * normal_cov) + normal_cov * (1 - TL))
#   
#   res <- apply(mat,2,function(TL_vec){
#     sum(sapply(hist$mids[peaks],function(peak_cov)  min(abs(TL_vec - peak_cov)) ))
#   })
# 
#   res_tab <- data.table(TL = test_TL_vec,score = res)
#   setorder(res_tab,score)
#   
#   return(res_tab)
# })
# 
# res2 <- lapply(res2,function(x) x[1:10])
# names(res2) <- unique(cov_tab$sample)


# select_sample = "211"
# 
# test <- densityMclust(test_cov_tab[sample == select_sample]$cov)
# plot(test, what = "BIC")

# test <- densityMclust(test_cov_tab[sample == select_sample]$cov,G = 2)
# 
# select_sample = "211"
# select_sample = "386"
# select_sample = "432"
# res <- lapply(unique(cov_tab$sample),function(select_sample) normalmixEM(test_cov_tab[sample == select_sample]$cov)$mu)


