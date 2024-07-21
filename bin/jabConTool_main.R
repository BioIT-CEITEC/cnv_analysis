library(data.table)
library(proxy)
library(fitdistrplus)
library(mclust)
library(tictoc)

# # develop and test
script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(paste0(script_dir,"/../../.."))
args <- readLines(con = "logs/all_samples/jabCoNtool/cnv_computation.log_Rargs")
args <- strsplit(args,split = " ")[[1]]
args <- sub("data/ceitec_cfg2","share/share",args)


#run as Rscript
# script_dir <- dirname(sub("--file=", "", commandArgs()[grep("--file=", commandArgs())]))


source(paste0(script_dir,"/jabConTool_func_load_inputs.R"))
source(paste0(script_dir,"/estimate_tumorload.R"))

#TODO sex chromosome estimation

#set_constants
sample_regex <<- ".*/(.*)/jabCoNtool.*"


min_corelation_threshold <<- 0.9
dist_transition_treshold <<- 50000


###########################
###########################
######               ######
######   SNP part    ######
######               ######
###########################
###########################


compute_snp_based_nloglike <- function(snp_tab,cn_het_var_count_table,complex_FP_probability = 0.005){

  # if(all(snp_tab$TL == 1)){
  #   cn_het_var_count_table <- cn_het_var_count_table[-1]
  # }


  snp_based_cn_nloglike <- data.table(id = rep(seq_along(snp_tab$chr),each = nrow(cn_het_var_count_table)),
                                      TL = rep(snp_tab$TL,each = nrow(cn_het_var_count_table)),
                                      pop_HET_probability = rep(snp_tab$pop_HET_probability,each = nrow(cn_het_var_count_table)),
                                      alt_count = rep(snp_tab$alt_count,each = nrow(cn_het_var_count_table)),
                                      ref_count = rep(snp_tab$ref_count,each = nrow(cn_het_var_count_table)),
                                      cn_category = rep(cn_het_var_count_table$cn_category,nrow(snp_tab)),
                                      het_var_count = rep(cn_het_var_count_table$het_var_count,nrow(snp_tab)),
                                      het_var_count_prob = rep(cn_het_var_count_table$het_var_count_prob,nrow(snp_tab)),
                                      cn_count = rep(cn_het_var_count_table$cn_count,nrow(snp_tab)))
  snp_based_cn_nloglike[,expected_het_ratio := ((1 - TL) + TL * het_var_count) / (2 * (1 - TL) + TL * cn_count)]
  snp_based_cn_nloglike[,snp_prob := dbeta(expected_het_ratio,alt_count + 1,ref_count + 1)]
  snp_based_cn_nloglike[,hom_0_prob := pbeta(0.1,alt_count + 1,ref_count + 1,lower.tail = T)]
  snp_based_cn_nloglike[,hom_1_prob := pbeta(0.9,alt_count + 1,ref_count + 1,lower.tail = F)]
  snp_based_cn_nloglike <- snp_based_cn_nloglike[,.(TL = TL[1],pop_HET_probability = pop_HET_probability[1],hom_prob = max(hom_0_prob,hom_1_prob),snp_prob = max(snp_prob * het_var_count_prob)),by = .(id,cn_category)]
  snp_based_cn_nloglike[,snp_prob := snp_prob / sum(snp_prob,na.rm = T),by = id]
  snp_based_cn_nloglike[,snp_prob := snp_prob + complex_FP_probability]
  snp_based_cn_nloglike[,snp_prob := snp_prob / sum(snp_prob,na.rm = T),by = id]
  snp_based_cn_nloglike[,snp_prob := snp_prob * pop_HET_probability + hom_prob * (1 - pop_HET_probability)]
  snp_based_cn_nloglike[,snp_prob := snp_prob / sum(snp_prob,na.rm = T),by = id]

  if(any(snp_based_cn_nloglike$TL == 1)){
    snp_based_cn_nloglike[TL == 1,snp_prob := (1 - 1 / length(snp_prob)) * snp_prob,by = id]
    snp_based_cn_nloglike[TL == 1 & cn_category == "0",snp_prob := 1 / length(unique(cn_het_var_count_table$cn_category))]
  }

  snp_based_cn_nloglike[,snp_nloglike := -log(snp_prob),by = id]

  snp_tab <- cbind(snp_tab,dcast.data.table(snp_based_cn_nloglike,formula = id ~ cn_category,value.var = "snp_nloglike"))
  snp_tab[,id := NULL]

  return(snp_tab)
}



###########################
###########################
######               ######
###### COVRAGE part  ######
######               ######
###########################
###########################


vector.weighted.mean <- function(value_vec,weight_matrix){
  return(colSums(value_vec * weight_matrix) / colSums(weight_matrix))
}

vector.weighted.sd <- function(value_vec,weight_matrix){
  mean_dif_mat <- matrix(value_vec - rep(vector.weighted.mean(value_vec,weight_matrix),each = length(value_vec)),ncol = length(value_vec))^2
  non_zero_weight_count <- rowSums(matrix(as.logical(weight_matrix),ncol = length(value_vec)))

  return(sqrt(colSums(mean_dif_mat * weight_matrix) / (((non_zero_weight_count - 1) / non_zero_weight_count)*colSums(weight_matrix))))
}

compute_distribution_pramaters <- function(cov_tab,library_type = "panel"){


  if(any(names(cov_tab) == "cov_norm_factor")){
    cov_tab[,norm_vec := TL * cov_norm_factor + (1 - TL)]
  } else {
    cov_tab[,norm_vec := 1]
  }

  if(library_type != "wgs"){
    mat <- dcast.data.table(cov_tab,formula = region_id ~ sample,value.var = "cov",fill = 0)
    mat_row_name <- mat$region_id
    mat <- as.matrix(mat[,-1,with = F])
    rownames(mat) <- mat_row_name
    cor_mat <- as.matrix(simil(mat,method = cor, by_rows = F,use = "complete.obs"))
    cor_mat <- cor_mat - min_corelation_threshold
    cor_mat[cor_mat < 0] <- 0
    cor_mat <- cor_mat / (1 - min_corelation_threshold)
    diag(cor_mat) <- 1

    cor_sums_tab <- rowSums(cor_mat)
    cor_sums_tab <- data.table(sample = names(cor_sums_tab),norm_weight = cor_sums_tab)

    cov_tab[,norm_dist_mean := vector.weighted.mean(cov,cor_mat),by = "region_id"]
    cov_tab[,norm_dist_sd := vector.weighted.sd(cov,cor_mat),by = "region_id"]
    cov_tab[is.na(norm_dist_sd),norm_dist_sd := 1]
    cov_tab <- merge(cov_tab,cor_sums_tab,by = "sample")

  } else {
    normalization_sample_type <- tail(sort(unique(cov_tab$type)),1)
    norm_cov_data <- cov_tab[type == normalization_sample_type]
    norm_cov_data[,cov := cov / norm_vec]
    norm_cov_data[,outlier_region := F]
    norm_filter_iter <- 6
    outlier_probability <- 0.1
    for(i in seq(1,by = 1,length.out = norm_filter_iter)){
      norm_nbinom_dist <- fitdist(as.integer(norm_cov_data[outlier_region == F]$cov),distr = "nbinom")
      outlier_low_value <- qnbinom(outlier_probability,size = norm_nbinom_dist$estimate["size"],mu = norm_nbinom_dist$estimate["mu"])
      outlier_high_value <- qnbinom(1 - outlier_probability,size = norm_nbinom_dist$estimate["size"],mu = norm_nbinom_dist$estimate["mu"])
      norm_cov_data[!(cov > outlier_low_value & cov < outlier_high_value),outlier_region := T]
    }
    norm_cov_data <- norm_cov_data[outlier_region == F]
    norm_nbinom_dist <- fitdist(as.integer(norm_cov_data$cov),distr = "nbinom")

    mu <- norm_nbinom_dist$estimate["mu"]
    size <- norm_nbinom_dist$estimate["size"]
    cov_tab[,nbinom_mean := mu]
    cov_tab[,nbinom_var := mu + mu^2/size]

  }

  cov_tab[,norm_vec := NULL]

  return(cov_tab)

}

#TODO implement differently with binom and byas
get_normal_negative_log_likelihoods <- function(cn_count,cov_tab,library_type){
  cn_count_to_normal <- cn_count / 2

  mean <- cov_tab$TL * (cn_count_to_normal * cov_tab$norm_dist_mean) + cov_tab$norm_dist_mean * (1 - cov_tab$TL)
  sd <- cov_tab$norm_dist_sd
  if(cn_count_to_normal == 0){
    mean[cov_tab$TL == 1] <- 0
    sd[cov_tab$TL == 1] <- cov_tab$norm_dist_sd
  }
  prob <- 1 - abs(pnorm(cov_tab$cov,mean = mean,sd = sd) - 0.5) * 2

  return(prob)
}

#TODO implement differently with binom and byas
get_nbinom_negative_log_likelihoods <- function(cn_count,cov_tab){
  cn_count_to_normal <- cn_count / 2

  mu <- cov_tab$TL * (cn_count_to_normal * cov_tab$nbinom_mean) + cov_tab$nbinom_mean * (1 - cov_tab$TL)
  var <- cov_tab$TL * (cn_count_to_normal * cov_tab$nbinom_var) + cov_tab$nbinom_var * (1 - cov_tab$TL)
  if(cn_count_to_normal == 0){
    mu[cov_tab$TL == 1] <- 1
    var[cov_tab$TL == 1] <- 50
  }

  size <- mu^2 / (var - mu)
  prob <- dnbinom(round(cov_tab$cov),mu = mu,size = size)
  return(prob)
}

compute_coverage_based_nloglike <- function(cov_tab,library_type,cn_categories_tab,complex_FP_probability = 0.005){

  cov_tab <- compute_distribution_pramaters(cov_tab,library_type)

  if(library_type != "wgs"){
    res_prob <- sapply(cn_categories_tab$cn_count,get_normal_negative_log_likelihoods,cov_tab = cov_tab)
    res_prob <- res_prob / rowSums(res_prob)
  } else {
    res_prob <- sapply(cn_categories_tab$cn_count,get_nbinom_negative_log_likelihoods,cov_tab = cov_tab)
    res_prob <- res_prob / rowSums(res_prob)
    res_prob[cov_tab$cov > cov_tab$nbinom_mean * max(cn_categories_tab$cov_norm_factor),which.max(cn_categories_tab$cn_count)] <- res_prob[cov_tab$cov > cov_tab$nbinom_mean * max(cn_categories_tab$cov_norm_factor),which.max(cn_categories_tab$cn_count)] + 10^-10
  }

  res_prob[is.na(res_prob)] <- 0
  res_prob <- res_prob + complex_FP_probability
  res_prob <- res_prob / rowSums(res_prob)
  res_nloglike <- -log(res_prob)
  res_nloglike[res_nloglike == Inf] <- 1000

  colnames(res_nloglike) <- cn_categories_tab$cn_category
  cov_tab <- cbind(cov_tab,res_nloglike)
  cov_tab <- cov_tab[type == "call"]


  setkey(cov_tab,chr,start,end,sample)

  return(cov_tab)
}


###########################
###########################
######               ######
###### bayes part   ######
######               ######
###########################
###########################


compute_byes_net_likelihood <-function(trans_mat_list,per_region_state_nloglike_matrix,only_model_nloglike = T,cn_categories_vec){

  # per_region_state_nloglike_matrix = combined_nloglike_matrix
  # only_model_nloglike = estimate_only

  n_regions <- length(trans_mat_list)
  prev = matrix(0, n_regions-1, length(cn_categories_vec))
  omega = matrix(0, n_regions, length(cn_categories_vec))


  omega[1,] = trans_mat_list[[1]][1,] + per_region_state_nloglike_matrix[1,]

  index_vec <- seq_along(trans_mat_list)[-1]

  trans_index <- 1:ncol(trans_mat_list[[1]]) - 1
  trans_mat_cols <- ncol(trans_mat_list[[1]])

  for(i in index_vec){
    probs_mat <- omega[i - 1, ] + trans_mat_list[[i]]
    prev[i - 1, ] <- apply(probs_mat,2,which.min)
    omega[i, ] <- probs_mat[prev[i - 1, ] + trans_index * trans_mat_cols] + per_region_state_nloglike_matrix[i,]
  }

  if(only_model_nloglike){
    return(min(omega[nrow(omega),]))
  }

  predict_state_vector <- integer(length = n_regions)
  predict_state_vector[length(predict_state_vector)] <- which.min(omega[nrow(omega),])


  for(i in (length(predict_state_vector)-1):1){
    predict_state_vector[i] <- prev[i,predict_state_vector[i + 1]]
  }

  return(list(min(omega[nrow(omega),]),predict_state_vector))
}

###########################
###########################
######               ######
###### wrapper part  ######
######               ######
###########################
###########################

predict_CNV_model <- function(sample_tab,trans_mat_list,cov_tab,snp_tab,library_type,categories_default_tabs,complex_FP_probability,TL_select = NULL,estimate_only = F){

  # sample_tab <- process_sample_tab
  sample_tab <- copy(sample_tab)

  if(!is.null(TL_select)){
    sample_tab[,TL := sample_tab[[TL_select]]]
  } else {
    if(!any(names(sample_tab) == "TL")){
      sample_tab[,TL := 1]
    }
  }

  tictoc::tic()
  print(paste0("cov nlog prob compute: "))

  if(any(names(cov_tab) == "cn_pred")){
    cov_tab <- merge(sample_tab[,.(sample,type,TL_new)],cov_tab,by = "sample")
    cov_tab[cn_pred == "2",TL := TL_new]
    cov_tab[,TL_new := NULL]
    cov_tab <- merge(cov_tab,categories_default_tabs$cn_categories_tab[,.(cn_pred = cn_category,cov_norm_factor)],by = "cn_pred")
  } else {
    cov_tab <- merge(sample_tab[,.(sample,type,TL)],cov_tab,by = "sample")
  }
  cov_nloglike_tab <- compute_coverage_based_nloglike(cov_tab,library_type,categories_default_tabs$cn_categories_tab,complex_FP_probability)

  tictoc::toc()


  if(!is.null(snp_tab)){
    tictoc::tic()
    print(paste0("snp nlog prob compute: "))
    if(any(names(snp_tab) == "cn_pred")){
      snp_tab <- merge(sample_tab[,.(sample,type,TL_new)],snp_tab,by = "sample")
      snp_tab[cn_pred == "2",TL := TL_new]
      snp_tab[,TL_new := NULL]

    } else {
      snp_tab <- merge(sample_tab[,.(sample,TL)],snp_tab,by = "sample")
    }

    snp_nloglike_tab <- compute_snp_based_nloglike(snp_tab,categories_default_tabs$cn_het_var_count_table,complex_FP_probability)
    tictoc::toc()
  }

  cn_categories_vec <- categories_default_tabs$cn_categories_tab$cn_category


  # Initializing a vector to keep track of samples that should be removed due to errors
  samples_to_remove <- character()
  
  # Using lapply for iteration and error handling
  per_sample_res <- lapply(sample_tab[type == "call"]$sample, function(sel_sample) {
    tryCatch({
      # Displaying the sample being processed
      print(paste0("sample: ", sel_sample))
      
      # Creating combined negative log-likelihood tables
      if (!is.null(snp_tab)) {
        tictoc::tic()
        combined_nloglike_tab <- rbind(cov_nloglike_tab[sample == sel_sample, c("region_id", cn_categories_vec), with = FALSE],
                                       snp_nloglike_tab[sample == sel_sample, c("region_id", cn_categories_vec), with = FALSE])
        setkey(combined_nloglike_tab, region_id)
        combined_nloglike_tab <- melt.data.table(combined_nloglike_tab, id.vars = "region_id")
        combined_nloglike_tab <- combined_nloglike_tab[, .(value = sum(value)), by = .(region_id, variable)]
        combined_nloglike_matrix <- matrix(data = combined_nloglike_tab$value, nrow = length(trans_mat_list))
        tictoc::toc()
      } else {
        combined_nloglike_matrix <- as.matrix(cov_nloglike_tab[sample == sel_sample, cn_categories_vec, with = FALSE])
      }
      
      # Computing Bayesian network likelihood
      tictoc::tic()
      res <- compute_byes_net_likelihood(trans_mat_list = trans_mat_list,
                                         per_region_state_nloglike_matrix = combined_nloglike_matrix,
                                         only_model_nloglike = estimate_only,
                                         cn_categories_vec = cn_categories_vec)
      tictoc::toc()
      
      # Return the result
      return(res)
    }, error = function(e) {
      # Print error message
      print(paste0("Error processing sample ", sel_sample, ": ", e$message))
      # Add to samples to be removed
      samples_to_remove <<- c(samples_to_remove, sel_sample)
      # Return NULL for this sample
      return(NULL)
    })
  })
  
  # Removing samples that caused errors from the sample_tab
  sample_tab <- sample_tab[!sample_tab$sample %in% samples_to_remove,]
  per_sample_res <- Filter(Negate(is.null), per_sample_res)
  
  

  if(estimate_only){
    TL_estimate_tab <- sample_tab[type == "call",c("sample",TL_select),with = F]
    TL_estimate_tab[,(paste0(TL_select,"_nloglike")) := unlist(per_sample_res)]

    return(TL_estimate_tab)
  } else {
    model_estimate_tab <- sample_tab[type == "call",.(sample,TL,nloglike = sapply(per_sample_res,function(x) x[[1]]))]
    cn_predict_tab <- lapply(per_sample_res,function(x) data.table(region_id = sort(unique(cov_tab$region_id)),cn_pred = x[[2]]))
    names(cn_predict_tab) <- sample_tab[type == "call"]$sample
    cn_predict_tab <- rbindlist(cn_predict_tab,use.names = T,idcol = "sample")

    if(!is.null(snp_tab)){
      if(library_type != "wgs"){
        cn_call_info_tab <- rbind(cov_nloglike_tab[,.(sample,region_id,cov,norm_dist_mean,norm_dist_sd,pos = NA,pop_HET_probability = NA,alt_count = NA,ref_count = NA)],
                                  snp_nloglike_tab[,.(sample,region_id,cov = NA,norm_dist_mean = NA,norm_dist_sd = NA,pos,pop_HET_probability,alt_count,ref_count)])
      } else {
        cn_call_info_tab <- rbind(cov_nloglike_tab[,.(sample,region_id,cov,nbinom_mean,nbinom_var,pos = NA,pop_HET_probability = NA,alt_count = NA,ref_count = NA)],
                                  snp_nloglike_tab[,.(sample,region_id,cov = NA,norm_dist_mean = NA,norm_dist_sd = NA,pos,pop_HET_probability,alt_count,ref_count)])
      }

      cn_call_info_tab <- cbind(cn_call_info_tab,rbind(cov_nloglike_tab[,cn_categories_vec,with = F],snp_nloglike_tab[,cn_categories_vec,with = F]))
      setorder(cn_call_info_tab,sample,region_id)
    } else {
      cov_nloglike_tab[,c("chr","start","end") := NULL]
      cn_call_info_tab <- cov_nloglike_tab
    }

    return(list(model_estimate_tab,cn_predict_tab,cn_call_info_tab))
  }
}


test_TL_from_N_longest_CNVs <- function(final_estimates,cn_categories_tab,N_longest_CNVs = 10,min_CNV_len = mean(final_estimates[,.(end - start)]$V1) * 3){
  CNV_tab <- final_estimates[cn_pred != 3]
  CNV_tab <- CNV_tab[,.(start = min(start),end = max(end),cov = median(cov),cn_pred_id = cn_pred[1]),by = .(sample,chr,cn_id)]
  CNV_tab[,cnv_length := end - start]
  # CNV_tab <- CNV_tab[cnv_length > min_CNV_len]
  setorder(CNV_tab,-cnv_length)

  CNV_tab <- CNV_tab[,.SD[1:min(N_longest_CNVs,nrow(.SD))],by = sample]
  CNV_tab <- merge(CNV_tab,final_estimates[,.(norm_cov = median(cov / cn_categories_tab$cov_norm_factor[cn_pred],na.rm = T)),by = sample],by = "sample")
  CNV_tab[,CNV_pred_TL := (cov / norm_cov - 1) / (cn_categories_tab$cn_count[cn_pred_id] / 2 - 1)]
  CNV_tab <- CNV_tab[!is.na(CNV_pred_TL)]
  TL_estimates_tab <- CNV_tab[,.(pred_TL = weighted.mean(x = CNV_pred_TL,w = cnv_length),pred_TL_sd = sqrt(weighted.mean( (CNV_pred_TL - weighted.mean(x = CNV_pred_TL,w = cnv_length))^2, cnv_length )),CNV_size_Mbp = sum(cnv_length) / 10^6),by =sample]
  TL_estimates_tab[,rel_CNV_size := round(CNV_size_Mbp / (sum(final_estimates[sample == sample[1],max(end),by = chr]$V1) / 10 ^ 6) * 100,2) ]
  TL_estimates_tab[pred_TL > 0.99,pred_TL := 0.99]
  return(TL_estimates_tab)
}

predict_CNVs <- function(sample_tab,cov_tab,snp_tab,library_type,trans_mat_list,categories_default_tabs,initial_TL,iterations,complex_FP_probability){

  CNV_pred_list <- list()
  TL_estimates_tab_list <- list()
  if(any(names(sample_tab) == "TL")){
    TL_estimates_tab_list[[1]] <- data.table(sample = sample_tab$sample,
                                             pred_TL = sample_tab$TL,
                                             pred_TL_sd = 0,
                                             CNV_size_kb = 0,
                                             rel_CNV_size = 0)
    sample_tab[,TL := NULL]
  } else {
    TL_estimates_tab_list[[1]] <- data.table(sample = sample_tab$sample,
                                             pred_TL = initial_TL,
                                             pred_TL_sd = 0,
                                             CNV_size_kb = 0,
                                             rel_CNV_size = 0)
  }


  for(i in seq_len(iterations)){
    # i <-1
    tictoc::tic()
    print(paste0("iter: ",i))
    process_sample_tab <- merge.data.table(sample_tab,TL_estimates_tab_list[[i]][,.(sample,TL = pred_TL)],by = "sample",all.x = T)
    process_sample_tab[is.na(TL),TL := 0.1]

    res <- predict_CNV_model(process_sample_tab,trans_mat_list,cov_tab,snp_tab,library_type,categories_default_tabs,complex_FP_probability)

    # save(res,file = "test_prediction.Rdata")
    # load("test_prediction.Rdata")
    TL_estimates_table <- res[[1]]
    final_estimates <- res[[2]]
    cn_call_info_tab <- res[[3]]

    final_estimates <- merge.data.table(cov_tab[,.(sample,region_id,chr,start,end,cov)],final_estimates,by = c("sample","region_id"))

    rle_res <- final_estimates[,rle(cn_pred),by = .(sample,chr)]
    rle_res[,cn_id := seq_along(values),by = .(sample,chr)]
    final_estimates[,cn_id := rep(rle_res$cn_id,rle_res$lengths)]

    if(initial_TL != 1){
      TL_estimates_tab_list[[i + 1]] <- test_TL_from_N_longest_CNVs(final_estimates,categories_default_tabs$cn_categories_tab,3)
      TL_estimates_tab_list[[i + 1]] <- merge.data.table(TL_estimates_tab_list[[i]][,.(sample)],TL_estimates_tab_list[[i + 1]],by = "sample",all.x = T)
      TL_estimates_tab_list[[i + 1]][is.na(CNV_size_Mbp),CNV_size_Mbp := 0]
      TL_estimates_tab_list[[i + 1]][is.na(rel_CNV_size),rel_CNV_size := 0]
    } else {
      TL_estimates_tab_list[[i + 1]] <- TL_estimates_tab_list[[1]]
    }

    final_estimates[,cov := NULL]
    final_estimates <- merge.data.table(final_estimates,cn_call_info_tab,by = c("sample","region_id"))
    CNV_pred_list[[i]] <- final_estimates

    tictoc::toc()
  }

  return(list(CNV_pred_list,TL_estimates_tab_list))

}


###########################
###########################
######               ######
######  postprocess  ######
######               ######
###########################
###########################



run_all <- function(args){
  out_filename <- args[1]
  panel_intervals_filename <- args[2]
  panel_snps_filename <- args[3] #filename or "no_use_snps"
  calling_type <- args[4] #tumor_only, tumor_normal, germline
  library_type <- args[5] #wgs, panel
  GC_normalization_file <- args[6] #filename or "no_GC_norm"
  cytoband_file <- args[7] #filename or "no_cytoband"
  cohort_data_filename <- args[8] #filename or "no_previous_cohort_data"
  prior_est_tumor_ratio <- as.logical(args[9])
  max_CNV_frequency_in_cohort <- as.numeric(args[10]) / 100

  dir.create(dirname(out_filename),recursive = T,showWarnings = F)

  #set defuault copy number and error probability if not set in params
  #TODO add to params (full vector or just non normal probability) for now is null
  default_cn_rel_prob_vec <- NULL
  if(is.null(default_cn_rel_prob_vec)){
    default_cn_rel_prob_vec <- c(1,2,500,2,1,0.5,2)
    default_cn_rel_prob_vec <- default_cn_rel_prob_vec / sum(default_cn_rel_prob_vec)
  }
  complex_FP_probability <- NULL
  if(is.null(complex_FP_probability)){
    complex_FP_probability <- (1 - max(default_cn_rel_prob_vec)) / 6
  }

  #create sample table
  if(calling_type == "tumor_normal"){
    tumor_cov_tab_filenames <- args[(which(args == "cov") + 1):(which(args == "norm_cov") - 1)]
    normal_cov_tab_filenames <- args[(which(args == "norm_cov") + 1):length(args)]
    sample_tab <- data.table(cov_tab_filenames = c(tumor_cov_tab_filenames,normal_cov_tab_filenames),
                             type = c(rep("call",length(tumor_cov_tab_filenames)),rep("normal",length(normal_cov_tab_filenames))))
    sample_tab[,sample := gsub(sample_regex,"\\1",cov_tab_filenames)]
    sample_tab[type == "normal",sample := paste0(sample,"_norm")]
  } else {
    sample_tab <- data.table(cov_tab_filenames = args[(which(args == "cov") + 1):length(args)],
                             type = c(rep("call",length(cov_tab_filenames))))
    sample_tab[,sample := gsub(sample_regex,"\\1",cov_tab_filenames)]
  }

  if(panel_snps_filename != "no_use_snps"){
    sample_tab[,snp_tab_filenames := gsub(".region_coverage.tsv",".snpAF.tsv",cov_tab_filenames)]
  }
  setcolorder(sample_tab,c("sample","type"))


  if(cohort_data_filename != "no_previous_cohort_data"){
    cohort_tab <- fread(cohort_data_filename)
    cohort_tab[,chr := as.character(chr)]
    cohort_tab <- cohort_tab[!(sample %in% sample_tab$sample)]
  } else {
    cohort_tab <- NULL
  }

  #load all data
  res <- load_and_prefilter_sample_data(sample_tab,
                                        panel_intervals_filename,
                                        panel_snps_filename,
                                        library_type,
                                        GC_normalization_file,
                                        cytoband_file,
                                        cohort_tab)
  cov_tab <- res[[1]]
  snp_tab <- res[[2]]
  sample_tab[,c("cov_tab_filenames","snp_tab_filenames") := NULL]

  if(cohort_data_filename != "no_previous_cohort_data"){
    #combine sample tab
    cohort_sample_tab <- data.table(sample = unique(cohort_tab$sample),type = "cohort")
    sample_tab <- rbind(sample_tab,cohort_sample_tab)
  }


  categories_default_tabs <- create_copy_number_categories_default_tabs(default_cn_rel_prob_vec,predict_LOH = !is.null(snp_tab))
  transition_matrix <- create_transition_matrix(categories_default_tabs$cn_categories_tab)
  trans_mat_list <- prepare_transition_matrix_list(transition_matrix,categories_default_tabs$cn_categories_tab,cov_tab,dist_transition_treshold,calling_type)

  if(calling_type != "germline" & library_type == "wgs" & prior_est_tumor_ratio == T){
    est_ratio_score_tab <- estimate_tumor_ratio_from_hist(cov_tab)
    est_tumor_ratio_tab <- est_ratio_score_tab[,.(TL = TL[1]),by = sample]
    sample_tab <- merge(sample_tab,est_tumor_ratio_tab,by = "sample")
  }

  # save(sample_tab,cov_tab,snp_tab,library_type,trans_mat_list,categories_default_tabs,file = "structural_varcalls/all_samples/jabCoNtool/test_run_objects.Rdata")
  # load("test_run_objects_260123.Rdata")
  # iterations = 1

  if(calling_type == "germline") {
    initial_TL = 1
    iterations = 1
  } else {
    initial_TL = 0.75
    iterations = 3
  }

  res <- predict_CNVs(sample_tab,cov_tab,snp_tab,library_type,trans_mat_list,categories_default_tabs,initial_TL,iterations,complex_FP_probability)

  # save(res,file = "test_CNV_call_all_samples_260123.Rdata")
  # load("test_optim_1_init_0_05.Rdata")

  final_cn_pred_info_table <- copy(res[[1]][[iterations]])
  final_cn_pred_info_table[,cn_pred := categories_default_tabs$cn_categories_tab$cn_category[cn_pred]]

  #store the cohort data
  if(cohort_data_filename != "no_previous_cohort_data"){
    calling_info_tab <- rbind(cohort_tab,final_cn_pred_info_table[,names(cohort_tab),with = F])
  } else {
    calling_info_tab <- final_cn_pred_info_table[,intersect(names(final_cn_pred_info_table),c("sample","chr","start","end","cn_pred","cov_raw","cov","pos","alt_count","ref_count")),with = F]
  }
  fwrite(calling_info_tab,file = paste0(dirname(out_filename),"/cohort_info_tab.tsv"),sep="\t")

  #remove and store too frequent CNVs in the cohort
  final_cn_pred_info_table[,CNV_in_cohort := length(unique(sample)),by = .(cn_pred,region_id)]

  if(cohort_data_filename != "no_previous_cohort_data"){

    CNV_in_previous_cohort_tab <- cohort_tab[,.(CNV_in_previous_cohort = length(unique(sample))),by = .(chr,start,end,cn_pred)]
    final_cn_pred_info_table <- merge(final_cn_pred_info_table,CNV_in_previous_cohort_tab,by = c("chr","start","end","cn_pred"),all.x = T)
    setcolorder(final_cn_pred_info_table,c("sample","region_id"))
    final_cn_pred_info_table[is.na(CNV_in_previous_cohort),CNV_in_previous_cohort := 0]
    final_cn_pred_info_table[,CNV_in_cohort := CNV_in_cohort + CNV_in_previous_cohort]
    final_cn_pred_info_table[,CNV_in_previous_cohort := NULL]
  }

  final_cn_pred_info_table[,samples_in_cohort := length(unique(cov_tab$sample))]
  final_cn_pred_info_table[,too_frequent_FP_CNVs := cn_pred != "2" & CNV_in_cohort / samples_in_cohort > max_CNV_frequency_in_cohort ]
  too_frequent_FP_CNVs_tab <- final_cn_pred_info_table[too_frequent_FP_CNVs == T,]
  fwrite(too_frequent_FP_CNVs_tab,file = paste0(dirname(out_filename),"/too_frequent_filtered_CNVs.tsv"),sep="\t")
  final_cn_pred_info_table[too_frequent_FP_CNVs == T,cn_pred := "2"]
  final_cn_pred_info_table[,too_frequent_FP_CNVs := NULL]

  fwrite(final_cn_pred_info_table,file = out_filename,sep="\t")


  if(calling_type != "germline") {
    tumor_cell_fraction_table <- res[[2]][[iterations]]
    fwrite(tumor_cell_fraction_table,file = paste0(dirname(out_filename),"/tc_fraction_prediction.tsv"),sep="\t")
  }
  save(res,too_frequent_FP_CNVs_tab,file = paste0(dirname(out_filename),"/support_data.Rdata"))
}

#run as Rscript

# args <- commandArgs(trailingOnly = T)
# print("start")
# timestamp()
# run_all(args)
# print("end")
# timestamp()




#### somatic_calling
# sample_tab[sample == "PC_3_WT",type := "norm"]
# cov_tab <-merge.data.table(sample_tab,cov_tab,by = "sample")
# median_norm_cov <- median(cov_tab[type == "norm"]$cov)
# cov_tab <- cov_tab[,.(chr = chr[1],
#            start = start[1],
#            end = end[1],
#            sample = sample[type == "call"],
#            cov_raw = cov_raw[type == "call"],
#            cov = cov[type == "call"] - mean(cov[type == "norm"]) + median_norm_cov),by = region_id]
# cov_tab[cov < 0,cov := 0]
# sample_tab <- sample_tab[type == "call"]



# combine snp tab
# if(!is.null(snp_tab)){
#
#   region_tab <- unique(snp_tab,by = c("region_id","pos"))
#   region_tab[,c("sample","alt_count","ref_count") := NULL]
#   cohort_snp_tab <- cohort_tab[!is.na(pos)]
#   cohort_snp_tab[,c("cov","start","end") := NULL]
#   cohort_snp_tab <- merge.data.table(cohort_snp_tab,region_tab,by = c("chr","pos"))
#   setcolorder(cohort_snp_tab,names(snp_tab))
#   snp_tab <- rbind(snp_tab,cohort_snp_tab)
# }


