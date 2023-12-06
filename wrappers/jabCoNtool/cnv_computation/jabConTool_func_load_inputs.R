

#
# functions for loading data
#
fread_vector_of_files <- function(file_list,sample_names){
  list_of_tabs <- lapply(file_list,fread)
  names(list_of_tabs) <- sample_names
  return(rbindlist(list_of_tabs,use.names = T,idcol = "sample"))
}

load_panel_intervals <- function(panel_intervals_filename,library_type){
  panel_intervals <- fread(panel_intervals_filename)
  if(length(panel_intervals) > 4){
    panel_intervals <- panel_intervals[,1:4,with = F]
  }
  if(length(panel_intervals) == 3){
    panel_intervals[,V4 := paste("reg",V1,V2,V3,sep = "_")]
  }
  setnames(panel_intervals,c("chr","start","end","region_name"))
  
  if(library_type == "wgs"){
    #remove regions smaller then 1/2 of size
    panel_intervals <- panel_intervals[end - start > (max(end - start) / 2) ,]
  }
  
  if(length(unique(panel_intervals$region_name)) < nrow(panel_intervals)){
    panel_intervals[,region_name_suffix := paste0("_",as.character(seq_along(start))),by = region_name]
    panel_intervals[,region_name_count := .N,by = region_name]
    panel_intervals[region_name_count == 1,region_name_suffix := ""]
    panel_intervals[,region_name := paste0(region_name,region_name_suffix)]
    panel_intervals[,c("region_name_count","region_name_suffix") := NULL]
  }
  
  panel_intervals[,region_id := seq_along(region_name)]
  setkeyv(panel_intervals, c("chr","start","end"))
}

get_cov_tab <- function(sample_tab,
                        panel_intervals,
                        cohort_tab,
                        library_type,
                        GC_normalization_file,
                        cytoband_file,
                        outlier_probability = 0,
                        wgs_outlier_filter_iter = 2,
                        usable_bases_ratio_threshold = 0.85,
                        join_intervals_distance = 0){
  
  cov_tab <- fread_vector_of_files(sample_tab$cov_tab_filenames,sample_tab$sample)
  cov_tab[,tail(names(cov_tab),3) := NULL]
  setnames(cov_tab,c("V1","V2","V3"),c("chr","start","end"))
  setnames(cov_tab,tail(names(cov_tab),1),c("cov_raw"))

  # chr X a Y pořešíme později, zatím analýza bez něj, aby nezkresloval výsledky
  #homsap hack TODO
  cov_tab <- cov_tab[chr %in% 1:22]

  #add panel_intervals grouping based on join
  if(join_intervals_distance > 0){
    current_interval_start <- panel_intervals$start[1]
    grouping_vec <- numeric(nrow(panel_intervals))
    current_group_id <- 1
    for(select_region_id in panel_intervals$region_id){
      if(panel_intervals$end[select_region_id] - current_interval_start > join_intervals_distance | panel_intervals$end[select_region_id] - current_interval_start < 0){
        current_group_id <- current_group_id + 1
        current_interval_start <- panel_intervals[select_region_id]$start
      }
      grouping_vec[select_region_id] <- current_group_id
    }

    panel_intervals[,grouping_vec := grouping_vec]
    cov_tab <- merge(cov_tab,panel_intervals[,.(chr,start,end,region_id,grouping_vec)],by = c("chr","start","end"))
    cov_tab <- cov_tab[,.(chr = chr[1],start = min(start),end = max(end),cov_raw = sum(cov_raw),region_id = grouping_vec[1]),by = .(sample,grouping_vec)]
    cov_tab[,grouping_vec := NULL]
  } else {
    cov_tab <- merge(cov_tab,panel_intervals[,.(chr,start,end,region_id)],by = c("chr","start","end"))
  }


  if(library_type != "wgs"){
  
    if(!is.null(cohort_tab)){
      # combine coverage tab
      region_tab <- unique(cov_tab,by = c("region_id","chr","start","end"))
      region_tab[,c("sample","cov_raw","cov") := NULL]
      cohort_cov_tab <- cohort_tab[!is.na(cov)]
      cohort_cov_tab <- merge.data.table(cohort_cov_tab,region_tab,by = c("chr","start","end"))
      cohort_cov_tab <- cohort_cov_tab[,names(cov_tab),with = F]
      cov_tab <- rbind(cov_tab,cohort_cov_tab)
    }
  
    #read coverage normalization
    overall_mean <- mean(cov_tab$cov_raw)
    cov_tab[,cov := cov_raw / mean(cov_raw) * overall_mean,by = sample]
  
  } else {
    
    if(cytoband_file != "no_cytoband"){
      centromere_tab <- fread(cytoband_file)
      
      setnames(centromere_tab,c("chr","start","end",as.character(seq_len(ncol(centromere_tab) - 4)),"band_type"))
      centromere_tab <- centromere_tab[band_type == "acen",.(start = min(start),end = max(end)),by = chr]
      centromere_tab[,acen_length := end - start]
      #set removed netromere area as centromere +/- 10% - #TODO set as parameter 
      centromere_tab[,start := start - (acen_length / 10)]
      centromere_tab[,end := end + (acen_length / 10)]
      setkey(centromere_tab,chr,start,end)
      
      cov_tab <- foverlaps(cov_tab,centromere_tab)
      cov_tab <- cov_tab[is.na(acen_length)]
      cov_tab[,c("start","end","acen_length") := NULL]
      setnames(cov_tab,c("i.start","i.end"),c("start","end"))
    }
    
    if(outlier_probability > 0){
      #test join cov_tab
      cov_tab[,median_cov := median(cov),by = region_id]
      cov_tab[,outlier_region := F]
      for(i in seq(1,by = 1,length.out = wgs_outlier_filter_iter)){
        norm_cov_tab <- cov_tab[sample %in% sample_tab[type == normalization_sample_type]$sample]
        norm_nbinom_dist <- fitdist(as.integer(norm_cov_tab[outlier_region == F]$cov),distr = "nbinom") 
        outlier_low_value <- qnbinom(outlier_probability,size = norm_nbinom_dist$estimate["size"],mu = norm_nbinom_dist$estimate["mu"])
        outlier_high_value <- qnbinom(1 - outlier_probability,size = norm_nbinom_dist$estimate["size"],mu = norm_nbinom_dist$estimate["mu"])
        cov_tab[!(median_cov > outlier_low_value & median_cov < outlier_high_value),outlier_region := T]
      }
      cov_tab <- cov_tab[outlier_region == F]
      cov_tab[,median_cov :=NULL]
      cov_tab[,outlier_region :=NULL]
    }
    
    if(nchar(GC_normalization_file) > 0 & GC_normalization_file != "no_GC_norm"){
      #load_GC_profile_file
      GC_bin_content_tab <- fread(GC_normalization_file)
      setnames(GC_bin_content_tab, c("chr","start","gc","usable_bases_ratio"))
      GC_bin_content_tab <- GC_bin_content_tab[usable_bases_ratio > usable_bases_ratio_threshold]
      GC_bin_content_tab[,usable_bases_ratio := NULL]
      
      #merge and normalize using linear model to predict correlation between coverage and GC
      #remove bins with too low gc content to be able to normalize #TODO
      min_predicted_value_ratio = 10
      
      cov_tab <- merge.data.table(cov_tab,GC_bin_content_tab,by = c("chr","start"))
      cov_tab[,c("intercept","a") := as.list(lm(cov_raw ~ gc,.SD)$coefficients),by = .(sample)]
      cov_tab[,lm_model_norm_cov := a * gc + intercept,by = .(sample)]
      cov_tab[,min_lm_model_norm_cov := (a[1] * mean(gc) + intercept[1]) / min_predicted_value_ratio,by = .(sample)]
      cov_tab <- cov_tab[gc > max(cov_tab[lm_model_norm_cov < min_lm_model_norm_cov]$gc)]
      cov_tab[,cov_norm_val_ratio := cov_raw / lm_model_norm_cov,by = .(sample)]
      cov_tab[,cov_norm := cov_norm_val_ratio * (a[1] * mean(gc) + intercept[1]),by = .(sample)]
      cov_tab[,cov := cov_norm]
      
      for(sa in unique(cov_tab$sample)){
        plot(cov_tab[sample == sa]$gc ,cov_tab[sample == sa]$cov_raw)
      }

      #filter out outlier regions
      cov_tab[,gc_hist_bin := cut(gc,seq(0,1,0.005))]
      cov_tab[,cov_norm_mean :=  mean(cov_norm),by = sample]
      cov_tab[,cov_norm_sd :=  sd(cov_norm),by = sample]
      cov_tab[,gc_hist_bin_count := sum(cov_norm > cov_norm_mean + 4 * cov_norm_sd),.(sample,gc_hist_bin)]
      #count threshold set to 15
      cov_tab[,sample_region_OK := F]
      cov_tab[gc_hist_bin_count < 15 | cov_norm < cov_norm_mean + 2 * cov_norm_sd,sample_region_OK := T]
      cov_tab[,sample_region_OK_count := .N,by = region_id]
      cov_tab <- cov_tab[sample_region_OK_count > length(unique(cov_tab$sample)) * 0.9]
      cov_tab[,c("gc",
                 "intercept",
                 "a",
                 "cov_norm",
                 "lm_model_norm_cov",
                 "min_lm_model_norm_cov",
                 "cov_norm_val_ratio",
                 "gc_hist_bin",
                 "cov_norm_mean",
                 "cov_norm_sd",
                 "gc_hist_bin_count",
                 "sample_region_OK",
                 "sample_region_OK_count") :=NULL]
      
    } else {
      cov_tab[,cov := cov_raw]
    } 
    
    # cov_tab[,cov_norm := log2(cov / cov_mode),by = sample]
    
  }
  
  setcolorder(cov_tab,c("sample","region_id","chr","start","end","cov_raw","cov"))
  return(cov_tab)
}

get_snp_tab <- function(sample_tab,panel_snps_filename,panel_intervals,individual_cov_threshold = 10,max_pos_cov_threshold = 20){
  snp_tab <- fread_vector_of_files(sample_tab$snp_tab_filenames,sample_tab$sample)
  setnames(snp_tab,c("sample","chr","pos","A","C","G","T","cov"))
  snp_tab <- snp_tab[cov > 0]
  panel_snps <- fread(panel_snps_filename)

  setnames(panel_snps,c("chr","pos","ref","alt","BAF"))
  #genetics law
  panel_snps[,pop_HET_probability := 2*BAF*(1-BAF)]

  panel_snps <- foverlaps(setkey(panel_snps[,.(chr,pos,pos2 = pos,alt,pop_HET_probability)]), panel_intervals, by.x=c("chr","pos","pos2"), by.y=c("chr","start","end"), nomatch = 0)

  panel_snps[,c("start","end","region_name","pos2") := NULL]
  setorder(panel_snps,region_id,-pop_HET_probability)
  panel_snps <- unique(panel_snps,by = c("chr","pos"))

  snp_tab <- merge(snp_tab,panel_snps,by = c("chr","pos"))
  row_vec <- (seq_along(snp_tab$chr) - 1L) * 4L + as.numeric(factor(snp_tab$alt,levels = c("A","C","G","T")))
  snp_tab[,alt_count := t(as.matrix(snp_tab[,.(A,C,G,`T`)]))[row_vec]]
  snp_tab <-snp_tab[,.(sample,region_id,chr,pos,pop_HET_probability,alt_count,ref_count = cov - alt_count,cov)]

  # chr X a Y pořešíme později, zatím analýza bez něj, aby nezkresloval výsledky
  #homsap hack TODO
  snp_tab <- snp_tab[chr %in% 1:22]

  snp_tab[,max_pos_cov := max(cov),by = .(chr,pos) ]

  snp_tab <- snp_tab[max_pos_cov > max_pos_cov_threshold]
  # snp_tab[cov < individual_cov_threshold,ref_count := NA]
  # snp_tab[cov < individual_cov_threshold,alt_count := NA]

  snp_tab[,cov := NULL]
  snp_tab[,max_pos_cov := NULL]

  #
  # snp_tab <- snp_tab[chr != "X"] # chr X pořešíme později, zatím analýza bez něj, aby nezkresloval výsledky

  setkey(snp_tab,chr,pos,sample)
  return(snp_tab)
}

load_and_prefilter_sample_data <- function(sample_tab,
                                           panel_intervals_filename,
                                           panel_snps_filename,
                                           library_type,
                                           GC_normalization_file,
                                           cytoband_file,
                                           cohort_tab){

  normalization_sample_type <- tail(sort(sample_tab$type),1)

  panel_intervals <- load_panel_intervals(panel_intervals_filename,library_type)
  cov_tab <- get_cov_tab(sample_tab,panel_intervals,cohort_tab,library_type,GC_normalization_file,cytoband_file)
  
  panel_intervals <- unique(cov_tab[,.(chr,start,end,region_id)])
  setkey(panel_intervals)
  
  
  if(panel_snps_filename != "no_use_snps"){
    snp_tab <- get_snp_tab(sample_tab,panel_snps_filename,panel_intervals)
  } else {
    snp_tab <- NULL
  }
  
  
  
  return(list(cov_tab,snp_tab))
}


###########################
###########################
######               ######
###### TRANS MATRIX  ######
######               ######
###########################
###########################



create_copy_number_categories_default_tabs <- function(default_cn_rel_prob_vec = c(1,2,500,2,1,0.5,2),cn_to_predict = "all",predict_LOH = T){
  

  if(cn_to_predict == "all"){
    cn_to_predict <- c("0","1","2","3","4","5","LOH")
  }
  
  if(!predict_LOH){
    cn_to_predict <- setdiff(cn_to_predict,"LOH")
  }
  
  cn_het_var_count_table <- data.table(
    cn_category = c("0","1","1","2","3","3","4","4","4","5","5","5","5","LOH","LOH"),
    cn_count = c(0,1,1,2,3,3,4,4,4,5,5,5,5,2,2),
    het_var_count = c(0,0,1,1,1,2,1,2,3,1,2,3,4,0,2),
    het_var_count_prob = c(1,1/2,1/2,1,1/2,1/2,1/4,1/2,1/4,1/8,3/8,3/8,1/8,1/2,1/2),
    cov_norm_factor = c(0,1/2,1/2,1,3/2,3/2,2,2,2,5/2,5/2,5/2,5/2,1,1)
  )

  cn_het_var_count_table <- cn_het_var_count_table[cn_category %in% cn_to_predict]
  cn_categories_tab <- unique(cn_het_var_count_table[,.(cn_category,cn_count,cov_norm_factor)])
  default_cn_rel_prob_vec <- default_cn_rel_prob_vec[match(cn_to_predict,cn_categories_tab$cn_category)]
  cn_categories_tab[,cn_rel_prob := default_cn_rel_prob_vec]
  cn_categories_tab[,cn_norm_prob := cn_rel_prob / sum(cn_rel_prob)]
  cn_categories_vec <- cn_categories_tab$cn_category
  
  return(list(cn_het_var_count_table = cn_het_var_count_table,
              cn_categories_tab = cn_categories_tab))
}
  
create_transition_matrix <- function(cn_categories_tab){
  
  transition_matrix <- matrix(1,nrow(cn_categories_tab),nrow(cn_categories_tab))
  diag(transition_matrix) <- 20000
  transition_matrix <- t(t(transition_matrix) / rowSums(t(transition_matrix)))
  
  # transition_matrix <- matrix(cn_categories_tab$cn_norm_prob,nrow(cn_categories_tab),nrow(cn_categories_tab))
  # transition_matrix[which(cn_categories_tab$cn_category == "2"),] <- 1/nrow(cn_categories_tab)
  # diag(transition_matrix) <- max(cn_categories_tab$cn_norm_prob)
  # transition_matrix <- t(t(transition_matrix) / rowSums(t(transition_matrix)))
  
  return(transition_matrix)
}

prepare_transition_matrix_list <- function(transition_matrix,cn_categories_tab,cov_tab,dist_transition_treshold,calling_type){
  region_list <- unique(cov_tab[,.(region_id,chr,start,end)])
  region_list[,dist := c(Inf,tail(start,-1) - head(end,-1)),by = chr]
  #compute trans_probability from distance and dist_transition_treshold
  region_list[,trans_prob := 1 - plogis(dist,dist_transition_treshold,dist_transition_treshold / 5)]

  #alternative simple yes no 
  # region_list[,trans_prob := as.numeric(dist < dist_transition_treshold)] # ---  simple yes no 
  
  #independent_prior_matrix == no information from predicted probs from previous region - distant non correlated regions
  independent_prior_matrix <- matrix(cn_categories_tab$cn_norm_prob,nrow(cn_categories_tab),nrow(cn_categories_tab))
  
  trans_mat_list <- lapply(region_list$trans_prob,function(trans_prob) {
    res <- t(transition_matrix * trans_prob + independent_prior_matrix * (1 - trans_prob))
    res <- -log(res)
    res[res == Inf] <- 1000
    return(res)
  })
  
  return(trans_mat_list)
}



