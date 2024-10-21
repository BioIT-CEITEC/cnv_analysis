library(data.table)


run_all <- function(args){
  input_file <- args[1]
  output_bed <- args[2]
  window <- as.integer(args[3])
  overlap <- as.integer(args[4])
  if(is.na(overlap)){
    overlap <- 0
  }

  tab <- fread(input_file,skip = "@SQ",select = c(2,3),header = F)
  tab[,V2 := gsub("SN:","",V2,fixed = T)]
  tab[,V3 := as.integer(gsub("LN:","",V3,fixed = T))]
  setnames(tab,c("chr","chr_size"))

  if(any(tab$chr_size > window)){
    out_bed <- tab[chr_size > window,.(end = c(seq(window,chr_size,by = window) - 1,chr_size)),by = chr]
    out_bed[,start := c(0,head(end,-1) + 1 - overlap),by = chr]
    setcolorder(out_bed,c("chr","start","end"))
  } else {
    out_bed <- data.table(chr = character(),start = numeric(),end = numeric())
  }



  out_bed <- rbind(out_bed,tab[chr_size <= window,.(chr,start = 0,end = chr_size - 1)])
  out_bed <- out_bed[order(match(chr, tab$chr))]

  dir.create(dirname(output_bed),recursive = T,showWarnings = F)
  options(scipen=99)
  fwrite(out_bed,file = output_bed,sep = "\t",col.names = F)
  options(scipen=0)
  #
  # shell("bedtools nuc -fi ../seq/TAIR10-31.fa -bed wgs_50kb.bed > GC_profile_50kb.bedtools.res")
  #
  # tab <- fread("/mnt/ssd/ssd_3/references/arabidopsis_thaliana/TAIR10-31/other/GC_content_profile/GC_profile_50kb.bedtools.res")
  # tab <- tab[,.(chr = `#1_usercol`,start = `2_usercol`,GC_content = `5_pct_gc`,coverage = round(`5_pct_gc` + `4_pct_at`,5))]
  # tab[coverage == 0,GC_content := -1]
  # fwrite(tab,file = "/mnt/ssd/ssd_3/references/arabidopsis_thaliana/TAIR10-31/other/GC_content_profile/GC_profile_50000.cnp",sep = "\t",col.names = F)

}


# develop and test
# args <- c("/mnt/ssd/ssd_3/references/homsap/GRCh38-p10/seq/GRCh38-p10.dict","/mnt/ssd/ssd_3/references/homsap/GRCh38-p10/other/GC_content_profile/intervals_50000.bed",50000,50)
# args <- c("/mnt/ssd/ssd_3/references/arabidopsis/TAIR10-31/seq/TAIR10-31.dict","/mnt/ssd/ssd_3/references/arabidopsis/TAIR10-31/other/GC_content_profile/wgs_50kb.bed",50000,0)

#run as Rscript
args <- commandArgs(trailingOnly = T)
run_all(args)


