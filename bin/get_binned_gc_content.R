library(data.table)


run_all <- function(args){
  input_file <- args[1]
  output_bed <- args[2]

  tab <- fread(input_file)
  tab <- tab[,.(chr = `#1_usercol`,start = `2_usercol`,GC_content = `5_pct_gc`,coverage = round(`5_pct_gc` + `4_pct_at`,5))]
  tab[coverage == 0,GC_content := -1]
  fwrite(tab,file = output_bed,sep = "\t",col.names = F)

}

# develop and test

#run as Rscript
args <- commandArgs(trailingOnly = T)
run_all(args)


