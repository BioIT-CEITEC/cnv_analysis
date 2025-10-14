suppressMessages(library(data.table))
suppressMessages(library(ExomeDepth))
suppressMessages(library(GenomicRanges))

Sys.setenv("R_ZIPCMD" = "zip")

run_all <- function(args){
  #read all params as variables
  BED <- args[1]
  FASTA <- args[2]
  Rdata <- args[3]
  BAMS <- args[4:length(args)]
  
  
  exons <- fread(BED,sep="\t")
  colnames(exons) <- c("chromosome","start","end","gene")
  exons <- exons[exons$end>exons$start,] # end > start
  exons$width <- exons$end-exons$start
  exons <- exons[exons$width>10,]
  exons <- subset(exons, select = -c(width))
  exons <- exons[exons$chromosome!="X",]
  exons <- exons[exons$chromosome!="Y",]
  exons <- unique(exons,by=c("chromosome", "start", "end"))
  exons <- as.data.table(exons)

  exons$No <- 1
  exons[, rowNo := cumsum(No), by = gene]
  exons[, name := paste0(gene,"_E",rowNo)]

  exons[,c("No","rowNo","gene"):=NULL]
  exons <- exons[,c("chromosome", "start", "end", "name")]
  setnames(exons,"name","gene")
  exons <- as.data.frame(exons)

  BAMcohort <- BAMS

  X_cohort <- ExomeDepth::getBamCounts(bed.frame = exons,
                           bam.files = BAMcohort,
                           include.chr = FALSE,
                           referenceFasta = FASTA)

  save(X_cohort,exons,file = Rdata)

}

script_dir <- dirname(sub("--file=", "", commandArgs()[grep("--file=", commandArgs())]))
args <- commandArgs(trailingOnly = T)
run_all(args)
