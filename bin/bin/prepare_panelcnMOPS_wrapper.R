suppressMessages(library(data.table))
suppressMessages(library(panelcn.mops))
suppressMessages(library(GenomicRanges))

Sys.setenv("R_ZIPCMD" = "zip")

run_all <- function(args){
  #read all params as variables
  BED <- args[1]
  Rdata <- args[2]
  BAMS <- args[3:length(args)]

  ## load exons
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
  exons[, name := paste(gene,paste0("E",rowNo),chromosome,start,end,sep = ".")]
  exons[, exon := paste0("E",rowNo)]

  exons[,c("No","rowNo"):=NULL]
  exons <- exons[,c("chromosome", "start", "end", "name", "gene", "exon")]
  exons <- as.data.frame(exons)

  BAMcohort <- BAMS

  X_cohort <- panelcn.mops::countBamListInGRanges(countWindows = exons,
                                    bam.files = BAMcohort, read.width = FALSE)

  save(X_cohort,exons,file=Rdata)


}

script_dir <- dirname(sub("--file=", "", commandArgs()[grep("--file=", commandArgs())]))
args <- commandArgs(trailingOnly = T)
run_all(args)