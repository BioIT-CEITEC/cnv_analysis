suppressMessages(library(data.table))
suppressMessages(library(cn.mops))
suppressMessages(library(GenomicRanges))

Sys.setenv("R_ZIPCMD" = "zip")

setMethod("$", "GRanges", function(x, name) { # {{{
  elementMetadata(x)[, name]
}) # }}}

setMethod("$<-", "GRanges", function(x, name, value) { # {{{
  elementMetadata(x)[[ name ]] <- value
  return(x)
}) # }}}

run_all <- function(args){
  #read all params as variables
  BED <- args[1]
  Rdata <- args[2]
  BAMS <- args[3:length(args)]

  exons <- fread(BED,sep="\t")
  colnames(exons) <- c("chromosome","start","end","gene")
  exons <- exons[exons$end>exons$start,] # end > start
  exons$width <- exons$end-exons$start
  exons <- exons[exons$width>10,]
  exons <- subset(exons, select = -c(width,gene))
  exons <- exons[exons$chromosome!="X",]
  exons <- exons[exons$chromosome!="Y",]
  exons <- unique(exons,by=c("chromosome", "start", "end"))
  exons <- as.data.frame(exons)

  exons <- GRanges(exons[,1],IRanges(exons[,2]-30,exons[,3]+30))
  exons <- reduce(exons)

  ## load COHORT BAMS
  BAMcohort <- BAMS

  X_cohort <- cn.mops::getSegmentReadCountsFromBAM(BAMcohort, GR=exons, parallel=8)

  save(X_cohort,exons,file=Rdata)

}

script_dir <- dirname(sub("--file=", "", commandArgs()[grep("--file=", commandArgs())]))
args <- commandArgs(trailingOnly = T)
run_all(args)
