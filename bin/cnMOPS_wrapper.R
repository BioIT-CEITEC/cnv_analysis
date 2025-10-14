suppressMessages(library(data.table))
suppressMessages(library(cn.mops))
suppressMessages(library(GenomicRanges))

Sys.setenv("R_ZIPCMD" = "zip")

setMethod("$", "GRanges", function(x, name) { # {{{
  elementMetadata(x)[, name]
}) 

setMethod("$<-", "GRanges", function(x, name, value) { # {{{
  elementMetadata(x)[[ name ]] <- value
  return(x)
})

run_all <- function(args){

  BAM_file <- args[1]
  cohort_Rdata <- args[2]
  sampleNameWildCard <- args[3]
  output_file <- args[4]

  if(!cohort_Rdata==""){
    load(cohort_Rdata)
    if(exists("exons.GRCh38")){exons <- exons.GRCh38}
  }else{error("Missing cohort Rdata file")}

  CNVs <- data.frame()
  tryCatch({


    X_sample <- cn.mops::getSegmentReadCountsFromBAM(BAM_file,sampleNames="SAMPLENAME", GR=exons)
    X_sampleDF <- data.frame(X_sample)

    X_cohort$SAMPLENAME <- X_sampleDF$SAMPLENAME

    resCNMOPS <- cn.mops::exomecn.mops(X_cohort)
    resCNMOPS <- cn.mops::calcIntegerCopyNumbers(resCNMOPS)
    CNVs <- as.data.frame(cnvs(resCNMOPS))
  },
  error = function(e) {

    message("Error detected: ", e$message)
    message("Running fallback instead...")
  })

  CNVs <- CNVs[CNVs$sampleName=="SAMPLENAME",]

  if(nrow(CNVs)>0){
    CNVs$sampleName <- sampleNameWildCard
    setnames(CNVs,"seqnames","chromosome")
    dir.create(file.path(dirname(output_file)), showWarnings = FALSE)
    fwrite(CNVs,file = output_file,sep = "\t")

  }else{
    CNVs <- data.frame(seqnames=character(0),start=character(0),end=character(0),width=character(0),strand=character(0),sampleName=character(0),
                       median=character(0),mean=character(0),CN=character(0))
    dir.create(file.path(dirname(output_file)), showWarnings = FALSE)
    fwrite(CNVs,file = output_file,sep = "\t")
  }
}

script_dir <- dirname(sub("--file=", "", commandArgs()[grep("--file=", commandArgs())]))
args <- commandArgs(trailingOnly = T)
run_all(args)