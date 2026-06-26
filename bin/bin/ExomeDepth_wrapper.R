suppressMessages(library(data.table))
suppressMessages(library(ExomeDepth))
suppressMessages(library(GenomicRanges))

Sys.setenv("R_ZIPCMD" = "zip")

run_all <- function(args){
  BAM_file <- args[1]
  cohort_Rdata <- args[2]
  sampleNameWildCard <- args[3]
  fastaFile <- args[4]
  output_file <- args[5]

  if(!cohort_Rdata==""){
    load(cohort_Rdata)
    if(exists("exons.GRCh38")){exons <- exons.GRCh38}
  }else{error("Missing cohort Rdata file")}

  DF <- data.frame()
  tryCatch({

    X_sample <- getBamCounts(bed.frame = exons,
                             bam.files = BAM_file,
                             include.chr = FALSE,
                             referenceFasta = fastaFile)

    idx <- grep(".bam$", colnames(X_sample))
    names(X_sample)[idx] <- "SAMPLENAME"

    X_cohort <- merge(X_cohort,X_sample,by=c("chromosome", "start", "end","exon","GC"))


    # Build the most appropriate reference set for selected sample
    all.samples <- colnames(X_cohort)[5:length(colnames(X_cohort))]
    my.ref.samples <- setdiff(all.samples,"SAMPLENAME")

    my.test <- X_cohort[,colnames(X_cohort)=="SAMPLENAME"]
    my.reference.set <- as.matrix(X_cohort[, my.ref.samples])
    my.choice <- select.reference.set (test.counts = my.test,
                                       reference.counts = my.reference.set,
                                       bin.length = (X_cohort$end - X_cohort$start)/1000,
                                       n.bins.reduced = 10000)

    my.matrix <- as.matrix( X_cohort[, my.choice$reference.choice, drop = FALSE])
    my.reference.selected <- apply(X = my.matrix,
                                   MAR = 1,
                                   FUN = sum)

    all.exons <- new('ExomeDepth',
                     test = my.test,
                     reference = my.reference.selected,
                     formula = 'cbind(test, reference) ~ 1')

    # CNV calling
    all.exons <- CallCNVs(x = all.exons,
                          transition.probability = 10^-4,
                          chromosome = X_cohort$chromosome,
                          start = X_cohort$start,
                          end = X_cohort$end,
                          name = X_cohort$exon)


    DF <- all.exons@CNV.calls
    DF$length <- DF$end-DF$start
    DF <- subset(DF, select=c("chromosome","start","end", "length",  "id", "type" ,        
                              "nexons","BF",
                              "reads.expected","reads.observed","reads.ratio"
    ))


  },
  error = function(e) {
    message("Error detected: ", e$message)
    message("Running fallback instead...")

  })
    if(nrow(DF)>0){
    dir.create(file.path(dirname(output_file)), showWarnings = FALSE)
    fwrite(DF,file = output_file,sep = "\t")
  }else{
    DF <- data.frame(chromosome=character(0),start=character(0),end=character(0),length=character(0),id=character(0),type=character(0),
                     nexons=character(0),BF=character(0),reads.expected=character(0),reads.observed=character(0),reads.ratio=character(0))
    dir.create(file.path(dirname(output_file)), showWarnings = FALSE)
    fwrite(DF,file = output_file,sep = "\t")
  }
}

script_dir <- dirname(sub("--file=", "", commandArgs()[grep("--file=", commandArgs())]))
args <- commandArgs(trailingOnly = T)
run_all(args)