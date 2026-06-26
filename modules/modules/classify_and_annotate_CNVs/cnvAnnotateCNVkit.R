
library(data.table)
library(writexl)

#nemelo by byt napevno... ale s detekci UCSC/Ensembl
listofCHR<-c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18",
             "chr19","chr20","chr21","chr22","chrY","chrX")

listofCHR2<-c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18",
              "19","20","21","22","X","Y")

# develop and test
# setwd("/media/rj/SSD_500GB/CNV_OVARIA/FFPE/")
# args <- c("variant_calls/1456-22_FFPE/cnvkit/cnvkit_CNV_preannot_1456-22_FFPE.tsv","variant_calls/1456-22_FFPE/cnvkit/1456-22_FFPE.classified.txt","processed_GFT.tsv","cnvkit_final/1456-22_FFPE_final_CNVs_annotated.tsv","cnvkit_final/1456-22_FFPE_final_CNVs_annotated.xlsx","")


run_all <- function(args){
  # arguments
  input_tsv <- args[1]
  classify_txt <- args[2]
  gtf_tsv <- args[3]
  output_tsv <- args[4]
  output_xlsx <- args[5]

  ##############################################################################
  inputDF <- fread(file=input_tsv, sep='\t', header = TRUE)
  names(inputDF)
  
  if(all(c("CHROM","POS","END","SVTYPE") %in% names(inputDF))){
    setnames(inputDF,c("CHROM","POS","END","SVTYPE"),c("CHR","START","STOP","CNVtype"))
    inputDF$VariantID <- paste(inputDF$CHR,inputDF$START,inputDF$STOP,inputDF$CNVtype,sep = "_") # chr1_970657_970704_DEL
  }
  if(all(c("CHROM","POS","END","CNVtype") %in% names(inputDF))){
    setnames(inputDF,c("CHROM","POS","END","CNVtype"),c("CHR","START","STOP","CNVtype"))
    inputDF$VariantID <- paste(inputDF$CHR,inputDF$START,inputDF$STOP,inputDF$CNVtype,sep = "_") # chr1_970657_970704_DEL
  }
  if(all(c("CHROM","START","END","consensus_type") %in% names(inputDF))){
    setnames(inputDF,c("CHROM","END","consensus_type"),c("CHR","STOP","CNVtype"))
    inputDF$VariantID <- paste(inputDF$CHR,inputDF$START,inputDF$STOP,inputDF$CNVtype,sep = "_") # chr1_970657_970704_DEL
  }
  
  ##############################################################################
  ## DETECT USCS or ENSEMBL
  if( grepl( "^chr",     inputDF[1,1], fixed = TRUE) ){
    UCSCorENSEMBL <- "ENS"
  }else{UCSCorENSEMBL <- "UCSC"}
  
  
  ##############################################################################
  # PROCESS GTF
  gtf_df_all <- fread(file=gtf_tsv, sep='\t', header = TRUE)
  
  if(UCSCorENSEMBL=="UCSC"){
    gtf_df_all$chromosome <- gsub("^chr","",gtf_df_all$chromosome) # carefull
    gtf_df_all <- gtf_df_all[chromosome %in% listofCHR2,]
  }else{
    gtf_df_all <- gtf_df_all[chromosome %in% listofCHR,]
  }
  
  setDT(gtf_df_all)



  ##############################################################################
  classifyCNV_res <- fread(file=classify_txt, sep='\t', header = TRUE)
  if("Chromosome" %in% names(classifyCNV_res)){setnames(classifyCNV_res,c("Chromosome"),c("CHR"))}
  if("Start" %in% names(classifyCNV_res)){setnames(classifyCNV_res,c("Start"),c("START"))}
  if("End" %in% names(classifyCNV_res)){setnames(classifyCNV_res,c("End"),c("STOP"))}
  if("Type" %in% names(classifyCNV_res)){setnames(classifyCNV_res,c("Type"),c("CNVtype"))}
  
  if(UCSCorENSEMBL=="UCSC"){
    classifyCNV_res$CHR <- gsub("^chr","",classifyCNV_res$CHR) # carefull
    classifyCNV_res <- classifyCNV_res[CHR %in% listofCHR2,]
    classifyCNV_res$VariantID <- gsub("^chr","",classifyCNV_res$VariantID)
  }else{
    classifyCNV_res <- classifyCNV_res[CHR %in% listofCHR,]
  }
  
  classifyCNV_res <- classifyCNV_res[,-c(8:42)] # delete categories
  classifyCNV_res[,c("CHR","START","STOP","CNVtype") :=NULL] # delete cause merge columns
  names(classifyCNV_res)
  ##############################################################################
  inputDF <- merge(inputDF,classifyCNV_res,by="VariantID")
  inputDF$IGV <- paste0(inputDF$CHR,":",inputDF$START,"-",inputDF$STOP) #chr1:144,874-969,268)
  
  # names(inputDF)
  # names(gtf_df_all)
  # nrow(inputDF)
  
  ##############################################################################
  # Put CHR, START, STOP at the front, keep others in their original order
  cols_to_front <- c("CHR", "START", "STOP")
  all_cols <- names(inputDF)
  remaining_cols <- setdiff(all_cols, cols_to_front)
  setcolorder(inputDF, c(cols_to_front, remaining_cols))
  


  ############################################################################################################################################################
  ## overlap inputDF and GTF
  inputDF$CHR <- as.character(inputDF$CHR)
  gtf_df_all$chromosome <- as.character(gtf_df_all$chromosome)
  
  setkey(inputDF, CHR,START,STOP)
  setkey(gtf_df_all, chromosome,start,end)
  
  overlapCNVxGTF <- foverlaps(inputDF, gtf_df_all, by.x = c("CHR","START","STOP"),
                              by.y = c("chromosome","start","end"), mult="all", type="any", nomatch=NULL, which=FALSE)
  
  setkey(overlapCNVxGTF, CHR, START,STOP)
  overlapCNVxGTF <- unique(overlapCNVxGTF, by = c("CHR","START","STOP","gene_name"))
  
  # remove those with ENSG in gene name if there are also duplicate TRUE
  # overlapCNVxGTF[, duplicate := .N > 1, by = key(overlapCNVxGTF)]
  # overlapCNVxGTF <- overlapCNVxGTF[!((gene_name %like% "ENSG") & (duplicate==TRUE)),]
  
  overlapCNVxGTF[, GFT_genes := paste(gene_name, collapse = ", "), by = VariantID]
  overlapCNVxGTF <- unique(overlapCNVxGTF, by = c("VariantID"))
  overlapCNVxGTF[,c("VariantID") :=NULL] # delete cause merge columns
  
  names(overlapCNVxGTF)
  
  desired_cols <- c("CHR","START","STOP","CNVtype","SVLEN","IGV",
                    "FOLD_CHANGE","FOLD_CHANGE_LOG","IMPRECISE","PROBES",
                    "GFT_genes",
                    "Classification","Total score","Known or predicted dosage-sensitive genes","All protein coding genes")
  present_cols <- intersect(desired_cols, names(overlapCNVxGTF))
  overlapCNVxGTF <- subset(overlapCNVxGTF, select=present_cols)
 
  overlapCNVxGTF$SVLEN <- overlapCNVxGTF$STOP-overlapCNVxGTF$START
  
  #######################################################################################################################
  dir.create(file.path(getwd(),dirname(output_tsv)), recursive = TRUE, showWarnings=FALSE)
  write.table(overlapCNVxGTF, file=output_tsv, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE )
  writexl::write_xlsx(overlapCNVxGTF,path=output_xlsx)
}


#run as Rscript
args <- commandArgs(trailingOnly = T)
run_all(args)