#!/usr/bin/env Rscript

# myfile = file("stdin")
# open(myfile, blocking=TRUE)
# myinput = readLines(myfile,n = 1)
# print(myinput)

result = tryCatch({
  dt <- data.table::fread('file:///dev/stdin')
  res <- apply(dt[,10:13,with=F],1,function(x){
    res <- fisher.test(matrix(x, nrow=2))
    return(c(round(res$p.value, 5),round(res$estimate, 5)))
  })
  
  write.table(data.frame(dt[,1:20], res[1,], res[2,], dt[,21:dim(dt)[2]]), file = "", quote = F, sep = "\t", eol = "\n", row.names=F, col.names=F)
  return("")
},error = function(e) {
  return("")
})


# myfile = file("stdin")
# open(myfile, blocking=TRUE)
# myinput = readLines(myfile) # read from stdin
# if (length(myinput) > 0 ){
#     mynumcols = sapply(gregexpr("\\t+", myinput[1]), length) + 1 # count num of tabs + 1
# }else{
#     mynumcols = 0
#     d = matrix(0,0,0)
# }
# if (mynumcols == 34 || mynumcols == 38){ # 34 columns for standard bed files, 38 for amplicon mode
#     d <- read.table( textConnection(myinput), sep = "\t", header = F, colClasses=c("character", NA, NA, NA, NA, "character", "character", NA, NA, NA, NA, NA, NA, "character", NA, NA, NA, NA, NA, NA, NA, NA), col.names=c(1:mynumcols))
# } else if (mynumcols > 0){
#     stop("Incorrect input detected in teststrandbias.R")
# }
# 
# if (nrow(d) > 0){
#     pvalues <- vector(mode="double", length=dim(d)[1])
#     oddratio <- vector(mode="double", length=dim(d)[1])
#     
#     for( i in 1:dim(d)[1] ) {
#         h <- fisher.test(matrix(c(d[i,10], d[i,11], d[i,12], d[i,13]), nrow=2))
#         pvalues[i] <- round(h$p.value, 5)
#         oddratio[i] <- round(h$estimate, 5)
#     }
#     write.table(data.frame(d[,1:20], pvalues, oddratio, d[,21:dim(d)[2]]), file = "", quote = F, sep = "\t", eol = "\n", row.names=F, col.names=F)
# }
