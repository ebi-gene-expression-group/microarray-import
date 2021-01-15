#!/usr/bin/env Rscript

library(affyio)

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("Please provide a directory with .CEL files to check", call.=FALSE)
}

base_dir<-args[1]
files<-list.files(path=base_dir, pattern = ".(CEL|cel)")

failure<-0
for(file in files) {
  out<-tryCatch(
    f<-read.celfile(paste(base_dir, file, sep="/"), intensity.means.only = TRUE),
    error = function(c) {
      return(list(value=1,message=conditionMessage(c)))
    }
  )
  if("value" %in% names(out)) {
    failure<-1
  }
}

quit(status=failure)
