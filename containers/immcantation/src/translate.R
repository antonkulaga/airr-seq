#!/usr/bin/env Rscript
debug <- FALSE

include = function(pkg){
  if(!suppressMessages(require(pkg, character.only = TRUE)))
  install.packages(pkg, character.only = TRUE)
  suppressMessages(library(pkg, pkg, character.only = TRUE))
}

include("docopt")
include("stringr")
include("alakazam")

translate_db = function(db_path){
    db <- readChangeoDb(db_path)
    without_gaps <- gsub("...", "", db$sequence_alignment, fixed=T)
    return(translateDNA(without_gaps))
}

with_translation = function(db_path){
    db <- readChangeoDb(db_path)
    without_gaps <- gsub("...", "", db$sequence_alignment, fixed=T)    
    return (cbind(db,Translation=translateDNA(without_gaps)))     
}

doc <- 'Usage:
  translate.R [--wd <wd>][--suffix <suffix>] <dbs> ...

  Options:   
   -w --wd <wd> [default: TRUE]
   -s --suffix <suffix> [default: _with_translation].
   -h --help     Show this screen.'

if(debug == TRUE) {
    ramos <- file.path("/data/samples/AIRR-Seq/ramos/test/merged")
    dbs <- file.path(ramos, c('heavy.tsv','light.tsv'))
    args <- union(c(TRUE, "_with_translation"), dbs)
    print(args)
    values <- docopt(doc, args = args, version="0.1")
} else {
    values <- docopt(doc, version="0.1")
}


wd <- values$wd
suffix <- values$suffix
dbs_pathes <- values$dbs
dbs_translated_pathes <- str_replace(dbs_pathes, ".tsv", paste0(suffix,".tsv"))
if(wd)
    dbs_translated_pathes <-file.path(getwd(), basename(dbs_translated_pathes))

print(str_interp("extending ${dbs_pathes} with Translation column!"))


for(i in 1:length(dbs_pathes)){
    path <- dbs_pathes[i]
    path_with_translation <- dbs_translated_pathes[i]
    db_with_translation <- with_translation(path)
    writeChangeoDb(db_with_translation, path_with_translation)
}

if(wd == TRUE){
    print("Saving files to working directory...")
}
print(str_interp("Execution successfully finished, file are saved as ${dbs_translated_pathes}"))