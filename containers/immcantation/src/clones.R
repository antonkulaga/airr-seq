include = function(pkg){
 if(!suppressMessages(require(pkg, character.only = TRUE)))
   install.packages(pkg, character.only = TRUE)
 suppressMessages(library(pkg, pkg, character.only = TRUE))
}
include("docopt")
include("stringr")
include("alakazam")
include("scoper")
include("dplyr")

doc <- 'Usage:
  clones.R [--name <name>] [--binwidth <binwidth>] [--wd <wd>][--suffix <suffix>] <tsv>

  Options:
   -w --wd <wd> [default: TRUE]
   -n --name <name>
   -b --binwidth <binwidth> [default: 0.02]
   -s --suffix <suffix> [default: _with_clones].
   -h --help     Show this screen.'

debug <- FALSE
if(debug == TRUE) {
 tsv <- file.path("/data", "samples", "AIRR-Seq", "OURS", "S3987Nr1", "S3987Nr1-PBMC_heavy", "changeo_igblast", "S3987Nr1-PBMC_heavy_f_parse-select_with_translation.tsv")# "clones", "changeo_clone",  "our_pbmc_germ-pass.tsv")
 args <- paste("--suffix _with_clones --name S3987Nr1-PBMC_heavy", tsv)
 print(args)
 values <- docopt(doc, args = args, version="0.1")
} else {
 values <- docopt(doc, version="0.1")
}


tsv <- values$tsv
name <- values$name
binwidth <-values$binwidth

db <- readChangeoDb(tsv)
# Clonal assignment using identical nucleotide sequences
results <- spectralClones(db, "vj")

write.table(results@vjl_groups, file=file.path(paste0(name,"_vjl_groups",'.tsv')), quote=FALSE, row.names = FALSE, sep='\t')


png(file=paste0(name+values$suffix,".png"), width = 800, height = 600)
plot(results, binwidth=0.02)
dev.off()