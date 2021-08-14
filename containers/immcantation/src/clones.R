#!/usr/bin/env Rscript
include <- function(pkg) {
  if (!suppressMessages(require(pkg, character.only = TRUE)))
    install.packages(pkg, character.only = TRUE)
  suppressMessages(library(pkg, pkg, character.only = TRUE))
}
include("docopt")
include("stringr")
include("alakazam")
include("scoper")
include("dplyr")

doc <- "Usage:
  clones.R [--name <name>] [--threads <threads>] [--binwidth <binwidth>] [--wd <wd>][--suffix <suffix>] <tsv>

  Options:
   -w --wd <wd> [default: .]
   -n --name <name> Name of the sample
   -t --threads <threads> Number of threads [type: int] [default: 8]
   -b --binwidth Width of the bin in plotting [type: num] [default: 0.02]
   -s --suffix <suffix> [default: _with_clones].
   -h --help     Show this screen."

debug <- TRUE
if (debug == TRUE) {
  tsv <- file.path("/data", "samples", "AIRR-Seq", "OURS", "S3987Nr1", "S3987Nr1-PBMC_heavy", "changeo_igblast", "S3987Nr1-PBMC_heavy_f_parse-select_with_translation.tsv")  # 'clones', 'changeo_clone',  'our_pbmc_germ-pass.tsv')
  args <- paste("--suffix _with_clones --name S3987Nr1-PBMC_heavy", tsv)
  print(args)
  values <- docopt(doc, args = args, version = "0.1")
  print(values)
} else {
  values <- docopt(doc, version = "0.1")
}

tsv <- values$tsv
name <- values$name
binwidth <- as.numeric(values$binwidth)
threads <- as.numeric(values$threads)

db <- readChangeoDb(tsv)

print(paste("starting spectral analyzis with ", threads, "threads"))

# Clonal assignment using identical nucleotide sequences
results <- spectralClones(db, "vj", nproc = threads)

print(paste("writing vjl groups", file.path(paste0(name, "_vjl_groups", ".tsv"))))
write.table(results@vjl_groups, file = file.path(paste0(name, "_vjl_groups", ".tsv")), quote = FALSE, row.names = FALSE, sep = "\t")
print(paste("writing spectral analyzes results", file.path(paste0(name, values$suffix, ".tsv"))))
writeChangeoDb(results@db, file.path(paste0(name, values$suffix, ".tsv")))

print(paste("writing spectral analyzes picture", file.path(paste0(name, values$suffix, ".png"))))
png(file = paste0(name, values$suffix, ".png"), width = 800, height = 600)
if (binwidth > 0) {
  print(paste("binwidth is", binwidth))
  plot(results, binwidth = binwidth)
} else plot(results)
if (debug == TRUE) {
  dev.off()
}

print("Computing clonotypes counts")
counts <- countClones(results@db, copy = "duplicate_count")
write.table(counts, file = file.path(paste0(name, "_clone_counts", ".tsv")), quote = FALSE, row.names = FALSE, sep = "\t")
print(paste("Clonotypes counts are written to",file.path(paste0(name, "_clone_counts", ".tsv"))))

print("Calculates abundancy with 95% confidence interval via 200 bootstrap realizations")
curve <- estimateAbundance(results@db, ci = 0.95, nboot = 200, clone = "clone_id", progress = TRUE)
write.table(curve@abundance, file = file.path(paste0(name, "_abundance_curve", ".tsv")), quote = FALSE, row.names = FALSE, sep = "\t")
print(paste("abundances curve values are written to",file.path(paste0(name, "_abundance_curve", ".tsv"))))

sample_colors <- c(`-1h` = "seagreen", `+7d` = "steelblue")
png(file = paste0(name, "_abundance_curve", ".png"), width = 800, height = 600)
# Plots a rank abundance curve of the relative clonal abundances
plot(curve, colors = sample_colors, legend_title = name)
if (debug == TRUE) {
  dev.off()
}
print(paste("abundances curve chart is written to",file.path(paste0(name, "_abundance_curve", ".png"))))