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
  clones.R analyze_clones [--name <name>] [--threads <threads>] [--binwidth <binwidth>] [--wd <wd>][--suffix <suffix>] [--spectral_method <spectral_method>] <tsv> 
  clones.R analyze_diversity [--name <name>] [--wd <wd>] <clones_tsv>  

  Options:   
   -w --wd <wd> [default: .]
   -n --name <name> Name of the sample
   -t --threads <threads> Number of threads [type: int] [default: 8]
   -b --binwidth Width of the bin in plotting [type: num] [default: 0.02]
   -s --suffix <suffix> [default: _with_clones].
   -sm --spectral_method <spectral_method> [default: novj]
   -h --help     Show this screen."

debug <- FALSE
debug_command <- "analyze_clones"  # {"analyze_clones", "analyze_diversity"}
if (debug == TRUE) {
    if (debug_command == "analyze_clones") {
        tsv <- file.path("/data", "samples", "AIRR-Seq", "OURS", "S3987Nr1", "S3987Nr1-PBMC_heavy", "changeo_igblast", "S3987Nr1-PBMC_heavy_f_parse-select_with_translation.tsv")  # 'clones', 'changeo_clone',  'our_pbmc_germ-pass.tsv')    
        args <- paste("analyze_clones --suffix _with_clones --name S3987Nr1-PBMC_heavy", tsv)
    } else if (debug_command == "analyze_diversity") {
        tsv <- file.path("/data", "samples", "AIRR-Seq", "OURS", "S3987Nr1", "S3987Nr1-PBMC_heavy", "subsamples", "sample_1000", "clones", "S3987Nr1-PBMC_heavy_with_clones.tsv")   
        args <- paste("analyze_diversity --name S3987Nr1-PBMC_heavy", tsv)   
    } else {
        stop(paste0("Debug command ", debug_command, " is not supported."))
    }
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
spectral_method <- values$spectral_method

clones_tsv <- values$clones_tsv
analyze_clones <- values$analyze_clones
analyze_diversity <- values$analyze_diversity

# Helper functions
build_filepath <- function(name, suffix, type="tsv") {
    if (type %in% c("tsv", "svg")) {
        return(
            file.path(paste0(name, "_", suffix, ".", type))
        )
    } else {
        stop(paste0("Filetype ", type, " not supported."))
    }
}

writeAnalysisTable <- function(table, filepath) {
    write.table(
        table, 
        file = filepath, 
        quote = FALSE, 
        row.names = FALSE, 
        sep = "\t")
    print(paste0("Output is written to ", filepath))
}

computeSpectralClones <- function(
    db, threads, spectralClones_method = "novj", verbose = TRUE) {
    
    results <- tryCatch({
            scoper::spectralClones(db, spectralClones_method, verbose = verbose, nproc = threads)
            },
            warning = function(cond) {
                iter_max_warning = 15000
                nstart_warning = 5000
                print(paste0(
                    "Warning issued. Retrying with iter_max=", iter_max_warning,
                    "; nstart=", nstart_warning
                ))
                scoper::spectralClones(
                    db,
                    spectralClones_method,
                    verbose = TRUE, 
                    iter_max = iter_max_warning,
                    nstart = nstart_warning
                )
        })
    return(results)
}

if (analyze_clones == TRUE) {
    db <- alakazam::readChangeoDb(tsv)

    print(paste("starting spectral analyzis with ", threads, "threads"))

    # Clonal assignment using identical nucleotide sequences
    results <- computeSpectralClones(db, threads, spectral_method, verbose=TRUE)
    
    groups_fp = build_filepath(name, paste0(spectral_method, "_groups"), "tsv")
    writeAnalysisTable(results@vjl_groups, groups_fp)
    
    changeo_fp = file.path(paste0(name, paste0("_", spectral_method, values$suffix), ".tsv"))
    print(paste("writing spectral analyzes results", changeo_fp))
    writeChangeoDb(results@db, changeo_fp)
    
    plot_fp = file.path(paste0(name, paste0("_", spectral_method, values$suffix), ".svg"))
    print(paste("writing spectral analyzes picture", plot_fp))
    svg(file = plot_fp, width = 800, height = 600)
    if (binwidth > 0) {
        print(paste("binwidth is", binwidth))
      plot(results, binwidth = binwidth)
    } else plot(results)
    if (debug == TRUE) {
      dev.off()
    }
}

computeCloneCounts <- function(clones_results, name) {
    print("Computing clonotypes counts")
    counts <- alakazam::countClones(clones_results, copy = "duplicate_count")
    writeAnalysisTable(counts, build_filepath(name, "clone_counts", "tsv"))
    return(counts)
}

computeCoverage <- function(counts) {
    
    computeCoverage_subroutine <- function(max_order, counts) {
        orders <- 1:max_order
        coverages <- sapply(
            orders,
            function(order) { alakazam::calcCoverage(counts$seq_count, order) }
        )
        return(coverages)
    }
    
    tryCatchComputeCoverage <- function(counts, env) {
        coverages <- tryCatch({
                order <- get("order", env=env)
                return(computeCoverage_subroutine(order, counts))
            },
            error=function(cond) {
                order <- get("order", env=env)
                assign("order", order-1, env=env)
                print(paste0("Attempted coverage computation with order=", order,
                             "; Encountered error: ", cond, "; Attempting with order=", order-1))
                return(tryCatchComputeCoverage(counts, env))
            }
        )
        return(coverages)
    }
    
    compute_coverage_env <- new.env()
    assign("order", 10, env=compute_coverage_env)
    coverages <- tryCatchComputeCoverage(counts, compute_coverage_env)

    order <- get("order", env=compute_coverage_env)
    df <- data.frame(1:order, coverages)
    writeAnalysisTable(df, build_filepath(name, "coverages", "tsv"))
}

plotAbundancy <- function(abundanceCurve, name, filepath, debug) {
    sample_colors <- c(`-1h` = "seagreen", `+7d` = "steelblue")
    svg(file = filepath, width = 800, height = 600)
    # Plots a rank abundance curve of the relative clonal abundances
    plot(abundanceCurve, colors = sample_colors, legend_title = name)
    if (debug == TRUE) {
      dev.off()
    }
    print(paste("abundances curve chart is written to", filepath))
    plot(abundanceCurve, colors = sample_colors, legend_title = name)  # to display in notebook
}
computeAbundancy <- function(clones_results, name, debug) {
    print("Calculates abundancy with 95% confidence interval via 200 bootstrap realizations")
    abundanceCurve <- alakazam::estimateAbundance(clones_results, ci = 0.95, nboot = 200, clone = "clone_id", progress = TRUE)
    writeAnalysisTable(abundanceCurve@abundance, build_filepath(name, "abundance_curve", "tsv"))
    plotAbundancy(abundanceCurve, name, build_filepath(name, "abundancy_curve", "svg"), debug)
    return(abundanceCurve)
}

plotDiversity <- function(diversity, filepath, debug) {
    svg(file = filepath, width = 800, height = 600)
    plot(diversity)
    if (debug == TRUE) {
      dev.off()
    }
    print(paste("diversity curve is written to", filepath))
    plot(diversity)
}
computeDiversity <- function(abundancyCurve, name, debug) {
    diversity <- alakazam::alphaDiversity(abundancyCurve)
    writeAnalysisTable(diversity@diversity, build_filepath(name, "diversity", "tsv"))
    plotDiversity(diversity, build_filepath(name, "diversity", "svg"), debug)
    return(diversity)
}

if (analyze_diversity == TRUE) {
    results <- alakazam::readChangeoDb(clones_tsv)
    counts <- computeCloneCounts(results, name)
    computeCoverage(counts)
    
    tryCatch({
        abundanceCurve <- computeAbundancy(results, name, debug)
        diversity <- computeDiversity(abundanceCurve, name, debug)
        },
        error = function(cond) {
            print(cond)
            print("Further execution halted.")
        }
    )
}
