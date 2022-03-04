#!/usr/bin/env Rscript

include <- function(pkg) {
  if (!suppressMessages(require(pkg, character.only = TRUE)))
    install.packages(pkg, character.only = TRUE)
  suppressMessages(library(pkg, pkg, character.only = TRUE))
}

include("alakazam")
include("shazam")
include("dplyr")
include("ggplot2")

# Helper
writeAnalysisTable <- function(table, filepath) {
    write.table(
        table, 
        file = filepath, 
        quote = FALSE, 
        row.names = FALSE, 
        sep = "\t")
    print(paste0("Output is written to ", filepath))
}

# Pipeline
clonal_analysis <- function(clones_path, name) {

    clones_filename = paste0(name, "_novj_with_clones.tsv")
    repertoire = read.csv(
        paste0(clones_path, "/", clones_filename), 
        sep='\t'
    )

    repertoire = repertoire[order(-repertoire[, "duplicate_count"]), ]
    
    # Collapse clones
    clonal_sequences = shazam::collapseClones(
        db = repertoire[, c("sequence_alignment", "germline_alignment", "clone_id")], 
        cloneColumn="clone_id", 
        sequenceColumn="sequence_alignment", 
        germlineColumn="germline_alignment", 
        regionDefinition=NULL, #shazam::IMGT_V, 
        method="thresholdedFreq", minimumFrequency=0.6,
        includeAmbiguous=FALSE, breakTiesStochastic=FALSE, 
        nproc=12
    )
    clonal_sequences_dt = data.table::data.table(clonal_sequences)
    clonal_sequences_dt$clone_id = as.integer(clonal_sequences_dt$clone_id)
    # Augment clones data.table with extra informative columns
    augment_clones = function(repertoire_dt, clonal_sequences_dt) {
        clones_dt = repertoire_dt[, 
                      list(
                          counts=sum(duplicate_count), 
                          num_seqs=length(unique((sequence_id)))
                      ), 
                      by="clone_id"][order(-counts)]

        clones_dt = merge(clones_dt, clonal_sequences_dt, by = "clone_id", all = FALSE)[order(-counts)]
        return(clones_dt)
    }

    repertoire_dt = data.table::data.table(repertoire)
    clones_dt = augment_clones(repertoire_dt, clonal_sequences_dt)
    writeAnalysisTable(clones_dt, paste0(clones_path, "/", name, "_collapse_clones.tsv"))

    # Observed mutations (sequences)
    repertoire_obs <- shazam::observedMutations(
        repertoire_dt, 
        sequenceColumn="sequence",
        germlineColumn="germline_alignment_d_mask",  # d_mask
#                             regionDefinition=shazam::IMGT_VDJ_BY_REGIONS,
        frequency=TRUE,
        combine=FALSE,
        nproc=12)
    writeAnalysisTable(repertoire_obs, paste0(clones_path, "/", name, "_novj_with_clones_and_muts.tsv"))
    
    # Selection pressure
    baseline <- shazam::calcBaseline(
        clones_dt, 
        testStatistic="focused", 
        regionDefinition=shazam::IMGT_V, 
        nproc=1, 
        calcStats = TRUE)
    writeAnalysisTable(baseline@stats, paste0(clones_path, "/", name, "_collapse_clones_with_selection_pressure.tsv"))
    
}


config = yaml::read_yaml("/data/sources/immune-repertoires-dash/config.yml")
sample_names = names(config$samples)


print(paste0("Number of samples: ", length(sample_names)))
for (i in 1:length(sample_names)) {
    sample_name = sample_names[i]
    sample_path = config$samples[[sample_name]]$sample_path
    clones_path = paste0(sample_path, "/clones")
    
    if (sample_name %in% c("S3987Nr2-PBMC1_heavy", "S3987Nr2-PBMC1_light", "S3987Nr2-RAMOS_heavy", "S3987Nr2-RAMOS_light")) {
        print(paste0("Sample in exclusion list, skipping ", sample_name))
        next
    }
    
    if (file.exists(paste0(clones_path, "/", sample_name, "_collapse_clones_with_selection_pressure.tsv"))) {
        print(paste0("Found baseline table, skipping ", sample_name))
        next
    }

    print(paste0("Clonal analysis for: ", i, " - ", sample_name))
    
    tryCatch(
        {
            clonal_analysis(clones_path, sample_name)            
        }, 
        error = function(cond) {
            message(paste0("Encountered error (message below) for: ", sample_name))
            message(cond)
        },
        warning = function(cond) {
            message(paste0("Encountered warning (message below) for: ", sample_name))
            message(cond)
        })
    
    print("\n")
}


