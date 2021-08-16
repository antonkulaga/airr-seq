# Helper functions
build_filepath <- function(name, suffix, type="tsv") {
    if (type %in% c("tsv", "png")) {
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