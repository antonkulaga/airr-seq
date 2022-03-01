version development

import "https://raw.githubusercontent.com/antonkulaga/bioworkflows/main/common/files.wdl" as files
import "https://raw.githubusercontent.com/antonkulaga/airr-seq/main/structural/absolut/absolut_rep.wdl" as absolut


workflow absolut_rep_many {
    input {
        String antigen
        File antigen_structures
        Array[File] sequences
        Int threads = 6
        String destination
    }

    scatter(sequence in sequences) {
        call absolut.absolut_rep{
        input:
            antigen = antigen,
            sequences = sequence,
            antigen_structures = antigen_structures,
            threads = threads,
            destination = destination + "/"+ basename(sequence, ".tsv")
        }
    }
}