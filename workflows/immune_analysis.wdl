version development

import "https://raw.githubusercontent.com/antonkulaga/airr-seq/main/workflows/immcantation.wdl" as imm

workflow immune_analysis{
    input {
        Array[File] sequences
        String species = "human"
        Int threads = 6
        Boolean ig = true
        Boolean only_functional = false
        Boolean partial_alignments = false
        String format = "airr"
        String destination
    }

    scatter(sequence in sequences){

        String prefix = basename(sequence)

        call imm.immcantation {
            input: sequence = sequence,
                species = species,
                name = prefix,
                ig = ig,
                only_functional = only_functional,
                partial_alignments = partial_alignments,
                format = format,
                destination = "/" + "immcantation"
        }
    }
}