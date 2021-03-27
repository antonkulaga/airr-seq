version development

import "https://raw.githubusercontent.com/antonkulaga/airr-seq/main/analysis/immcantation.wdl" as imm
import "https://raw.githubusercontent.com/antonkulaga/bioworkflows/main/common/files.wdl" as files


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
        Boolean include_tigger = true
        String name = ""
    }


    scatter(sequence in sequences){

        String prefix = basename(sequence)

        call imm.changeo_igblast as igblast {
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

    Array[File] airrs = select_all(igblast.airr_tsv)

    call files.copy {
        input: files = airrs, destination = destination / "airrs"

    }

    String merged_name = if(name=="") then "merged.tsv" else name + ".tsv"
    call merge_airr {
        input: files =  airrs, filename = merged_name
    }

    call files.copy {
        input: files = [merge_airr.out], destination = destination
    }


}


task merge_airr {
    input {
        Array[File] files
        String filename = ".tsv"
    }

    command {
        ParseDb.py merge -d ~{sep=' ' files} -o ~{filename}
    }

    runtime {
        docker: "immcantation/suite:devel"
    }

    output {
        File out = filename
    }
}