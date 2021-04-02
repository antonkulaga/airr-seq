version development

import "https://raw.githubusercontent.com/antonkulaga/bioworkflows/main/common/files.wdl" as files

#production version
import "https://raw.githubusercontent.com/antonkulaga/airr-seq/main/analysis/tasks.wdl" as imm
import "https://raw.githubusercontent.com/antonkulaga/airr-seq/main/analysis/clonal_analysis.wdl" as clonal

#local debug version (uncomment for debugging and comment the production version)
#import  "immcantation.wdl" as imm


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
        #Boolean include_tigger = true
        String name = "samples"
    }


    scatter(sequence in sequences){

        String prefix = basename(sequence)

        call imm.changeo_igblast as igblast {
            input: sequence = sequence,
                species = species,
                prefix = prefix,
                ig = ig,
                only_functional = only_functional,
                partial_alignments = partial_alignments,
                format = format
                #destination = "/" + "immcantation"
        }
    }


    Array[File] airrs = select_all(igblast.airr_tsv)
    Array[File] airrs_failed = select_all(igblast.airr_fail_tsv)

    call files.copy as copy_airrs {
        input: files = airrs, destination = destination + "/" + "airrs"
    }

    call files.copy as copy_failed {
        input: files = airrs_failed, destination = destination + "/" + "airrs" + "/" + "failed"
    }


    String merged_name = if(name=="") then "merged.tsv" else name + ".tsv"

    call merge_airr as merge_passed {
        input: files =  airrs, filename = merged_name
    }


    String failed_merged_name = if(name=="") then "failed_merged.tsv" else "failed_"+name + ".tsv"

    call merge_airr as merge_failed {
        input: files =  airrs_failed, filename = failed_merged_name
    }

    Array[File] airrs_merged = select_all([merge_passed.out,merge_failed.out])

    call files.copy as copy_merged{
        input: files = airrs_merged, destination = destination
    }



    #call files.copy as copy_chains{
    #    input: files = translate.out, destination = destination
    #}

    #output {
    #    String airrs_folder =  destination + "/" + "airrs"
    #    Array[File] merged = copy_merged.out
    #    Array[File] chains = copy_chains.out
    #}

}


task merge_airr {
    input {
        Array[File] files
        String filename = "merged.tsv"
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