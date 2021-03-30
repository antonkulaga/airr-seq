version development

import "https://raw.githubusercontent.com/antonkulaga/bioworkflows/main/common/files.wdl" as files

#production version
import "https://raw.githubusercontent.com/antonkulaga/airr-seq/main/analysis/immcantation.wdl" as imm

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

    call files.copy as copy_airrs {
        input: files = airrs, destination = destination + "/" + "airrs"
    }


    String merged_name = if(name=="") then "merged.tsv" else name + ".tsv"

    call merge_airr {
        input: files =  airrs, filename = merged_name
    }

    call files.copy as copy_merged{
        input: files = [merge_airr.out], destination = destination
    }

    call split_chains {
        input: airr =  copy_merged.out[0]
    }
    call translate {
        input: files = [split_chains.heavy, split_chains.light]
    }

    call files.copy as copy_chains{
        input: files = [translate[1], translate[0]], destination = destination
    }

    output {
        String airrs_folder =  destination + "/" + "airrs"
        File merged = copy_merged.out[0]
        File light =  copy_chains.out[0]
        File heavy =  copy_chains.out[1]
    }

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

task split_chains {
    input {
        File airr
    }

    String link = basename(airr) # to avoid parsedb creating stuff in the input folder

    command {
        ln -s ~{airr} ~{link}
        ParseDb.py select -d ~{link} -f v_call -u "IGH" --logic all --regex --outname heavy
        mv heavy_parse-select.tsv heavy.tsv
        ParseDb.py select -d ~{link} -f v_call -u "IG[LK]" --logic all --regex --outname light
        mv light_parse-select.tsv light.tsv
    }

    runtime {
        docker: "immcantation/suite:devel"
    }

    output {
        File light = "light.tsv"
        File heavy = "heavy.tsv"
    }
}

task translate {
    input {
        Array[File] files
        String suffix = "_with_translation"
    }
    String gl = "*"+suffix+".tsv"
    command {
        translate.R --wd TRUE --suffix ~{suffix} ~{sep=" " files}
    }

    runtime {
        docker: "quay.io/comp-bio-aging/immcantation"
    }
    output {
        Array[File] out = glob(gl)
    }
}