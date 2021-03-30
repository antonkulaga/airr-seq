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
    Array[File] airrs_failed = select_all(igblast.airr_fail_tsv)

    call files.copy as copy_airrs {
        input: files = airrs, destination = destination + "/" + "airrs"
    }

    call files.copy as copy_failed {
        input: files = airrs_failed, destination = destination + "/" + "airrs" + "/" + "failed"
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
        input: files = [split_chains.heavy, split_chains.heavy_productive, split_chains.heavy_not_productive, split_chains.light, split_chains.light_productive, split_chains.light_not_productive],suffix = "_with_translation"
    }

    call files.copy as copy_chains{
        input: files = translate.out, destination = destination
    }

    output {
        String airrs_folder =  destination + "/" + "airrs"
        File merged = copy_merged.out[0]
        Array[File] chains = copy_chains.out
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
        ParseDb.py split -d ~{link} heavy.tsv -f productive
        mv productive-T.tsv heavy_productive.tsv
        mv productive-F.tsv heavy_not_productive.tsv
        ParseDb.py select -d ~{link} -f v_call -u "IG[LK]" --logic all --regex --outname light -f productive -u T
        mv light_parse-select.tsv light.tsv
        ParseDb.py split -d ~{link} light.tsv -f productive
        mv productive-T.tsv light_productive.tsv
        mv productive-F.tsv light_not_productive.tsv
    }

    runtime {
        docker: "immcantation/suite:devel"
    }

    output {
        File light = "light.tsv"
        File heavy = "heavy.tsv"
        File light_productive = "light_productive"
        File light_not_productive = "light_not_productive"
        File heavy_productive = "heavy_productive"
        File heavy_not_productive = "heavy_not_productive"
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