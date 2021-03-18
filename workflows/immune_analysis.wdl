version development

import "https://raw.githubusercontent.com/antonkulaga/bioworkflows/main/workflows/extract_run.wdl" as extractor

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
        call changeo_igblast {
            input: sequence = sequence,
                species = species,
                prefix = prefix,
                ig = ig,
                only_functional = only_functional,
                partial_alignments = partial_alignments,
                format = format
        }

        call extractor.copy as copy_changeo_igblast {
            input:
                files = [changeo_igblast.out], destination = destination + "/" + "changeo"
        }
    }



}


task changeo_igblast {
    input {
        File sequence
        String species
        Boolean ig
        String prefix
        Boolean only_functional
        Boolean partial_alignments
        String format
    }

    #changeo-igblast -h
    #Usage: changeo-igblast [OPTIONS]
    #-s  FASTA or FASTQ sequence file.
    #-r  Directory containing IMGT-gapped reference germlines.
    #Defaults to /usr/local/share/germlines/imgt/human/vdj when species is human.
    #Defaults to /usr/local/share/germlines/imgt/mouse/vdj when species is mouse.
    #-g  Species name. One of human or mouse. Defaults to human.
    #-t  Receptor type. One of ig or tr. Defaults to ig.
    #-b  IgBLAST IGDATA directory, which contains the IgBLAST database, optional_file
    #and auxillary_data directories. Defaults to /usr/local/share/igblast.
    #-n  Sample identifier which will be used as the output file prefix.
    #Defaults to a truncated version of the sequence filename.
    #-o  Output directory. Will be created if it does not exist.
    #Defaults to a directory matching the sample identifier in the current working directory.
    #-f  Output format. One of airr (default) or changeo. Defaults to airr.
    #-p  Number of subprocesses for multiprocessing tools.
    #Defaults to the available cores.
    #-k  Specify to filter the output to only productive/functional sequences.
    #-i  Specify to allow partial alignments.
    #-h  This message.
    command {
        changeo-igblast -s ~{sequence} -g ~{species} -t ~{if(ig) then "ig" else "tr"}  -n ~{prefix} -f ~{format}  ~{if(only_functional) then "-k" else ""} ~{if(partial_alignments) then "-i" else ""}
    }

    runtime {
        docker: "immcantation/suite:devel"
    }

    output {
        File out = prefix
        File? airr_tsv = prefix + "/" + prefix+"_db-pass.tsv"
        File? fmt7 = prefix + "/" + prefix+"_igblast.fmt7"
    }
}