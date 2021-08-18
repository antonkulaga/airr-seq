version development

import "https://raw.githubusercontent.com/antonkulaga/bioworkflows/main/common/files.wdl" as files

workflow anarci {
    input {
        String species = "human"
        String scheme = "martin"
        File? fasta
        String? sequence #either sequence or fasta
        String name
        Int threads
        String restrict
        Boolean assign_germline
        String destination
    }

    call anarci{
        input: species = species, scheme = scheme, fasta = fasta,
            name = name,  threads = threads, sequence = sequence,
            restrict = restrict, assign_germline  = assign_germline
    }

    call files.copy as copy {
        input: destination = destination, files = [anarci.out]
    }

    output {
        File out = copy.destination_folder
    }
}

task anarci {
    input {
        String species = "human"
        String scheme = "martin"
        File? fasta
        String? sequence #either sequence or fasta
        String name
        Int threads = 1
        String restrict
        Boolean assign_germline = false
        String output_prefix = "anarci_results"
    }

    String output_dir =  output_prefix + "_"  + name

    command {
        # Make output directory and empty log files
        mkdir -p ~{output_dir}; cd ~{output_dir}

        ANARCI --sequence ~{fasta} ~{sequence} --outfile ~{name} --scheme ~{scheme} \
        --csv --outfile_hits ~{name}_hit.txt --ncpu ~{threads} \
        ~{if(assign_germline) then "--assign_germline" else ""} -r ~{restrict}
    }

    runtime {
        docker: "quay.io/biocontainers/anarci:2021.02.04--pyhdfd78af_0"

    }

    output {
        File out = output_dir
    }
}