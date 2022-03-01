version development

import "https://raw.githubusercontent.com/antonkulaga/bioworkflows/main/common/files.wdl" as files


workflow absolut_rep {
    input {
        String antigen
        File antigen_structures
        File sequences
        Int threads = 12
        String destination
    }

    call absolut_repertoire {
        input: antigen = antigen, sequences = sequences, antigen_structures = antigen_structures, threads = threads
    }
    call absolut_features{
        input: antigen = antigen, binding = absolut_repertoire.out, output_name = "features_" + antigen + "_"+ basename(sequences)
    }

    call files.copy as copy{
        input: destination = destination, files = [absolut_repertoire.out, absolut_features.out]
     }

    output {
        File binding = copy.out[0]
        File features = copy.out[1]
    }
}

task absolut_repertoire {
    input {
        String antigen
        File antigen_structures
        File sequences
        Int threads
    }
    String structures_name = basename(antigen_structures)
    String rep_name = basename(sequences)

    command {
        ln -s ~{antigen_structures} ~{structures_name}
        ln -s ~{sequences} ~{rep_name}
        Absolut repertoire ~{antigen} ~{rep_name} ~{threads}
    }

    runtime {
        docker: "quay.io/comp-bio-aging/absolut:latest"
    }

    output {
        File out = antigen+"FinalBindings_Process_1_Of_1.txt"
    }
}

task absolut_features {
    input {
        File binding
        String antigen
        String output_name
    }


    command {
        Absolut getFeatures ~{antigen} ~{binding} ~{output_name}
    }

    runtime {
        docker: "quay.io/comp-bio-aging/absolut:latest"
    }

    output {
        File out = output_name
    }
}