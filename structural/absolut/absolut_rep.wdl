version development

workflow absolut_rep {
    input {
        String antigen
        File antigen_structures
        File sequences
        Int threads = 12
    }

    call absolut_repertoire {
        input: antigen = antigen, sequences = sequences, antigen_structures = antigen_structures, threads = threads
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

    }
}