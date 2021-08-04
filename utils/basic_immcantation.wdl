version 1.0

workflow basic_immcantation {
    input {
        File fasta
        Int threads = 6
    }

    call igblast {
        input: fasta = fasta, threads = threads
    }

    call changeo{
        input: fasta = fasta, fmt7 = igblast.fmt7
    }
}

task igblast {
    input {
        File fasta
        Int threads
    }
    command {
        mkdir -p igblast
        AssignGenes.py igblast -s ~{fasta} \
            -b /usr/local/share/igblast --organism human --loci ig \
            --format blast --outdir igblast --nproc ~{threads}
    }
    runtime {
        docker: "immcantation/suite:devel"
    }

    output {
        File results = "igblast"
        File fmt7 = "igblast/input_igblast.fmt7"
    }
}

task changeo {
    input {
        File fasta
        File fmt7
    }

    command {
        mkdir -p changeo
        MakeDb.py igblast \
        -s ~{fasta} -i ~{fmt7} \
        -r /usr/local/share/germlines/imgt/human/vdj/ --outdir changeo \
        --regions --scores --outname data

        ParseDb.py select -d changeo/data_db-pass.tab -f FUNCTIONAL -u T TRUE --outname data_f
        ParseDb.py select -d changeo/data_f_parse-select.tab -f V_CALL -u IGHV --regex --outname data_fh
    }

    runtime {
        docker: "immcantation/suite:devel"
    }

    output {
        File results = "changeo"
        File? heavy_chains = "changeo/data_f"
        File? light_chains = "changeo/data_fh"
    }

}