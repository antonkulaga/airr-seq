version 1.0

workflow changeo_igblast {
    input {
        File fastq
        Int threads
        String name
    }

    call fastq_conversion {
        input: fastq = fastq
    }


    call igblast {
        input: fasta = fastq_conversion.out, threads = threads
    }

    call changeo{
        input: fasta = fastq_conversion.out, fmt7 = igblast.fmt7, name = name
    }
    output {
        File out = changeo.out
    }
}

task fastq_conversion {
    input{
        File fastq
    }
    command {
        fastq2fasta.py ~{fastq}
    }
    runtime {
        docker: "immcantation/suite:devel"
    }
    output {
        File out = basename(fastq, ".fastq") + ".fasta"
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
        File out = "igblast"
        File fmt7 = "igblast/"+basename(fasta, "_igblast.fasta")+".fmt7"
    }
}


task changeo {
    input {
        File fasta
        File fmt7
        String name
    }

    command {
        mkdir -p changeo
        MakeDb.py igblast \
        -s ~{fasta} -i ~{fmt7} \
        -r /usr/local/share/germlines/imgt/human/vdj/ --outdir changeo \
        --regions --scores --outname ~{name}  --extended

        ParseDb.py select -d changeo/~{name}_db-pass.tab -f FUNCTIONAL -u T TRUE --outname ~{name}_f
    }

    runtime {
        docker: "immcantation/suite:devel"
    }

    output {
        File out = "changeo"
    }

}