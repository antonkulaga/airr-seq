version development

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
        File fmt7 = "igblast" + "/" + basename(fasta, ".fasta")+"_igblast.fmt7"
    }
}


task changeo {
    input {
        File fasta
        File fmt7
        String name
        String format = "airr"
        String outdir = "changeo_igblast"
    }

    command {
        mkdir -p changeo_igblast
        MakeDb.py igblast \
        -s ~{fasta} -i ~{fmt7} \
        -r /usr/local/share/germlines/imgt/human/vdj/ --outdir ~{outdir}--outname ~{name} --extended --format ~{format}

        ParseDb.py select -d  ~{outdir}/~{name}_db-pass.tab -f FUNCTIONAL -u T TRUE --outname ~{name}_f
    }

    runtime {
        docker: "immcantation/suite:devel"
    }

    output {
        File out = outdir
    }

}