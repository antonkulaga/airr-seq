version development

import "https://raw.githubusercontent.com/antonkulaga/bioworkflows/main/common/files.wdl" as files


workflow changeo_igblast {
    input {
        File fastq
        Int threads
        String name
        String destination
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

    call files.copy as copy {
        input: destination = destination, files = [changeo.out]
    }

    output {
        File out = copy.out[0]
        File airr_tsv = changeo.airr_tsv
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
        -r /usr/local/share/germlines/imgt/human/vdj/ --outdir ~{outdir} --outname ~{name} --extended --format ~{format}

        ParseDb.py select -d  ~{outdir}/~{name}_db-pass.tsv -f productive -u T TRUE --outname ~{name}_f
    }

    runtime {
        docker: "immcantation/suite:devel"
    }

    output {
        File out = outdir
        File airr_tsv = outdir + "/" + name + "_db-pass.tsv"
        File functional = outdir + "/" + name + "_f_parse-select.tsv"
    }

}