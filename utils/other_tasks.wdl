version development

task download_imgt{
    input {

    }

    command {
        fetch_imgtdb.sh
    }

    runtime {
        docker: "immcantation/suite:devel"
    }

    output {

    }
}


task convert_headers {
    input {
        String format = "imgt"
    }

    command {
        ConvertHeaders ~{format} -s /data/imgt/human/irep/vdj/imgt_human_IGLV.fasta
    }

    runtime {
        docker: "immcantation/suite:devel"
    }

    output {

    }
}

task clear_umi {

    input {
        Array[File] reads
        String pattern = "XXXNNNNNN"
    }

    String read_1_name =  basename(reads[0], ".fastq")
    String read_2_name =  basename(reads[1], ".fastq")

    command {
        umi_tools extract --stdin="~{reads[0]}" --read2-in="~{reads[1]}" --bc-pattern=~{pattern} \
        --stdout=~{read_1_name + "_processed.fastq"} --read2-out=~{read_2_name + "_processed.fastq"}
    }

    runtime {
        docker: "quay.io/biocontainers/umi_tools:1.1.2--py38h4a8c8d9_0"
    }

    output {
        Array[File] out = [read_1_name + "_processed.fastq", read_2_name + "_processed.fastq"]
    }
}