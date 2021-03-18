version development

workflow cell_ranger_index {
    input {
        String index_name = "results"
        File fasta
        File gtf
        String destination
    }

    call mkref {
        input: output_name = index_name, fasta = fasta, genes = gtf
    }

    call copy{
       input: files = [mkref.out], destination = destination
    }

}

task mkref {
    input {
        String output_name
        File fasta
        File genes
    }

    command {
        cellranger mkref --genome=~{output_name} --fasta=~{fasta} --genes=~{genes}
    }

    runtime {
        docker: "docker.io/cumulusprod/cellranger:6.0.0"
    }

    output {
        File out = output_name
    }

}



task copy {
    input {
        Array[File] files
        String destination
    }

    String where = sub(destination, ";", "_")

    command {
        mkdir -p ~{where}
        cp -L -R -u ~{sep=' ' files} ~{where}
        declare -a files=(~{sep=' ' files})
        for i in ~{"$"+"{files[@]}"};
        do
        value=$(basename ~{"$"}i)
        echo ~{where}/~{"$"}value
        done
    }

    output {
        Array[File] out = read_lines(stdout())
    }
}