version development

workflow cellranger{
    input {
        File transcriptome
        String name
        File fastq_folder
        String prefix
        Int threads = 12
        String destination


    }

    call count{
    input:
           name = name,
           transcriptome = transcriptome,
           fastq_folder= fastq_folder,
           threads = threads,
            prefix = prefix
    }
    call copy{
        input: files = [count.out], destination = destination
    }

    output {
        Array[File] out = copy.out
    }

}

task count{

    input {
        String name
        File transcriptome
        File fastq_folder
        String prefix
        String? description
        Int threads
    }

    command {
        cellranger count --id ~{name} --transcriptome ~{transcriptome} ~{"--description" + description} --fastqs ~{fastq_folder} --localcores ~{threads}
    }

    runtime {
        docker: "docker.io/cumulusprod/cellranger:6.0.0"
    }

    output {
        File out = name
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

