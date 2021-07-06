version development
import "https://raw.githubusercontent.com/antonkulaga/bioworkflows/main/common/files.wdl" as files

workflow trust_index {
    input {
        File genome
        File gtf
        File gene_names
        String filename = "bcrtcr.fa"
        String destination
    }
    call trust_index{
        input: fasta = genome, gtf = gtf, gene_names = gene_names, filename = filename
    }


}

task trust_index {
    input {
        File fasta
        File gtf
        File gene_names
        String filename
    }

    #To generate the file specified by "-f", you need the reference genome of the species you are interested in and corresponding genome annotation GTF file. Then you can use command
    command {
        perl /usr/local/bin/BuildDatabaseFa.pl ~{fasta} ~{gtf} ~{gene_names} > bcrtcr.fa
    }

    output {
        File out = "bcrtcr.fa"
    }

    runtime {
        docker: "quay.io/biocontainers/trust4:1.0.4--h2e03b76_0"
    }
}