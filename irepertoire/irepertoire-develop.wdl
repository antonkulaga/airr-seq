version development

import "https://raw.githubusercontent.com/antonkulaga/bioworkflows/main/common/files.wdl" as files

#production version
import "https://raw.githubusercontent.com/antonkulaga/airr-seq/main/analysis/changeo_igblast.wdl" as changeo
import "https://raw.githubusercontent.com/ursueugen/airr-seq/main/analysis/clonal_analysis/clonal_analysis.wdl" as clonal
import "https://raw.githubusercontent.com/ursueugen/airr-seq/main/analysis/preprocess/presto/presto.wdl" as presto


workflow irepertoire{

    input {
        Array[File] reads

        String destination
        String name

        File IGHC_fasta = ""

        String species = "human"

        String coordinates = "illumina" #sra
        Int threads = 12
        Int min_reads_per_ig = 2
        Int min_length = 196
        Int min_quality = 20
        Int start_v = 10
        Int min_dupcount = 2
        Float clones_bin_width = 0.02
        Int max_memory_gb = 96
        String spectral_method = "novj"
    }

    call fastqc {
        input: output_dir = "fastqc_results", reads = reads
    }

    call files.copy as copy_fastqc {
        input: destination = destination, files = [fastqc.results]
    }

    call presto.presto as presto {
        input: name = name, output_dir = "presto",
        reads = reads, NPROC = threads, min_quality = min_quality, min_length = min_length, start_v = start_v, dupcount = min_dupcount,
        coordinates = coordinates, max_memory = max_memory_gb, IGHC_fasta = IGHC_fasta
    }

    call files.copy as copy_presto { 
        input: destination = destination, files = [presto.results]
    }
    
    call changeo.changeo_igblast as igblast {
        input: fastq = presto.out, threads = threads, name = name, destination = destination
    }

    call clonal.clonal_analysis as clonal_analysis {
        input: airr_tsv = igblast.airr_tsv_translated_functional, 
        destination = destination, name = name, method = spectral_method,
        threads = threads, clones_bin_width = clones_bin_width, 
        max_memory_gb = max_memory_gb
    }

    output {
        File presto_results = copy_presto.out[0]
        File changeo_results = igblast.out
        File clones = clonal_analysis.clones
        File diversity = clonal_analysis.diversity
    }

}


task fastqc {
    input {
        String output_dir = "fastqc_results"
        Array[File] reads
    }

    command {
        mkdir -p ~{output_dir}
        fastqc --outdir ~{output_dir} ~{reads[0]} ~{reads[1]}
    }

    runtime {
        docker: "biocontainers/fastqc:v0.11.9_cv8"
    }

    output {
        File results = output_dir
    }
}

task tigger_genotype {
    input {
        File airr_tsv #/data/changeo/sample/sample_db-pass.tab
        String name
        Int minSeq = 25 #50 #5 #x
        Int minGerm = 100 #200 #20 #y
        String outdir = "tigger_genotype"
    }

    command {
        mkdir tigger_genotype
        tigger-genotype -d ~{airr_tsv} -n ~{name} -o ~{outdir} -x ~{minSeq} -y ~{minGerm}
    }

    runtime {
        docker: "immcantation/suite:devel"
    }

    output {
        File out = outdir
        File genotype_fasta = outdir + "/" + name + "_genotype.fasta"
        File genotype = outdir + "/" + name + "_genotype.pdf"
        File genotype_tsv = outdir + "/" + name + "_genotyped.tsv"
    }
}
