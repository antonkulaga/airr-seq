version development

import "https://raw.githubusercontent.com/antonkulaga/bioworkflows/main/common/files.wdl" as files
import "https://raw.githubusercontent.com/antonkulaga/rna-seq/master/pipelines/rna-alignment/align_rna_run.wdl" as star
import "https://raw.githubusercontent.com/antonkulaga/airr-seq/main/reconstruction/trust.wdl" as trust_me

# workflow that uses star for rna-seq alignment and trust4 for reconstruction

workflow in_star_we_trust_sc {
    input {
        String title = ""
        File imgt_reference
        File star_index_dir
        Array[String] runs
        String experiment_folder
        String key = "0a1d74f32382b8a154acacc3a024bdce3709"
        Int trust_threads = 12
        Int extract_threads = 12
        Int max_memory_gb = 42
        Int align_threads = 12
        Boolean copy_extracted = true
        Boolean copy_cleaned = true
        Boolean aspera_download = true
        Boolean skip_technical = true
        Boolean original_names = false
        Boolean deep_folder_structure = true
    }


    scatter(run in runs) {
        call download_sc {
            input: sra = run, aspera_download = aspera_download
        }
        call star.star_align as star_align{
          input: run = run,
                reads = download_sc.out,
                index_dir = star_index_dir
        }

        AlignedRun aligned =  star_align.out
        String copy_to = experiment_folder + "/" + run  + "/" + "aligned"
        Pair[AlignedRun, String] mapping = (aligned, copy_to)

        call files.copy as copy_star {
            input:
                files = [aligned.sorted, aligned.to_transcriptome, aligned.summary, aligned.log, aligned.progress, aligned.reads_per_gene, aligned.junctions],
                destination =  copy_to,
        }

        call trust_me.trust as trust {
            input:
                bam = mapping.left.sorted,
                threads = trust_threads,
                imgt_reference = imgt_reference,
                destination = mapping.right + "/" + "reconstruction"
        }
    }

}

task download_sc {
    input {
        String sra
        Boolean aspera_download
    }
    #prefetch --ascp-path "/root/.aspera/connect/bin/ascp|/root/.aspera/connect/etc/asperaweb_id_dsa.openssh" --force yes -O results ~{sra}
    command {
        ~{if(aspera_download) then "download_sra_aspera.sh " else "prefetch -X 9999999999999 --force yes -O results -t http "} ~{sra}
    }

    #https://github.com/antonkulaga/biocontainers/tree/master/downloaders/sra

    runtime {
        docker: "quay.io/comp-bio-aging/download_sra:latest"
        maxRetries: 1
    }

    output {
        Array[File] out = glob(sra + "/*")
    }
}
