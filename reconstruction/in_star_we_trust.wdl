version development

import "https://raw.githubusercontent.com/antonkulaga/rna-seq/master/pipelines/rna-alignment/align_rna_run.wdl" as star
import "https://raw.githubusercontent.com/antonkulaga/airr-seq/main/reconstruction/trust.wdl" as trust_me

# workflow that uses star for rna-seq alignment and trust4 for reconstruction

workflow in_star_we_trust {
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

    call star.star_align_run as star_align{
        input: title = title, runs = runs, experiment_folder = experiment_folder,
            index_dir = star_index_dir, key = key, extract_threads = extract_threads,
            max_memory_gb = max_memory_gb, align_threads = align_threads, copy_extracted = copy_extracted,
            copy_cleaned = copy_cleaned, aspera_download = aspera_download, skip_technical = skip_technical,
            original_names = original_names,
            deep_folder_structure = deep_folder_structure
    }

    scatter(mapping in  star_align.mappings) {
        call trust_me.trust as trust {
            input:
                bam = mapping.left.sorted,
                threads = trust_threads,
                imgt_reference = imgt_reference,
                destination = mapping.right + "/" + "reconstruction"
        }
    }
}