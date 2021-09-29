version development


import "https://raw.githubusercontent.com/antonkulaga/bioworkflows/main/download/download_runs.wdl" as downloader

import "https://raw.githubusercontent.com/antonkulaga/bioworkflows/main/common/files.wdl" as files

import "https://raw.githubusercontent.com/antonkulaga/airr-seq/main/analysis/changeo_igblast/changeo_igblast.wdl" as changeo
import "https://raw.githubusercontent.com/antonkulaga/airr-seq/main/analysis/clonal_analysis/clonal_analysis.wdl" as clonal


workflow cowan_airr {
    input {
        String title = ""
        Array[String] runs
        String experiment_folder
        String key = "0a1d74f32382b8a154acacc3a024bdce3709"

        File V_primers #/data/samples/AIRR-Seq/RA/PRJNA561156/primers/C_primers.fasta
        File C_primers #/data/samples/AIRR-Seq/RA/PRJNA561156/primers/V_primers.fasta

        Int extract_threads = 12
        Int max_memory_gb = 42
        Int align_threads = 12
        #Int sort_threads = 12
        Boolean copy_extracted = true
        Boolean copy_cleaned = true
        Boolean aspera_download = true
        Boolean skip_technical = true
        Boolean original_names = false
        String coordinates = "sra"
        Int threads = 12
        Int min_reads_per_ig = 2
        Int min_length = 196
        Int min_quality = 20
        Int start_v = 6
        Int min_dupcount = 2
        Float clones_bin_width = 0.02

    }

    call downloader.download_runs as download_runs{
        input:
            title = title,
            runs = runs,
            experiment_folder = experiment_folder,
            key = key,
            extract_threads = extract_threads,
            copy_cleaned = copy_cleaned,
            aspera_download = aspera_download,
            skip_technical = skip_technical,
            original_names = original_names,
            copy_extracted = copy_extracted,
    }
    Array[CleanedRun] cleaned_runs =  download_runs.out

    scatter(run in cleaned_runs) {
        String name = run.run
        call presto{
            input: name = name, output_dir = "presto",
                reads = run.original_reads, NPROC = threads, min_quality = min_quality, min_length = min_length, start_v = start_v, dupcount = min_dupcount,
                coordinates = coordinates, max_memory = max_memory_gb, C_primers = C_primers, V_primers = V_primers
        }

    }
}

task presto {
    input {
        String output_dir = "presto"
        String name
        Array[File] reads
        String coordinates = "illumina"
        File C_primers
        File V_primers

        Int NPROC = 12
        Int collapse_max_missing = 20

        Int dupcount = 2
        Int min_quality = 20

        Int const_length = 27
        Int variable_length = 21
        Int min_length

        Int start_v = 0
        Int max_memory
    }

    Array[String] basenames = [basename(reads[0], ".fastq.gz"),basename(reads[1],".fastq.gz")]

    command {
        # Make output directory and empty log files
        mkdir -p ~{output_dir}; cd ~{output_dir}

        # Start
        echo "OUTPUT DIRECTORY: ~{output_dir}"
        echo -e "START"
        STEP=0

        ln -s ~{reads[0]} ~{basenames[0]}.fastq.gz
        ln -s ~{reads[1]} ~{basenames[1]}.fastq.gz

        gunzip -f ~{basenames[0]}.fastq.gz ~{basenames[1]}.fastq.gz

        # Assemble paired ends via mate-pair alignment
        printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "AssemblePairs align"
        AssemblePairs.py align -1 ~{basenames[1]}.fastq -2 ~{basenames[0]}.fastq --coord sra --rc tail --outname ~{name}  --outdir . --log AP.log --nproc ~{NPROC}

        # Remove low quality reads
        printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "FilterSeq quality"
        FilterSeq.py quality -s "~{name}_assemble-pass.fastq" -q 20 --outname ~{name} --log FS.log --nproc ~{NPROC}

        # Identify forward and reverse primers
        printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "MaskPrimers score"
        MaskPrimers.py score -s "~{name}_quality-pass.fastq" -p ~{C_primers} \
        --mode mask --start ~{start_v} --pf VPRIMER \
        --outname "~{name}-FWD" --log MPV.log --nproc ~{NPROC}

        MaskPrimers.py score -s "~{name}-FWD_primers-pass.fastq" -p ~{V_primers} \
        --mode cut --start 7 --revpr --pf CPRIMER \
        --outname "~{name}-REV" --log MPC.log --nproc ~{NPROC}
    }

    runtime {
        docker: "immcantation/suite:devel"
        docker_memory: "~{max_memory}G"
        docker_cpu: "~{NPROC+1}"
    }

    output {
        File results = output_dir
        File AP_log = output_dir + "/" + "AP.log"
        File assembly = output_dir + "/" + name+"_assemble-pass.fastq"
        #File AP_table = output_dir + "/" + "AP_table.tab"
        #File out = output_dir + "/" + name + "_atleast-" + dupcount + ".fastq"
        #File headers = output_dir + "/" + name + "_atleast-" + dupcount + "_headers.tab"
    }
}

