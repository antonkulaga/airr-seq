version development


import "https://raw.githubusercontent.com/antonkulaga/bioworkflows/main/download/download_runs.wdl" as downloader

import "https://raw.githubusercontent.com/antonkulaga/bioworkflows/main/common/files.wdl" as files

import "https://raw.githubusercontent.com/antonkulaga/airr-seq/main/analysis/changeo_igblast/changeo_igblast.wdl" as changeo
import "https://raw.githubusercontent.com/antonkulaga/airr-seq/main/analysis/clonal_analysis/clonal_analysis.wdl" as clonal


workflow cowan_airr {
    input {
        String title = ""
        Array[String] runs
        String destination
        String key = "0a1d74f32382b8a154acacc3a024bdce3709"

        #File V_primers #/data/samples/AIRR-Seq/RA/PRJNA561156/primers/C_primers.fasta
        #File C_primers #/data/samples/AIRR-Seq/RA/PRJNA561156/primers/V_primers.fasta

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
        #Int min_reads_per_ig = 2
        Int min_length = 196
        Int min_quality = 20
        #Int start_v = 6
        #Int start_c = 18
        #Int barcode_length = 6
        #Int constant_length = 12
        Int min_dupcount = 2
        Float clones_bin_width = 0.02

    }

    call downloader.download_runs as download_runs{
        input:
            title = title,
            runs = runs,
            experiment_folder = destination,
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
        call presto_no_primers as presto{
            input: name = name, output_dir = "presto",
                    reads = run.original_reads, NPROC = threads,
                    min_quality = min_quality, min_length = min_length,
                    dupcount = min_dupcount, coordinates = coordinates, max_memory = max_memory_gb
        }
        #call presto{
        #    input: name = name, output_dir = "presto",
        #        reads = run.original_reads, NPROC = threads,
        #        min_quality = min_quality, min_length = min_length,
        #        start_v = start_v, start_c = start_c, dupcount = min_dupcount,
        #        coordinates = coordinates, max_memory = max_memory_gb, C_primers = C_primers, V_primers = V_primers,
        #        barcode_length = barcode_length,
        #        constant_length = constant_length
        #}

        String run_folder = destination + "/" + name

        call files.copy as copy_presto {
            input: destination = run_folder, files = [presto.results]
        }

        call changeo.changeo_igblast as igblast {
            input: fastq = presto.out, threads = threads, name = name, destination = destination
        }

        call clonal.clonal_analysis as clonal_analysis {
            input:
                airr_tsv = igblast.airr_tsv,
                destination = run_folder,
                name = name,
                threads = threads,
                clones_bin_width = clones_bin_width,
                max_memory_gb = max_memory_gb
        }

    }
}

task presto_no_primers{
    input {
        String output_dir = "presto"
        String name
        Array[File] reads
        String coordinates = "illumina"
        Int NPROC = 12


        Int dupcount = 2
        Int min_quality = 20
        Int collapse_max_missing = 20

        Int min_length

        Int max_memory
    }

    Array[String] basenames = [basename(reads[0], ".fastq"),basename(reads[1],".fastq")]


    command {
        # Make output directory and empty log files
        mkdir -p ~{output_dir}; cd ~{output_dir}

        # Start
        echo "OUTPUT DIRECTORY: ~{output_dir}"
        echo -e "START"
        STEP=0

        ln -s ~{reads[0]} ~{basenames[0]}.fastq
        ln -s ~{reads[1]} ~{basenames[1]}.fastq

        # Assemble paired ends via mate-pair alignment
        printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "AssemblePairs align"
        AssemblePairs.py align -1 ~{basenames[1]}.fastq -2 ~{basenames[0]}.fastq --coord ~{coordinates} --rc tail --outname ~{name}  --outdir . --log AP.log --nproc ~{NPROC}

        # Remove low quality reads
        printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "FilterSeq quality"
        FilterSeq.py quality -s "~{name}_assemble-pass.fastq" -q ~{min_quality} --outname ~{name} --log FS_quality.log --nproc ~{NPROC}
        printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "FilterSeq length"
        FilterSeq.py length -s "~{name}_quality-pass.fastq" -n ~{min_length} --outname ~{name} --log FS_length.log --nproc ~{NPROC}

        ParseLog.py -l AP.log -f ID LENGTH OVERLAP ERROR PVALUE FIELDS1 FIELDS2
        CollapseSeq.py -s "~{name}_length-pass.fastq" -n ~{collapse_max_missing} --inner --outname ~{name}
        SplitSeq.py group -s ~{name}_collapse-unique.fastq -f DUPCOUNT --num ~{dupcount} --outname ~{name}
        ParseHeaders.py table -s "~{name}_atleast-~{dupcount}.fastq" -f ID DUPCOUNT
        rm -rf ~{basenames[0]}.fastq ~{basenames[1]}.fastq
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
        File out = output_dir + "/" + name + "_atleast-" + dupcount + ".fastq"
        File headers = output_dir + "/" + name + "_atleast-" + dupcount + "_headers.tab"
    }
}

task presto_primers { #todo, rewrite
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

        Int constant_length
        #Int variable_length = 21
        Int barcode_length = 6
        Int min_length

        Int start_v = 0
        Int start_c = 0
        Int max_memory
    }

    Array[String] basenames = [basename(reads[0], ".fastq"),basename(reads[1],".fastq")]

    command {
        # Make output directory and empty log files
        mkdir -p ~{output_dir}; cd ~{output_dir}

        # Start
        echo "OUTPUT DIRECTORY: ~{output_dir}"
        echo -e "START"
        STEP=0

        ln -s ~{reads[0]} ~{basenames[0]}.fastq
        ln -s ~{reads[1]} ~{basenames[1]}.fastq

        # Assemble paired ends via mate-pair alignment
        printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "AssemblePairs align"
        AssemblePairs.py align -1 ~{basenames[1]}_consensus-pass.fastq -2 ~{basenames[0]}_consensus-pass.fastq --coord sra --rc tail --outname ~{name}  --outdir . --log AP.log --nproc ~{NPROC}

        # Remove low quality reads
        printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "FilterSeq quality"
        FilterSeq.py quality -s "~{name}_assemble-pass.fastq" -q 20 --outname ~{name} --log FS.log --nproc ~{NPROC}
        ParseHeaders.py collapse -s ~{name}_assemble-pass.fastq -f CONSCOUNT --act min
        CollapseSeq.py -s ~{name}*reheader.fastq -n 20 --inner --uf PRCONS \
        --cf CONSCOUNT --act sum --outname ~{name}
        SplitSeq.py group -s ~{name}_collapse-unique.fastq -f CONSCOUNT --num 2 --outname ~{name}
        ParseHeaders.py table -s ~{name}_atleast-2.fastq -f ID PRCONS CONSCOUNT DUPCOUNT
        ParseLog.py -l FS1.log FS2.log -f ID QUALITY
        ParseLog.py -l MP1.log MP2.log -f ID PRIMER BARCODE ERROR
        ParseLog.py -l BC1.log BC2.log -f BARCODE SEQCOUNT CONSCOUNT PRIMER PRCONS PRCOUNT \
        PRFREQ ERROR
        ParseLog.py -l AP.log -f ID LENGTH OVERLAP ERROR PVALUE FIELDS1 FIELDS2
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

