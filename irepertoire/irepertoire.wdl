## # Immune repertoire analysis 
## This workflow performs end-to-end immune repertoire analysis from MiSeq 2x250 bp runs.

version development

# Imports tasks for files: copy, merge
import "https://raw.githubusercontent.com/antonkulaga/bioworkflows/main/common/files.wdl" as files

#production version
import "https://raw.githubusercontent.com/antonkulaga/airr-seq/main/analysis/changeo_igblast/changeo_igblast.wdl" as changeo
import "https://raw.githubusercontent.com/antonkulaga/airr-seq/main/analysis/clonal_analysis/clonal_analysis.wdl" as clonal


workflow irepertoire {

    input {
        Array[File] reads

        String destination
        String name

        # String species = "human"
        File? IGHC_fasta
        String coordinates = "illumina" #sra
        # Int min_reads_per_ig = 2
        Int min_length = 196
        Int min_quality = 20
        Int start_v = 10
        Int min_dupcount = 2
        Float clones_bin_width = 0.02
        Int threads = 12
        Int max_memory_gb = 96

        #"cowan_airr.C_primers": "/data/samples/AIRR-Seq/RA/PRJNA561156/primers/C_primers.fasta",
        #"cowan_airr.V_primers": "/data/samples/AIRR-Seq/RA/PRJNA561156/primers/V_primers.fasta",
    }

    meta {
        description: "End-to-end workflow for immune repertoire analysis"
    }
    
    parameter_meta {
        reads: "Paths 2-array to the paired fastq files"
        destination: "Path to directory to store inputs, intermediary and output files"
        name: "Analysis name, used throughout the workflow for file naming"
        # species: "Not used"
        coordinates: "See `presto` task"
        # min_reads_per_ig: "Not used"
        min_length: "See `presto` task"
        min_quality: "See `presto` task"
        start_v: "See `presto` task"
        min_dupcount: "See `presto` dupcount"
        clones_bin_width: "See `clonal.clonal_analysis` task"
        threads: "Constraints resource consumption, can crash containers if insufficient resources."
        max_memory_gb: "Constraints resource consumption, can crash containers if insufficient resources."
    }



    call fastqc {
        input: 
            output_dir = "fastqc_results", 
            reads = reads
    }

    call files.copy as copy_fastqc {
        input: 
            destination = destination, 
            files = [fastqc.results]
    }

    call presto as presto {
        input: 
            name = name, 
            output_dir = "presto",
            reads = reads,
            IGHC_fasta = IGHC_fasta,
            NPROC = threads, 
            min_quality = min_quality, 
            min_length = min_length, 
            start_v = start_v, 
            dupcount = min_dupcount,
            coordinates = coordinates, 
            max_memory = max_memory_gb
    }

    call files.copy as copy_presto { 
        input: 
            destination = destination, 
            files = [presto.results]
    }
    
    call changeo.changeo_igblast as igblast {
        input: 
            fastq = presto.out, 
            threads = threads, 
            name = name, 
            destination = destination
    }

    call clonal.clonal_analysis as clonal_analysis {
            input:
                # airr_tsv = igblast.out,
                airr_tsv = igblast.airr_tsv_translated_functional,
                destination = destination,
                name = name,
                threads = threads,
                clones_bin_width = clones_bin_width,
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

    meta {
        description: "Perform quality control"
    }

    parameter_meta {
        output_dir: "Path to results subdirectory with fastqc results"
        reads: "Paths to the paired fastqc files"
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

task presto {
    input {
        String output_dir = "results"
        String name
        Array[File] reads
        String coordinates = "illumina"

        File? IGHC_fasta

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

        # Mask primers and infer constant region
        MaskPrimers.py extract -s ~{basenames[0]}.fastq --len ~{const_length} --mode mask --pf CPRIMER
        MaskPrimers.py extract -s ~{basenames[1]}.fastq --start ~{start_v} ~{if(start_v>0) then "--barcode --bf V_barcode" else ""} --len ~{variable_length} --mode cut --pf VPRIMER
        PairSeq.py -1 ~{basenames[1]}_primers-pass.fastq -2 ~{basenames[0]}_primers-pass.fastq  --1f VPRIMER ~{if(start_v>0) then "V_barcode" else ""} --2f CPRIMER

        # Assemble paired ends via mate-pair alignment
        printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "AssemblePairs align"

        AssemblePairs.py align -1 ~{basenames[1]}_primers-pass_pair-pass.fastq -2 ~{basenames[0]}_primers-pass_pair-pass.fastq \
          --coord ~{coordinates} --rc tail --2f CPRIMER VPRIMER ~{if(start_v>0) then "V_barcode" else ""}  \
          --outname ~{name}  --outdir . --log AP.log --nproc ~{NPROC}

        ParseLog.py -l AP.log -f ID LENGTH OVERLAP ERROR PVALUE

        # Remove low quality reads
        printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "FilterSeq quality"
        FilterSeq.py quality -s "~{name}_assemble-pass.fastq" -q ~{min_quality} --outname ~{name} --log FS_quality.log --nproc ~{NPROC}
        FilterSeq.py length -s "~{name}_quality-pass.fastq" -n ~{min_length} --outname ~{name} --log FS_length.log --nproc ~{NPROC}

        ParseLog.py -l FS_quality.log -f ID QUALITY

        # Call constant region (for IGH) and furter processing
        if ~{defined(IGHC_fasta)}
        then

            # C_CALL for IGH
            printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "C_CALL"
            MaskPrimers.py align -s "~{name}_length-pass.fastq" -p ~{IGHC_fasta} --revpr --skiprc --maxerror 0.3 --mode tag --outname "~{name}_c_call" --pf C_CALL
            mv "~{name}_c_call_primers-pass.fastq" "~{name}_c_called.fastq"
            
             # Remove duplicate sequences
            printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "CollapseSeq"
            CollapseSeq.py -s "~{name}_c_called.fastq" -n ~{collapse_max_missing} --uf CPRIMER C_CALL --cf VPRIMER --act set --inner --outname ~{name}
            SplitSeq.py group -s ~{name}_collapse-unique.fastq -f DUPCOUNT --num ~{dupcount} --outname ~{name}
            ParseHeaders.py table -s "~{name}_atleast-~{dupcount}.fastq" -f ID DUPCOUNT CPRIMER VPRIMER C_CALL ~{if(start_v>0) then "V_barcode" else ""}
        
        else
            # Remove duplicate sequences
            printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "CollapseSeq"
            CollapseSeq.py -s "~{name}_length-pass.fastq" -n ~{collapse_max_missing} --uf CPRIMER --cf VPRIMER --act set --inner --outname ~{name}
            SplitSeq.py group -s ~{name}_collapse-unique.fastq -f DUPCOUNT --num ~{dupcount} --outname ~{name}
            ParseHeaders.py table -s "~{name}_atleast-~{dupcount}.fastq" -f ID DUPCOUNT CPRIMER VPRIMER ~{if(start_v>0) then "V_barcode" else ""}
        fi

        rm -rf ~{basenames[0]}.fastq ~{basenames[1]}.fastq
    }

    runtime {
        # docker: "immcantation/suite:devel"
        docker: "immcantation/suite:4.2.0"
        docker_memory: "~{max_memory}G"
        docker_cpu: "~{NPROC+1}"
    }

    output {
        File results = output_dir
        File AP_table = output_dir + "/" + "AP_table.tab"
        File out = output_dir + "/" + name + "_atleast-" + dupcount + ".fastq"
        File under_dupcount = output_dir + "/" + name + "_under-" + dupcount + ".fastq"
        File headers = output_dir + "/" + name + "_atleast-" + dupcount + "_headers.tab"
        # File primer_pass_1 = output_dir + "/" + name + ".1_primers-pass.fastq"
        # File primer_pass_2 = output_dir + "/" + name + ".2_primers-pass.fastq"
        # File pair_pass_1 = output_dir + "/" + name + ".1_primers-pass_pair-pass.fastq"
        # File pair_pass_2 = output_dir + "/" + name + ".2_primers-pass_pair-pass.fastq"
    }
}