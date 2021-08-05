version development

import "https://raw.githubusercontent.com/antonkulaga/bioworkflows/main/common/files.wdl" as files

workflow presto_irep{
    input {
        Array[File] reads
        String destination
        String name
        String coord = "illumina" #sra
        Int threads = 12
        Int min_reads_per_ig = 2
        Int min_length = 196
        Int min_quality = 20
        #File constant
        #File variable

        #Float maxerror = 0.2
        #Float maxgap = 0.5

    }

    call presto{
        input: name = name, output_dir = "presto",
        reads = reads, NPROC = threads, min_quality  = min_quality, min_length = min_length
    }

    call files.copy{ input: destination = destination, files = [presto.out] }

}

task presto {
    input {
        String output_dir = "results"
        String name
        Array[File] reads
        String coordinates = "illumina"
        #File constant
        #File variable

        Int NPROC = 12
        Int collapse_max_missing = 20

        Int dupcount = 2
        Int min_quality = 20

        Int const_length = 27
        Int variable_length = 21
        Int min_length
    }

    Array[String] basenames = [basename(reads[0], ".fastq"),basename(reads[1],".fastq")]

    #String constant_name = basename(constant)
    #String variable_name = basename(variable)

    command {

        # Make output directory and empty log files
        mkdir -p ~{output_dir}; cd ~{output_dir}

        # Start
        echo "OUTPUT DIRECTORY: ~{output_dir}"
        echo -e "START"
        STEP=0

        ln -s ~{reads[0]} ~{basenames[0]}.fastq
        ln -s ~{reads[1]} ~{basenames[1]}.fastq

        # Mask primers

        MaskPrimers.py extract -s ~{basenames[0]}.fastq --len ~{const_length} --mode mask --pf CPRIMER
        MaskPrimers.py extract -s ~{basenames[0]}_primers-pass.fastq --revpr --len ~{variable_length} --mode mask --pf VPRIMER
        PairSeq.py -2 ~{basenames[0]}_primers-pass_primers-pass.fastq -1 ~{basenames[1]}.fastq --2f CPRIMER VPRIMER

        # Assemble paired ends via mate-pair alignment
        printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "AssemblePairs align"

        AssemblePairs.py align -1 ~{basenames[1]}_pair-pass.fastq -2 ~{basenames[0]}_primers-pass_primers-pass_pair-pass.fastq \
          --coord ~{coordinates} --rc tail --2f CPRIMER VPRIMER  \
          --outname ~{name}  --outdir . --log AP.log --nproc ~{NPROC}

        # Remove low quality reads
        printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "FilterSeq quality"
        FilterSeq.py quality -s "~{name}_assemble-pass.fastq" -q ~{min_quality} --outname ~{name} --log FS_quality.log --nproc ~{NPROC}
        FilterSeq.py length -s "~{name}_quality-pass.fastq" -n ~{min_length} --outname ~{name} --log FS_length.log --nproc ~{NPROC}

        # Remove duplicate sequences
        printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "CollapseSeq"
        CollapseSeq.py -s "~{name}_length-pass.fastq" -n ~{collapse_max_missing} --uf CPRIMER --cf VPRIMER --act set --inner --outname ~{name}
        SplitSeq.py group -s ~{name}_collapse-unique.fastq -f DUPCOUNT --num ~{dupcount} --outname ~{name}
        rm -rf ~{basenames[0]}.fastq ~{basenames[1]}.fastq

    }

    runtime {
        docker: "immcantation/suite:devel"
    }

    output {
        File out = output_dir
        #File ap_lof = output_dir + "/" + "AP.log"
        #File ap_table = output_dir + "/" + "AP_table.tab"
        #File fs_log = output_dir + "/" + "FS.log"
        #File fs_table = output_dir + "/" + "FS_table.tab"
        #File mpv_log = output_dir + "/" + "MPV.log"
        #File mpv_table = output_dir + "/" + "MPV_table.tab"
        #File assembly_pass =  output_dir + "/" + name + "_assemble-pass.fastq"
        #File quality_pass =  output_dir + "/" + name + "_quality-pass.fastq"
    }
}