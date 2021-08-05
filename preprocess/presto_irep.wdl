version development

import "https://raw.githubusercontent.com/antonkulaga/bioworkflows/main/common/files.wdl" as files

workflow presto_irep{
    input {
        Array[File] reads
        String destination
        String name
        String coord = "illumina" #sra
        Int threads = 4
        Boolean flip_reads = false
        Int min_reads_per_ig = 2
        File primer_constant
        File primer_variable

        Boolean rc = true
        Float maxerror = 0.2
        Float maxgap = 0.5

    }

    call presto{
        input: name = name, output_dir = "presto",
        reads = if(flip_reads) then [reads[0], reads[1]] else reads, NPROC = threads, primer_constant = primer_constant, primer_variable = primer_variable
    }

    call files.copy{ input: destination = destination, files = [presto.out] }

}

task presto {
    input {
        String output_dir = "results"
        String name
        Array[File] reads
        String coordinates = "illumina"
        File primer_constant
        File primer_variable

        Int NPROC = 4
        Boolean rc = true
        Boolean revpr = true
        Int collapse_max_missing = 20
        Int dupcount = 2
        Int min_quality = 20
    }

    Array[String] basenames = [basename(reads[0]),basename(reads[1])]

    String primer_constant_name = basename(primer_constant)
    String primer_variable_name = basename(primer_variable)

    command {

        # Make output directory and empty log files
        mkdir -p ~{output_dir}; cd ~{output_dir}

        # Start
        echo "OUTPUT DIRECTORY: ~{output_dir}"
        echo -e "START"
        STEP=0

        ln -s ~{reads[0]} ~{basenames[0]}
        ln -s ~{reads[1]} ~{basenames[1]}
        ln -s ~{primer_constant} ~{primer_constant_name}
        ln -s ~{primer_variable} ~{primer_variable_name}

        # Assemble paired ends via mate-pair alignment
        printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "AssemblePairs align"
        AssemblePairs.py align -1 ~{basenames[1]} -2 ~{basenames[0]} --coord ~{coordinates} ~{if(rc) then "--rc tail" else ""} --outname ~{name}  --outdir . --log AP.log --nproc ~{NPROC}

        # Remove low quality reads
        printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "FilterSeq quality"
        FilterSeq.py quality -s "~{name}_assemble-pass.fastq" -q ~{min_quality} --outname ~{name} --log FS.log --nproc ~{NPROC}

        # Identify forward and reverse primers
        printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "MaskPrimers score"

        MaskPrimers.py score -s "~{name}_quality-pass.fastq" -p ~{primer_variable_name} \
        --mode mask --pf VPRIMER \
        --outname "~{name}-FWD" --log MPV.log --nproc ~{NPROC}

        MaskPrimers.py score -s "~{name}-FWD_primers-pass.fastq" -p ~{primer_constant_name} \
        --mode mask ~{if(revpr) then "--revpr" else ""} --pf CPRIMER \
        --outname "~{name}-REV" --log MPC.log --nproc ~{NPROC}

        # Remove duplicate sequences
        printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "CollapseSeq"
        CollapseSeq.py -s "~{name}-REV_primers-pass.fastq" -n ~{collapse_max_missing} \
        --uf CPRIMER --cf VPRIMER --act set --inner \
        --outname ~{name}

        # Subset to sequences observed at least twice
        printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "SplitSeq group"
        SplitSeq.py group -s "~{name}_collapse-unique.fastq" -f DUPCOUNT --num ~{dupcount} \
        --outname ~{name}

        # Create annotation table of final unique sequences
        printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "ParseHeaders table"
        ParseHeaders.py table -s "~{name}_atleast-~{dupcount}.fastq" \
        -f ID DUPCOUNT CPRIMER VPRIMER

        # Process log files
        printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "ParseLog"
        ParseLog.py -l AP.log -f ID LENGTH OVERLAP ERROR PVALUE > /dev/null &
        ParseLog.py -l FS.log -f ID QUALITY > /dev/null &
        ParseLog.py -l MP[VC].log -f ID PRIMER ERROR > /dev/null &

        rm ~{basenames[0]} ~{basenames[1]} ~{primer_constant_name} ~{primer_variable_name}
        printf "DONE\n\n"
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