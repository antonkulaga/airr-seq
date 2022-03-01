version development

import "https://raw.githubusercontent.com/antonkulaga/bioworkflows/main/common/files.wdl" as files

workflow presto {
    input {
        Array[File] reads
        String destination
        String name
        Array[File] primers
        Int threads = 4
    }

    call presto{
        input:
            R1_FILE = reads[0],
            R2_FILE = reads[1],
            R1_PRIMERS = primers[0],
            R2_PRIMERS = primers[1],
            NPROC = threads,
            name = name

    }

    call files.copy as copy {
        input: destination = destination + "/" + name, files = [presto.out]
    }


}

task presto {
    input {
        String output_dir = "results"
        String name = "M4"

        File R1_FILE #ERR346600_1.fastq
        File R2_FILE #ERR346600_2.fastq
        File R1_PRIMERS #Greiff2014_CPrimers.fasta
        File R2_PRIMERS #Greiff2014_VPrimers.fasta
        Int NPROC = 4
    }


    command {

        # Make output directory and empty log files
        mkdir -p ~{output_dir}; cd ~{output_dir}

        # Start
        echo "OUTPUT DIRECTORY: ~{output_dir}"
        echo -e "START"
        STEP=0

        # Assemble paired ends via mate-pair alignment
        printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "AssemblePairs align"
        AssemblePairs.py align -1 ~{R2_FILE} -2 ~{R1_FILE} --coord sra --rc tail --outname ~{name}  --outdir . --log AP.log --nproc ~{NPROC}

        # Remove low quality reads
        printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "FilterSeq quality"
        FilterSeq.py quality -s "~{name}_assemble-pass.fastq" -q 20 --outname ~{name} --log FS.log --nproc ~{NPROC}

        # Identify forward and reverse primers
        printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "MaskPrimers score"
        MaskPrimers.py score -s "~{name}_quality-pass.fastq" -p ~{R2_PRIMERS} \
        --mode mask --start 4 --pf VPRIMER \
        --outname "~{name}-FWD" --log MPV.log --nproc ~{NPROC}
        MaskPrimers.py score -s "~{name}-FWD_primers-pass.fastq" -p ~{R1_PRIMERS} \
        --mode cut --start 4 --revpr --pf CPRIMER \
        --outname "~{name}-REV" --log MPC.log --nproc ~{NPROC}

        # Remove duplicate sequences
        printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "CollapseSeq"
        CollapseSeq.py -s "~{name}-REV_primers-pass.fastq" -n 20 \
        --uf CPRIMER --cf VPRIMER --act set --inner \
        --outname ~{name}

        # Subset to sequences observed at least twice
        printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "SplitSeq group"
        SplitSeq.py group -s "~{name}_collapse-unique.fastq" -f DUPCOUNT --num 2 \
        --outname ~{name}

        # Create annotation table of final unique sequences
        printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "ParseHeaders table"
        ParseHeaders.py table -s "~{name}_atleast-2.fastq" -f ID DUPCOUNT CPRIMER VPRIMER

        # Process log files
        printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "ParseLog"
        ParseLog.py -l AP.log -f ID LENGTH OVERLAP ERROR PVALUE > /dev/null &
        ParseLog.py -l FS.log -f ID QUALITY > /dev/null &
        ParseLog.py -l MP[VC].log -f ID PRIMER ERROR > /dev/null &
        printf "DONE\n\n"
    }

    runtime {
        docker: "immcantation/suite:devel"
    }

    output {
        File out = output_dir
        File ap_lof = output_dir + "/" + "AP.log"
        File ap_table = output_dir + "/" + "AP_table.tab"
        File fs_log = output_dir + "/" + "FS.log"
        File fs_table = output_dir + "/" + "FS_table.tab"
        File mpv_log = output_dir + "/" + "MPV.log"
        File mpv_table = output_dir + "/" + "MPV_table.tab"
        File assembly_pass =  output_dir + "/" + name + "_assemble-pass.fastq"
        File quality_pass =  output_dir + "/" + name + "_quality-pass.fastq"
    }
}