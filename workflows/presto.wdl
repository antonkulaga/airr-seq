version development

workflow quant_index {
    input {

        Array[File] reads = ["/data/sources/airr-seq/workflows/inputs/Greiff2014/ERR346600_1.fastq", "/data/sources/airr-seq/workflows/inputs/Greiff2014/ERR346600_2.fastq"]
        Array[File] primers = ["/data/sources/airr-seq/workflows/inputs/Greiff2014/Greiff2014_CPrimers.fasta", "/data/sources/airr-seq/workflows/inputs/Greiff2014/Greiff2014_VPrimers.fasta"]
        #String destination

    }

    call presto{
        input: OUTDIR = "test_input", OUTNAME = "test_output", R1_FILE = reads[0], R2_FILE = reads[1], R1_PRIMERS = primers[0], R2_PRIMERS = primers[1], NPROC = 4
    }

    #call download { input: sra = run, aspera_download = aspera_download }
    #call extract {input: sra = download.out, is_paired = is_paired, threads = extract_threads}


}


task presto {
    input {
        String OUTDIR = "output"
        String OUTNAME = "M4"

        File R1_FILE #ERR346600_1.fastq
        File R2_FILE #ERR346600_2.fastq
        File R1_PRIMERS #Greiff2014_CPrimers.fasta
        File R2_PRIMERS #Greiff2014_VPrimers.fasta
        Int NPROC = 4
    }


    command {

        # Make output directory and empty log files
        mkdir -p ~{OUTDIR}; cd ~{OUTDIR}

        # Start
        echo "OUTPUT DIRECTORY: ~{OUTDIR}"
        echo -e "START"
        STEP=0

        # Assemble paired ends via mate-pair alignment
        printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "AssemblePairs align"
        AssemblePairs.py align -1 ~{R1_FILE} -2 ~{R2_FILE} --coord sra --rc tail --outname ~{OUTNAME} --outdir . --log AP.log --nproc ~{NPROC}

        # Remove low quality reads
        printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "FilterSeq quality"
        FilterSeq.py quality -s "~{OUTNAME}_assemble-pass.fastq" -q 20 --outname ~{OUTNAME} --log FS.log --nproc ~{NPROC}

        # Identify forward and reverse primers
        printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "MaskPrimers score"
        MaskPrimers.py score -s "~{OUTNAME}_quality-pass.fastq" -p ~{R2_PRIMERS} \
        --mode mask --start 4 --pf VPRIMER \
        --outname "~{OUTNAME}-FWD" --log MPV.log --nproc ~{NPROC}
        MaskPrimers.py score -s "~{OUTNAME}-FWD_primers-pass.fastq" -p ~{R1_PRIMERS} \
        --mode cut --start 4 --revpr --pf CPRIMER \
        --outname "~{OUTNAME}-REV" --log MPC.log --nproc ~{NPROC}

        # Remove duplicate sequences
        printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "CollapseSeq"
        CollapseSeq.py -s "~{OUTNAME}-REV_primers-pass.fastq" -n 20 \
        --uf CPRIMER --cf VPRIMER --act set --inner \
        --outname ~{OUTNAME}

        # Subset to sequences observed at least twice
        printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "SplitSeq group"
        SplitSeq.py group -s "~{OUTNAME}_collapse-unique.fastq" -f DUPCOUNT --num 2 \
        --outname ~{OUTNAME}

        # Create annotation table of final unique sequences
        printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "ParseHeaders table"
        ParseHeaders.py table -s "~{OUTNAME}_atleast-2.fastq" \
        -f ID DUPCOUNT CPRIMER VPRIMER

        # Process log files
        printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "ParseLog"
        ParseLog.py -l AP.log -f ID LENGTH OVERLAP ERROR PVALUE > /dev/null &
        ParseLog.py -l FS.log -f ID QUALITY > /dev/null &
        ParseLog.py -l MP[VC].log -f ID PRIMER ERROR > /dev/null &
    }

    runtime {
        docker: "immcantation/suite:devel"
    }

    output {
        File ap_log = "AP.log"
        File fs_log = "FS.log"
        File out = OUTDIR

    }
}



task copy {
    input {
        Array[File] files
        String destination
    }

    String where = sub(destination, ";", "_")

    command {
        mkdir -p ~{where}
        cp -L -R -u ~{sep=' ' files} ~{where}
        declare -a files=(~{sep=' ' files})
        for i in ~{"$"+"{files[@]}"};
        do
        value=$(basename ~{"$"}i)
        echo ~{where}/~{"$"}value
        done
    }

    output {
        Array[File] out = read_lines(stdout())
    }
}
