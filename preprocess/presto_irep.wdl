version development

import "https://raw.githubusercontent.com/antonkulaga/bioworkflows/main/common/files.wdl" as files

workflow presto_irep{
    input {
        Array[File] reads
        String destination
        String name
        String coord = "illumina" #sra
        Int umi_start = 3
        Int umi_length = 6
        Int threads = 4
        Float maxerror = 0.2
        Float maxgap = 0.5
        Boolean flip_reads = false
    }
    call extract_umi {
        input:
            reads = reads, start = 3, len = 6, coord = coord, name = name, maxerror = maxerror, maxgap = maxgap
    }

    call assemble_pairs {
        input:
            reads = if(flip_reads) then [extract_umi.consensus_reads[1], extract_umi.consensus_reads[0]] else extract_umi.consensus_reads, coord = "presto",
            name = name,
            NPROC = threads

    }

    #call files.copy as copy {
    #    input: destination = destination + "/" + name, files = [presto.out]
    #}

}

task extract_umi {
    input {
        Int start
        Int len
        String name
        Array[File] reads
        String coord
        String output_dir = "results"
        Float maxerror = 0.2
        Float maxgap = 0.5

                                    }
    # MaskPrimers.py extract -s in.fastq --start 0 --len 15 --pf UMI -o out.fastq

    String read_1_name = basename(reads[0])
    String read_2_name = basename(reads[1])

    command {
        # Make output directory and empty log files
        mkdir -p ~{output_dir}; cd ~{output_dir}

        ln -s ~{reads[0]} ~{read_1_name}
        ln -s ~{reads[1]} ~{read_2_name}

        # Start
        echo "OUTPUT DIRECTORY: ~{output_dir}"
        echo -e "START"

        echo "FilterSeq"

        FilterSeq.py quality -s ~{read_1_name} -q 20 --outname ~{name}_R1_umi --log FS1.log
        FilterSeq.py quality -s ~{read_2_name} -q 20 --outname ~{name}_R2 --log FS2.log

        echo "Extracting UMI"

        MaskPrimers.py extract -s ~{name}_R1_umi_quality-pass.fastq --start ~{start} --len ~{len} --pf BARCODE -o ~{name}_R1_quality-pass.fastq

        echo "Move UMI annotations to second read"

        PairSeq.py -1  ~{name}_R1_quality-pass.fastq -2  ~{name}_R2_quality-pass.fastq --1f BARCODE --coord ~{coord}

        echo "Build consensus"

        BuildConsensus.py -s ~{name}_R1_quality-pass_pair-pass.fastq --bf BARCODE --maxerror ~{maxerror} --maxgap ~{maxgap} --outname ~{name}_R1 --log BC1.log
        BuildConsensus.py -s ~{name}_R1_quality-pass_pair-pass.fastq --bf BARCODE --maxerror ~{maxerror} --maxgap ~{maxgap} --outname ~{name}_R2 --log BC2.log

        echo "pairing consensus sequences together and converting coordinates to presto"

        PairSeq.py -1 ~{name}_R1_consensus-pass.fastq -2 ~{name}_R2_consensus-pass.fastq --coord presto
    }

    output {
        File out = output_dir
        Array[File] consensus_reads = [output_dir + "/" + name + "_R1_consensus-pass_pair-pass.fastq", output_dir + "/" + name + "_R2_consensus-pass_pair-pass.fastq" ]
    }

    runtime {
        docker: "immcantation/suite:devel"
    }
}



task assemble_pairs {
    input {
        String output_dir = "results"
        String name

        Array[File] reads
        Int NPROC
        String coord = "presto"
    }

    String assembly_pass = name + "_assemble-pass.fastq"


    command {

        # Make output directory and empty log files
        mkdir -p ~{output_dir}; cd ~{output_dir}

        # Start
        echo "OUTPUT DIRECTORY: ~{output_dir}"
        echo -e "START"
        STEP=0

        # Assemble paired ends via mate-pair alignment
        printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "AssemblePairs align"
        AssemblePairs.py align -1 ~{reads[0]} -2 ~{reads[1]} --coord ~{coord } --rc tail --outname ~{name}  --outdir . --log AP.log --nproc ~{NPROC}

        # Remove low quality reads
        printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "FilterSeq quality"
        FilterSeq.py quality -s "~{assembly_pass}" -q 20 --outname ~{name} --log FS.log --nproc ~{NPROC}

        # Remove duplicate sequences
        printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "CollapseSeq"
        CollapseSeq.py -s "~{assembly_pass}" -n 20 --act set --inner --outname ~{name}
        printf "DONE\n\n"
    }

    runtime {
        docker: "immcantation/suite:devel"
    }

    output {
        File out = output_dir
    }
}