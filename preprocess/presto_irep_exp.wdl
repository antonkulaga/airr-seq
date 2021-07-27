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
        Boolean flip_reads = false
        Int min_reads_per_ig = 2
        File constant_region

    }
    call extract_umi {
        input:
            reads = reads, start = 3, len = 6, coord = coord, name = name,  constant_region = constant_region
    }

    call files.copy as copy_umi_extraction {
        input: destination = destination + "/" + name + "/" + "presto" , files = [extract_umi.out]
    }
}

task extract_umi {
    input {
        Int start
        Int len
        String name
        File constant_region
        Array[File] reads
        String coord
        String output_dir = "extract_umi_results"

        Float maxerror = 0.3

        String mode = "cut"

        Int maxlen = 120

    }
    # echo "Extracting UMI"
    # MaskPrimers.py extract -s ~{name}_R1_umi_quality-pass.fastq --start ~{start} --len ~{len} --pf BARCODE -o ~{name}_R1_quality-pass.fastq

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
        
        echo "Mask constant regions R1 REVP=false"
        MaskPrimers.py align -s ~{name}_R1_quality-pass_pair-pass.fastq -p ~{constant_region} --maxlen ~{maxlen} --maxerror ~{maxerror} --mode ~{mode}  --pf C_CALL

        echo "Mask constant regions R2 REVP=false"
        MaskPrimers.py align -s ~{name}_R2_quality-pass_pair-pass.fastq -p ~{constant_region} --maxlen ~{maxlen} --maxerror ~{maxerror} --mode ~{mode}  --pf C_CALL

        echo "Mask constant regions R1 REVP=true"
        MaskPrimers.py align -s ~{name}_R1_quality-pass_pair-pass.fastq -p ~{constant_region} --maxlen ~{maxlen} --maxerror ~{maxerror} --mode ~{mode} --revpr --pf C_CALL

        echo "Mask constant regions R2 REVP=true"
        MaskPrimers.py align -s ~{name}_R2_quality-pass_pair-pass.fastq -p ~{constant_region} --maxlen ~{maxlen} --maxerror ~{maxerror} --mode ~{mode} --revpr --pf C_CALL


    }

    output {
        File out = output_dir
    }

    runtime {
        docker: "immcantation/suite:devel"
    }
}