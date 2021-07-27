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

        Boolean revpr = true
        Float maxerror = 0.2
        Float maxgap = 0.5

    }
    call extract_umi {
        input:
            reads = reads, start = 3, len = 6, coord = coord, name = name, maxerror = maxerror, maxgap = maxgap, constant_region = constant_region,
            revpr = revpr
    }

    call files.copy as copy_umi_extraction {
        input: destination = destination + "/" + name + "/" + "presto" , files = [extract_umi.out]
    }

    call assemble_pairs {
        input:
            reads = if(flip_reads) then [extract_umi.consensus_reads[1], extract_umi.consensus_reads[0]] else extract_umi.consensus_reads, coord = "presto",
            name = name,
            NPROC = threads

    }

    call files.copy as copy_assemble_pairs {
        input: destination = destination + "/" + name + "/" + "presto" , files = [assemble_pairs.out]
    }

    call collapse {
        input :
            assembly = assemble_pairs.assembly, name = name, num  = min_reads_per_ig, n = 20
    }

    call files.copy as copy_collapse {
        input: destination = destination + "/" + name + "/" + "presto" , files = [collapse.out]
    }


    output {
        File out = destination + "/" + name
        #File collapse_results = copy_collapse.out
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
        Float maxgap = 0.5

        Boolean revpr
        String mode = "cut"

        Int maxlen = 100

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


        echo "Mask constant regions"
        MaskPrimers.py align -s ~{name}_R1_quality-pass_pair-pass.fastq -p ~{constant_region} --maxlen ~{maxlen} --maxerror ~{maxerror} --mode ~{mode} ~{if(revpr) then "--revpr" else ""} --pf C_CALL

        echo "Build consensus"

        BuildConsensus.py -s ~{name}_R1_quality-pass_pair-pass.fastq --bf BARCODE --maxerror ~{maxerror} --maxgap ~{maxgap} --outname ~{name}_R1 --log BC1.log
        BuildConsensus.py -s ~{name}_R2_quality-pass_pair-pass.fastq --bf BARCODE --maxerror ~{maxerror} --maxgap ~{maxgap} --outname ~{name}_R2 --log BC2.log

        echo "pairing consensus sequences together and converting coordinates to presto"

        PairSeq.py -1 ~{name}_R1_consensus-pass.fastq -2 ~{name}_R2_consensus-pass.fastq --coord presto
        rm  ~{read_1_name}
        rm  ~{read_2_name}
    }

    output {
        File out = output_dir
        Array[File] consensus_reads = [
                                      output_dir + "/" + name + "_R1_consensus-pass_pair-pass.fastq",
                                      output_dir + "/" + name + "_R2_consensus-pass_pair-pass.fastq"
                                      ]
        Array[File] logs_fs = [output_dir + "/" + "FS1.log", output_dir + "/" + "FS2.log"]
        Array[File] logs_bc = [output_dir + "/" + "BC1.log", output_dir + "/" + "BC2.log"]
    }

    runtime {
        docker: "immcantation/suite:devel"
    }
}



task assemble_pairs {
    input {
        String output_dir = "assemble_pairs_results"
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
        File assembly = output_dir + "/" + name + "_assemble-pass.fastq"
        File log_AP = output_dir + "/" + "AP.log"
        File log_FS = output_dir + "/" + "FS.log"
    }
}


task collapse {
    input {
        File assembly
        String name
        String output_dir = "collapse_results"
        Int n = 20
        Int num = 2
        String uf = "PRCONS"
    }

    command {
        # Make output directory and empty log files
        mkdir -p ~{output_dir}; cd ~{output_dir}

        ln -s ~{assembly} ~{basename(assembly)}

        # Start
        echo "OUTPUT DIRECTORY: ~{output_dir}"
        echo -e "START"
        STEP=0

        ParseHeaders.py collapse -s ~{basename(assembly)} -f CONSCOUNT --act min
        CollapseSeq.py -s ~{name}_assemble-pass_reheader.fastq -n ~{n} --inner --uf ~{uf} --cf CONSCOUNT --act sum --outname ~{name}
        SplitSeq.py group -s ~{name}_collapse-unique.fastq -f CONSCOUNT --num ~{num} --outname ~{name}
        rm ~{basename(assembly)}
    }

    runtime {
        docker: "immcantation/suite:devel"
    }

    output {
        File out = output_dir
        File unique = output_dir + "/" + name + "_collapse-unique.fastq"
    }
}