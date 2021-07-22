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
            reads = reads,
            primers = primers,
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
        String name
        Array[File] reads
        Array[File] primers #usually C than V
        Int NPROC = 4
    }


    command {
        # Make output directory and empty log files
        mkdir -p ~{output_dir}; cd ~{output_dir}

        ln -s ~{reads[0]} ~{basename(reads[0])}
        ln -s ~{reads[1]} ~{basename(reads[1])}

        # Start
        echo "OUTPUT DIRECTORY: ~{output_dir}"
        echo -e "START"
        
        FilterSeq.py quality -s ~{basename(reads[0])} -q 20 --outname ~{name}_R1 --log FS1.log
        FilterSeq.py quality -s ~{basename(reads[1])} -q 20 --outname ~{name}_R2 --log FS2.log
        MaskPrimers.py score -s ~{name}_R1_quality-pass.fastq -p ~{primers[0]} \
        --start 15 --mode cut --barcode --outname ~{name}_R1 --log MP1.log
        MaskPrimers.py score -s ~{name}_R2_quality-pass.fastq -p ~{primers[1]} \
        --start 0 --mode mask --outname ~{name}_R2 --log MP2.log
        PairSeq.py -1 ~{name}_R1_primers-pass.fastq -2 ~{name}_R2_primers-pass.fastq \
        --1f BARCODE --coord sra
        BuildConsensus.py -s ~{name}_R1_primers-pass_pair-pass.fastq --bf BARCODE --pf PRIMER \
        --prcons 0.6 --maxerror 0.1 --maxgap 0.5 --outname ~{name}_R1 --log BC1.log
        BuildConsensus.py -s ~{name}_R2_primers-pass_pair-pass.fastq --bf BARCODE --pf PRIMER \
        --maxerror 0.1 --maxgap 0.5 --outname ~{name}_R2 --log BC2.log
        PairSeq.py -1 ~{name}_R1_consensus-pass.fastq -2 ~{name}_R2_consensus-pass.fastq \
        --coord presto
        AssemblePairs.py align -1 ~{name}_R2_consensus-pass_pair-pass.fastq \
        -2 ~{name}_R1_consensus-pass_pair-pass.fastq --coord presto --rc tail \
        --1f CONSCOUNT --2f CONSCOUNT PRCONS --outname ~{name} --log AP.log
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
        rm ~{basename(reads[0])}
        rm ~{basename(reads[1])}
    }

    runtime {
        docker: "immcantation/suite:devel"
    }

    output {
        File out = output_dir
    }
}