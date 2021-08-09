version development

import "https://raw.githubusercontent.com/antonkulaga/bioworkflows/main/common/files.wdl" as files

#production version
import "https://raw.githubusercontent.com/antonkulaga/airr-seq/main/analysis/changeo_igblast.wdl" as changeo
import "https://raw.githubusercontent.com/antonkulaga/airr-seq/main/analysis/clonal_analysis.wdl" as clonal


workflow irepertoire{

    input {
        Array[File] reads
        String destination
        String name

        String species = "human"

        String coordinates = "illumina" #sra
        Int threads = 12
        Int min_reads_per_ig = 2
        Int min_length = 196
        Int min_quality = 20
        Int start_v = 10
        Int min_dupcount = 2
        Boolean ig = true
        Boolean only_functional = false
        Boolean partial_alignments = true
        String format = "airr"
        #File constant
        #File variable

        #Float maxerror = 0.2
        #Float maxgap = 0.5

    }

    call presto{
        input: name = name, output_dir = "presto",
        reads = reads, NPROC = threads, min_quality  = min_quality, min_length = min_length, start_v = start_v, dupcount = min_dupcount,
        coordinates = coordinates
    }

    call files.copy as copy_presto{ input: destination = destination, files = [presto.results] }
    call changeo.changeo_igblast as igblast{
        input: fastq = presto.out, threads = threads, name = name, destination = destination
    }

    call clonal.clonal_analysis as clonal_analysis{
        input:
          airr_tsv = igblast.airr_tsv, name = name, distance_model = "ham", format = "airr",
          threshold_model = "gamma-gamma", shazam_method = "density", destination = destination, only_functional = only_functional
    }
    #call imm.translate as translate {
    #    input: files = [igblast.airr_tsv, clonal_analysis.clones],suffix = "_with_translation"
    #}

    output {
        File presto_results = copy_presto.out[0]
        File changeo_results = igblast.out
        #File clonal = clonal_analysis.out
    }
    #File out = prefix
    #File? airr_tsv = prefix + "/" + prefix+"_db-pass.tsv"
    #File? airr_fail_tsv = prefix + "/" + prefix+"_db-fail.tsv"
    #File? fmt7 = prefix + "/" + prefix+"_igblast.fmt7"
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
        Int start_v = 0
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
        MaskPrimers.py extract -s ~{basenames[1]}.fastq --start ~{start_v} ~{if(start_v>0) then "--barcode --bf V_barcode" else ""} --len ~{variable_length} --mode cut --pf VPRIMER
        PairSeq.py -1 ~{basenames[1]}_primers-pass.fastq -2 ~{basenames[0]}_primers-pass.fastq  --1f VPRIMER ~{if(start_v>0) then "V_barcode" else ""} --2f CPRIMER

        # Assemble paired ends via mate-pair alignment
        printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "AssemblePairs align"

        AssemblePairs.py align -1 ~{basenames[1]}_primers-pass_pair-pass.fastq -2 ~{basenames[0]}_primers-pass_pair-pass.fastq \
          --coord ~{coordinates} --rc tail --2f CPRIMER VPRIMER ~{if(start_v>0) then "V_barcode" else ""}  \
          --outname ~{name}  --outdir . --log AP.log --nproc ~{NPROC}

        # Remove low quality reads
        printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "FilterSeq quality"
        FilterSeq.py quality -s "~{name}_assemble-pass.fastq" -q ~{min_quality} --outname ~{name} --log FS_quality.log --nproc ~{NPROC}
        FilterSeq.py length -s "~{name}_quality-pass.fastq" -n ~{min_length} --outname ~{name} --log FS_length.log --nproc ~{NPROC}

        # Remove duplicate sequences
        printf "  %2d: %-*s $(date +'%H:%M %D')\n" $((++STEP)) 24 "CollapseSeq"
        CollapseSeq.py -s "~{name}_length-pass.fastq" -n ~{collapse_max_missing} --uf CPRIMER --cf VPRIMER --act set --inner --outname ~{name}
        SplitSeq.py group -s ~{name}_collapse-unique.fastq -f DUPCOUNT --num ~{dupcount} --outname ~{name}
        ParseHeaders.py table -s "~{name}_atleast-~{dupcount}.fastq" -f ID DUPCOUNT CPRIMER VPRIMER ~{if(start_v>0) then "V_barcode" else ""}
        rm -rf ~{basenames[0]}.fastq ~{basenames[1]}.fastq
    }

    runtime {
        docker: "immcantation/suite:devel"
    }

    output {
        File results = output_dir
        File out = output_dir + "/" + name + "_atleast-" + dupcount + ".fastq"
        File headers = output_dir + "/" + name + "_atleast-" + dupcount + "_headers.tab"
    }
}