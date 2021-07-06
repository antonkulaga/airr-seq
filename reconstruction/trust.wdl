version development
import "https://raw.githubusercontent.com/antonkulaga/bioworkflows/main/common/files.wdl" as files

workflow trust {
    input {
        File? bam
        File? fastq_1
        File? fastq_2
        Int threads = 12
        File imgt_reference
        String destination
    }

    call trust4 {
        input: bam = bam, fastq_1 = fastq_1, fastq_2 = fastq_2, threads = threads, fasta = imgt_reference
    }


    call files.copy as copy_trust {
        input:
            files = [trust4.out],
            destination =  destination,
    }

    output {
        File out = copy_trust.out[0]
    }
}

task trust4 {

    input {
        File fasta # -f STRING: path to the fasta file coordinate and sequence of V/D/J/C genes
        File? bam
        File? fastq_1
        File? fastq_2
        Int threads
        Boolean abnormalUnmapFlag = false
    }

    #
    #Required:
    #    -b STRING: path to bam file
    #    -1 STRING -2 STRING: path to paired-end read files
    #    -u STRING: path to single-end read file
    #    -f STRING: path to the fasta file coordinate and sequence of V/D/J/C genes
    #
    # Optional:
    #        --ref STRING: path to detailed V/D/J/C gene reference file from IMGT database. (default: not used but recommended)
    #        -o STRING: prefix of output files. (default: inferred from file prefix)
    #        --od STRING: the directory for output files. (default: ./)
    #        -t INT: number of threads (default: 1)
    #        --barcode STRING: if -b, bam field for barcode; if -1 -2/-u, file containing barcodes (default: not used)
    #        --barcodeRange INT INT CHAR: start, end(-1 for length-1), strand in a barcode is the true barcode (default: 0 -1 +)
    #        --barcodeWhitelist STRING: path to the barcode whitelist (default: not used)
    #        --read1Range INT INT: start, end(-1 for length-1) in -1/-u files for genomic sequence (default: 0 -1)
    #        --read2Range INT INT: start, end(-1 for length-1) in -2 files for genomic sequence (default: 0 -1)
    #        --mateIdSuffixLen INT: the suffix length in read id for mate. (default: not used)
    #        --skipMateExtension: do not extend assemblies with mate information, useful for SMART-seq (default: not used)
    #        --abnormalUnmapFlag: the flag in BAM for the unmapped read-pair is nonconcordant (default: not set)
    #        --noExtraction: directly use the files from provided -1 -2/-u to assemble (default: extraction first)
    #        --repseq: the data is from TCR-seq or BCR-seq (default: not set)
    #        --stage INT: start TRUST4 on specified stage (default: 0):
    #                0: start from beginning (candidate read extraction)
    #                1: start from assembly
    #                2: start from annotation
    #                3: start from generating the report table
    #


    command {
        run-trust4 --od trust_results  -f ~{fasta} -t ~{threads} ~{"-b " + bam} ~{if(defined(fastq_2)) then "-1 " + fastq_1 else "-u " + fastq_1}  ~{"-2 " + fastq_2} ~{if(abnormalUnmapFlag) then "--abnormalUnmapFlag" else ""}
    }

    runtime {
        docker: "quay.io/biocontainers/trust4:1.0.4--h2e03b76_0"
    }

    output {
        File out = "trust_results"
        #File raw = "TRUST_" +"_raw.out"
        #File final = "TRUST_" +"_final.out"
        #File fasta_annotated = "TRUST_" +"_annot.fa"
        #File cdr3 = "TRUST_" +"_cdr3.out"
        #File out = "TRUST_" +"_report.tsv"

        #TRUST4 outputs several files. trust_raw.out, trust_final.out are the contigs and corresponding nucleotide weight.
        #trust_annot.fa is in fasta format for the annotation of the consensus assembly.
        #trust_cdr3.out reports the CDR1,2,3 and gene information for each consensus assemblies.
        #And trust_report.tsv is a report file focusing on CDR3 and is compatible with other repertoire analysis tool such as VDJTools.
    }

}
