version development
import "https://raw.githubusercontent.com/antonkulaga/bioworkflows/main/common/files.wdl" as files

import "https://raw.githubusercontent.com/antonkulaga/airr-seq/main/structural/anarci/anarci.wdl" as anarci

workflow anarci_airr {
    input {
        String destination
        File airr
        String sequence_id_field = "sequence_id"
        String sequence_field = "Translation"
        String species = "human"
        String scheme = "martin"

        String name
        Int threads = 12
        String? chains_restriction
        Array[String] meta_fields = ["duplicate_count", "clone_id"]
        Boolean assign_germline = true
    }

    call airr_to_fasta {
        input: airr_tsv = airr, sequence_id_field = sequence_id_field, sequence_field = sequence_field, name = name,  meta_fields =  meta_fields
    }

    call files.copy as copy_fasta{
        input: files=[airr_to_fasta.out], destination = destination
          }

    File created_fasta = copy_fasta.out[0] + "/" + basename(airr, ".tsv") + "_sequences.fasta"

    call anarci.anarci as anarci{
        input: destination = destination, fasta = created_fasta,
            species = species, scheme = scheme,
           name = name, threads = threads,
           restrict = chains_restriction,
           assign_germline = assign_germline
    }

    output{
        File fasta = copy_fasta.out[0]
        File out = anarci.out
    }

}


task airr_to_fasta {
    input {
        String name
        File airr_tsv
        String sequence_field
        String sequence_id_field = "sequence_id"
        String output_prefix = "clonotypes_as_fasta"
        Array[String] meta_fields
    }
    String output_dir =  output_prefix + "_"  + name

    command {
        # Make output directory
        mkdir -p ~{output_dir}
        ConvertDb.py fasta -d ~{airr_tsv} --outdir ~{output_dir} --if ~{sequence_id_field} --sf ~{sequence_field} ~{if(length(meta_fields)>0) then "--mf" else ""} ~{sep=" " meta_fields}
    }

    runtime {
        docker: "immcantation/suite:devel"
    }

    output {
        File out = output_dir
    }
}