version development

import "https://raw.githubusercontent.com/antonkulaga/bioworkflows/main/workflows/extract_run.wdl" as extractor

workflow immcantation {
    input {
        File sequence #fasta or fastq
        String species = "human"
        Int threads = 6
        String name
        Boolean ig = true
        Boolean only_functional = false
        Boolean partial_alignments = false
        Boolean include_tigger = false
        String format = "airr"
        String destination

    }

    call changeo_igblast {
        input: sequence = sequence,
            species = species,
            prefix = name,
            ig = ig,
            only_functional = only_functional,
            partial_alignments = partial_alignments,
            format = format
    }

    call extractor.copy as copy_changeo_igblast {
        input:
            files = [changeo_igblast.out], destination = destination + "/" + "changeo_igblast"
    }

    if (include_tigger && defined(changeo_igblast.airr_tsv)){
        File airr_tsv = select_first([changeo_igblast.airr_tsv])
         call shazam_threshold {
            input:
                airr_tsv = airr_tsv,
                name = name
        }

        call tigger_genotype {
            input:
                airr_tsv = airr_tsv,
                name = name
        }

        call extractor.copy as copy_tigger {
            input:
                files = [tigger_genotype.out], destination = destination + "/" + "tigger"
        }

        call changeo_clone {
            input:
                airr_tsv = airr_tsv,
                name = name,
                threshold_tsv = shazam_threshold.threshold_tsv
        }

        call extractor.copy as copy_changeo_clone {
            input:
                files = [changeo_clone.out], destination = destination + "/" + "changeo_clone"
        }

    }

    output {
        Array[File] changeo = copy_changeo_igblast.out
        Array[File] clones = select_first([copy_changeo_clone.out])
    }


}


task changeo_igblast {
    input {
            File sequence
            String species
            Boolean ig
            String prefix
            Boolean only_functional
            Boolean partial_alignments
            String format
        }

    #changeo-igblast -h
    #Usage: changeo-igblast [OPTIONS]
    #-s  FASTA or FASTQ sequence file.
    #-r  Directory containing IMGT-gapped reference germlines.
    #Defaults to /usr/local/share/germlines/imgt/human/vdj when species is human.
    #Defaults to /usr/local/share/germlines/imgt/mouse/vdj when species is mouse.
    #-g  Species name. One of human or mouse. Defaults to human.
    #-t  Receptor type. One of ig or tr. Defaults to ig.
    #-b  IgBLAST IGDATA directory, which contains the IgBLAST database, optional_file
    #and auxillary_data directories. Defaults to /usr/local/share/igblast.
    #-n  Sample identifier which will be used as the output file prefix.
    #Defaults to a truncated version of the sequence filename.
    #-o  Output directory. Will be created if it does not exist.
    #Defaults to a directory matching the sample identifier in the current working directory.
    #-f  Output format. One of airr (default) or changeo. Defaults to airr.
    #-p  Number of subprocesses for multiprocessing tools.
    #Defaults to the available cores.
    #-k  Specify to filter the output to only productive/functional sequences.
    #-i  Specify to allow partial alignments.
    #-h  This message.
    command {
        changeo-igblast -s ~{sequence} -g ~{species} -t ~{if(ig) then "ig" else "tr"}  -n ~{prefix} -f ~{format}  ~{if(only_functional) then "-k" else ""} ~{if(partial_alignments) then "-i" else ""}
    }

    runtime {
        docker: "immcantation/suite:devel"
    }

    output {
        File out = prefix
        File? airr_tsv = prefix + "/" + prefix+"_db-pass.tsv"
        File? fmt7 = prefix + "/" + prefix+"_igblast.fmt7"
    }
}

task shazam_threshold {
    input {
        File airr_tsv
        String name
        String method = "density"
        String outdir = "shazam"
        String format = "airr"
        String model = "gamma-gamma"
    }

    #    -d DB, --db=DB
    #    Tabulated data file, in Change-O (TAB) or AIRR format (TSV).
    #    -m METHOD, --method=METHOD
    #    Threshold inferrence to use. One of gmm, density, or none.
    #    If none, the distance-to-nearest distribution is plotted without threshold detection.
    #    Defaults to density.
    #    -n NAME, --name=NAME
    #    Sample name or run identifier which will be used as the output file prefix.
    #    Defaults to a truncated version of the input filename.
    #    -o OUTDIR, --outdir=OUTDIR
    #    Output directory. Will be created if it does not exist.
    #    Defaults to the current working directory.
    #    -f FORMAT, --format=FORMAT
    #    File format. One of 'airr' (default) or 'changeo'.
    #    -p NPROC, --nproc=NPROC
    #    Number of subprocesses for multiprocessing tools.
    #    Defaults to the available processing units.
    #    --model=MODEL
    #    Model to use for the gmm model.
    #    One of gamma-gamma, gamma-norm, norm-norm or norm-gamma.
    #    Defaults to gamma-gamma.
    #    --subsample=SUBSAMPLE
    #    Number of distances to downsample the data to before threshold calculation.
    #    By default, subsampling is not performed.
    #    --repeats=REPEATS
    #    Number of times to recalculate.
    #    Defaults to 1.
    command {
        shazam-threshold -d ~{airr_tsv} --method ~{method} --name ~{name} --format ~{format} --model ~{model} -o ~{outdir}
    }

    runtime {
        docker: "immcantation/suite:devel"
    }

    output {
        File out = outdir
        Array[Array[String]] threshold_tsv = read_tsv(outdir + "/" + name+ "_threshold-values.tab")
    }
}

task changeo_clone {

    input {
        File airr_tsv
        String name
        Array[Array[String]] threshold_tsv
        #Int distance_threshold
        String distance_model = "ham"
        String format = "airr"
        Boolean only_functional = false
        Boolean full_dataset = true
    }
    String distance_threshold = threshold_tsv[1][1]


        #        Usage: changeo-clone [OPTIONS]
        #        -d  Change-O formatted TSV (TAB) file.
        #        -x  Distance threshold for clonal assignment.
        #        -m  Distance model for clonal assignment.
        #        Defaults to the nucleotide Hamming distance model (ham).
        #        -r  Directory containing IMGT-gapped reference germlines.
        #        Defaults to /usr/local/share/germlines/imgt/human/vdj.
        #        -n  Sample identifier which will be used as the output file prefix.
        #        Defaults to a truncated version of the input filename.
        #        -o  Output directory. Will be created if it does not exist.
        #        Defaults to a directory matching the sample identifier in the current working directory.
        #        -f  Output format. One of airr (default) or changeo.
        #        -p  Number of subprocesses for multiprocessing tools.
        #        Defaults to the available cores.
        #        -a  Specify to clone the full data set.
        #        By default the data will be filtering to only productive/functional sequences.
        #        -h  This message.
    # ~{if(full_dataset) then "-a" else ""}
    command {
        changeo-clone -d ~{airr_tsv} -x ~{distance_threshold} -m ~{distance_model} -n ~{name} -f ~{format}
    }

    runtime {
        docker: "immcantation/suite:devel"
    }

    output {
        File out = name
        File germ_pass_tsv =  name + "/"+ name + "_germ-pass.tsv"
    }
}


task tigger_genotype {
    input {
        File airr_tsv #/data/changeo/sample/sample_db-pass.tab
        String name
    }

    #         Usage: /usr/local/bin/tigger-genotype [options]
    #
    #
    #        Options:
    #        -d DB, --db=DB
    #        Change-O formatted TSV (TAB) file.

    #        -r REF, --ref=REF
    #        FASTA file containing IMGT-gapped V segment reference germlines.
    #        Defaults to /usr/local/share/germlines/imgt/human/vdj/imgt_human_IGHV.fasta.

    #        -v VFIELD, --vfield=VFIELD
    #        Name of the output field containing genotyped V assignments.
    #        Defaults to V_CALL_GENOTYPED.

    #        -x MINSEQ, --minseq=MINSEQ
    #        Minimum number of sequences in the mutation/coordinate range.
    #        Samples with insufficient sequences will be excluded.
    #        Defaults to 50.

    #       -y MINGERM, --mingerm=MINGERM
    #        Minimum number of sequences required to analyze a germline allele.
    #        Defaults to 200.

    #        -n NAME, --name=NAME
    #        Sample name or run identifier which will be used as the output file prefix.
    #        Defaults to a truncated version of the input filename.

    #        -o OUTDIR, --outdir=OUTDIR
    #        Output directory. Will be created if it does not exist.
    #        Defaults to the current working directory.

    #        -f FORMAT, --format=FORMAT
    #        File format. One of 'airr' (default) or 'changeo'.

    #        -p NPROC, --nproc=NPROC
    #        Number of subprocesses for multiprocessing tools.
    #        Defaults to the available processing units.


    command {
        mkdir tigger_genotype
        tigger-genotype -d ~{airr_tsv} -n ~{name} -o tigger_genotype
    }

    runtime {
        docker: "immcantation/suite:devel"
    }
    output {
        File out = "tigger_genotype"
    }
}
