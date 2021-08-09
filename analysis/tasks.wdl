version development

# tasks that did not get into specific pipelines

task changeo_igblast {
    input {
        File sequence
        String species
        Boolean ig
        String prefix
        Boolean only_functional
        Boolean partial_alignments
        String format
        Boolean avoid_temp = true
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
        changeo-igblast -s ~{sequence} -g ~{species} \
        -r /usr/local/share/germlines/imgt/human/vdj -b /usr/local/share/igblast \
        -t ~{if(ig) then "ig" else "tr"}  -n ~{prefix} ~{if(avoid_temp) then "-z" else ""} \
        -f ~{format}  ~{if(only_functional) then "-k" else ""} ~{if(partial_alignments) then "-i" else ""}
    }

    runtime {
        docker: "immcantation/suite:devel"
    }

    output {
        File out = prefix
        File airr_tsv = prefix + "/" + prefix+"_db-pass.tsv"
        File? airr_fail_tsv = prefix + "/" + prefix+"_db-fail.tsv"
        File? fmt7 = prefix + "/" + prefix+"_igblast.fmt7"
    }
}

task split_chains {
    input {
        File airr
    }

    String link = basename(airr) # to avoid parsedb creating stuff in the input folder

    command {
        ln -s ~{airr} ~{link}
        ParseDb.py select -d ~{link} -f v_call -u "IGH" --logic all --regex --outname heavy
        mv heavy_parse-select.tsv heavy.tsv
        ParseDb.py split -d heavy.tsv -f productive
        ParseDb.py select -d ~{link} -f v_call -u "IG[LK]" --logic all --regex --outname light
        mv light_parse-select.tsv light.tsv
        ParseDb.py split -d light.tsv -f productive
    }

    runtime {
        docker: "immcantation/suite:devel"
    }

    output {
        File light = "light.tsv"
        File heavy = "heavy.tsv"
        File light_productive = "light_productive-T.tsv"
        File light_not_productive = "light_productive-F.tsv"
        File heavy_productive = "heavy_productive-T.tsv"
        File heavy_not_productive = "heavy_productive-F.tsv"
    }
}

task translate {
    input {
        Array[File] files
        String suffix = "_with_translation"
    }
    String gl = "*"+suffix+".tsv"

    command {
        /home/magus/notebooks/translate.R --wd TRUE --suffix ~{suffix} ~{sep=" " files}
    }

    runtime {
        docker: "quay.io/comp-bio-aging/immcantation"
    }
    output {
        Array[File] out = glob(gl)
    }
}


# not really used yet, as there is a separate immcantation pipeline that has both germlines and define_clones
task define_clones {
    input {
        File airr
        String model = "ham" #--model {ham,aa,hh_s1f,hh_s5f,mk_rs1nf,mk_rs5nf,hs1f_compat,m1n_compat}
        String format = "airr"
        Float distance_threshold = 0.0
    }

    String link = basename(airr) # to avoid parsedb creating stuff in the input folder
    String prefix = basename(airr, ".tsv")+"_"

    command {
        ln -s ~{airr} ~{link}
        DefineClones.py -d ~{link} --dist ~{distance_threshold} --format ~{format} --model ~{model}
    }

    runtime {
        docker: "immcantation/suite:devel"
    }

    output
    {
        File out = prefix + "clone-pass.tsv"
        #File failed = prefix + "clone-fail.tsv"
    }
}