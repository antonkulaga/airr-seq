version development

import "https://raw.githubusercontent.com/antonkulaga/bioworkflows/main/common/files.wdl" as files


workflow sample {
    input {
        Array[File] reads
        String coord = "illumina"
        Array[Int] numbers
        String destination
        String name
    }

    call sample_pair {
        input: reads = reads, coord = coord, numbers = numbers, name = name
    }

    call files.copy as copy{
        input: destination = destination, files = [sample_pair.out]
    }
    output {
        File out = copy.out[0]
    }
}

task sample_pair {
    input {
        Array[File] reads
        String coord = "illumina"
        Array[Int] numbers
        String name
        String output_dir = "subsamples"
    }

    Array[String] basenames = [basename(reads[0]),basename(reads[1])]
    command {

        # Make output directory and empty log files
        mkdir -p ~{output_dir}; cd ~{output_dir}

        # Start
        echo "OUTPUT DIRECTORY: ~{output_dir}"
        echo -e "START"

        ln -s ~{reads[0]} ~{basenames[0]}
        ln -s ~{reads[1]} ~{basenames[1]}

        SplitSeq.py samplepair -1 ~{basenames[0]} -2  ~{basenames[1]} --coord ~{coord} -n ~{sep=" " numbers} --outname ~{name}
        rm ~{basenames[0]} ~{basenames[1]}
        echo -e "DONE"
    }

    runtime {
        docker: "immcantation/suite:devel"
    }

    output {
        File out = output_dir
    }
}
