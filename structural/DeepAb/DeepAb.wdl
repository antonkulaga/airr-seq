version development
#work in progress
workflow DeepAb {
    input {

    }
}

task predict {
    input {

    }


    runtime {
        python predict.py data/sample_files/4h0h.fasta --decoys 5 --renumber
    }

    output {

    }
}