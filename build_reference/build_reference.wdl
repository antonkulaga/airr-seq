version development

import "https://raw.githubusercontent.com/antonkulaga/bioworkflows/main/common/files.wdl" as files

#see https://changeo.readthedocs.io/en/stable/examples/igblast.html for more info

workflow build_igblast {
    input {
        String destination
    }

    call download_reference{
        input:
    }
    call files.copy as copy_reference{
        input:
        files = [download_reference.igblast, download_reference.imgt], destination = destination
    }

    output {
        File igblast = copy_reference.out[0]
        File imgt = copy_reference.out[1]
   }
}

task download_reference {
    input {

    }
    command {
        # Download reference databases
        fetch_igblastdb.sh -o igblast
        fetch_imgtdb.sh -o imgt
        # Build IgBLAST database from IMGT reference sequences
        imgt2igblast.sh -i imgt -o igblast
    }
    runtime {
        docker: "immcantation/suite:devel"
    }

    output {
        File igblast = "igblast"
        File imgt = "imgt"
    }

}
