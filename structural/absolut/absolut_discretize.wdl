version development

import "https://raw.githubusercontent.com/antonkulaga/bioworkflows/main/common/files.wdl" as files


workflow absolut_discretize {
    input {
        String pdb
        String chains
        Float resolution = 5.25
        String pos = "FuC"
        String destination
    }

    call absolut_discretize {
        input: pdb = pdb, chains = chains, resolution = resolution, pos = pos
    }

    call files.copy as copy{
        input: destination = destination,
            files = [absolut_discretize.pdb_de_ins,absolut_discretize.dpb_prepared, absolut_discretize.discretize, absolut_discretize.lattice]
    }

    output {
        File out = copy.destination_folder
    }
}


task absolut_discretize {
    input {
        String pdb
        String chains
        Float resolution = 5.25
        String pos = "FuC"
    }


    command {
        Absolut discretize ~{pdb} ~{chains} ~{resolution} ~{pos}
    }

    runtime {
        docker: "quay.io/comp-bio-aging/absolut:latest"
    }

    output {
      File pdb_de_ins = pdb+ "deIns.pdb" # PDB with removed insertions (using pdb-tools)
      File dpb_prepared = pdb + "_VWprepared.pdb" # new PDB with only the chains of interest
      File discretize = pdb + "discretized" + resolution +pos + ".pdb"
      File lattice = pdb +"_VWInLattice.txt"  #Description of the discretized (lattice) antigen [Each chain is described as a starting position in the lattice (6-digits number) and a list of moves in space (straight S, up U, down D, left L, right R). See ‘info_position’ to convert lattice positions.
    }
}