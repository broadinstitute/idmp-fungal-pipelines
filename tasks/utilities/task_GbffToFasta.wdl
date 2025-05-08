task GbffToFasta {
    File ref_gbff
    String ref_gbff_basename = basename(ref_gbff, ".fasta")

    Int disk_size = 50
    Int mem_size_gb = 16
    String docker = "python:3.11-slim"

    command <<<
        set -euo pipefail
        pip install biopython
        python <<CODE
        from Bio import SeqIO
        input_file = "${ref_gbff}"
        output_file = "${ref_gbff_basename}.fasta"
        with open(output_file, "w") as out_fasta:
            for record in SeqIO.parse(input_file, "genbank"):
                 SeqIO.write(record, out_fasta, "fasta")
        CODE
        # strip the .1 in the reference sequence name
        sed -i -E 's/^(>[^ ]+)\.[0-9]+/\1/' "${ref_gbff_basename}.fasta"
    >>>
    output {
        File reference_fasta =  "${ref_gbff_basename}.fasta"
    }
    runtime {
        docker: docker
        memory: mem_size_gb + " GB"
        disks: "local-disk " + disk_size + " HDD"

    }
}