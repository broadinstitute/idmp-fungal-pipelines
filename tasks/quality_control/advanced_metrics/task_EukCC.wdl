version 1.0

task EukCC {

    input {
        String docker = "us.gcr.io/broad-gotc-prod/eukcc:2.1.3"
        String memory = "8"
        String disk_size = "50"
        String cpu = "2"
        File assembly

    }
    command <<<

        mkdir eukccdb
        cd eukccdb
        wget http://ftp.ebi.ac.uk/pub/databases/metagenomics/eukcc/eukcc2_db_ver_1.1.tar.gz
        tar -xzvf eukcc2_db_ver_1.1.tar.gz
        pwd
        cd ..

        eukcc single --out outfolder --threads 8 ~{assembly}

    >>>
    output {
        File outfolder = "outfolder"

    }
    runtime {
        docker: "~{docker}"
        memory: "~{memory} GB"
        cpu: cpu
        disks: "local-disk " + disk_size + " SSD"
        disk: disk_size + " GB"
        preemptible: 1
    }
}