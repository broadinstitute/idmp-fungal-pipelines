version 1.0

task EukCC {

    input {
        String docker = "us.gcr.io/broad-gotc-prod/eukcc:2.1.3"
        String memory = "16"
        String disk_size = "50"
        String cpu = "8"
        File assembly

    }
    command <<<

        mkdir eukccdb
        cd eukccdb
        wget http://ftp.ebi.ac.uk/pub/databases/metagenomics/eukcc/eukcc2_db_ver_1.1.tar.gz
        tar -xzvf eukcc2_db_ver_1.1.tar.gz
        export EUKCC2_DB=`pwd`/eukcc2_db_ver_1.1
        pwd
        cd ..

        eukcc single --out outfolder --threads 8 ~{assembly}

        tar -czvf outfolder.tar.gz outfolder

    >>>
    output {
        File outfolder = "outfolder.tar.gz"

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