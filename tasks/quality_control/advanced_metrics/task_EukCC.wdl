version 1.0

task EukCC {

    input {
        String docker = "us.gcr.io/broad-gotc-prod/eukcc:2.1.3"
        String memory = "16"
        String disk_size = "50"
        String cpu = "8"
        File assembly
        Float contamination_percent_threshold

    }
    command <<<
        set -euo pipefail

        # Load the EukCC2 database
        mkdir eukccdb
        cd eukccdb
        wget http://ftp.ebi.ac.uk/pub/databases/metagenomics/eukcc/eukcc2_db_ver_1.1.tar.gz
        tar -xzvf eukcc2_db_ver_1.1.tar.gz
        export EUKCC2_DB=`pwd`/eukcc2_db_ver_1.1
        cd ..

        # Run EukCC2
        eukcc single --out outfolder --threads ~{cpu} ~{assembly}

        # Extract contamination percentage from 3rd column, skipping header
        contamination=$(awk -F'\t' 'NR==2 {print $3}' outfolder/eukcc.csv)

        # Fail if contamination > 5%
        if awk "BEGIN {exit !($contamination > ~{contamination_percent_threshold})}"; then
            echo "Contamination level too high: ${contamination}"
            exit 1
        fi

        tar -czvf outfolder.tar.gz outfolder

    >>>
    output {
        File outfolder = "outfolder.tar.gz"

    }
    runtime {
        docker: "~{docker}"
        memory: "~{memory} GB"
        cpu: cpu
        disks: "local-disk " + disk_size + " HDD"
        disk: disk_size + " GB"
        preemptible: 1
    }
}