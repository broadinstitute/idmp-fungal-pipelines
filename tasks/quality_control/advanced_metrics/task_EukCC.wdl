version 1.0

task EukCC {

    input {
        String docker = "us.gcr.io/broad-gotc-prod/eukcc/eukcc_1.0.0:latest"
        String memory = "16"
        String disk_size = "200"
        String cpu = "8"
        String? eukcc_db_path
        File assembly
        Float contamination_percent_threshold

    }
    command <<<
        set -euo pipefail

        eukcc_db_path=~{eukcc_db_path}

        # Determine EukCC DB path
        # If no eukcc_db_path is provided, then use the db baked into the docker (http://ftp.ebi.ac.uk/pub/databases/metagenomics/eukcc/eukcc2_db_ver_1.1.tar.gz)
        if [ -z "$eukcc_db_path" ]; then
            echo "eukcc_db_path is empty...Using the db baked into the docker"
        else
            # If eukcc_db_path is provided, use it
            echo "eukcc_db_path is not empty..."
            echo "Downloading EukCC database..."
            mkdir -p /app/db && \
            wget -O /app/db/eukcc_db.tar.gz $eukcc_db_path
            tar -C /app/db/ -xzvf /app/db/eukcc_db.tar.gz && \
            rm /app/db/eukcc_db.tar.gz
            # Find the top-level directory that was extracted and set EUKCC2_DB to its absolute path
            export EUKCC2_DB=$(find /app/db -mindepth 1 -maxdepth 1 -type d | head -n 1)
            echo "EUKCC2_DB set to $EUKCC2_DB"
        fi

        # Run EukCC2
        eukcc single --out outfolder --threads ~{cpu} ~{assembly}

        # Extract contamination percentage from 3rd column, skipping header
        contamination=$(awk -F'\t' 'NR==2 {print $3}' outfolder/eukcc.csv)

        # Determine if contamination is too high (≥ threshold)
        contamination_too_high=$(awk "BEGIN {print ($contamination >= ~{contamination_percent_threshold})}")

        if [[ "$contamination_too_high" == "1" ]]; then
            echo "The contamination level too high: ${contamination}"
            exit 1
        fi


        tar -czvf EukCC_output_folder.tar.gz outfolder
    >>>
    output {
        File EukCC_output_folder = "EukCC_output_folder.tar.gz"

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