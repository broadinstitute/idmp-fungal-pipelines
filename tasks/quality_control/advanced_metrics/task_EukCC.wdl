version 1.0

task EukCC {

    input {
        String docker = "us.gcr.io/broad-gotc-prod/eukcc/eukcc_1.0.0:latest"
        String memory = "16"
        String disk_size = "50"
        String cpu = "8"
        String? eukcc_db_path
        File assembly
        Float contamination_percent_threshold

    }
    command <<<
        set -euo pipefail

        # Determine which EukCC DB to use.
        # If eukcc_db_path is provided, use it; otherwise use the baked-in one (http://ftp.ebi.ac.uk/pub/databases/metagenomics/eukcc/eukcc2_db_ver_1.1.tar.gz)
        DB_PATH="~{if defined(eukcc_db_path) then eukcc_db_path else "/app/db/"}"

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