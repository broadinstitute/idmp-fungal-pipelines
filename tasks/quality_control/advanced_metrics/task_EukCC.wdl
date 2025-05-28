version 1.0

task EukCC {

    input {
        String docker = "us-central1-docker.pkg.dev/gcid-bacterial/gcid-bacterial/eukcc:1.0.0"
        String memory = "16"
        String disk_size = "50"
        String cpu = "8"
        String? eukcc_db_path
        File assembly
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
            mkdir -p ./db
            wget -O ./db/eukcc_db.tar.gz $eukcc_db_path
            tar -C ./db/ -xzvf ./db/eukcc_db.tar.gz
            rm ./db/eukcc_db.tar.gz

            # Find the top-level directory that was extracted and set EUKCC2_DB to its absolute path
            export EUKCC2_DB=$(find ./db -mindepth 1 -maxdepth 1 -type d | head -n 1)
            echo "EUKCC2_DB set to $EUKCC2_DB"
        fi

        # Run EukCC2
        eukcc single --out outfolder --threads ~{cpu} ~{assembly}

        # Extract completeness and contamination from eukcc.csv (columns 2 and 3)
        awk -F'\t' 'NR==2 {print $2}' outfolder/eukcc.csv > COMPLETENESS
        awk -F'\t' 'NR==2 {print $3}' outfolder/eukcc.csv > CONTAMINATION

    >>>
    output {
        File eukcc_csv = "outfolder/eukcc.csv"
        Float completeness = read_float("COMPLETENESS")
        Float contamination = read_float("CONTAMINATION")
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