version 1.0

task kraken2 {

    input {
        File read1
        File read2
        String samplename
        String? kraken2_db_path
        Int cpu = 4
        Int memory = 32
        String docker = "us-central1-docker.pkg.dev/gcid-bacterial/gcid-bacterial/kraken2:1.0.0"
    }
    command <<<
        set -euo pipefail
        kraken2_db_path=~{kraken2_db_path}

        # Determine Kraken2 DB path
        # If no kraken2_db_path is provided, then use the db baked into the docker (https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_16gb_20250402.tar.gz)
        if [ -z "$kraken2_db_path" ]; then
          echo "kraken2_db_path is empty...Using the db baked into the docker"
          DB_PATH="/app/db"
        else
          # If kraken2_db_path is provided, download and extract it into the working directory
          echo "kraken2_db_path is not empty..."
          echo "Downloading Kraken2 database into working directory..."
          mkdir -p ./db
          wget -O ./db/kraken2_db.tar.gz "$kraken2_db_path"
          tar -C ./db/ -xzvf ./db/kraken2_db.tar.gz
          rm ./db/kraken2_db.tar.gz
          DB_PATH="./db"
          echo "DB_PATH set to $DB_PATH"
        fi

        # Run Kraken2
        kraken2 --paired \
        --db $DB_PATH \
        --threads ~{cpu} \
        --report-zero-counts \
        --report ~{samplename}.kraken2.report.txt \
        --output ~{samplename}.kraken2.classified_reads.txt \
        ~{read1} ~{read2}

        cat ~{samplename}.kraken2.report.txt

        # Extract top hit from report - taking species with highest abundance
        # Format of report: percent abundance,  clade reads  #reads_taxon  rank  taxid  sci_name
        # Filters for lines where the rank is "S", i.e., species level, sorts numerically in reverse (% abundance), take the species with highest percent abundance

        awk '
          $4 == "S" {
            percent = $1
            name = ""
            for (i = 6; i <= NF; i++) {
              name = name $i " "
            }
            gsub(/ $/, "", name)
            if (percent + 0 > max_percent + 0) {
              max_percent = percent
              top_name = name
            }
          }
          END {
            if (top_name != "") print top_name
          }
        ' ~{samplename}.kraken2.report.txt > TOP_TAXON_NAME
      >>>
    output {
        File kraken2_report = "~{samplename}.kraken2.report.txt"
        File classified_reads = "~{samplename}.kraken2.classified_reads.txt"
        String kraken2_report_taxon_name = read_string("TOP_TAXON_NAME")
    }
    runtime {
        docker: docker
        memory: "~{memory} GB"
        cpu: cpu
        disks: "local-disk 50 HDD"
        preemptible: 0
    }
}
