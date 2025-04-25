version 1.0

task kraken2 {

    input {
        File read1
        File read2
        String samplename
        String? kraken2_db_path
        #String kraken2_db_path = "https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_16gb_20250402.tar.gz"
        Int cpu = 4
        Int memory = 32
        String docker = "us.gcr.io/broad-gotc-prod/kraken2/kraken2_1.0.0"
    }
    command <<<
        kraken2_db_path=~{kraken2_db_path}

        # Determine Kraken2 DB path
        # If no kraken2_db_path is provided, then use the db baked into the docker (https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_16gb_20250402.tar.gz)
        if [ -z "$kraken2_db_path" ]; then
            echo "kraken2_db_path is empty...Using the db baked into the docker"
            DB_PATH="/app/db"
        else
            # If kraken2_db_path is provided, use it
            echo "kraken2_db_path is not empty..."
            echo "Downloading Kraken2 database..."
            mkdir -p /app/db && \
            wget -O /app/db/kraken2_db.tar.gz $kraken2_db_path
            tar -C /app/db/ -xzvf /app/db/kraken2_db.tar.gz && \
            rm /app/db/kraken2_db.tar.gz
            DB_PATH="/app/db"
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
        grep -P "^\s+\d+\.\d+\s+\d+\s+\d+\s+S\s+\d+\s+" ~{samplename}.kraken2.report.txt | sort -k1,1nr | head -n 1 | awk '{for(i=6;i<=NF;i++) printf "%s ", $i; print ""}' | sed 's/ $//' > TOP_TAXON_NAME

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
        disks: "local-disk 250 HDD"
        preemptible: 0
    }
}
