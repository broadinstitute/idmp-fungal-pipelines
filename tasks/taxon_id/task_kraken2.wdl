version 1.0

task kraken2 {

    input {
        File read1
        File read2
        String samplename
        String kraken2_db_path = "https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_16gb_20250402.tar.gz"
        Int cpu = 4
        Int memory = 32
        String docker = "marcoteix/bracken:1.0.0"
    }
    command <<<
        # Create database directory and download/extract database
        mkdir -p db
        wget -O kraken2_db.tar.gz ~{kraken2_db_path}
        tar -C ./db/ -xzvf kraken2_db.tar.gz

        # Run Kraken2
        kraken2 --paired \
        --db ./db/ \
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
        disks: "local-disk 50 HDD"
        preemptible: 0
    }
}
