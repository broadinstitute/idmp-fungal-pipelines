version 1.0

task qc_flags {
  input {
    String? raw_read_screen
    String? clean_read_screen
    Float? est_coverage_clean
    String? kraken2_top_taxon_name
    String? gambit_predicted_taxon
    Float? eukcc_completeness
    Float? eukcc_contamination
    Float min_coverage = 40.0
    Float max_contamination = 5.0
    Float min_completeness = 80.0
    String docker = "marcoteix/gemstone-qc:1.0.0"
  }
  command <<<

    python /scripts/bin/cauris_row.py \
      -rs "~{raw_read_screen}" \
      -cs "~{clean_read_screen}" \
      -x ~{est_coverage_clean} \
      -g "~{gambit_predicted_taxon}" \
      -k "~{kraken2_top_taxon_name}" \
      -c ~{eukcc_contamination} \
      -C ~{eukcc_completeness} \
      -mx ~{min_coverage} \
      -Mc ~{max_contamination} \
      -mc ~{min_completeness} \
      -o "home/qc"

  >>>
  output {
    String qc_check = read_string("home/qc/qc_check")
    String qc_note = read_string("home/qc/qc_note")
  }
  runtime {
    docker: "~{docker}"
    memory: "8 GB"
    cpu: 1
    disks: "local-disk 10 SSD"
    disk: "10 GB"
    maxRetries: 0
    preemptible: 0
  }
}