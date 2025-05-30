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

        python3 <<CODE

        from pathlib import Path

        coverage = (
        float("~{est_coverage_clean}")
        if "~{est_coverage_clean}" != ""
        else 0
        )

        completeness = (
        float("~{eukcc_completeness}")
        if "~{eukcc_completeness}" != ""
        else 0.0
        )

        contamination = (
        float("~{eukcc_contamination}")
        if "~{eukcc_contamination}" != ""
        else 100.0
        )

        kraken_taxon = "~{kraken2_top_taxon_name}"
        gambit_taxon = "~{gambit_predicted_taxon}"

        min_coverage = float("~{min_coverage}")
        max_contamination = float("~{max_contamination}")
        min_completeness = float("~{min_completeness}")

        qc_check, qc_note = "PASS", ""
        if "~{raw_read_screen}" != "PASS" or "~{clean_read_screen}" != "PASS":
            # If the isolate fails the raw or clean read QC, it should fail global QC
            qc_check = "FAIL"
            qc_note = "Low yield/quality"
        elif coverage < min_coverage:
            # If the isolate does not meet the minimum coverage, it should fail global QC
            qc_check = "FAIL"
            qc_note = "Low coverage"
        else:
            if contamination > max_contamination:
                qc_check = "FAIL"
                qc_note = "Contamination"
            elif completeness < min_completeness:
                # Set incomplete samples to "FAIL"
                qc_check = "FAIL"
                qc_note = "Low completeness"    
            elif (
                gambit_taxon != "Candidozyma auris" or \
                kraken2_taxon != "Candidozyma auris"
            ): 
                # If it is not Candidozyma auris, fail QC
                qc_check = "ALERT"
                qc_note = "Taxonomic mismatch"

        # Write outputs
        print("QC check: ", qc_check)
        print("QC note: ", qc_note)

        with open("qc_check", "w") as f:
        f.write(qc_check)
        with open("qc_note", "w") as f:
        f.write(qc_note)

        CODE

    >>>
    output {
        String qc_check = read_string("qc_check")
        String qc_note = read_string("qc_note")
    }
    runtime {
        docker: docker
        memory: "8 GB"
        cpu: 1
        disks: "local-disk 10 SSD"
        disk: "10 GB"
        maxRetries: 1
        preemptible: 0
    }
}