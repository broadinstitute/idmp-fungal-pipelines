version 1.0

task bam_filter_fixmates {
    input {
        File input_bam
        String output_prefix
        Int cpu = 2
        Int memory = 8
        String docker = "us-central1-docker.pkg.dev/gcid-bacterial/gcid-bacterial/samtools:1.19"
    }

    command <<<
        set -euo pipefail

        # Copy original BAM
        cp "~{input_bam}" "~{output_prefix}_original.bam"
        samtools index "~{output_prefix}_original.bam"

        # Filter BAM: keep reads that are properly paired
        samtools view -h -f 2 -b "~{output_prefix}_original.bam" -o "~{output_prefix}_filtered.bam"

        # Index the filtered BAM
        samtools index "~{output_prefix}_filtered.bam"

    >>>

    output {
        File filtered_bam = "~{output_prefix}_filtered.bam"
        File filtered_bai = "~{output_prefix}_filtered.bam.bai"
    }

    runtime {
        docker: docker
        memory: "~{memory} GB"
        cpu: cpu
        disks: "local-disk 50 HDD"
    }
}
