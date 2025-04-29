version 1.0

task bam_filter_fixmates {
    input {
        File input_bam
        String output_prefix
        Int cpu = 2
        Int memory = 8
        String docker = "staphb/samtools:1.19"
    }

    command <<<
        set -euo pipefail

        # Copy original BAM
        cp "~{input_bam}" "~{output_prefix}_original.bam"
        samtools index "~{output_prefix}_original.bam"

        # Filter BAM: keep reads that are properly paired, OR unpaired reads with no first/second flags
        samtools view -h "~{input_bam}" | \
        samtools view -b -o "~{output_prefix}_filtered.bam" -h - \
        -F 192

        # Index the filtered BAM
        samtools index "~{output_prefix}_filtered.bam"

    >>>

    output {
        File filtered_bam = "~{output_prefix}_filtered.bam"
        File filtered_bai = "~{output_prefix}_filtered.bam.bai"
        File original_bam = "~{output_prefix}_original.bam"
        File original_bai = "~{output_prefix}_original.bam.bai"
    }

    runtime {
        docker: docker
        memory: "~{memory} GB"
        cpu: cpu
        disks: "local-disk 50 HDD"
    }
}
