version 1.0

task MarkDuplicates {

    input {
    File sorted_bam
    String sample_name

    Int mem_size_gb = 16
    String docker = "us-central1-docker.pkg.dev/gcid-bacterial/gcid-bacterial/fungi-gatk3:v1.0"
    Int cmd_mem_size_gb = mem_size_gb - 1
    Int disk_size = 50
    }

    command {
        set -euo pipefail

        java -Xmx${mem_size_gb}G -jar /opt/picard.jar MarkDuplicates \
            I=${sorted_bam} \
            O=${sample_name}.marked_duplicates.bam \
            M=${sample_name}.marked_duplicates.metrics
    }

    output {
        File bam = "${sample_name}.marked_duplicates.bam"
    }

    runtime {
        preemptible: 3
        docker: docker
        memory: mem_size_gb + " GB"
        disks: "local-disk " + disk_size + " HDD"
    }
}
