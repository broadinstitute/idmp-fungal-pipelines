version 1.0

task GenotypeGVCFs {

    input {
    File ref
    File ref_dict
    File ref_index
    File vcf_file
    File vcf_index_file
    String vcf_basename = basename(vcf_file, ".vcf.gz")

    String docker = "us-central1-docker.pkg.dev/gcid-bacterial/gcid-bacterial/fungi-gatk3:v1.0"
    Int disk_size_gb = ceil(size(vcf_file, "GiB") * 2) + 20
    Int mem_size_gb = ceil(size(vcf_file, "GiB") * 2.5) + 10
    }

    command {
        java -Xmx~{mem_size_gb - 1}G -jar /opt/GenomeAnalysisTK.jar \
            -T GenotypeGVCFs \
            -R ~{ref} \
            -o ~{vcf_basename}.vcf.gz \
            --variant ~{vcf_file}
    }
    output {
        File output_vcf_name = "~{vcf_basename}.vcf.gz"
        File output_vcf_index_name = "~{vcf_basename}.vcf.gz.tbi"
    }

    runtime {
        preemptible: 4
        docker: docker
        memory: mem_size_gb + " GB"
        disks: "local-disk " + disk_size_gb + " HDD"
    }
}
