task GenotypeGVCFs {
    File ref
    File ref_dict
    File ref_index
    File vcf_file
    File vcf_index_file
    String vcf_basename = basename(vcf_file, ".vcf.gz")


    Int disk_size = 100
    Int mem_size_gb = 30
    String docker = "us-central1-docker.pkg.dev/gcid-bacterial/gcid-bacterial/fungi-gatk3:v1.0"
    Int cmd_mem_size_gb = mem_size_gb - 1

    command {
        java -Xmx${cmd_mem_size_gb}G -jar /opt/GenomeAnalysisTK.jar \
            -T GenotypeGVCFs \
            -R ${ref} \
            -o ${vcf_basename}.vcf.gz \
            --variant ${vcf_file}
    }
    output {
        File output_vcf_name = "${vcf_basename}.vcf.gz"
        File output_vcf_index_name = "${vcf_basename}.vcf.gz.tbi"
    }

    runtime {
        preemptible: 4
        docker: docker
        memory: mem_size_gb + " GB"
        disks: "local-disk " + disk_size + " HDD"
    }
}
