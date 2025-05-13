task GenotypeGVCFs {
    File ref
    File ref_dict
    File ref_index
    File vcf_file
    File vcf_index_file
    String vcf_basename = basename(vcf_file, ".vcf.gz")


    #Int disk_size = 100
    #Int mem_size_gb = 30
    Int disk_size_gb = ceil(size(vcf_file, "GiB") * 2) + 10
    Int memory_mb = ceil(size(vcf_file, "MiB") * 2.5) + 40000
    String docker = "xiaoli2020/fungi-gatk3:v1.0"
    Int cmd_mem_size_mb = memory_mb - 1000

    command {
        java -Xmx${cmd_mem_size_mb}M -jar /opt/GenomeAnalysisTK.jar \
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
        memory: memory_mb + " MiB"
        disks: "local-disk " + disk_size_gb + " HDD"
    }
}
