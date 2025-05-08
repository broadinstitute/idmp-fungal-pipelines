task CombineGVCFs {
    File ref
    File ref_dict
    File ref_index
    Array[File] vcf_files
    Array[File] vcf_index_files

    String gvcf_out = "combined_gvcfs.vcf.gz"
    String gvcf_out_index = "combined_gvcfs.vcf.gz.tbi"

    Int disk_size = 100
    Int mem_size_gb = 30
    String docker = "xiaoli2020/fungi-gatk3:v1.0"
    Int cmd_mem_size_gb = mem_size_gb - 1

    command {
        java -Xmx${cmd_mem_size_gb}G -jar /opt/GenomeAnalysisTK.jar \
        -T CombineGVCFs \
        -R ${ref} \
        -o ${gvcf_out} \
        --variant ${sep=" --variant " vcf_files}
    }
    output {
        File out = gvcf_out
        File out_index = gvcf_out_index
    }

    runtime {
        preemptible: 4
        docker: docker
        memory: mem_size_gb + " GB"
        disks: "local-disk " + disk_size + " HDD"
    }
}
