version 1.0

task HardFiltration {

    input {
    # hard-filter a vcf, if vqsr not available
    # http://gatkforums.broadinstitute.org/gatk/discussion/2806/howto-apply-hard-filters-to-a-call-set
    File ref
    File vcf
    File vcf_index
    File ref_dict
    File ref_index
    String output_filename

    String snp_filter_expr
    String indel_filter_expr

    Int disk_size_gb = ceil(size(vcf, "GiB") * 2) + 20
    Int mem_size_gb = ceil(size(vcf, "GiB") * 2.5) + 10
    String docker = "us-central1-docker.pkg.dev/gcid-bacterial/gcid-bacterial/fungi-gatk3:v1.0"

    }

    command {
        # select snps
        java -Xmx~{mem_size_gb - 1}G -jar /opt/GenomeAnalysisTK.jar \
            -T SelectVariants \
            -R ~{ref} \
            -V ~{vcf} \
            -selectType SNP \
            -o raw_snps.g.vcf

        # filter snps
        java -Xmx~{mem_size_gb}G -jar /opt/GenomeAnalysisTK.jar \
            -T VariantFiltration \
            -R ~{ref} \
            -V ~{vcf} \
            --filterExpression "~{snp_filter_expr}" \
            --filterName "snp_filter" \
            -o filtered_snps.g.vcf

        # select indels
        java -Xmx~{mem_size_gb}G -jar /opt/GenomeAnalysisTK.jar \
           -T SelectVariants \
           -R ~{ref} \
           -V ~{vcf} \
           -selectType INDEL \
           -o raw_indels.g.vcf

        # filter indels
        java -Xmx~{mem_size_gb}G -jar /opt/GenomeAnalysisTK.jar \
            -T VariantFiltration \
            -R ~{ref} \
            -V ~{vcf} \
            --filterExpression "~{indel_filter_expr}" \
            --filterName "indel_filter" \
            -o filtered_indels.g.vcf

        # combine variants
        java -Xmx~{mem_size_gb}G -jar /opt/GenomeAnalysisTK.jar\
            -T CombineVariants \
            -R ~{ref} \
            --variant filtered_snps.g.vcf \
            --variant filtered_indels.g.vcf \
            -o ~{output_filename} \
            --genotypemergeoption UNSORTED
    }

    output {
        File out = "~{output_filename}"
    }

    runtime {
        preemptible: 3
        docker: docker
        memory: mem_size_gb + " GB"
        disks: "local-disk " + disk_size_gb + " HDD"
    }
}
