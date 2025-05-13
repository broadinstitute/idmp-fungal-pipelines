
task HardFiltration {
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

    Int memory_mb = ceil(size(vcf, "MiB") * 2.5) + 4000
    Int disk_gb = ceil(size(vcf, "GiB") * 2) + 5
    String docker = "us-central1-docker.pkg.dev/gcid-bacterial/gcid-bacterial/fungi-gatk3:v1.0"


    Int cmd_mem_size_mb = memory_mb - 1000

    command {
        # select snps
        java -Xmx${cmd_mem_size_mb}M -jar /opt/GenomeAnalysisTK.jar \
            -T SelectVariants \
            -R ${ref} \
            -V ${vcf} \
            -selectType SNP \
            -o raw_snps.g.vcf

        # filter snps
        java -Xmx${cmd_mem_size_mb}M -jar /opt/GenomeAnalysisTK.jar \
            -T VariantFiltration \
            -R ${ref} \
            -V ${vcf} \
            --filterExpression "${snp_filter_expr}" \
            --filterName "snp_filter" \
            -o filtered_snps.g.vcf

        # select indels
        java -Xmx${cmd_mem_size_mb}M -jar /opt/GenomeAnalysisTK.jar \
           -T SelectVariants \
           -R ${ref} \
           -V ${vcf} \
           -selectType INDEL \
           -o raw_indels.g.vcf

        # filter indels
        java -Xmx${cmd_mem_size_mb}M -jar /opt/GenomeAnalysisTK.jar \
            -T VariantFiltration \
            -R ${ref} \
            -V ${vcf} \
            --filterExpression "${indel_filter_expr}" \
            --filterName "indel_filter" \
            -o filtered_indels.g.vcf

        # combine variants
        java -Xmx${cmd_mem_size_mb}M -jar /opt/GenomeAnalysisTK.jar\
            -T CombineVariants \
            -R ${ref} \
            --variant filtered_snps.g.vcf \
            --variant filtered_indels.g.vcf \
            -o ${output_filename} \
            --genotypemergeoption UNSORTED
    }

    output {
        File out = "${output_filename}"
    }

    runtime {
        preemptible: 3
        docker: docker
        memory: memory_mb + " MiB"
        disks: "local-disk " + disk_gb + " HDD"
    }
}
