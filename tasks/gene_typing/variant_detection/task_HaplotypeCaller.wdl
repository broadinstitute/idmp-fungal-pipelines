version 1.0

task HaplotypeCaller {

    input {
    File input_bam
    File input_bam_index

    String gvcf_name
    String gvcf_index

    File ref_dict
    File ref
    File ref_index
    String sample_name

    Int disk_size = 250
    Int mem_size_gb = 30
    Int cmd_mem_size_gb = mem_size_gb - 1

    String docker = "us-central1-docker.pkg.dev/gcid-bacterial/gcid-bacterial/fungi-gatk3:v1.0"
    }
    command {
        java -Xmx${cmd_mem_size_gb}G -jar /opt/GenomeAnalysisTK.jar \
            -T HaplotypeCaller \
            -R ${ref} \
            -I ${input_bam} \
            -o ${gvcf_name} \
            -ERC "GVCF" \
            -ploidy 1 \
            -variant_index_type LINEAR \
            -variant_index_parameter 128000 \
            --read_filter OverclippedRead
    }

    output {
        File output_gvcf = "${gvcf_name}"
        File output_gvcf_index = "${gvcf_index}"
    }

    runtime {
        task_name: "HaplotypeCaller"
        preemptible: 3
        docker: docker
        memory: mem_size_gb + " GB"
        disks: "local-disk " + disk_size + " HDD"
    }

    parameter_meta {
        ref: "fasta file of reference genome"
        sample_name: "The name of the sample as indicated by the 1st column of the gatk.samples_file json input."
        input_bam: "The bam file to call HaplotypeCaller on."
        output_gvcf: "VCF file produced by haplotype caller."
    }
}
