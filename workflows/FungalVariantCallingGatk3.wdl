## Copyright Broad Institute 2020
##
## Variant Calling pipeline for fungal haploid genomes
## Developed by Xiao Li (xiaoli@broadinstitute.org)
## Fungal Genomics Group, Infectious Disease and Microbiome Program.
## The Broad Institute of MIT and Harvard
##
## Verions of the software used in this WDL:
##   PICARD_VER=1.782
##   GATK37_VER=3.7-93-ge9d8068
##   SAMTOOLS_VER=1.3.1
##   BWA_VER=0.7.12
##   TABIX_VER=0.2.5_r1005
##   BGZIP_VER=1.3
##
## Cromwell version support
## - Successfully tested on v28, v30 and v36
## - Does not work on versions < v23 due to output syntax
##
## LICENSING :
## This script is released under the WDL source code license (BSD-3) (see LICENSE in
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the dockers
## for detailed licensing information pertaining to the included programs.


workflow FungalVariantCallingGatk3 {

    meta {
        description: ""
        allowNestedInputs: "true"
    }

    ## config params
    # input data
    String analysis_name
    Array[String] input_samples
    Array[File] input_bams
    File reference_fasta
    File ref_dict
    File ref_index

    # hard filtering params: both of these params are required
    String snp_filter_expr
    String indel_filter_expr



    # run pipeline on each sample, in parallel
    scatter(i in range(length(input_samples))) {
        String sample_name = input_samples[i]
        String input_bam = input_bams[i]

        call MarkDuplicates {
            input:
            sample_name = sample_name,
            sorted_bam = input_bam
        }

        call ReorderBam {
            input:
            bam = MarkDuplicates.bam,
            ref = reference_fasta,
            dict = ref_dict
        }

        call HaplotypeCaller {
            input:
            input_bam = ReorderBam.out,
            input_bam_index = ReorderBam.out_index,
            sample_name = sample_name,
            gvcf_name = "${sample_name}.g.vcf",
            gvcf_index = "${sample_name}.g.vcf.idx",
            ref = reference_fasta,
            ref_dict = ref_dict,
            ref_index = ref_index
        }
    }

    call CombineGVCFs {
        input:
        vcf_files = HaplotypeCaller.output_gvcf,
        vcf_index_files = HaplotypeCaller.output_gvcf_index,
        ref = reference_fasta,
        ref_dict = ref_dict,
        ref_index = ref_index
    }

    call GenotypeGVCFs {
        input:
        vcf_file = CombineGVCFs.out,
        vcf_index_file = CombineGVCFs.out_index,
        ref = reference_fasta,
        ref_dict = ref_dict,
        ref_index = ref_index
    }

    call HardFiltration {
        input:
        vcf = GenotypeGVCFs.output_vcf_name,
        vcf_index = GenotypeGVCFs.output_vcf_index_name,
        snp_filter_expr = snp_filter_expr,
        indel_filter_expr = indel_filter_expr,
        ref = reference_fasta,
        ref_dict = ref_dict,
        ref_index = ref_index,
        output_filename = "${analysis_name}.hard_filtered.vcf.gz"
    }

    output {
        File hard_filtered_gvcf = HardFiltration.out

    }
}


# mark duplicate reads in bam
task MarkDuplicates {
    File sorted_bam
    String sample_name

    Int mem_size_gb = 16
    String docker = "xiaoli2020/fungi-gatk3:v1.0"
    Int cmd_mem_size_gb = mem_size_gb - 1
    Int disk_size = 50

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


# reorder and index a bam
task ReorderBam {
    File ref
    File dict
    File bam
    String bam_prefix = basename(bam, '.bam')

    Int disk_size = 50
    Int mem_size_gb = 30
    String docker = "xiaoli2020/fungi-gatk3:v1.0"

    Int cmd_mem_size_gb = mem_size_gb - 1

    command {
        # reorder bam
        java -Xmx${cmd_mem_size_gb}G -jar /opt/picard.jar ReorderSam \
            I=${bam} \
            O=${bam_prefix}.reordered.bam \
            R=${ref}

        # then index
        java -Xmx${cmd_mem_size_gb}G -jar /opt/picard.jar BuildBamIndex \
            I=${bam_prefix}.reordered.bam
    }

    output {
        File out = "${bam_prefix}.reordered.bam"
        File out_index = "${bam_prefix}.reordered.bai"
    }

    runtime {
        preemptible: 3
        docker: docker
        memory: mem_size_gb + " GB"
        disks: "local-disk " + disk_size + " HDD"
    }
}


# merge vcfs before genotyping
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


# genotype gvcfs
task GenotypeGVCFs {
    File ref
    File ref_dict
    File ref_index
    File vcf_file
    File vcf_index_file
    String vcf_basename = basename(vcf_file, ".vcf.gz")


    Int disk_size = 100
    Int mem_size_gb = 30
    String docker = "xiaoli2020/fungi-gatk3:v1.0"
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

    Int disk_size = 200
    Int mem_size_gb = 60
    String docker = "xiaoli2020/fungi-gatk3:v1.0"

    Int cmd_mem_size_gb = mem_size_gb - 1

    command {
        # select snps
        java -Xmx${cmd_mem_size_gb}G -jar /opt/GenomeAnalysisTK.jar \
            -T SelectVariants \
            -R ${ref} \
            -V ${vcf} \
            -selectType SNP \
            -o raw_snps.g.vcf

        # filter snps
        java -Xmx${cmd_mem_size_gb}G -jar /opt/GenomeAnalysisTK.jar \
            -T VariantFiltration \
            -R ${ref} \
            -V ${vcf} \
            --filterExpression "${snp_filter_expr}" \
            --filterName "snp_filter" \
            -o filtered_snps.g.vcf

        # select indels
        java -Xmx${cmd_mem_size_gb}G -jar /opt/GenomeAnalysisTK.jar \
            -T SelectVariants \
            -R ${ref} \
            -V ${vcf} \
            -selectType INDEL \
            -o raw_indels.g.vcf

        # filter indels
        java -Xmx${cmd_mem_size_gb}G -jar /opt/GenomeAnalysisTK.jar \
            -T VariantFiltration \
            -R ${ref} \
            -V ${vcf} \
            --filterExpression "${indel_filter_expr}" \
            --filterName "indel_filter" \
            -o filtered_indels.g.vcf

        # combine variants
        java -Xmx${cmd_mem_size_gb}G -jar /opt/GenomeAnalysisTK.jar\
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
        memory: mem_size_gb + " GB"
        disks: "local-disk " + disk_size + " HDD"
    }
}


# run haplotype caller
task HaplotypeCaller {
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

  String docker = "xiaoli2020/fungi-gatk3:v1.0"

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
