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
version 1.0

import "../tasks/alignment_preprocessing/task_MarkDuplicates.wdl" as MarkDuplicates
import "../tasks/alignment_preprocessing/task_ReorderBam.wdl" as ReorderBam
import "../tasks/gene_typing/variant_detection/task_HaplotypeCaller.wdl" as HaplotypeCaller
import "../tasks/gene_typing/variant_detection/task_CombineGVCFs.wdl" as CombineGVCFs
import "../tasks/gene_typing/variant_detection/task_GenotypeGVCFs.wdl" as GenotypeGVCFs
import "../tasks/gene_typing/variant_detection/task_HardFiltration.wdl" as HardFiltration


workflow FungalVariantCallingGatk3 {

    meta {
        allowNestedInputs: true
    }

    input {
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
    }

    # run pipeline on each sample, in parallel
    scatter(i in range(length(input_samples))) {
        String sample_name = input_samples[i]
        String input_bam = input_bams[i]

        call MarkDuplicates.MarkDuplicates as MarkDuplicates {
            input:
            sample_name = sample_name,
            sorted_bam = input_bam
        }

        call ReorderBam.ReorderBam as ReorderBam {
            input:
            bam = MarkDuplicates.bam,
            ref = reference_fasta,
            dict = ref_dict
        }

        call HaplotypeCaller.HaplotypeCaller as HaplotypeCaller {
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

    call CombineGVCFs.CombineGVCFs as CombineGVCFs {
        input:
        vcf_files = HaplotypeCaller.output_gvcf,
        vcf_index_files = HaplotypeCaller.output_gvcf_index,
        ref = reference_fasta,
        ref_dict = ref_dict,
        ref_index = ref_index
    }

    call GenotypeGVCFs.GenotypeGVCFs as GenotypeGVCFs {
        input:
        vcf_file = CombineGVCFs.out,
        vcf_index_file = CombineGVCFs.out_index,
        ref = reference_fasta,
        ref_dict = ref_dict,
        ref_index = ref_index
    }

    call HardFiltration.HardFiltration as HardFiltration {
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

