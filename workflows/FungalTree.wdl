version 1.0

import "../workflows/FungalVariantCallingGatk3.wdl" as FungalVariantCallingGatk3
import "../tasks/utilities/task_GbffToFasta.wdl" as GbffToFasta
import "../tasks/utilities/task_GenerateRefFiles.wdl" as GenerateRefFiles
import "../tasks/utilities/task_VCFToFasta.wdl" as VCFToFasta
import "../tasks/species_typing/task_IqTree2.wdl" as IqTree2


workflow FungalTree {

    meta {
        description: "FungalTree is a WDL-based pipeline for variant calling and phylogenetic analysis in fungal haploid genomes, using a reference GenBank file and aligned BAM files to generate filtered variant calls and a maximum-likelihood phylogenetic tree via IQ-TREE2."
        allowNestedInputs: "true"
    }

    ## config params
    # input data
    input {

    String analysis_name
    File ref_gbff
    Array[String] input_samples
    Array[File] input_bams

    # hard filtering params: both of these params are required
    String snp_filter_expr
    String indel_filter_expr

    # IQ-TREE2 parameters
    String? iqtree2_model
    Int iqtree2_bootstraps = 1000
    Int alrt = 1000
    String? iqtree2_opts
    }

    call GbffToFasta.GbffToFasta as GbffToFasta {
        input:
            ref_gbff = ref_gbff
    }

    call GenerateRefFiles.GenerateRefFiles as GenerateRefFiles{
        input:
            ref_fasta = GbffToFasta.reference_fasta,
            input_bam = input_bams[0]
    }

    call FungalVariantCallingGatk3.FungalVariantCallingGatk3 as FungalVariantCallingGatk3 {
    input:
           analysis_name = analysis_name,
           input_samples = input_samples,
           input_bams = input_bams,
           reference_fasta = GenerateRefFiles.reference_fasta,
           ref_dict = GenerateRefFiles.ref_dict,
           ref_index = GenerateRefFiles.ref_index,
           snp_filter_expr = snp_filter_expr ,
           indel_filter_expr = indel_filter_expr
    }


    call VCFToFasta.VCFToFasta as VCFToFasta {
        input:
        vcf_file = FungalVariantCallingGatk3.hard_filtered_gvcf
    }

    call IqTree2.IqTree2 as IqTree2 {
        input:
        alignment = VCFToFasta.alignment_fasta,
        iqtree2_model = iqtree2_model,
        iqtree2_bootstraps = iqtree2_bootstraps,
        alrt = alrt,
        iqtree2_opts = iqtree2_opts,
        cluster_name = analysis_name
    }

    output {
        File hard_filtered_gvcf = FungalVariantCallingGatk3.hard_filtered_gvcf
        File IqTree2_ml_tree = IqTree2.ml_tree
        String IqTree2_run_date = IqTree2.date
        String IqTree2_model_used = IqTree2.iqtree2_model_used
    }
}

